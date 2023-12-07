// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

#define MIN(a, b) ((a < (b)) ? (a) : (b))

// Structure for a thread's parameters.
typedef struct _params {
    int tid;                            // thread ID
    int P;                              // thread count
    pthread_barrier_t *barrier;         // pointer to the barrier

    ppm_image **image, **scaled_image;  // pointers to the images
    ppm_image ***contour_map;           // pointer to the contour map
    unsigned char ***grid;              // pointer to the grid map
} params;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(unsigned char **grid, ppm_image *image, int step_x, int step_y, unsigned char sigma, int thread_count, int thread_id, pthread_barrier_t *barrier) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    // Parallelization by indexes: start and end hold the limits for each thread,
    // values dependent on the thread count, thread ID and the max value.

    int start, end;
    start = thread_id * p / thread_count;
    end = MIN((thread_id + 1) * p / thread_count, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    pthread_barrier_wait(barrier);

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    pthread_barrier_wait(barrier);

    start = thread_id * q / thread_count;
    end = MIN((thread_id + 1) * q / thread_count, q);   

    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    pthread_barrier_wait(barrier);

    return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y, int thread_count, int thread_id) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start, end;
    start = (thread_id * p) / thread_count;
    end = MIN(((thread_id + 1) * p) / thread_count, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

ppm_image *rescale_image(ppm_image *new_image, ppm_image *image, int thread_count, int thread_id, pthread_barrier_t *barrier) {
    uint8_t sample[3];

    // we only rescale downwards
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    int start, end;
    start = (thread_id * new_image->x) / thread_count;
    end = MIN(((thread_id + 1) * new_image->x) / thread_count, new_image->x);

    // use bicubic interpolation for scaling
    for (int i = start; i < end; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }

    pthread_barrier_wait(barrier);

    return new_image;
}

void *f(void *args) {
    params *p = (params *) args;

    // 1. Rescale the image
    (*(p->scaled_image)) = rescale_image((*(p->scaled_image)), (*(p->image)), p->P, p->tid, p->barrier);

    pthread_barrier_wait(p->barrier);

    // 2. Sample the grid
    (*(p->grid)) = sample_grid((*(p->grid)), (*(p->scaled_image)), STEP, STEP, SIGMA, p->P, p->tid, p->barrier);

    pthread_barrier_wait(p->barrier);

    // 3. March the squares
    march((*(p->scaled_image)), (*(p->grid)), (*(p->contour_map)), STEP, STEP, p->P, p->tid);
    
    pthread_barrier_wait(p->barrier);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    // vvvvvvvvvvvvvvvv
    // This code block has been moved from their original position
    // So that each thread refers to the same memory zone.
    // vvvvvvvvvvvvvvvv

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // alloc memory for image
    ppm_image *scaled_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!scaled_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    scaled_image->x = RESCALE_X;
    scaled_image->y = RESCALE_Y;

    scaled_image->data = (ppm_pixel*)malloc(scaled_image->x * scaled_image->y * sizeof(ppm_pixel));
    if (!scaled_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    unsigned char **grid = (unsigned char **)malloc((scaled_image->x / step_x + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= scaled_image->x / step_x; i++) {
        grid[i] = (unsigned char *)malloc((scaled_image->y / step_y + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // ^^^^^^^^^^^^^^^^
    // This code block has been moved from their original position
    // So that each thread refers to the same memory zone.
    // ^^^^^^^^^^^^^^^^

    int r;
    int P = atoi(argv[3]);

    pthread_t threads[P];
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, P);

    params p[P];

    for (int i = 0; i < P; i++) {
        p[i].image = &image;
        p[i].scaled_image = &scaled_image;
        p[i].contour_map = &contour_map;
        p[i].grid = &grid;
        p[i].tid = i;
        p[i].P = P;
        p[i].barrier = &barrier;

        r = pthread_create(&threads[i], NULL, f, &p[i]);
        if (r) {
            printf("Eroare la crearea thread-ului %d\n", i);
            exit(-1);
        }
    }

    for (int i = 0; i < P; i++) {
        r = pthread_join(threads[i], NULL);
        if (r) {
            printf("Eroare la asteptarea thread-ului %d\n", i);
            exit(-1);
        }
    }
    pthread_barrier_destroy(&barrier);

    // 4. Write output
    write_ppm(scaled_image, argv[2]);

    free_resources(scaled_image, contour_map, grid, step_x);

    return 0;
}
