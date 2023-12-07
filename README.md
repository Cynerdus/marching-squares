# HW 1 APD / PDA

First Homework of 'Parallel and Distributed Algorithms'.
Featuring a performance boost for the Marching Squares Algorithm.

## Optimization Approach

The task foresaw the usage of 2 to 4 threads for a parallel curvature mapping,
which led to the conclusion that barriers must be used.
Other means of 'traffic' flow controller were not needed (mutex, semaphores etc).

For starters, the ```main``` function of the program has been fragmented,
and separated as such on the much philosophical criterion of
		*'Will this piece still work if I assign it to minions?'*

Indeed a question of all times it is, which led to the parallelization of
three out of five functions, respectively ```rescale_image```, ```sample_grid``` and ```march```.

These functions have been carefully moved (totally not leaking in terms of memory),
to an auxiliary function handler, ```f```. Now, in order to perform this,
a structure was implemented, which would cover all the necessary values and pointers
needed in the computation.

The structure used, ```params```, covers as such:
- Thread ID (*tid*);
- Thread count (*P*);
- Pointer to barrier (*barrier*);
- Pointers to the original image map (*image*), the scaled image map (*scaled_image*),
the contour map (*contour_map*) and the grid map (*grid*).

Pointers have been used, so that the threads will know where to refer
when passing the needed arguments to the computing functions. A worthy mention
would be my first attempt at this task, which included separate memory
management for each thread, resulting in extra steps of joining the
result, upon calling back and joining the threads.
Dreadful experience it was, leaving me unwilling of life.

## Feedback

Calendaristic mapping: 2 days.
Actual work time: 5 hours, distributed evenly between sunset, midnight and early morning.

My thoughts: An enjoyable homework, allowing me to enhance my understanding of the
concepts discussed in the first laboratory sessions. 10/10 would recommend.

##### Copyright Popescu Iulia Cleopatra 333CA
