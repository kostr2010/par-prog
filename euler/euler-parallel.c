#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

long double factorial(int n, int depth) {
    if (n <= 0) {
        return 1;
    }

    long double fac = 1;

    for (int i = n; i > ((depth != 0) ? (n - depth) : (0)); i--) {
        fac *= i;
    }

    return fac;
}

int main(int argc, char** argv) {
    clock_t t_begin = clock();

    assert(argc >= 2);
    long double res = 0.0;

    MPI_Init(&argc, &argv);

    int steps = atoi(argv[argc - 1]);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int lower_bound = 0;
    int upper_bound = 0;

    lower_bound = 1 + my_rank * (steps / world_size);
    upper_bound = (my_rank == world_size - 1) ? (steps) : ((my_rank + 1) * (steps / world_size));

    printf("process %d calculates [%d, %d]\n", my_rank, lower_bound, upper_bound);

    long double part_res = 0;
    for (int i = lower_bound; i <= upper_bound; i++) {
        part_res += (long double)factorial(steps, steps - i);
    }

    MPI_Reduce(&part_res, &res, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    res /= (long double)factorial(steps, 0);

    clock_t t_end = clock();
    double elapsed_time = (double)(t_end - t_begin) / CLOCKS_PER_SEC;

    if (my_rank == 0) {
        printf("euler:      %.20llf, steps: %d, processes: %d\n", res, steps, world_size);
        printf("true euler: 2.71828182845904523536028747135266249775724709369995...\n");
        printf("\nelapsed time: %.10f s\n", elapsed_time);
    }

    return 0;
}