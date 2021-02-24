#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    assert(argc == 2 /* argument STEPS needed */);

    clock_t t_begin = clock();

    int steps = atoi(argv[1]);
    long double euler = 0, fac = 1;

    for (int i = steps; i >= 1; i--) {
        fac *= i;
        euler += fac;
    }

    euler /= fac;

    clock_t t_end = clock();
    double elapsed_time = (double)(t_end - t_begin) / CLOCKS_PER_SEC;

    printf("euler:      %.20llf\n", euler);
    printf("true euler: 2.71828182845904523536028747135266249775724709369995...\n");
    printf("\nelapsed time: %.10f s\n", elapsed_time);

    return 0;
}