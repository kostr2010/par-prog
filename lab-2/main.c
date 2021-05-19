#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "calc.h"

const size_t THREADS_N = 10;
const double CALCULATE_FROM = 1e-3;
const double CALCULATE_TO = 1.0;

int main(int argc, char** argv)
{
    calc_space_t space = {.intervals = (pair_t*)calloc(sizeof(pair_t), INTERVALS_N),
                          .intervals_cur = 0};
    double intervals_step = (CALCULATE_TO - CALCULATE_FROM) / INTERVALS_N;

    for (size_t i = 0; i < INTERVALS_N; i++)
    {
        space.intervals[i].first = CALCULATE_FROM + i * intervals_step;

        if (i == INTERVALS_N - 1)
        {
            space.intervals[i].second = CALCULATE_TO;
        }
        else
        {
            space.intervals[i].second = CALCULATE_FROM + (i + 1) * intervals_step;
        }
    }

    pthread_t* threads = calloc(sizeof(pthread_t), THREADS_N);

    pthread_mutex_init(&(space.mux), NULL);

    for (size_t i = 0; i < THREADS_N; i++)
    {
        pthread_create(threads + i, NULL, Routine, &space);
    }

    double integral = 0.0;

    for (size_t i = 0; i < THREADS_N; ++i)
    {
        double* res = NULL;
        pthread_join(threads[i], (void**)&res);
        integral += *(res);
        free(res);
    }

    pthread_mutex_destroy(&(space.mux));
    printf("result:         %.10f\n", integral);
    printf("correct result: 0.504066\n");
    return 0;
}
