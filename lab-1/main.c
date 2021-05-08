#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG
// #undef DEBUG

typedef double point_t;

typedef point_t (*func_bound_t)(point_t);
typedef point_t (*func_t)(point_t, point_t);

const point_t X_FROM = 0.0;
const point_t X_TO = 1.0;
const size_t X_STEPS = 10;

const point_t T_FROM = 0.0;
const point_t T_TO = 1.0;
const size_t T_STEPS = 10;

const int ROOT_PROC = 0;
int PROC_RANK = 0;
int N_PROC = 0;

// returns value of u_(m)^(k + 1)
point_t sheme_cross(const size_t k, const point_t tau, const size_t m, const point_t h,
                    const point_t* u, const point_t f_k_m, const size_t segment_size) {
    return (f_k_m - (u[segment_size * k + m + 1] - u[segment_size * k + m - 1]) / (2 * h)) *
               (2 * tau) +
           u[segment_size * (k - 1) + m];
}

// returns value of u_(m)^(k + 1)
point_t sheme_langle(const size_t k, const point_t tau, const size_t m, const point_t h,
                     const point_t* u, const point_t f_k_m, const size_t segment_size) {
    return (f_k_m - (u[segment_size * k + m] - u[segment_size * k + m - 1]) / h) * tau +
           u[segment_size * k + m];
}

point_t f(point_t x, point_t t) {
    return x - t;
}

point_t phi(point_t x) {
    return x;
}

point_t ksi(point_t t) {
    return t;
}

struct pair_t {
    size_t begin;
    size_t end;
};
typedef struct pair_t pair_t;

void calculate(point_t* result, pair_t* map);

int main(int argc, char** argv) {
    int err = 0;
    err = MPI_Init(&argc, &argv);
    assert(err == 0);

    err = MPI_Comm_rank(MPI_COMM_WORLD, &PROC_RANK);
    assert(err == 0);
    err = MPI_Comm_size(MPI_COMM_WORLD, &N_PROC);
    assert(err == 0);

    // indexes of first and last t points for each process
    pair_t* map = (pair_t*)calloc(N_PROC, sizeof(pair_t));

    if (N_PROC > X_STEPS) {
        for (size_t i = 0; i < N_PROC; i++) {
            if (i < X_STEPS) {
                map[i].begin = i;
                map[i].end = i + 1;
            } else {
                map[i].begin = -1;
                map[i].end = -1;
            }
        }
    } else {
        size_t steps_per_proc = X_STEPS / N_PROC;
        size_t underhead = X_STEPS - steps_per_proc * N_PROC;

        for (size_t i = 0; i < N_PROC; i++) {
            map[i].begin = steps_per_proc * i;

            if (i == N_PROC - 1) {
                map[i].end = X_STEPS;
            } else {
                map[i].end = map[i].begin + steps_per_proc;
            }

            // balance excess jobs
            if (N_PROC - i <= underhead) {
                map[i].begin += i + underhead - N_PROC;
                map[i].end = map[i].begin + steps_per_proc + 1;
            }
        }
    }

#ifdef DEBUG
    if (PROC_RANK == ROOT_PROC) {
        for (size_t i = 0; i < N_PROC; i++) {
            printf("proc %2d shall calculate nodes [%2d, %2d) of x (%d in total)\n",
                   i,
                   map[i].begin,
                   map[i].end,
                   map[i].end - map[i].begin);
        }
    }
#endif

    // allocate memory to store results
    point_t* result = NULL;

    if (PROC_RANK == ROOT_PROC) {
        // master stores the entire grid
        size_t n_elements = X_STEPS * T_STEPS;
        result = (point_t*)calloc(n_elements, sizeof(point_t));
#ifdef DEBUG
        printf("process %d, allocated %d bytes\n", PROC_RANK, n_elements * sizeof(point_t));
#endif
    } else if (PROC_RANK != N_PROC - 1) {
        // slave stores only his span, and two nodes extra - one left and one right
        // except for the lat one, that doesn't need one extra on the right
        size_t n_elements = T_STEPS * (map[PROC_RANK].end - map[PROC_RANK].begin + 2);
        result = (point_t*)calloc(n_elements, sizeof(point_t));
#ifdef DEBUG
        printf("process %d, allocated %d bytes\n", PROC_RANK, n_elements * sizeof(point_t));
#endif
    } else {
        size_t n_elements = T_STEPS * (map[PROC_RANK].end - map[PROC_RANK].begin + 1);
        result = (point_t*)calloc(n_elements, sizeof(point_t));
#ifdef DEBUG
        printf("process %d, allocated %d bytes\n", PROC_RANK, n_elements * sizeof(point_t));
#endif
    }

    // syncronize
    assert(result != NULL);
    MPI_Barrier(MPI_COMM_WORLD);

    // calculate results
    calculate(result, map);

    err = MPI_Finalize();
    assert(err == 0);

    return 0;
}

void calculate(point_t* result, pair_t* map) {
    // pre calculating some constants
    const point_t step_x = (X_TO - X_FROM) / X_STEPS;
    const point_t step_t = (T_TO - T_FROM) / T_STEPS;

    const size_t region_x_begin = map[PROC_RANK].begin;
    const size_t region_x_end = map[PROC_RANK].end;

    size_t row_length = 0;
    if (PROC_RANK == ROOT_PROC) {
        row_length = X_STEPS;
    } else if (PROC_RANK != N_PROC - 1) {
        row_length = region_x_end - region_x_begin + 2;
    } else {
        row_length = region_x_end - region_x_begin + 1;
    }

    /**
     * memmaps:
     * ROOT:
     * [0][1] ... [N]  (entire grid)
     * OTHER:
     * [-1][0] ... [n][n+1] (two extra: one on the left, one on the right)
     */

    // iterating upon t for specific x region, then syncronizig results with other nodes.
    // repeat
    for (size_t t = 0; t < T_STEPS; t++) {
#ifdef DEBUG
        printf("process %2d is calculating for t = %2d, x = [%2d, %2d)\n",
               PROC_RANK,
               t,
               region_x_begin,
               region_x_end);
#endif
        for (size_t x = region_x_begin; x < region_x_end; x++) {
            // for iterating upon local array
            size_t x_local = x - region_x_begin;
            if (PROC_RANK != ROOT_PROC) {
                x_local += 1;
            }

            if (t == 0) {
                result[row_length * t + x_local] = phi(X_FROM + x * step_x);
                continue;
            }

            if (x == 0) {
                result[row_length * t + x_local] = ksi(T_FROM + t * step_t);
                continue;
            }

            if (x == 1 || t == T_STEPS - 1) {
                result[row_length * t + x_local] = 1; // FIXME: langle
                continue;
            }

            result[row_length * t + x_local] = 2.0; // FIXME: cross
#ifdef DEBUG
            printf("process %2d, offset %d\n", PROC_RANK, row_length * t + x_local);
#endif
        }

        // syncronize master and slaves. first, master recieves data from each process. then, each
        // process recieves needed nodes
        MPI_Status status;

#ifdef DEBUG
        printf("process %2d finished calculating for t = %2d, x = [%2d, %2d)\n",
               PROC_RANK,
               t,
               region_x_begin,
               region_x_end);
#endif

        // recieve directly calculated values from each process
        if (PROC_RANK == ROOT_PROC) {
            size_t offset = map[0].end - map[0].begin;
            for (size_t i = 1; i < N_PROC; i++) {
                size_t n_elements = map[i].end - map[i].begin;

                MPI_Recv(&result[X_STEPS * t + offset],
                         n_elements,
                         MPI_DOUBLE,
                         i,
                         0,
                         MPI_COMM_WORLD,
                         &status);

                offset += n_elements;

#ifdef DEBUG
                printf(
                    "ROOT process recieved %d elements from process %d, starting from offset %d\n",
                    n_elements,
                    i,
                    X_STEPS * t + offset - n_elements);
#endif
            }
        } else {
            MPI_Send(&result[row_length * t + 1], // +1 to offset the -1th element
                     region_x_end - region_x_begin,
                     MPI_DOUBLE,
                     0,
                     0,
                     MPI_COMM_WORLD);
#ifdef DEBUG
            printf("process %d sent %d elements to ROOT, starting from offset %d\n",
                   PROC_RANK,
                   region_x_end - region_x_begin,
                   row_length * t + 1);
#endif
        }

        // wait untill each message has been sent MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // send -1th elements
        if (PROC_RANK == ROOT_PROC) {
            size_t offset = map[0].end - map[0].begin;
            for (size_t i = 1; i < N_PROC; i++) {
                size_t n_elements = map[i].end - map[i].begin;

                MPI_Send(&result[X_STEPS * t + offset - 1], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

                offset += n_elements;

#ifdef DEBUG
                printf("ROOT process sent %d th element to process %d as -1th\n",
                       X_STEPS * t + offset - 1 - n_elements,
                       i);
#endif
            }
        } else {
            MPI_Recv(&result[row_length * t], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG
            printf(
                "process %d recieved %d th element from ROOT as -1th\n", PROC_RANK, row_length * t);
#endif
        }

        // wait untill each message has been sent
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        // send +1th elements
        if (PROC_RANK == ROOT_PROC) {
            size_t offset = map[0].end - map[0].begin;
            for (size_t i = 1; i < N_PROC - 1; i++) {
                size_t n_elements = map[i].end - map[i].begin;

                MPI_Send(&result[X_STEPS * t + offset + n_elements],
                         1,
                         MPI_DOUBLE,
                         i,
                         0,
                         MPI_COMM_WORLD);

                offset += n_elements;

#ifdef DEBUG
                printf("ROOT process sent %d th element to process %d as +1th\n",
                       X_STEPS * t + offset,
                       i);
#endif
            }
        } else if (PROC_RANK != N_PROC - 1) {
            MPI_Recv(
                &result[row_length * (t + 1) - 1], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG
            printf("process %d recieved %d th element from ROOT as +1th\n",
                   PROC_RANK,
                   row_length * (t + 1) - 1);
#endif
        }

        // wait untill each message has been sent
        MPI_Barrier(MPI_COMM_WORLD);
    }
}