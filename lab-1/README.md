# DESCRIPTION

solving convection equation using cross-scheme (and left-angle scheme on the borders)

# BUILD & RUN

to compile, run `mpicc -o main.out main.c`

to run, execute `mpirun -np <N_PROC>`

to see results, execute `pyhton3 plot.py`

last command will generate plot of the solution of given equation
execution time will be output in the cout

# BENCHMARKS

file `benchmark.txt` contains respective execution times for parallel and single-thread programs. to see for yourself def / undef word `PARALLEL` in `main.c`
