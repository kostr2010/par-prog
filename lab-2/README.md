# DESCRIPTION

calculating integral of sin(1 / x) from 1e-3 to 1.0 using pthreads, with adjusting step and set precision

# BUILD & RUN

to compile, run `make all`

to run, execute `./main.out`

program will print result to stdout

# NOTES

as sin(1/x) oscillates very fast in zero proximity, it's impossible to use traditional integration methods properly. due to the lack of techniques available, the best working approachh is to calculate step based on change of argument for 2pi
