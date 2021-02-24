MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

all: build euler-l euler-p euler-p-la

euler-l: euler-linear.c
	g++ -o bin/euler-l.out euler-linear.c

euler-p: euler-parallel.c
	mpicc $(MPI_COMPILE_FLAGS) euler-parallel.c $(MPI_LINK_FLAGS) -o bin/euler-parallel.out

euler-p-la: euler-parallel-la.c
	mpicc $(MPI_COMPILE_FLAGS) euler-parallel-la.c $(MPI_LINK_FLAGS) -o bin/euler-parallel-la.out

build:
	mkdir -p bin/

clean:
	rm -rf bin/ *.out