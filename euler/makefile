MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

all: build euler-l euler-p euler-p-la

euler-l: euler-linear.cpp
	g++ -o bin/euler-l.out euler-linear.cpp

euler-p: euler-parallel.cpp
	mpic++ $(MPI_COMPILE_FLAGS) euler-parallel.cpp $(MPI_LINK_FLAGS) -o bin/euler-parallel.out

euler-p-la: euler-parallel-la.cpp
	mpic++ $(MPI_COMPILE_FLAGS) euler-parallel-la.cpp $(MPI_LINK_FLAGS) -o bin/euler-parallel-la.out

build:
	mkdir -p bin/

clean:
	rm -rf bin/ *.out