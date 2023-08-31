#FLAGS=-I_MPI_WAIT_MODE=0 -I_MPI_THREAD_YIELD=3 -I_MPI_THREAD_SLEEP=10


all: knn.cpp
	mpic++ -Wno-narrowing point.cpp knn.cpp -o knn

run:
	mpirun -n 4 --oversubscribe knn environment.txt challenge.txt 20 result.txt

generator:
	g++ generator.cpp -o generator

generate_env:
	./generator 20 1000 10 environment.txt

generate_challenge:
	./generator 500 10 challenge.txt