all: clcg4.h clcg4.c simulation.c
	gcc -I. -Wall -O3 -c clcg4.c -o clcg4.o
	mpicc -I. -Wall -O3 simulation.c clcg4.o -o simulation -lpthread
