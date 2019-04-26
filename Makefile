home: clcg4.h clcg4.c simulation.c
	gcc -I. -Wall -O3 -c clcg4.c -o clcg4.o
	mpicc -I. -Wall -O3 sim_simple.c clcg4.o -o simulation -lpthread

onedim: clcg4.h clcg4.c sim_interstate.c
	gcc -I. -Wall -O3 -c clcg4.c -o clcg4.o
	mpicc -I. -Wall -O3 sim_interstate.c clcg4.o -o sim -lpthread

all:
	mpixlc -O5 mpi-hello.c -o mpi-hello.xl

debug:
	mpixlc -g mpi-hello.c -o mpi-hello.xl
