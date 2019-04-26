#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<stddef.h>

#include "clcg4.h"
#include<mpi.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif

#ifdef BGQ
#define SIZE 32768
#else
#define SIZE 1024
#endif

#define NUM_TICKS 256
#define CAR 2

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/


double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

int mpi_myrank;
int mpi_commsize;

MPI_Request requests[4];
MPI_Status status[4];

unsigned int miles_per_rank;
unsigned int global_row;
unsigned int *interstate_west; // index 0 on eastern side, index miles_per_rank-1 on western side
unsigned int *interstate_east; // index 0 on western side, index miles_per_rank-1 on western side
unsigned int global_reached[NUM_TICKS]; // number of cars that have reached their destination.
unsigned int local_reached[NUM_TICKS];
unsigned int global_total_cars;
unsigned int local_miles_travelled;
unsigned int global_miles_travelled;

/***************************************************************************/
/* Functions ***********************************************************/
/***************************************************************************/

unsigned int mk_cars(float proportion);
int update_roads(unsigned int miles_per_rank, unsigned int *interstate_west,
        unsigned int *interstate_east, int can_move_west, int can_move_east);
int reached_dest(unsigned int *interstate, int i, int east);
int run_simulation(unsigned int *interstate_west, unsigned int *interstate_east);

void print_interstate(unsigned int *interstate_west, unsigned int *interstate_east);



int main(int argc, char *argv[]) {
    // Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    miles_per_rank = SIZE/mpi_commsize;
    global_row = CAR*mpi_myrank*miles_per_rank;

    InitDefault();

    if (0 == mpi_myrank) g_start_cycles = GetTimeBase();

    interstate_west = calloc(CAR*miles_per_rank, sizeof(int));
    interstate_east = calloc(CAR*miles_per_rank, sizeof(int));

    for (size_t i = 0; i < CAR*miles_per_rank; ++i) {
        interstate_west[i] = 0;
        interstate_east[i] = 0;
    }

    double proportion = 0.05;
    int total_cars = mk_cars(proportion);
    printf("total cars: %d\n", total_cars);
    print_interstate(interstate_west, interstate_east);
    local_miles_travelled = run_simulation(interstate_west, interstate_east);

    MPI_Reduce(&total_cars, &global_total_cars, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_miles_travelled, &global_miles_travelled, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_reached, global_reached, NUM_TICKS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (0 == mpi_myrank) {
        g_end_cycles = GetTimeBase();
        g_time_in_secs = ((double) (g_end_cycles - g_start_cycles)) /
                    g_processor_frequency;
    }

    if (0 == mpi_myrank) {
        for (size_t i = 0; i < NUM_TICKS; ++i) {
            printf("global_state: %d\n", global_reached[i]);
        }
        printf("total cars: %d\n", global_total_cars);
        printf("total miles travelled: %d\n", global_miles_travelled);
        printf("Completed simulation of interstate traffic in %lf seconds.\n", g_time_in_secs);
    }

    free(interstate_west);
    free(interstate_east);
    MPI_Finalize();
    return 0;
}


// Make cars at random locations on the highway. A car is marked by a 1 (as opposed to a 0)
// at an even index, and the next entry contains the destination of the car.
unsigned int mk_cars(float proportion){ 
    int row, total_cars = 0;
    // do I double check here that the car is facing the correct direction?
    for (int i = 0; i < CAR*miles_per_rank; i += CAR) {
        row = i + mpi_myrank*miles_per_rank;
        interstate_west[i] = (int) (GenVal(row) < proportion);
        if (interstate_west[i]) {
            interstate_west[i+1] = (int) (GenVal(row)*SIZE); // should be in range i+2...2*miles_per_rank
            total_cars += 1;
        }
        else interstate_west[i+1] = 0;

        interstate_east[i] = (int) (GenVal(row+1) < proportion);
        if (interstate_east[i]) {
            interstate_east[i+1] = (int) (GenVal(row+1)*SIZE); // should be in range i+2...2*miles_per_rank
            total_cars += 1;
        }
        else interstate_east[i+1] = 0;
    }
    return total_cars;
}

// Update roads by moving cars to the east or west if there is space in front of them.
int update_roads(unsigned int miles_per_rank, unsigned int *interstate_west,
        unsigned int *interstate_east, int can_move_west, int can_move_east) {
    unsigned int miles_travelled = 0;
    if (can_move_west) { // erase car from eastern side of road
        interstate_west[CAR*miles_per_rank-CAR] = 0;
        interstate_west[CAR*miles_per_rank-CAR+1] = 0;
    }
    if (can_move_east) { // erase car from eastern side of road
        interstate_east[CAR*miles_per_rank-CAR] = 0;
        interstate_east[CAR*miles_per_rank-CAR+1] = 0;
    }
    for (int i = CAR*miles_per_rank-CAR; i >= CAR; i -= CAR) {
        if (!interstate_west[i] && interstate_west[i-CAR]) {
            interstate_west[i] = interstate_west[i-CAR];
            interstate_west[i+1] = interstate_west[i-CAR+1];
            interstate_west[i-CAR] = 0;
            miles_travelled += 1; // can I just add the total number when the car has reached its destination?
        }

        if (!interstate_east[i] && interstate_east[i-CAR]) {
            interstate_east[i] = interstate_east[i-CAR];
            interstate_east[i+1] = interstate_east[i-CAR+1];
            interstate_east[i-CAR] = 0;
            miles_travelled += 1;
        }
    }
    interstate_west[0] = 0;
    interstate_west[1] = 0;
    interstate_east[0] = 0;
    interstate_east[1] = 0;

    return miles_travelled;
}

int reached_dest(unsigned int *interstate, int i, int east) {
    if ((east && interstate[i+1] == global_row+i) ||
            (!east && interstate[CAR*miles_per_rank-CAR-i+1] == global_row+i)) { // global_row-2-i+1. is plus 1 right? is it global_row-2-(i+1)?
        printf("%d %d\n", global_row+i, global_row+CAR*miles_per_rank-CAR-i);
        // printf("%d %d %d\n", interstate[i+1], interstate[CAR*miles_per_rank-CAR-i+1], global_row+i);
        return 1;
    }
    return 0;
}

int run_simulation(unsigned int *interstate_west, unsigned int *interstate_east) {
    int c1[1], c2[1], m1[1], m2[2];
    int total_reached=0, null_val=-1, empty_val=0, full_val=1, can_move_west=0, can_move_east=0, miles_travelled=0;
    int *null = &null_val;
    int *empty = &empty_val;
    int *full = &full_val;
    for (int i = 0; i < NUM_TICKS; ++i) {
        // one rank is sending a message to another rank to see if it can send a car. If sent, it will be set.
        // How does it get the message back that it can move?
        // get neighboring cars
        // if car is at beginning of highway, can't receive car... need to send a message saying no
        if (mpi_myrank > 0) { // receive from west
            MPI_Irecv(c1, 1, MPI_INT, mpi_myrank-1, mpi_myrank-1, MPI_COMM_WORLD, &requests[0]); // receive from sender, which may or may not be empty
            MPI_Irecv(m1, 1, MPI_INT, mpi_myrank-1, (mpi_myrank-1)+mpi_commsize, MPI_COMM_WORLD, &requests[2]); // receive from lower rank about whether it can send
        }
        if (mpi_myrank < mpi_commsize-1) { // receive from east
            MPI_Irecv(c2, 1, MPI_INT, mpi_myrank+1, mpi_myrank+1, MPI_COMM_WORLD, &requests[1]); // receive from sender, which may or may not be empty
            MPI_Irecv(m2, 1, MPI_INT, mpi_myrank-1, mpi_myrank+1+mpi_commsize, MPI_COMM_WORLD, &requests[3]); // receive from lower rank about whether it can send
        }

        if (mpi_myrank > 0) { // send to west, say whether or not the space is free
            // empty or full to western east
            if (interstate_west[0])
                MPI_Isend(empty, 1, MPI_INT, mpi_myrank-1, mpi_myrank+mpi_commsize, MPI_COMM_WORLD, &requests[3]);
            else
                MPI_Isend(full, 1, MPI_INT, mpi_myrank-1, mpi_myrank+mpi_commsize, MPI_COMM_WORLD, &requests[3]);
            // send car or null to western west
            if (interstate_west[CAR*miles_per_rank-1])
                MPI_Isend(interstate_west+CAR*miles_per_rank-1, 1, MPI_INT, mpi_myrank-1, mpi_myrank, MPI_COMM_WORLD, &requests[1]);
            else
                MPI_Isend(null, 1, MPI_INT, mpi_myrank-1, mpi_myrank, MPI_COMM_WORLD, &requests[1]);
        }
        if (mpi_myrank < mpi_commsize-1) { // send to east, say whether or not the space is free
            // empty or full to eastern west
            if (interstate_east[0])
                MPI_Isend(empty, 1, MPI_INT, mpi_myrank+1, mpi_myrank+mpi_commsize, MPI_COMM_WORLD, &requests[2]);
            else
                MPI_Isend(full, 1, MPI_INT, mpi_myrank+1, mpi_myrank+mpi_commsize, MPI_COMM_WORLD, &requests[2]);
            // send car or null to eastern east
            if (interstate_east[CAR*miles_per_rank-1])
                MPI_Isend(interstate_east+CAR*miles_per_rank-CAR+1, 1, MPI_INT, mpi_myrank+1, mpi_myrank, MPI_COMM_WORLD, &requests[0]);
            else
                MPI_Isend(null, 1, MPI_INT, mpi_myrank+1, mpi_myrank, MPI_COMM_WORLD, &requests[0]);
        }

        if (mpi_myrank > 0) { // wait on west
            MPI_Wait(&requests[0], &status[0]); // receive car or null
            MPI_Wait(&requests[2], &status[2]); // receive something about whether or not space is free
            if (c1[0] >= 0) {
                interstate_east[0] = 1; // car exists
                interstate_east[1] = *c1; // destination
                // printf("c1: %d %d\n", c1[0], c1[1]);
            }
            if (m1[0]>=0) can_move_west = 1;
        }
        if (mpi_myrank < mpi_commsize-1) { // wait on east
            MPI_Wait(&requests[1], &status[1]); // receive car or null
            MPI_Wait(&requests[3], &status[3]); // receive something about whether or not space is free
            if (c2[0] >= 0) {
                interstate_west[0] = 1; // car exists
                interstate_west[1] = *c2; // destination
                // printf("c2: %d %d\n", c2[0], c2[1]);
            }
            if (m2[0]>=0) can_move_east = 1;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // update roads
        miles_travelled += update_roads(miles_per_rank, interstate_west, interstate_east, can_move_west, can_move_east);

        // print_interstate(interstate_west, interstate_east);

        // check if reached destination
        for (int j = 0; j < CAR*miles_per_rank; j += CAR) {
            if (interstate_west[j] && reached_dest(interstate_west, j, 0)) { // if '1' and at correct position
                interstate_west[j] = 0; // remove from interstate, they're taking their exit!
                interstate_west[j+1] = 0;
                total_reached += 1;
                printf("reached destination at %d! total_reached = %d\n", j, total_reached);
            }
            if (interstate_east[j] && reached_dest(interstate_east, j, 1)) {
                interstate_east[j] = 0;
                interstate_east[j+1] = 0;
                total_reached += 1;
                printf("reached destination at %d! total_reached = %d\n", j, total_reached);
            }
        }

        local_reached[i] = total_reached;

        MPI_Barrier(MPI_COMM_WORLD);

        printf("MPI_rank %d: total reached: %d\n", mpi_myrank, total_reached);
    }

    return miles_travelled;
}

void print_interstate(unsigned int *interstate_west, unsigned int *interstate_east) {
    printf("Rank %d:\n", mpi_myrank);
    for (int i = CAR*miles_per_rank-CAR; i >= 0; i -= CAR) {
        printf("%d ", interstate_west[i]);
    }
    printf("\n");
    for (int i = 0; i < CAR*miles_per_rank; i += CAR) {
        printf("%d ", interstate_east[i]);
    }
    printf("\n");
}

// make sure count is accurate
// try to turn cars around if facing wrong direction
// add distances the cars travel
// add weather events