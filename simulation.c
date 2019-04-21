/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include <stddef.h>

#include "clcg4.h"
#include<mpi.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif
#define DEBUG 1

#include <assert.h>

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#ifdef BGQ
#define SIDE_LENGTH 32768
#else
#define SIDE_LENGTH 8
#endif


#define ROAD_CAP 4

/***************************************************************************/
/* Structs *****************************************************************/
/***************************************************************************/


// denotes a place on the grid with column #, row # and which part of the street
// typedef struct location {
//     unsigned long col;
//     unsigned long row;
//     unsigned short idx;
// } loc;

// denotes a car
// s stands for start e stands for end
typedef struct vehicle {
    unsigned long s_col;
    unsigned long s_row;
    unsigned short s_idx;
    unsigned long e_col;
    unsigned long e_row;
    unsigned short e_idx;
} car;

typedef struct streets {
    car* go_es[ROAD_CAP];
    car* go_wn[ROAD_CAP];
} street;

typedef struct intrsctns {
    street* nrth;
    street* soth;
    street* west;
    street* east;
} intrsctn;


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// this checks the grid for validity
void check_grid(unsigned int rpr, intrsctn* is, 
        street* s_ew, street* s_ns, street* g_n_ns, street* g_s_ns);

void check_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col);

void mk_grid(unsigned int rpr, intrsctn* is, street* s_ns, street* s_ew, 
        street* g_ns_n, street* g_ns_s, int mpi_myrank, int mpi_commsize);

void generate_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col, float proportion);

int compare_sloc(car* c, unsigned long col, unsigned long row, unsigned short idx) {
    return c->s_col == col && c->s_row == row && c->s_idx == idx;
}
int compare_eloc(car* c, unsigned long col, unsigned long row, unsigned short idx) {
    return c->e_col == col && c->e_row == row && c->e_idx == idx;
}

MPI_Datatype mkcartype();

void pack_transfer(car* tt_n, street* gs_n, 
        car* tt_s, street* gs_s);
void unpack_transfer(car* tr_n, street* gs_n,
        car* tr_s, street* gs_s);

void transfer(car* tt_n, car* tr_n, car* tt_s, car* tr_s, int mpi_myrank, int mpi_commsize, 
        MPI_Datatype t_type);

void update_intersections(unsigned int rpr, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt);

void update_streets(street *ew_now, street *ns_now, street *ew_nxt, street *ns_nxt,
        street *ghost_nrth_now, street *ghost_soth_now, street *ghost_nrth_nxt,
            street *ghost_soth_nxt);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/
int main(int argc, char *argv[])
{
    int mpi_myrank;
    int mpi_commsize;
    // Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

    // constants 
    // Rows per rank
    unsigned int rpr = SIDE_LENGTH/mpi_commsize;
    float proportion = 1;
    unsigned long long capacity = 2*(SIDE_LENGTH-1)*(SIDE_LENGTH)*ROAD_CAP;
    unsigned long long n_cars = (unsigned long long) (((long double) proportion)*capacity);
    unsigned int glbl_index = 2*rpr*mpi_myrank;
    InitDefault();

    // To save space even ticks will be computed in now variables and odd
    // ticks will be computed in in nxt variables.
    // The pointers will switch after every tick
    
    // array of streets
    street* streets_ew_now = calloc((rpr)*(SIDE_LENGTH-1), sizeof(street));
    street* streets_ew_nxt = calloc((rpr)*(SIDE_LENGTH-1), sizeof(street));
    street* streets_ns_now = calloc((rpr-1)*SIDE_LENGTH, sizeof(street));
    street* streets_ns_nxt = calloc((rpr-1)*SIDE_LENGTH, sizeof(street));
    // ghost streets
    // these streets are the ones in the middle of the rectangle
    street* ghost_ns_nrth_now = calloc(SIDE_LENGTH, sizeof(street));
    street* ghost_ns_soth_now = calloc(SIDE_LENGTH, sizeof(street));
    street* ghost_ns_nrth_nxt = calloc(SIDE_LENGTH, sizeof(street));
    street* ghost_ns_soth_nxt = calloc(SIDE_LENGTH, sizeof(street));
    // array of intersections
    intrsctn* intrsctn_now = calloc(rpr*SIDE_LENGTH, sizeof(intrsctn));
    intrsctn* intrsctn_nxt = calloc(rpr*SIDE_LENGTH, sizeof(intrsctn));
    // ghost cars to transfer via mpi
    car* to_transfer_nrth = calloc(SIDE_LENGTH, sizeof(car));
    car* to_transfer_soth = calloc(SIDE_LENGTH, sizeof(car));
    car* to_receive_nrth = calloc(SIDE_LENGTH, sizeof(car));
    car* to_receive_soth = calloc(SIDE_LENGTH, sizeof(car));

    MPI_Datatype mpi_car_type = mkcartype();


    mk_grid(rpr, intrsctn_now, streets_ns_now, streets_ew_now,
            ghost_ns_nrth_now, ghost_ns_soth_now, mpi_myrank, mpi_commsize);
    mk_grid(rpr, intrsctn_nxt, streets_ns_nxt, streets_ew_nxt,
            ghost_ns_nrth_nxt, ghost_ns_soth_nxt, mpi_myrank, mpi_commsize);


// check if grids are valid
#ifdef DEBUG
    check_grid(rpr, intrsctn_now, streets_ew_now, streets_ns_now, 
            ghost_ns_nrth_now, ghost_ns_soth_nxt);
    check_grid(rpr, intrsctn_nxt, streets_ew_nxt, streets_ns_nxt, 
            ghost_ns_nrth_nxt, ghost_ns_soth_nxt);
#endif
    // do e/w
    generate_cars(rpr*(SIDE_LENGTH-1), glbl_index, streets_ew_now, 0, proportion);
    // do n/s
    generate_cars((rpr-1)*(SIDE_LENGTH), glbl_index+1, streets_ns_now, 1, proportion);    
    // do ghost this will be a duplicate just for easy
    if(mpi_myrank != mpi_commsize-1) {
        generate_cars(SIDE_LENGTH, glbl_index+2*rpr-1, ghost_ns_soth_now, 1, proportion);
    }
    if(mpi_myrank != 0) {
        generate_cars(SIDE_LENGTH, glbl_index-1, ghost_ns_nrth_now, 1, proportion);
    }
#ifdef DEBUG
    check_cars(rpr*(SIDE_LENGTH-1), glbl_index, streets_ew_now, 0); 
    check_cars((rpr-1)*(SIDE_LENGTH), glbl_index+1, streets_ns_now, 1);
    if(mpi_myrank != 0) {
        check_cars(SIDE_LENGTH, glbl_index-1, ghost_ns_nrth_now, 1);
    }
    if(mpi_myrank != mpi_commsize-1) {
        check_cars(SIDE_LENGTH, glbl_index+2*rpr-1, ghost_ns_soth_now, 1);
    }
#endif
    
    MPI_Barrier( MPI_COMM_WORLD );

    // for loop of ticks
    // send and receive (blocking)
    // NEED TO SEND THE FIRST SLOT GOING SOUTH TO THE SOUTH
    // NEED TO SEND THE FIRST SLOT GOING NORTH TO THE NORTH
    pack_transfer(to_transfer_nrth, ghost_ns_nrth_now,
            to_transfer_soth, ghost_ns_soth_now);
    transfer(to_transfer_nrth, to_receive_nrth, to_transfer_soth, to_receive_soth,
            mpi_myrank, mpi_commsize, mpi_car_type);
    // NEED TO RECEIVE THE FIRST SLOT GOING SOUTH TO THE NORTH
    // NEED TO RECEIVE THE FIRST SLOT GOING NORTH TO THE SOUTH
    unpack_transfer(to_receive_nrth, ghost_ns_nrth_now,
            to_receive_soth, ghost_ns_soth_now);
    // update streets
    // run intersections

    // frees
    free(streets_ew_now);
    free(streets_ew_nxt);
    free(streets_ns_now);
    free(streets_ns_nxt);
    free(ghost_ns_nrth_now);
    free(ghost_ns_soth_now);
    free(ghost_ns_nrth_nxt);
    free(ghost_ns_soth_nxt);
    free(intrsctn_now);
    free(intrsctn_nxt);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

void check_row(intrsctn* is, int n_i) {
    for(size_t i = 0; i < n_i-1; i++)
    {
        assert(is[i].east == is[i+1].west);
    }
}
void check_col(intrsctn* is, int n_i) {
    for(size_t i = 0; i < n_i-1; i++)
    {
        assert(is[i].soth == is[i+SIDE_LENGTH].nrth);
    }
}

// FOR TESTING ONLY
// CHECKS IF GRID IS BUILT PROPERLY
void check_grid(unsigned int rpr, intrsctn* is, 
        street* s_ew, street* s_ns, street* g_n_ns, street* g_s_ns){
    for(size_t i = 0; i < rpr; i++)
    {
        check_row(&is[i*SIDE_LENGTH], SIDE_LENGTH);
    }

    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        check_col(&is[i], rpr);
    }
    
    
}

void check_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col){
    unsigned long row, col;
    for(size_t i = 0; i < n; i++)
    {
        // find current global positions
        row = glbl_row_idx + 2*(i/(SIDE_LENGTH-(!row_or_col)));
        col = 2*(i%(SIDE_LENGTH-!row_or_col)) + !row_or_col;
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            // only one side of street should be full            
            assert(!strts[i].go_es[j] || !strts[i].go_wn[j]);
            // if a street has a car it should be at the starting location
            assert(!strts[i].go_es[j] || (strts[i].go_es[j] && 
                    compare_sloc(strts[i].go_es[j], 
                    col, row, j)));
            assert(!strts[i].go_wn[j] || (strts[i].go_wn[j] && 
                    compare_sloc(strts[i].go_wn[j], 
                    col, row, j)));
            
        }
        
    }
    
}

void mk_grid(unsigned int rpr, intrsctn* is, street* s_ns, street* s_ew, 
        street* g_ns_n, street* g_ns_s, int mpi_myrank, int mpi_commsize){
    // connect up the intersections
    for(size_t i = 0; i < rpr*SIDE_LENGTH; i++)
    {
        // if (mpi_myrank == 0) {
        //     printf("%lu\n", i);
        // }
        
        if(i/SIDE_LENGTH == 0) {
            // if touching north side
            if(mpi_myrank != 0) {
                is[i].nrth = &g_ns_n[i];
            } else {
                is[i].nrth = NULL;
            }
        } else {
            is[i].nrth = &s_ns[i-SIDE_LENGTH];
        }
        if(i/SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching south side
            if(mpi_myrank != mpi_commsize-1) {
                is[i].soth = &g_ns_s[i%SIDE_LENGTH];
            } else {
                is[i].soth = NULL;
            }
        } else {
            is[i].soth = &s_ns[i];
        }
        if(i%SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching east side
            is[i].east = NULL;
        } else {
            is[i].east = &s_ew[i-i/SIDE_LENGTH];
        }
        if(i%SIDE_LENGTH == 0) {
            // if touching west side
            is[i].west = NULL;
        } else {
            is[i].west = &s_ew[i-1-i/SIDE_LENGTH];
        }
    }
}

void generate_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col, float proportion){
    // initial cars
    // each car will have a starting point and an end point
    // each row will use its corresponding rng number generator
    // ew rows use generators 0, 2, 4, etc
    // ns rows use generators 1, 3, 5, etc
    // each slot will have a threshold chance to have a starting car
    // then it will pick a random slot on the grid to finish at
    // DUE to the very large nature of the possible number of slots
    // 4 random numbers will be used to define a finished location
    //  1 for on e/w or n/s
    //  1 for row
    //  1 for column
    //  1 for street position

    
    // for row_or_col, row=0, col=1
    unsigned long row, col;
    for(size_t i =0; i < n; i++) {
        // find current global positions
        row = glbl_row_idx + 2*(i/(SIDE_LENGTH-(!row_or_col)));
        col = 2*(i%(SIDE_LENGTH-!row_or_col)) + !row_or_col;
        for(size_t j = 0; j < ROAD_CAP; j++) {
            if(GenVal(row) < proportion) {
                // start location
                car* c = calloc(1, sizeof(car));
                c->s_col = col;
                c->s_row = row;
                c->s_idx = j;

                // end location (non-unique)
                c->e_col = (unsigned long) GenVal(row)*(2*SIDE_LENGTH-1);
                c->e_row = (unsigned long) GenVal(row)*(2*SIDE_LENGTH-1);
                c->e_idx = (int) GenVal(row)*(ROAD_CAP);

                if ((c->e_col >= c->s_col && !row_or_col) ||
                        (c->e_row >= c->s_row && row_or_col)) {
                    strts[i].go_es[j] = c;
                }else{
                    strts[i].go_wn[j] = c;
                }
            }
        }
    }
}

MPI_Datatype mkcartype(){
        // mpi struct
    const int n_items = 6;
    int blocklengths[6] = {1,1,1,1,1,1};
    MPI_Datatype types[6] = {MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, 
                            MPI_UNSIGNED_SHORT, MPI_UNSIGNED_LONG,
                            MPI_UNSIGNED_LONG, MPI_UNSIGNED_SHORT};
    MPI_Datatype mpi_car_type;
    MPI_Aint offsets[6];
    offsets[0] = offsetof(car, s_col);
    offsets[1] = offsetof(car, s_row);
    offsets[2] = offsetof(car, s_idx);
    offsets[3] = offsetof(car, e_col);
    offsets[4] = offsetof(car, e_row);
    offsets[5] = offsetof(car, e_idx);
    MPI_Type_create_struct(n_items, blocklengths, offsets, types, &mpi_car_type);
    MPI_Type_commit(&mpi_car_type);
    return mpi_car_type;
}

void pack_transfer(car* tt_n, street* gs_n, 
        car* tt_s, street* gs_s){
    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        memcpy(&tt_n[i], &gs_n[i].go_wn[i], sizeof(car));
        memcpy(&tt_s[i], &gs_s[i].go_es[i], sizeof(car));
    }
    

}
void unpack_transfer(car* tr_n, street* gs_n,
        car* tr_s, street* gs_s){
    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        memcpy(&tr_n[i], &gs_n[i].go_es[i], sizeof(car));
        memcpy(&tr_s[i], &gs_s[i].go_wn[i], sizeof(car));
    }
}

void transfer(car* tt_n, car* tr_n, car* tt_s, car* tr_s, int mpi_myrank, int mpi_commsize,
        MPI_Datatype t_type){
    MPI_Request nrth;
    MPI_Request soth;
    MPI_Request request;
    if(mpi_myrank != 0) {
        MPI_Irecv(tr_n, SIDE_LENGTH, t_type, mpi_myrank-1, 0, MPI_COMM_WORLD, &nrth);
    }
    if(mpi_myrank != mpi_commsize -1) {
        MPI_Irecv(tr_s, SIDE_LENGTH, t_type, mpi_myrank+1, 1, MPI_COMM_WORLD, &soth);
    }
    if(mpi_myrank != mpi_commsize -1) {
        MPI_Isend(tt_s, SIDE_LENGTH, t_type, mpi_myrank+1, 0, MPI_COMM_WORLD, &request);
    }
    if(mpi_myrank != 0) {
        MPI_Isend(tt_n, SIDE_LENGTH, t_type, mpi_myrank-1, 1, MPI_COMM_WORLD, &request);
    }

    MPI_Wait(&nrth, MPI_STATUS_IGNORE);
    MPI_Wait(&soth, MPI_STATUS_IGNORE);
}

// move everything down the streets
void update_streets(street *ew_now, street *ns_now, street *ew_nxt, street *ns_nxt,
        street *ghost_nrth_now, street *ghost_soth_now, street *ghost_nrth_nxt,
            street *ghost_soth_nxt) {

}

// make sure to switch which are passed...
// first give priority to 
void update_intersections(unsigned int rpr, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt) {
    // update nxt based on now
    for(size_t i = 0; i < rpr*SIDE_LENGTH; i++) {

        if(i/SIDE_LENGTH == 0) {
            // if touching north side
            // if(mpi_myrank != 0) {
            //     is[i].nrth = &g_ns_n[i];
            // } else {
            //     is[i].nrth = NULL;
            // }
        } else {
            // is[i].nrth = &s_ns[i-SIDE_LENGTH];
        }
        if(i/SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching south side
            // if(mpi_myrank != mpi_commsize-1) {
            //     is[i].soth = &g_ns_s[i%SIDE_LENGTH];
            // } else {
            //     is[i].soth = NULL;
            // }
        } else {
            // is[i].soth = &s_ns[i];
        }
        if(i%SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching east side
            // is[i].east = NULL;
        } else {
            // is[i].east = &s_ew[i-i/SIDE_LENGTH];
        }
        if(i%SIDE_LENGTH == 0) {
            // if touching west side
            // is[i].west = NULL;
        } else {
            // is[i].west = &s_ew[i-1-i/SIDE_LENGTH];
        }
    }
}