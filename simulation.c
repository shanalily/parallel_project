/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

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
typedef struct location {
    unsigned long col;
    unsigned long row;
    unsigned short idx;
} loc;

// denotes a car
typedef struct vehicle {
    loc* start;
    loc* end;
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

void mk_grid(unsigned int rpr, intrsctn* is, street* s_ns, street* s_ew, 
        street* g_ns_n, street* g_ns_s, int mpi_myrank, int mpi_commsize);

void generate_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col, float proportion);

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
    float proportion = 0.25;
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
    street* ghost_ns_nrth = calloc(SIDE_LENGTH, sizeof(street));
    street* ghost_ns_soth = calloc(SIDE_LENGTH, sizeof(street));
    // array of intersections
    intrsctn* intrsctn_now = calloc(rpr*SIDE_LENGTH, sizeof(intrsctn));
    intrsctn* intrsctn_nxt = calloc(rpr*SIDE_LENGTH, sizeof(intrsctn));

    mk_grid(rpr, intrsctn_now, streets_ns_now, streets_ew_now,
            ghost_ns_nrth, ghost_ns_soth, mpi_myrank, mpi_commsize);
    mk_grid(rpr, intrsctn_nxt, streets_ns_nxt, streets_ew_nxt,
            ghost_ns_nrth, ghost_ns_soth, mpi_myrank, mpi_commsize);


// check if grids are valid
#ifdef DEBUG
    check_grid(rpr, intrsctn_now, streets_ew_now, streets_ns_now, 
            ghost_ns_nrth, ghost_ns_soth);
    check_grid(rpr, intrsctn_nxt, streets_ew_nxt, streets_ns_nxt, 
            ghost_ns_nrth, ghost_ns_soth);
#endif
    // do e/w
    generate_cars(rpr*(SIDE_LENGTH-1), mpi_myrank*rpr, streets_ew_now, 0, proportion);
    // do n/s
    generate_cars((rpr-1)*(SIDE_LENGTH), mpi_myrank*rpr, streets_ns_now, 1, proportion);    

    
    MPI_Barrier( MPI_COMM_WORLD );

    // frees
    free(streets_ew_now);
    free(streets_ew_nxt);
    free(streets_ns_now);
    free(streets_ns_nxt);
    free(ghost_ns_nrth);
    free(ghost_ns_soth);
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
    unsigned long idx, row, col;
    for(size_t i =0; i < n; i++) {
        // find current global positions
        idx = glbl_row_idx + 2*(i/(SIDE_LENGTH-(!row_or_col)));
        row = 2*(i/(SIDE_LENGTH-!row_or_col)) + row_or_col;
        col = 2*(i%(SIDE_LENGTH-!row_or_col)) + !row_or_col;
        for(size_t j = 0; j < ROAD_CAP; j++) {
            if(GenVal(idx) < proportion) {
                // start location
                loc* strt = calloc(1, sizeof(loc));
                strt->col = col;
                strt->row = row;
                strt->idx = j;

                // end location (non-unique)
                loc* _end = calloc(1, sizeof(loc));
                _end->col = (unsigned long) GenVal(idx)*(2*SIDE_LENGTH-1);
                _end->row = (unsigned long) GenVal(idx)*(2*SIDE_LENGTH-1);
                _end->row = (int) GenVal(idx)*(ROAD_CAP);
                car* c = calloc(1, sizeof(car));
                c->start = strt;
                c->end = _end;

                if ((_end->col >= strt->col && !row_or_col) ||
                        (_end->row >= strt->row && row_or_col)) {
                    strts[i].go_es[j] = c;
                }else{
                    strts[i].go_wn[j] = c;
                }
            }
        }
    }
}
