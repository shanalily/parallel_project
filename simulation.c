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
#define SIDE_LENGTH 128
#endif


#define ROAD_CAP 4

/***************************************************************************/
/* Structs *****************************************************************/
/***************************************************************************/


// denotes a place on the grid with column #, row # and which part of the street
typedef struct location {
    unsigned int col;
    unsigned int row;
    unsigned int idx;
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

void print_grid();

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
    unsigned int glbl_index = rpr*mpi_myrank;
    InitDefault();

    // To save space even ticks will be computed in now variables and odd
    // ticks will be computed in in nxt variables.
    // The pointers will switch after every tick
    
    // array of streets
    street* streets_ew_now = malloc((rpr)*(SIDE_LENGTH-1)*sizeof(street));
    street* streets_ew_nxt = malloc((rpr)*(SIDE_LENGTH-1)*sizeof(street));
    street* streets_ns_now = malloc((rpr-1)*SIDE_LENGTH*sizeof(street));
    street* streets_ns_nxt = malloc((rpr-1)*SIDE_LENGTH*sizeof(street));
    // ghost streets
    // these streets are the ones in the middle of the rectangle
    street* ghost_ns_nrth = malloc(SIDE_LENGTH*sizeof(street));
    street* ghost_ns_soth = malloc(SIDE_LENGTH*sizeof(street));
    // array of intersections
    intrsctn* intrsctn_now = malloc(rpr*SIDE_LENGTH*sizeof(intrsctn));
    intrsctn* intrsctn_nxt = malloc(rpr*SIDE_LENGTH*sizeof(intrsctn));

    // connect up the intersections
    // TO DO: check to see if it works and repeat for next ones
    for(size_t i = 0; i < rpr*SIDE_LENGTH; i++)
    {
        // if (mpi_myrank == 0) {
        //     printf("%lu\n", i);
        // }
        
        if(i/SIDE_LENGTH == 0) {
            // if touching north side
            if(mpi_myrank != 0) {
                intrsctn_now[i].nrth = &ghost_ns_nrth[i];
                intrsctn_nxt[i].nrth = &ghost_ns_nrth[i];
            } else {
                intrsctn_now[i].nrth = NULL;
                intrsctn_nxt[i].nrth = NULL;
            }
        } else {
            intrsctn_now[i].nrth = &streets_ns_now[i-SIDE_LENGTH];
            intrsctn_nxt[i].nrth = &streets_ns_nxt[i-SIDE_LENGTH];
        }
        if(i/SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching south side
            if(mpi_myrank != mpi_commsize-1) {
                intrsctn_now[i].soth = &ghost_ns_soth[i%SIDE_LENGTH];
                intrsctn_nxt[i].soth = &ghost_ns_soth[i%SIDE_LENGTH];
            } else {
                intrsctn_now[i].soth = NULL;
                intrsctn_nxt[i].soth = NULL;
            }
        } else {
            intrsctn_now[i].soth = &streets_ns_now[i];
            intrsctn_nxt[i].soth = &streets_ns_nxt[i];
        }
        if(i%SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching east side
            intrsctn_now[i].east = NULL;
            intrsctn_nxt[i].east = NULL;
        } else {
            intrsctn_now[i].east = &streets_ew_now[i-i/SIDE_LENGTH];
            intrsctn_nxt[i].east = &streets_ew_nxt[i-i/SIDE_LENGTH];
        }
        if(i%SIDE_LENGTH == 0) {
            // if touching west side
            intrsctn_now[i].west = NULL;
            intrsctn_nxt[i].west = NULL;
        } else {
            intrsctn_now[i].west = &streets_ew_now[i-1-i/SIDE_LENGTH];
            intrsctn_nxt[i].west = &streets_ew_nxt[i-1-i/SIDE_LENGTH];
        }
    }

// check if grids are valid
#ifdef DEBUG
    check_grid(rpr, intrsctn_now, streets_ew_now, streets_ns_now, 
            ghost_ns_nrth, ghost_ns_soth);
    check_grid(rpr, intrsctn_nxt, streets_ew_nxt, streets_ns_nxt, 
            ghost_ns_nrth, ghost_ns_soth);
#endif

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

    // generate east_west cars
    int idx;
    int row;
    int col;
    for(size_t i = 0; i < rpr*(SIDE_LENGTH-1); i++)
    {
        idx = glbl_index + i/(SIDE_LENGTH-1);
        row = i/(SIDE_LENGTH-1);
        col = i%(SIDE_LENGTH-1);
        for(size_t j = 0; j < ROAD_CAP; j++) {
            if(GenVal(idx) < proportion) {
                // start location
                loc* strt = malloc(sizeof(loc));
                strt->col = col;
                strt->row = row;
                strt->idx = j;

                // end location (non-unique)
                loc* _end = malloc(sizeof(loc));
                _end->col = (int) GenVal(idx)*(SIDE_LENGTH-1);
                _end->row = (int) GenVal(idx)*(SIDE_LENGTH-1);
                _end->row = (int) GenVal(idx)*(ROAD_CAP);
                car* c = malloc(sizeof(car));
                c->start = strt;
                c->end = _end;

                // if it needs to go right or left
                if (_end->col >= strt->col) {
                    streets_ew_now[i].go_es[j] = c;
                }else{
                    streets_ew_now[i].go_wn[j] = c;
                }
                
            } 
        }
    }

    // generate north_south cars
    for(size_t i = 0; i < (rpr-1)*(SIDE_LENGTH); i++)
    {
        idx = glbl_index + i/(SIDE_LENGTH);
        row = i/(SIDE_LENGTH);
        col = i%(SIDE_LENGTH);
        for(size_t j = 0; j < ROAD_CAP; j++) {
            if(GenVal(idx) < proportion) {
                // start location
                loc* strt = malloc(sizeof(loc));
                strt->col = col;
                strt->row = row;
                strt->idx = j;

                // end location (non-unique)
                loc* _end = malloc(sizeof(loc));
                _end->col = (int) GenVal(idx)*(SIDE_LENGTH-1);
                _end->row = (int) GenVal(idx)*(SIDE_LENGTH-1);
                _end->row = (int) GenVal(idx)*(ROAD_CAP);
                car* c = malloc(sizeof(car));
                c->start = strt;
                c->end = _end;

                // if it needs to go right or left
                if (_end->row >= strt->row) {
                    streets_ns_now[i].go_es[j] = c;
                }else{
                    streets_ns_now[i].go_wn[j] = c;
                }
                
            } 
        }
    }
    
    

    
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
