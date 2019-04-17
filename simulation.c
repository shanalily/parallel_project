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

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#ifdef BGQ
#define SIDE_LENGTH 32768
#else
#define SIDE_LENGTH 128
#endif

// #define DEBUG 1

#define ROAD_CAP 4

/***************************************************************************/
/* Structs *****************************************************************/
/***************************************************************************/

typedef struct streets {
    short frwd[ROAD_CAP];
    short bckw[ROAD_CAP];
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
check_grid(intrsctn* is, street* s_ew, street* s_ns, street* g_n_ns, street* g_s_ns, street* g_s_ew);

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

    // To save space even ticks will be computed in now variables and odd
    // ticks will be computed in in nxt variables.
    // The pointers will switch after every tick
    
    // array of streets
    street* streets_ew_now = malloc((rpr)*(SIDE_LENGTH-1)*sizeof(street));
    street* streets_ew_nxt = malloc((rpr)*(SIDE_LENGTH-1)*sizeof(street));
    street* streets_ns_now = malloc((rpr-1)*SIDE_LENGTH*sizeof(street));
    street* streets_ns_nxt = malloc((rpr-1)*SIDE_LENGTH*sizeof(street));
    // ghost streets
    // need EW streets from rank North and NS from rank both North and South
    // BOTTOM ROW MUST keep track of GHOST_ew_soth
    street* ghost_ns_nrth = malloc(SIDE_LENGTH*sizeof(street));
    street* ghost_ns_soth = malloc(SIDE_LENGTH*sizeof(street));
    street* ghost_ew_soth = malloc((SIDE_LENGTH-1)*sizeof(street));
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
            intrsctn_now[i].nrth = &streets_ns_now[i];
            intrsctn_nxt[i].nrth = &streets_ns_nxt[i];
        }
        if(i/SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching south side
            if(mpi_myrank != mpi_commsize-1) {
                intrsctn_now[i].soth = &ghost_ns_soth[i];
                intrsctn_nxt[i].soth = &ghost_ns_soth[i];
            } else {
                intrsctn_now[i].soth = NULL;
                intrsctn_nxt[i].soth = NULL;
            }
        } else {
            intrsctn_now[i].soth = &streets_ns_now[i-SIDE_LENGTH];
            intrsctn_nxt[i].soth = &streets_ns_nxt[i-SIDE_LENGTH];
        }
        if(i%SIDE_LENGTH == 0) {
            // if touching east side
            intrsctn_now[i].east = NULL;
            intrsctn_nxt[i].east = NULL;
        } else {
            if(i/SIDE_LENGTH == SIDE_LENGTH-1) {
                // if touching south side
                intrsctn_now[i].east = &ghost_ew_soth[i];
                intrsctn_nxt[i].east = &ghost_ew_soth[i];
            } else {
                intrsctn_now[i].east = &streets_ew_now[i-i/SIDE_LENGTH];
                intrsctn_nxt[i].east = &streets_ew_nxt[i-i/SIDE_LENGTH];
            }
        }
        if(i%SIDE_LENGTH == SIDE_LENGTH-1) {
            // if touching west side
            intrsctn_now[i].west = NULL;
            intrsctn_nxt[i].west = NULL;
        } else {
            if(mpi_myrank != mpi_commsize-1) {
                intrsctn_now[i].west = &ghost_ew_soth[i-1];
                intrsctn_nxt[i].west = &ghost_ew_soth[i-1];
            } else {
                intrsctn_now[i].west = &streets_ew_now[i-i/SIDE_LENGTH];
                intrsctn_nxt[i].west = &streets_ew_now[i-i/SIDE_LENGTH];
            }
        }
    }

// check if grids are valid
#ifdef DEBUG
    check_grid(intrsctn_now, streets_ew_now, streets_ns_now, 
            ghost_ns_nrth, ghost_ns_soth, ghost_ew_soth);
    check_grid(intrsctn_nxt, streets_ew_nxt, streets_ns_nxt, 
            ghost_ns_nrth, ghost_ns_soth, ghost_ew_soth);
#endif

    // initial cars
    // each car will have a starting point and an end point
    // each row will use its corresponding rng number generator
    // ew rows use generators 0, 2, 4, etc
    // ns rows use generators 1, 3, 5, etc
    // each slot will have a threshold chance to have a starting car
    // then it will pick a random slot on the grid to finish at
    

    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
