/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include"clcg4.h"

#include<mpi.h>

#define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime            
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

// Each rank is responsible for a SIDE_LENGTH by SIDE_LENGTH square at the 
// maximun number of ranks. Size of region responsible grows by 4x
// number or ranks/nodes must be a square. I.E. 4, 16, 64, 256

#ifdef BGQ
#define MAX_RANKS 16384
#define MAX_SIDE_RANKS 128
#define MIN_REGION_LENGTH 64
#else
#define MAX_RANKS 16
#define MAX_SIDE_RANKS 4
#define MIN_REGION_LENGTH 16
#endif


#define SIDE_LENGTH MIN_REGION_LENGTH*MAX_SIDE_RANKS
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
    unsigned int sqr_size = MIN_REGION_LENGTH*((int) sqrt((float)(MAX_RANKS/mpi_commsize)));
    unsigned int n_intrsctn = (sqr_size+1)*(sqr_size+1);
    short westmost = mpi_myrank%MAX_SIDE_RANKS == MAX_SIDE_RANKS-1;
    short sothmost = mpi_myrank/MAX_SIDE_RANKS == MAX_SIDE_RANKS-1;


    // To save space even ticks will be computed in now variables and odd
    // ticks will be computed in in nxt variables.
    // The pointers will switch after every tick
    
    // array of streets
    street* streets_ew_now = malloc(sqr_size*sqr_size*sizeof(street));
    street* streets_ew_nxt = malloc(sqr_size*sqr_size*sizeof(street));
    street* streets_ns_now = malloc(sqr_size*sqr_size*sizeof(street));
    street* streets_ns_nxt = malloc(sqr_size*sqr_size*sizeof(street));
    // array of intersections
    intrsctn* intrsctn_now = malloc(n_intrsctn*sizeof(intrsctn));;
    intrsctn* intrsctn_nxt = malloc(n_intrsctn*sizeof(intrsctn));;
    // ghost streets
    // direction from square followed by direction of street
    street* ghost_east_soth = malloc((sqr_size+1)*sizeof(street));
    street* ghost_east_east = malloc((sqr_size+1)*sizeof(street));
    street* ghost_soth_east = malloc((sqr_size+1)*sizeof(street));
    street* ghost_soth_soth = malloc((sqr_size+1)*sizeof(street));
    street* ghost_west_west = malloc((sqr_size+1)*sizeof(street));
    street* ghost_nrth_nrth = malloc((sqr_size+1)*sizeof(street));
    /* NOTE: 
    for streets to the west at the sw corner it is from the sw square
    for streets to the north at the ne corner it is from the ne square
    Thus must communicate with at most 6 neighbor squares
    */

    // only squares of type 2, 3 and 4 will use the last intersection
    // tie together interesections
    for(size_t i = 0; i < n_intrsctn; i++)
    {
        if (i == 0) {
            // NW corner
            intrsctn_now[i].nrth = ghost_nrth_nrth[i];
            intrsctn_now[i].east = strets_ew_now[i];
            intrsctn_now[i].soth = strets_ns_now[i];
            intrsctn_now[i].east = strets_now[i];
        }
        else if(i == sqr_size) {
            // NE corner
        }
        else if(i == (sqr_size+1)*sqr_size) {
            // SW corner
        }
        else if(i == n_intrsction-1) {
            // SE corner
        }
        else if (i%(sqr_size+1) == 0) {
            // W side
        }
        else if (i/(sqr_size+1) == 0) {
            // N side
        }
        else if (i%(sqr_size+1) == sqr_size) {
            // E side
        }
        else if (i/(sqr_size+1) == sqr_size) {
            // S side
        }
        
    }
    


    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
