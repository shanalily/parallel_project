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
#define DEBUG_IS 1

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

#define EMPTY NULL
#define TURN_OPTS 3
#define RIGHT 0
#define STRGHT 1
#define LEFT 2

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

// each car at an intersection has three possibilities for how to turn.
// A value of 1 means they can turn, 0 means they can't.
// index 0 -> can turn right
// index 1 -> can go straight
// index 2 -> can turn left
typedef struct right_of_way {
    unsigned short nrth[TURN_OPTS];
    unsigned short soth[TURN_OPTS];
    unsigned short west[TURN_OPTS];
    unsigned short east[TURN_OPTS];
} intrsctn_rules;


/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// this checks the grid for validity
void check_grid(unsigned int rpr, intrsctn* is, 
        street* s_ew, street* s_ns, street* g_n_ns, street* g_s_ns);

void check_cars_start(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col);

int check_car_count(unsigned int rpr, intrsctn *is);

void mk_grid(unsigned int rpr, intrsctn* is, street* s_ns, street* s_ew, 
        street* g_ns_n, street* g_ns_s, int mpi_myrank, int mpi_commsize);

void generate_cars(unsigned long n, unsigned long glbl_row_idx, street* strts, 
        int row_or_col, float proportion);

int compare_sloc(car* c, unsigned long col, unsigned long row, unsigned short idx) {
    return c->s_col == col && c->s_row == row && c->s_idx == idx;
}
unsigned long dist_eloc(car* _c, unsigned long col, unsigned long row) {
    return abs(_c->e_col - col) + abs(_c->e_row - row);
}

MPI_Datatype mkcartype();

void pack_transfer(car* tt_n, street* gs_n, 
        car* tt_s, street* gs_s);
void unpack_transfer(car* tr_n, street* gs_n,
        car* tr_s, street* gs_s);

void transfer(car* tt_n, car* tr_n, car* tt_s, car* tr_s, int mpi_myrank, int mpi_commsize, 
        MPI_Datatype t_type);

void update_streets(unsigned int n, unsigned long glbl_row_idx, street *streets_now, street *streets_nxt,
        int row_or_col);
void update_ghost_streets(unsigned int n, street* ghost_now, street* ghost_nxt, int n_or_s);

int move_nrth(intrsctn *is, car *c, intrsctn_rules *r);
int move_soth(intrsctn *is, car *c, intrsctn_rules *r);
int move_west(intrsctn *is, car *c, intrsctn_rules *r);
int move_east(intrsctn *is, car *c, intrsctn_rules *r);
void nrth_to_soth_rules(intrsctn_rules *r);
void east_to_west_rules(intrsctn_rules *r);
void soth_to_nrth_rules(intrsctn_rules *r);
void west_to_east_rules(intrsctn_rules *r);
int reached_dest(unsigned long glbl_row_idx, unsigned long glbl_col_idx, unsigned short idx, car *c, unsigned short es);
void streets_check_dest(unsigned int n, unsigned long glbl_row_idx, street *streets, int row_or_col, int chck_nw, int check_se);

void update_intersections(unsigned int rpr, unsigned long glbl_index, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt);
void update_intersections2(unsigned int rpr, unsigned long glbl_index, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt, intrsctn_rules *r);

unsigned long total_grid_dist_to_travel(unsigned long glbl_row_idx, street* sts, unsigned int n);
unsigned long total_grid_dist_to_travel_ghost(unsigned long glbl_row_idx, street* sts, unsigned int n, int n_or_s);




unsigned long total_street_dist_to_travel();

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
    unsigned int glbl_index = 2*rpr*mpi_myrank;
    unsigned long num_ticks = 2;
    InitDefault();

    // To save space even ticks will be computed in now variables and odd
    // ticks will be computed in in nxt variables.
    // The pointers will switch after every tick
    
    // array of streets
    street* streets_ew_now = calloc(rpr*(SIDE_LENGTH-1), sizeof(street));
    street* streets_ew_nxt = calloc(rpr*(SIDE_LENGTH-1), sizeof(street));
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
    // struct to contain state of right of way
    intrsctn_rules r;

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
    check_cars_start(rpr*(SIDE_LENGTH-1), glbl_index, streets_ew_now, 0); 
    check_cars_start((rpr-1)*(SIDE_LENGTH), glbl_index+1, streets_ns_now, 1);
    if(mpi_myrank != 0)
        check_cars_start(SIDE_LENGTH, glbl_index-1, ghost_ns_nrth_now, 1);
    if(mpi_myrank != mpi_commsize-1)
        check_cars_start(SIDE_LENGTH, glbl_index+2*rpr-1, ghost_ns_soth_now, 1);
    int car_count = check_car_count(rpr, intrsctn_now);
#endif
    
    MPI_Barrier( MPI_COMM_WORLD );

    unsigned long dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
    unsigned long dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
    unsigned long dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
    unsigned long dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
    unsigned long dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
    printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, 0, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);

    // intrsctn *cur_intrsctn_now, *cur_intrsctn_nxt;
    // street *cur_streets_ew_now, *cur_streets_ns_now, *cur_streets_ew_nxt, *cur_streets_ns_nxt;
    // for loop of ticks
    for (size_t i = 1; i < num_ticks; ++i) {
        

        streets_check_dest(rpr*(SIDE_LENGTH-1), glbl_index, streets_ew_now, 0, 1, 1);
        streets_check_dest((rpr-1)*SIDE_LENGTH, glbl_index, streets_ns_now, 1, 1, 1);
        streets_check_dest(SIDE_LENGTH, glbl_index-1, ghost_ns_nrth_now, 0, 0, 1);
        streets_check_dest(SIDE_LENGTH, glbl_index+2*rpr-1, ghost_ns_soth_now, 0, 1, 0);
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_nxt, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_nxt, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_nxt, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_nxt, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left nxt is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        
        
        // send and receive (blocking)
        pack_transfer(to_transfer_nrth, ghost_ns_nrth_now,
                to_transfer_soth, ghost_ns_soth_now);
        transfer(to_transfer_nrth, to_receive_nrth, to_transfer_soth, to_receive_soth,
                mpi_myrank, mpi_commsize, mpi_car_type);
        unpack_transfer(to_receive_nrth, ghost_ns_nrth_now,
                to_receive_soth, ghost_ns_soth_now);
        
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_nxt, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_nxt, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_nxt, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_nxt, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left nxt is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // update streets
        // update east/west
        // update_streets(rpr*(SIDE_LENGTH-1), glbl_index, streets_ew_now, streets_ew_nxt, 0); // something in this function is wrong
        // // update north/south
        // update_streets((rpr-1)*SIDE_LENGTH, glbl_index, streets_ns_now, streets_ns_nxt, 1);
        // // printf("Rank %d f\n", mpi_myrank);
        // // update ghost rows
        // update_ghost_streets(SIDE_LENGTH, ghost_ns_nrth_now, ghost_ns_nrth_nxt, 0);
        // // printf("Rank %d g\n", mpi_myrank);
        // update_ghost_streets(SIDE_LENGTH, ghost_ns_soth_now, ghost_ns_soth_nxt, 1);

        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_nxt, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_nxt, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_nxt, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_nxt, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left nxt is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // printf("Rank %d h\n", mpi_myrank);
        // run intersections
        update_intersections(rpr, glbl_index, intrsctn_now, intrsctn_nxt);

        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_nxt, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_nxt, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_nxt, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_nxt, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left nxt is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);
        // printf("Rank %d i\n", mpi_myrank);

        street* tmp;
        intrsctn* tmp2;

        // do the exchange
        tmp =  streets_ew_now;
        streets_ew_now = streets_ew_nxt;
        streets_ew_nxt = tmp;
        tmp =  streets_ns_now;
        streets_ns_now = streets_ns_nxt;
        streets_ns_nxt = tmp;
        tmp = ghost_ns_nrth_now;
        ghost_ns_nrth_now = ghost_ns_nrth_nxt;
        ghost_ns_nrth_nxt = tmp;
        tmp = ghost_ns_soth_now;
        ghost_ns_soth_now = ghost_ns_soth_nxt;
        ghost_ns_soth_nxt = tmp;
        // array of intersections
        tmp2 = intrsctn_now;
        intrsctn_now = intrsctn_nxt;
        intrsctn_nxt = tmp2;
        


        // dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
        // dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
        // dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
        // dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
        // dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
        // printf("Rank %d: Tick %d: Cars left is %lu %lu %lu %lu %lu\n", mpi_myrank, i, dist_left_ew, dist_left_ns, dist_left_gn, dist_left_gs, dist_left);

    }
#ifdef DEBUG
    check_grid(rpr, intrsctn_now, streets_ew_now, streets_ns_now, 
            ghost_ns_nrth_now, ghost_ns_soth_nxt);
    check_grid(rpr, intrsctn_nxt, streets_ew_nxt, streets_ns_nxt, 
            ghost_ns_nrth_nxt, ghost_ns_soth_nxt);
    assert(car_count == check_car_count(rpr, intrsctn_nxt));
#endif

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

void check_cars_start(unsigned long n, unsigned long glbl_row_idx, street* strts, 
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
            assert(!strts[i].go_es[ROAD_CAP-1-j] || !strts[i].go_wn[j]);
            // if a street has a car it should be at the starting location
            assert(!strts[i].go_es[j] || (strts[i].go_es[j] && 
                    compare_sloc(strts[i].go_es[j], 
                    col, row, j)));
            assert(!strts[i].go_wn[j] || (strts[i].go_wn[j] && 
                    compare_sloc(strts[i].go_wn[j], 
                    col, row, ROAD_CAP-1-j)));
            
        }
        
    }
    
}

int check_row_count(intrsctn *is, int n_i) {
    unsigned long count = 0;
    for(size_t i = 0; i < n_i; i++) {
        // assert(is[i].east == is[i+1].west);
        if (is[i].east != EMPTY) count += 1;
    }
    return count;
}

int check_col_count(intrsctn *is, int n_i) {
    unsigned long count = 0;
    for(size_t i = 0; i < n_i; i++) {
        // assert(is[i].soth == is[i+SIDE_LENGTH].nrth);
        if (is[i].soth != EMPTY) count += 1;
    }
    return count;
}

// for testing. but what if cars leave?
int check_car_count(unsigned int rpr, intrsctn *is) {
    // TODO: implement
    unsigned long car_count = 0;
    for(size_t i = 0; i < rpr; i++)
    {
        car_count += check_row_count(&is[i*SIDE_LENGTH], SIDE_LENGTH);
    }

    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        car_count += check_col_count(&is[i], rpr);
    }
    return car_count; // can use this for assertion, should be same before and after update_intrsctns
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
                c->s_idx = (unsigned short) j;

                // end location (non-unique)
                c->e_col = (unsigned long) (GenVal(row)*(2*SIDE_LENGTH-1));
                c->e_row = 2*((unsigned long) (GenVal(row)*(SIDE_LENGTH-1)))+!(c->e_col%2);
                c->e_idx = (int) (GenVal(row)*(ROAD_CAP));
                #ifdef DEBUG
                assert(c->e_col%2 != c->e_row%2);
                #endif
                // printf("%d %d %d %d %d %d\n", c->e_col, c->e_row, c->s_idx, c->s_col, c->s_row, c->s_idx);
                if(c->e_row == c->s_row && c->e_col == c->s_col && c->e_idx >= c->s_idx) {
                    strts[i].go_es[j] = c;
                }
                else if(c->e_row == c->s_row && c->e_col == c->s_col && c->e_idx < c->s_idx) {
                    strts[i].go_wn[ROAD_CAP-1-j] = c;
                }
                else if ((c->e_row > c->s_row && !row_or_col) ||
                        (c->e_col > c->s_col && row_or_col)) {
                    strts[i].go_es[j] = c;
                }else{
                    strts[i].go_wn[ROAD_CAP-1-j] = c;
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
    car sample;
    sample.s_col = 0;
    sample.s_row = 0;
    sample.s_idx = 0;
    sample.e_col = 0;
    sample.e_row = 0;
    sample.e_idx = 0;
    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        if (gs_n[i].go_wn[0]) {
            memcpy(&(tt_n[i]), gs_n[i].go_wn[0], sizeof(car));
            free(gs_n[i].go_wn[0]);
        }
        else{
            memcpy(&(tt_n[i]), &sample, sizeof(car));
        }
        if(gs_s[i].go_es[0]){
            memcpy(&(tt_s[i]), gs_s[i].go_es[0], sizeof(car));
            free(gs_s[i].go_es[0]);
        }
        else {
            memcpy(&(tt_s[i]), &sample, sizeof(car));
        }
    }
    

}
void unpack_transfer(car* tr_n, street* gs_n,
        car* tr_s, street* gs_s){
    for(size_t i = 0; i < SIDE_LENGTH; i++)
    {
        if (tr_n[i].s_col != 0 || tr_n[i].s_row != 0) {
            gs_n[i].go_es[0] = calloc(1, sizeof(car));
            memcpy(gs_n[i].go_es[0], &tr_n[i], sizeof(car));
        }
        else {
            gs_n[i].go_es[0] = NULL;
        }
        if (tr_s[i].s_col != 0 || tr_s[i].s_row != 0) {
            gs_s[i].go_wn[0] = calloc(1, sizeof(car));
            memcpy(gs_s[i].go_wn[0], &tr_s[i], sizeof(car));
        }
        else {
            gs_s[i].go_wn[0] = NULL;
        }
    }
}

void transfer(car* tt_n, car* tr_n, car* tt_s, car* tr_s, int mpi_myrank, int mpi_commsize,
        MPI_Datatype t_type){
    MPI_Request nrth;
    MPI_Request soth;
    MPI_Request request;

    // printf("Rank %d 1\n", mpi_myrank);
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
    // printf("Rank %d 2\n", mpi_myrank);

    MPI_Type_commit(&mpi_car_type);
    // printf("Rank %d 3\n", mpi_myrank);

    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    if(mpi_myrank != 0) {
        MPI_Irecv(&tr_n[0], SIDE_LENGTH, mpi_car_type, mpi_myrank-1, 0, MPI_COMM_WORLD, &nrth);
    }
    if(mpi_myrank != mpi_commsize -1) {
        MPI_Irecv(&tr_s[0], SIDE_LENGTH, mpi_car_type, mpi_myrank+1, 1, MPI_COMM_WORLD, &soth);
    }
    if(mpi_myrank != mpi_commsize -1) {
        MPI_Isend(&tt_s[0], SIDE_LENGTH, mpi_car_type, mpi_myrank+1, 0, MPI_COMM_WORLD, &request);
    }
    if(mpi_myrank != 0) {
        MPI_Isend(&tt_n[0], SIDE_LENGTH, mpi_car_type, mpi_myrank-1, 1, MPI_COMM_WORLD, &request);
    }
    // printf("Rank %d 4\n", mpi_myrank);

    if(mpi_myrank != 0) MPI_Wait(&nrth, MPI_STATUS_IGNORE);
    if(mpi_myrank != mpi_commsize -1) MPI_Wait(&soth, MPI_STATUS_IGNORE);
    // printf("Rank %d 5\n", mpi_myrank);
    MPI_Type_free(&mpi_car_type);
    // printf("Rank %d 6\n", mpi_myrank);
}

void streets_check_dest(unsigned int n, unsigned long glbl_row_idx, street *streets, int row_or_col, int check_nw, int check_se) {
    unsigned long row_idx, col_idx;
    for (size_t i =0; i<n; ++i) {
        row_idx = glbl_row_idx + 2*(i/(SIDE_LENGTH-(!row_or_col))); // global row of block of street. -glbl_col_idx?
        col_idx = 2*(i%(SIDE_LENGTH-!row_or_col)) + !row_or_col; // global column of block of street. -glbl_col_idx?
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            // if (streets[i].go_es[j]) {
            //     printf("%d %d %d %d %d %d\n", row_idx, col_idx, j, streets[i].go_es[j]->e_col, streets[i].go_es[j]->e_row, streets[i].go_es[j]->e_idx);
            // }
            // if (streets[i].go_wn[j]) {
            //     printf("%d %d %d %d %d %d\n", row_idx, col_idx, j, streets[i].go_wn[j]->e_col, streets[i].go_wn[j]->e_row, streets[i].go_wn[j]->e_idx);
            // }
            
            if (check_se && streets[i].go_es[j] && reached_dest(row_idx, col_idx, j, streets[i].go_es[j], 1)) {
                printf("reached destination! %lu %lu %lu\n", row_idx, col_idx, j);
                free(streets[i].go_es[j]);
                streets[i].go_es[j] = NULL;
            }
            if (check_nw && streets[i].go_wn[j] && reached_dest(row_idx, col_idx, j, streets[i].go_wn[j], 0)) {
                printf("reached destination! %lu %lu %lu\n", row_idx, col_idx, j);
                free(streets[i].go_wn[j]);
                streets[i].go_wn[j] = NULL;
            }
            
        }
        
    }
}

// n is either rpr*(SIDE_LENGTH-1) or (rpr-1)*SIDE_LENGTH, which is the number of cells in either the n/s or e/w streets
// Move everything down the streets. Slot ROAD_CAP-1 is the location on the street closest
// to the intersection at the end of the block.
// 0           ----> ROAD_CAP-1
// ROAD_CAP-1 <----  0
void update_streets(unsigned int n, unsigned long glbl_row_idx, street *streets_now, street *streets_nxt,
        int row_or_col) { // row_or_col is 0 if e/w and 1 if n/s but I don't know why...
    // unsigned long row_idx, col_idx;
    for (size_t i = 0; i < n; ++i) {
        // if location on street is empty, move up the next car (if it exists) from previous location
        // set previous location to empty
        // row_idx = glbl_row_idx + 2*(i/(SIDE_LENGTH-(!row_or_col))); // global row of block of street. -glbl_col_idx?
        // col_idx = 2*(i%(SIDE_LENGTH-!row_or_col)) + !row_or_col; // global column of block of street. -glbl_col_idx?

        for (unsigned int j = ROAD_CAP-1; j > 0; --j) {
            if (streets_now[i].go_es[j] != EMPTY) {
                if(j < ROAD_CAP-1 && streets_now[i].go_es[j+1] == EMPTY){
                    streets_nxt[i].go_es[j+1] = streets_now[i].go_es[j];
                    streets_now[i].go_es[j] = EMPTY;
                }
                else {
                    streets_nxt[i].go_es[j] = streets_now[i].go_es[j];
                    streets_now[i].go_es[j] = EMPTY;
                }
            }
            if (streets_now[i].go_wn[j] != EMPTY) {
                if(j < ROAD_CAP-1 && streets_now[i].go_wn[j+1] == EMPTY){
                    streets_nxt[i].go_wn[j+1] = streets_now[i].go_wn[j];
                    streets_now[i].go_wn[j] = EMPTY;
                }
                else {
                    streets_nxt[i].go_wn[j] = streets_now[i].go_wn[j];
                    streets_now[i].go_wn[j] = EMPTY;
                }
            }
        }
    }
}

// can I reach my destination here, or can I always assume another rank will handle it?
void update_ghost_streets(unsigned int n, street* ghost_now, street* ghost_nxt, int n_or_s){
    // n_or_s == 0 -> n 1->s
    for (size_t i = 0; i < n; ++i) {
        // if location on street is empty, move up the next car (if it exists) from previous location
        // set previous location to empty
        for (size_t j = ROAD_CAP-1; j >= 1; --j) {
            if (n_or_s && ghost_now[i].go_es[j] == EMPTY) {
                ghost_nxt[i].go_es[j] = ghost_now[i].go_es[j-1];
                ghost_now[i].go_es[j-1] = EMPTY;
            }
            if (!n_or_s && ghost_now[i].go_wn[j] == EMPTY) {
                ghost_nxt[i].go_wn[j] = ghost_now[i].go_wn[j-1];
                ghost_now[i].go_wn[j-1] = EMPTY;
            }
        }
        // last slot of intrsctn_now will have it's previous value still, can be used for update_intersections.
    }
}

// If car has reached it's end point, return 1. Otherwise, return 0.
// If row is even, EW/WE. If row is odd, SN/NS. It doesn't matter which side of street the car is on, but index needs to be adjusted accordingly.
// input should be the even glbl_row_idx and even glbl_col_idx from update_intersections. Assume e_idx if is pointing west or north.
int reached_dest(unsigned long glbl_row_idx, unsigned long glbl_col_idx, unsigned short idx, car *c, unsigned short es) {
    if (glbl_row_idx == c->e_row && glbl_col_idx == c->e_col
            && ((es && idx == c->e_idx) || (!es && idx == ROAD_CAP-1-c->e_idx))) {
        return 1;
    }
    return 0;
}

// Set first location on northern street to the car being moved
int move_nrth(intrsctn *is, car *c, intrsctn_rules *r) {
    is->nrth->go_wn[0] = c;
    r->west[LEFT] = 0; // car from west can't go left
    r->east[RIGHT] = 0; // car from east can't go right
    r->soth[STRGHT] = 0;
    return 1;
}

// Set first location on southern street to the car being moved
int move_soth(intrsctn *is, car *c, intrsctn_rules *r) {
    is->soth->go_es[0] = c;
    r->west[RIGHT] = 0;
    r->east[LEFT] = 0;
    r->nrth[STRGHT] = 0;
    return 1;
}

// Set first location on western street to the car being moved
int move_west(intrsctn *is, car *c, intrsctn_rules *r) {
    is->west->go_wn[0] = c;
    r->nrth[RIGHT] = 0;
    r->soth[LEFT] = 0;
    r->east[STRGHT] = 0;
    return 1;
}

// Set first location on eastern street to the car being moved
int move_east(intrsctn *is, car *c, intrsctn_rules *r) {
    is->east->go_es[0] = c;
    r->nrth[LEFT] = 0;
    r->soth[RIGHT] = 0;
    r->west[STRGHT] = 0;
    return 1;
}

// Set and reset intersection rules.
void reset_intrsctn(intrsctn intrsctn_now, intrsctn intrsctn_nxt, intrsctn_rules *r) {
    unsigned short road_start = 0, road_end = ROAD_CAP-1;
    if (intrsctn_nxt.nrth) {
        intrsctn_nxt.nrth->go_wn[road_start] = intrsctn_now.nrth->go_wn[road_start]; // copy over cars that don't change. Do these even exist?
        intrsctn_nxt.nrth->go_es[road_end] = EMPTY;
    }
    if (intrsctn_nxt.soth) {
        intrsctn_nxt.soth->go_es[road_start] = intrsctn_now.soth->go_es[road_start];
        intrsctn_nxt.soth->go_wn[road_end] = EMPTY;
    }
    if (intrsctn_nxt.west) {
        intrsctn_nxt.west->go_wn[road_start] = intrsctn_now.west->go_wn[road_start];
        intrsctn_nxt.west->go_es[road_end] = EMPTY;
    }
    if (intrsctn_nxt.east) {
        intrsctn_nxt.east->go_es[road_start] = intrsctn_now.east->go_es[road_start];
        intrsctn_nxt.east->go_wn[road_end] = EMPTY;
    }
    // shouldn't be valid if there are four cars on road already. Look at current intersection.
    // CAN TURN RIGHT
    r->nrth[RIGHT] = (intrsctn_now.west && intrsctn_now.west->go_wn[0] == EMPTY) ? 1 : 0;
    r->soth[RIGHT] = (intrsctn_now.east && intrsctn_now.east->go_es[0] == EMPTY) ? 1 : 0;
    r->west[RIGHT] = (intrsctn_now.soth && intrsctn_now.soth->go_es[0] == EMPTY) ? 1 : 0;
    r->east[RIGHT] = (intrsctn_now.nrth && intrsctn_now.nrth->go_wn[0] == EMPTY) ? 1 : 0;
    // CAN GO STRAIGHT
    r->nrth[STRGHT] = (intrsctn_now.soth && intrsctn_now.soth->go_es[0] == EMPTY) ? 1 : 0;
    r->soth[STRGHT] = (intrsctn_now.nrth && intrsctn_now.nrth->go_wn[0] == EMPTY) ? 1 : 0;
    r->west[STRGHT] = (intrsctn_now.east && intrsctn_now.east->go_es[0] == EMPTY) ? 1 : 0;
    r->east[STRGHT] = (intrsctn_now.west && intrsctn_now.west->go_wn[0] == EMPTY) ? 1 : 0;
    // CAN GO LEFT
    r->nrth[LEFT] = (intrsctn_now.east && intrsctn_now.east->go_es[0] == EMPTY) ? 1 : 0;
    r->soth[LEFT] = (intrsctn_now.west && intrsctn_now.west->go_wn[0] == EMPTY) ? 1 : 0;
    r->west[LEFT] = (intrsctn_now.nrth && intrsctn_now.nrth->go_wn[0] == EMPTY) ? 1 : 0;
    r->east[LEFT] = (intrsctn_now.soth && intrsctn_now.soth->go_es[0] == EMPTY) ? 1 : 0;
}

// If a car from north goes straight south, other cars are restricted in the following ways.
void nrth_to_soth_rules(intrsctn_rules *r) {
    r->east[STRGHT] = 0; // east can't go straight
    r->soth[LEFT] = 0; // south can't go left
    r->west[STRGHT] = 0; // west can't go straight or left
    r->west[LEFT] = 0;
}

// If a car from east goes straight west, other cars are restricted in the following ways.
void east_to_west_rules(intrsctn_rules *r) {
    r->soth[STRGHT] = 0; // south can't go straight
    r->west[LEFT] = 0; // west can't go left
    r->nrth[STRGHT] = 0; // north can't go straight or left
    r->nrth[LEFT] = 0;
}

// If a car from south goes straight north, other cars are restricted in the following ways.
void soth_to_nrth_rules(intrsctn_rules *r) {
    r->west[STRGHT] = 0; // west can't go straight
    r->nrth[LEFT] = 0; // north can't go left
    r->east[STRGHT] = 0; // east can't go straight or left
    r->east[LEFT] = 0;
}

// If a car from west goes straight east, other cars are restricted in the following ways.
void west_to_east_rules(intrsctn_rules *r) {
    r->nrth[STRGHT] = 0; // north can't go straight
    r->east[LEFT] = 0; // east can't go left
    r->soth[STRGHT] = 0; // south can't go straight or left
    r->soth[LEFT] = 0;
}

// For debugging. This might not be helpful, idk. Will delete later.
// I think I have go_wn and go_es mixed up ugh.
void print_intersection(intrsctn intersection) {
    printf("-------\n");
    if (intersection.nrth && intersection.nrth->go_wn[0]) { // going north
        if (intersection.nrth->go_es[ROAD_CAP-1]) printf("  sn   \n"); // heading south
        else printf("   n  \n");
    }
    else if (intersection.nrth && intersection.nrth->go_es[ROAD_CAP-1]) printf("  s    \n");
    else printf("      \n");

    if (intersection.west && intersection.west->go_wn[0]) { // heading west
        if (intersection.east && intersection.east->go_wn[ROAD_CAP-1]) printf("w    w\n"); // coming from east, going west
        else printf("w     \n");
    }
    else if (intersection.east && intersection.east->go_wn[ROAD_CAP-1]) printf("     w\n");
    else printf("      \n");

    if (intersection.east && intersection.east->go_es[0]) { // coming from west, going east
        if (intersection.west && intersection.west->go_es[ROAD_CAP-1]) printf("e    e\n");
        else printf("     e\n");
    }
    else if (intersection.west && intersection.west->go_es[ROAD_CAP-1]) printf("e     \n");
    else printf("      \n");

    if (intersection.soth && intersection.soth->go_wn[ROAD_CAP-1]) { // coming from south, heading north
        if (intersection.soth->go_es[0])
            printf("  sn  \n");
        else
            printf("   n  \n");
    } else if (intersection.soth && intersection.soth->go_es[0]) { // coming from north, going south
        printf("  s   \n");
    }
    else printf("      \n");
    printf("-------\n");
}

void update_intersections(unsigned int rpr, unsigned long glbl_row_idx, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt) {
    unsigned long row_idx, col_idx;
    unsigned int block_end = ROAD_CAP-1;
    // int nrth, soth, east, west;
    int strt_ew = 0, strt_ns = 0;
    int left_ne_sw = 0, left_wn_es = 0;
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    for (size_t i = 0; i < rpr*SIDE_LENGTH; ++i) {
        row_idx = glbl_row_idx + 2*(i/SIDE_LENGTH);
        col_idx = 2*(i%SIDE_LENGTH);
        strt_ew = 0;
        strt_ns = 0;
        left_ne_sw = 0;
        left_wn_es = 0;
    #ifdef DEBUG_IS
        if(rank == 0){
            printf("before %d %d %d\n", i, row_idx, col_idx);
            print_intersection(intrsctn_now[i]);
        }
    #endif
        // right turns
        // north going to the west
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_es[block_end] && 
                intrsctn_nxt[i].west && !intrsctn_nxt[i].west->go_wn[0] && 
                intrsctn_now[i].nrth->go_es[block_end]->e_col < col_idx){
            intrsctn_nxt[i].west->go_wn[0] = intrsctn_now[i].nrth->go_es[block_end];
            intrsctn_now[i].nrth->go_es[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("North taking a right turn to the west\n");
            }
        #endif
        }
        // south going to the east
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_wn[block_end] && 
                intrsctn_nxt[i].east && !intrsctn_nxt[i].east->go_es[0] && 
                intrsctn_now[i].soth->go_wn[block_end]->e_col > col_idx){
            intrsctn_nxt[i].east->go_es[0] = intrsctn_now[i].soth->go_wn[block_end];
            intrsctn_now[i].soth->go_wn[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("South taking a right turn to the east\n");
            }
        #endif
        }
        // west going to the south
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_es[block_end] && 
                intrsctn_nxt[i].soth && !intrsctn_nxt[i].soth->go_es[0] && 
                intrsctn_now[i].west->go_es[block_end]->e_row > row_idx){
            intrsctn_nxt[i].soth->go_es[0] = intrsctn_now[i].west->go_es[block_end];
            intrsctn_now[i].west->go_es[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("West taking a right turn to the south\n");
            }
        #endif
        }
        // east going to the north
        // east going south
        // if(i == 8 && rank == 0){
        //     printf("%d %d %d %d %d\n", intrsctn_now[i].east, intrsctn_now[i].east->go_wn[block_end], 
        //             intrsctn_nxt[i].nrth, !intrsctn_nxt[i].nrth->go_wn[0], intrsctn_now[i].east->go_wn[block_end]->e_row < row_idx);
        // }
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_wn[block_end] && 
                intrsctn_nxt[i].nrth && !intrsctn_nxt[i].nrth->go_wn[0] && 
                intrsctn_now[i].east->go_wn[block_end]->e_row < row_idx){
            intrsctn_nxt[i].nrth->go_wn[0] = intrsctn_now[i].east->go_wn[block_end];
            intrsctn_now[i].east->go_wn[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("East taking a right turn to the north\n");
            }
        #endif
        }

        // straight turns
        // north going south
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_es[block_end] && 
                intrsctn_nxt[i].soth && !intrsctn_nxt[i].soth->go_es[0] && 
                !strt_ew && intrsctn_now[i].nrth->go_es[block_end]->e_row > row_idx){
            intrsctn_nxt[i].soth->go_es[0] = intrsctn_now[i].nrth->go_es[block_end];
            intrsctn_now[i].nrth->go_es[block_end] = NULL;
            strt_ns = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("north taking a straight turn to the south\n");
            }
        #endif
        }
        // south going north
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_wn[block_end] && 
                intrsctn_nxt[i].nrth && !intrsctn_nxt[i].nrth->go_wn[0] && 
                !strt_ew && intrsctn_now[i].soth->go_wn[block_end]->e_row < row_idx){
            intrsctn_nxt[i].nrth->go_wn[0] = intrsctn_now[i].soth->go_wn[block_end];
            intrsctn_now[i].soth->go_wn[block_end] = NULL;
            strt_ns = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("south taking a straight turn to the north\n");
            }
        #endif
        }
        // west going east
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_es[block_end] && 
                intrsctn_nxt[i].east && !intrsctn_nxt[i].east->go_es[0] && 
                !strt_ns && intrsctn_now[i].west->go_es[block_end]->e_col > col_idx){
            intrsctn_nxt[i].east->go_es[0] = intrsctn_now[i].west->go_es[block_end];
            intrsctn_now[i].west->go_es[block_end] = NULL;
            strt_ew = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("west taking a straight turn to the east\n");
            }
        #endif
        }
        // east going west
        // if(i == 30 && rank == 0){
        //     printf("%d %d %d %d %d %d\n", intrsctn_now[i].east, intrsctn_now[i].east->go_wn[block_end], 
        //             intrsctn_nxt[i].west, !intrsctn_nxt[i].west->go_wn[0], !strt_ns, intrsctn_now[i].east->go_wn[block_end]->e_col < col_idx);
        // }
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_wn[block_end] && 
                intrsctn_nxt[i].west && !intrsctn_nxt[i].west->go_wn[0] && 
                !strt_ns && intrsctn_now[i].east->go_wn[block_end]->e_col < col_idx){
            intrsctn_nxt[i].west->go_wn[0] = intrsctn_now[i].east->go_wn[block_end];
            intrsctn_now[i].east->go_wn[block_end] = NULL;
            strt_ew = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("east taking a straight turn to the west\n");
            }
        #endif
        }

        // left turns
        // north going east
        // if(i == 9 && rank == 0){
        //     printf("%d %d %d %d %d %d %d\n", intrsctn_now[i].nrth, intrsctn_now[i].nrth->go_es[block_end], 
        //             intrsctn_nxt[i].east, !intrsctn_nxt[i].east->go_es[0], !left_wn_es, !strt_ns, !strt_ew, intrsctn_now[i].nrth->go_es[block_end]->e_col > col_idx);
        // }
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_es[block_end] && 
                intrsctn_nxt[i].east && !intrsctn_nxt[i].east->go_es[0] && 
                !left_wn_es && !strt_ns && !strt_ew && intrsctn_now[i].nrth->go_es[block_end]->e_col > col_idx){
            intrsctn_nxt[i].east->go_es[0] = intrsctn_now[i].nrth->go_es[block_end];
            intrsctn_now[i].nrth->go_es[block_end] = NULL;
            left_ne_sw = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("north taking a left turn to the east\n");
            }
        #endif
        }
        // south going west
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_wn[block_end] && 
                intrsctn_nxt[i].west && !intrsctn_nxt[i].west->go_wn[0] && 
                !left_wn_es && !strt_ns && !strt_ew && intrsctn_now[i].soth->go_wn[block_end]->e_col < col_idx){
            intrsctn_nxt[i].west->go_wn[0] = intrsctn_now[i].soth->go_wn[block_end];
            intrsctn_now[i].soth->go_wn[block_end] = NULL;
            left_ne_sw = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("south taking a left turn to the west\n");
            }
        #endif
        }
        // west going north
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_es[block_end] && 
                intrsctn_nxt[i].nrth && !intrsctn_nxt[i].nrth->go_wn[0] && 
                !left_ne_sw && !strt_ns && !strt_ew && intrsctn_now[i].west->go_es[block_end]->e_row < row_idx){
            intrsctn_nxt[i].nrth->go_wn[0] = intrsctn_now[i].west->go_es[block_end];
            intrsctn_now[i].west->go_es[block_end] = NULL;
            left_wn_es = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("west taking a left turn to the north\n");
            }
        #endif
        }
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_wn[block_end] && 
                intrsctn_nxt[i].soth && !intrsctn_nxt[i].soth->go_es[0] && 
                !left_ne_sw && !strt_ns && !strt_ew && intrsctn_now[i].east->go_wn[block_end]->e_row > row_idx){
            intrsctn_nxt[i].soth->go_es[0] = intrsctn_now[i].east->go_wn[block_end];
            intrsctn_now[i].east->go_wn[block_end] = NULL;
            left_wn_es = 1;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("east taking a left turn to the south\n");
            }
        #endif
        }

        // staying still
        // north
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_es[block_end]){
            intrsctn_nxt[i].nrth->go_es[block_end] = intrsctn_now[i].nrth->go_es[block_end];
            intrsctn_now[i].nrth->go_es[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("north staying still, end dest %d %d\n", 
                        intrsctn_nxt[i].nrth->go_es[block_end]->e_row, intrsctn_nxt[i].nrth->go_es[block_end]->e_col);
            }
        #endif
        }
        // south
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_wn[block_end]){
            intrsctn_nxt[i].soth->go_wn[block_end] = intrsctn_now[i].soth->go_wn[block_end];
            intrsctn_now[i].soth->go_wn[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("south staying still, end dest %d %d\n", 
                        intrsctn_nxt[i].soth->go_wn[block_end]->e_row, intrsctn_nxt[i].soth->go_wn[block_end]->e_col);
            }
        #endif
        }
        // west
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_es[block_end]){
            intrsctn_nxt[i].west->go_es[block_end] = intrsctn_now[i].west->go_es[block_end];
            intrsctn_now[i].west->go_es[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("west staying still, end dest %d %d\n",
                        intrsctn_nxt[i].west->go_es[block_end]->e_row, intrsctn_nxt[i].west->go_es[block_end]->e_col);
            }
        #endif
        }
        // east
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_wn[block_end]){
            intrsctn_nxt[i].east->go_wn[block_end] = intrsctn_now[i].east->go_wn[block_end];
            intrsctn_now[i].east->go_wn[block_end] = NULL;
        #ifdef DEBUG_IS
            if(rank==0){
                printf("east staying still, end dest %d %d\n",
                    intrsctn_nxt[i].east->go_wn[block_end]->e_row, intrsctn_nxt[i].east->go_wn[block_end]->e_col);
            }
        #endif
        }

    #ifdef DEBUG_IS
        if(rank==0){
            printf("after\n");
            print_intersection(intrsctn_nxt[i]);
        }
    #endif
    }
}

// make sure to switch which intersections are passed between now and nxt.
// EW streets have even rows numbers and odd columns numbers
// NS streets have odd rows numbers and even column numbers
void update_intersections2(unsigned int rpr, unsigned long glbl_row_idx, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt, intrsctn_rules *r) {
    // update nxt based on now. Determine who can move first, corresponding to traffic rules.
    unsigned int block_end = ROAD_CAP-1;
    int left_turn = 0;
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);

    for (size_t i = 0; i < rpr*SIDE_LENGTH; ++i) {
        // printf("Rank %d: At %d\n", rank, i);
        unsigned long row_idx = glbl_row_idx + 2*(i/SIDE_LENGTH); // index is for a horizontal street
        unsigned long col_idx = 2*(i%SIDE_LENGTH);
        reset_intrsctn(intrsctn_now[i], intrsctn_nxt[i], r);
        
        if(rank == 0){
            printf("before\n");
            print_intersection(intrsctn_now[i]);
        }
        
        // MOVE CAR ON NORTHERN SIDE OF INTERSECTION
        // Check where car is heading.
        // prioritize going straight, then right, then left, depending on direction of destination
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_es[block_end]) {

            car *c = intrsctn_now[i].nrth->go_es[block_end];
            // coming from north, should they turn straight, right, or left?
            // check | c->e_row == row_idx+1? row_idx+1 is NS
            // left or right? or arrived? arrived should be checked at end of this function though
            // sooo what is row AND column equal? then I either need to cross to east or south, depending on if odd or even?
            // north: If even, turn left. If odd, stay straight.
            if ((c->e_row == row_idx || c->e_row == row_idx+1) && c->e_col < col_idx && r->nrth[RIGHT]) { // go west/right?.   Is my problem that streets aren't correctly set for nxt?
                move_west(&intrsctn_nxt[i], c, r);
                // now other cars can't take that spot, make sure it's marked!
                if(rank == 0) printf("north moved right\n");
            } else if ((c->e_row == row_idx || c->e_row == row_idx+1) && r->nrth[LEFT]) { // go east/left? could be at destination street.
                move_east(&intrsctn_nxt[i], c, r); // should I be yielding the right of way?
                left_turn = 1; // no other cars can move
                if(rank == 0) printf("north moved left\n");
            }
            else if (r->nrth[STRGHT]) {
                move_soth(&intrsctn_nxt[i], c, r); // go straight?
                nrth_to_soth_rules(r);
                if(rank == 0) printf("north moved straight\n");
            }
            else {
                intrsctn_nxt[i].nrth->go_es[block_end] = c; // stay in place
            }

        }

        // MOVE CAR ON EASTERN SIDE OF INTERSECTION
        // sooo are my row and column checks actually going to work
        // when row and column are equal:
        // east: If row is even, doesn't make sense, gone too far. If row is odd, turn left?
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_wn[block_end]) {
            car *c = intrsctn_now[i].east->go_wn[block_end];
            // north or south?
            if (!left_turn && c->e_col == col_idx && c->e_row < row_idx && r->east[RIGHT]) { // north/right?
                move_nrth(&intrsctn_nxt[i], c, r);
                if(rank == 0) printf("east moved right\n");
            } else if (!left_turn && c->e_col == col_idx && c->e_row > row_idx && r->east[LEFT]) { // shouldn't turn if rows are equal, on street already
                move_soth(&intrsctn_nxt[i], c, r);
                left_turn = 1;
                if(rank == 0) printf("east moved left\n");
            }
            else if (!left_turn && r->east[STRGHT]) {
                if(rank == 0) printf("east moved straight\n");
                move_west(&intrsctn_nxt[i], c, r);
                east_to_west_rules(r);
            }
            else {
                intrsctn_nxt[i].east->go_wn[block_end] = c;
            }

        }

        // MOVE CAR ON SOUTH SIDE OF INTERSECTION
        // If rows and columns equal:
        // south: If row is even, turn right. If row is odd, doesn't make sense, gone too far.
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_wn[block_end]) {

            car *c = intrsctn_now[i].soth->go_wn[block_end];
            // east or west?
            if (!left_turn && (c->e_row == row_idx || c->e_row == row_idx+1) && c->e_col <= col_idx && r->soth[RIGHT]) { // go east/right?
                move_east(&intrsctn_nxt[i], c, r);
                if(rank == 0) printf("south moved right\n");
            } else if (!left_turn && (c->e_row == row_idx || c->e_row == row_idx+1) && r->soth[LEFT]) {
                move_west(&intrsctn_nxt[i], c, r);
                left_turn = 1;
                if(rank == 0) printf("south moved left\n");
            }
            else if (!left_turn && r->soth[STRGHT]) {
                if(rank == 0)printf("south moved straight\n");
                move_nrth(&intrsctn_nxt[i], c, r);
                soth_to_nrth_rules(r);
            }
            else {
                intrsctn_nxt[i].soth->go_wn[block_end] = c;
            }

        }

        // MOVE CAR ON WESTERN SIDE OF INTERSECTION
        // When row and column equal:
        // west: If row is even, go straight. If row is odd, turn right.
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_es[block_end]) { // coming from west
            // I can't go
            car *c = intrsctn_now[i].west->go_es[block_end];
            // up or down?
            if (!left_turn && c->e_col == col_idx && c->e_row > row_idx && r->west[RIGHT]) { // south?
                move_soth(&intrsctn_nxt[i], c, r);
                if(rank == 0) printf("west moved right\n");
            } else if (!left_turn && c->e_col == col_idx && r->west[LEFT]) { // up?
                move_nrth(&intrsctn_nxt[i], c, r);
                if(rank == 0) printf("west moved left\n");
            }
            else if (!left_turn && r->east[STRGHT]) { // try to go straight
                if(rank == 0) printf("west moved straight\n");
                move_east(&intrsctn_nxt[i], c, r);
                intrsctn_nxt[i].east->go_es[0] = c;
                west_to_east_rules(r);
            }
            else {
                intrsctn_nxt[i].west->go_es[block_end] = c;
            }

        }
        // printf("east: %d\nsouth: %d\nwest: %d\nleft turn: %d\n\n", east, south, west, left_turn);
        if(rank==0){
            printf("after\n");
            print_intersection(intrsctn_nxt[i]);
        }
        // printf("\n");
        left_turn = 0;
    }

}

unsigned long total_grid_dist_to_travel(unsigned long glbl_row_idx, street* sts, unsigned int n){
    unsigned long sum = 0;
    // unsigned long r, c;
    // unsigned long glbl_col_idx = !(glbl_row_idx%2);
    for(size_t i = 0; i < n; i++)
    {
        // r = glbl_row_idx + 2*(i/(SIDE_LENGTH-glbl_col_idx));
        // c = glbl_col_idx + 2*(i%(SIDE_LENGTH-glbl_col_idx));
        // printf("%lu %lu\n", r, c);
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            // sum += sts[i].go_es[j] ? dist_eloc(sts[i].go_es[j], c, r) : 0;
            // sum += sts[i].go_wn[j] ? dist_eloc(sts[i].go_wn[j], c, r) : 0;
            sum += sts[i].go_es[j] ? 1 : 0;
            sum += sts[i].go_wn[j] ? 1 : 0;
        }
        
    }
    return sum;
}

unsigned long total_grid_dist_to_travel_ghost(unsigned long glbl_row_idx, street* sts, unsigned int n, int n_or_s){
    // n_or_s 0 -> only count north 1-> only count south
    unsigned long sum = 0;
    // unsigned long r, c;
    // unsigned long glbl_col_idx = !(glbl_row_idx%2);
    for(size_t i = 0; i < n; i++)
    {
        // r = glbl_row_idx + 2*(i/SIDE_LENGTH);
        // c = glbl_col_idx + 2*(i%SIDE_LENGTH);
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            // if (n_or_s) {sum += sts[i].go_es[j] ? dist_eloc(sts[i].go_es[j], c, r) : 0;}
            // else {sum += sts[i].go_wn[j] ? dist_eloc(sts[i].go_wn[j], c, r) : 0;}
            if (n_or_s) {sum += sts[i].go_es[j] ? 1 : 0;}
            else {sum += sts[i].go_wn[j] ? 1 : 0;}

        }
        
    }
    return sum;
}