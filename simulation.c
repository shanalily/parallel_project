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

void update_streets(unsigned int n, street *streets_now, street *streets_nxt);

void update_ghost_streets(unsigned int n, street* ghost_now, street* ghost_nxt, int n_or_s);

int move_nrth(intrsctn *is, car *c, intrsctn_rules *r);
int move_soth(intrsctn *is, car *c, intrsctn_rules *r);
int move_west(intrsctn *is, car *c, intrsctn_rules *r);
int move_east(intrsctn *is, car *c, intrsctn_rules *r);
void move_from_wn(street *nxt_str);
void move_from_es(street *nxt_str);
void nrth_to_soth_rules(intrsctn_rules *r);
void east_to_west_rules(intrsctn_rules *r);
void soth_to_nrth_rules(intrsctn_rules *r);
void west_to_east_rules(intrsctn_rules *r);
int reached_dest(unsigned long glbl_row_idx, unsigned long glbl_col_idx, unsigned short idx, car *c, unsigned short wn);

void update_intersections(unsigned int rpr, unsigned long glbl_index, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt, intrsctn_rules *r);

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
    unsigned long long capacity = 2*(SIDE_LENGTH-1)*(SIDE_LENGTH)*ROAD_CAP;
    unsigned long long n_cars = (unsigned long long) (((long double) proportion)*capacity);
    unsigned int glbl_index = 2*rpr*mpi_myrank;
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
#endif
    
    MPI_Barrier( MPI_COMM_WORLD );

    unsigned long dist_left_ew = total_grid_dist_to_travel(glbl_index, streets_ew_now, (SIDE_LENGTH-1)*rpr);
    unsigned long dist_left_ns = total_grid_dist_to_travel(glbl_index+1, streets_ns_now, SIDE_LENGTH*(rpr-1));
    unsigned long dist_left_gn = mpi_myrank != 0 ? total_grid_dist_to_travel_ghost(glbl_index-1, ghost_ns_nrth_now, SIDE_LENGTH, 0) : 0;
    unsigned long dist_left_gs = mpi_myrank != mpi_commsize ? total_grid_dist_to_travel_ghost(glbl_index+2*rpr-1, ghost_ns_soth_now, SIDE_LENGTH, 1) : 0;
    unsigned long dist_left = dist_left_ew+dist_left_ns+dist_left_gn+dist_left_gs;
    printf("Rank %d: Tick %d: Distance left is %lu\n", mpi_myrank, 0, dist_left);

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
    // update east/west
    update_streets(rpr*(SIDE_LENGTH-1), streets_ew_now, streets_ew_nxt);
    // update north/south
    update_streets((rpr-1)*SIDE_LENGTH, streets_ns_now, streets_ns_nxt);
    // update ghost rows
    update_ghost_streets(SIDE_LENGTH, ghost_ns_nrth_now, ghost_ns_nrth_nxt, 0);
    update_ghost_streets(SIDE_LENGTH, ghost_ns_soth_now, ghost_ns_soth_nxt, 1);
    // run intersections. how do I make sure this gets updated ghost row streets?
#ifdef DEBUG
    int car_count = check_car_count(rpr, intrsctn_now);
#endif
    update_intersections(rpr, glbl_index, intrsctn_now, intrsctn_nxt, &r);
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
                c->e_row = (unsigned long) (GenVal(row)*(2*SIDE_LENGTH-1));
                c->e_idx = (int) (GenVal(row)*(ROAD_CAP));

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

    if(mpi_myrank != 0) MPI_Wait(&nrth, MPI_STATUS_IGNORE);
    if(mpi_myrank != mpi_commsize -1) MPI_Wait(&soth, MPI_STATUS_IGNORE);

    MPI_Type_free(&mpi_car_type);
}

// Move everything down the streets. Slot ROAD_CAP-1 is the location on the street closest
// to the intersection at the end of the block.
// 0           ----> ROAD_CAP-1
// ROAD_CAP-1 <----  0
void update_streets(unsigned int n, street *streets_now, street *streets_nxt) {
    for (size_t i = 0; i < n; ++i) {
        // if location on street is empty, move up the next car (if it exists) from previous location
        // set previous location to empty
        for (size_t j = ROAD_CAP-1; j >= 1; --j) {
            if (streets_now[i].go_es[j] == EMPTY) {
                streets_nxt[i].go_es[j] = streets_now[i].go_es[j-1];
                streets_now[i].go_es[j-1] = EMPTY;
            }
            if (streets_now[i].go_wn[j] == EMPTY) {
                streets_nxt[i].go_wn[j] = streets_now[i].go_wn[j-1];
                streets_now[i].go_wn[j-1] = EMPTY;
            }
        }
        // last slot of now will have it's previous value still.
    }
}

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

// Set first location on northern street to the car being moved
int move_nrth(intrsctn *is, car *c, intrsctn_rules *r) {
    is->nrth->go_wn[0] = c;
    r->west[2] = 0;
    r->east[0] = 0;
    return 1;
}

// Set first location on southern street to the car being moved
int move_soth(intrsctn *is, car *c, intrsctn_rules *r) {
    is->soth->go_es[0] = c;
    r->west[0] = 0;
    r->east[2] = 0;
    return 1;
}

// Set first location on western street to the car being moved
int move_west(intrsctn *is, car *c, intrsctn_rules *r) {
    is->west->go_wn[0] = c;
    r->nrth[0] = 0;
    r->soth[2] = 0;
    return 1;
}

// Set first location on eastern street to the car being moved
int move_east(intrsctn *is, car *c, intrsctn_rules *r) {
    is->east->go_es[0] = c;
    r->nrth[2] = 0;
    r->soth[0] = 0;
    return 1;
}

// Set last location on street, which is at intersection on street going west or north, to empty
void move_from_wn(street *nxt_str) {
    nxt_str->go_wn[ROAD_CAP-1] = EMPTY;
}

// Set last location on street, which is at intersection on street going east or south, to empty
void move_from_es(street *nxt_str) {
    nxt_str->go_es[ROAD_CAP-1] = EMPTY;
}

// Set and reset intersection rules.
void reset_intrsctn(intrsctn intrsctn_nxt, intrsctn_rules *r) {
    // intrsctn_nxt->nrth = EMPTY;
    // intrsctn_nxt->soth = EMPTY;
    // intrsctn_nxt->west = EMPTY;
    // intrsctn_nxt->east = EMPTY;
    for (size_t i = 0; i < TURN_OPTS; ++i) {
        r->nrth[i] = 1;
        r->soth[i] = 1;
        r->west[i] = 1;
        r->east[i] = 1;
    }
}

void nrth_to_soth_rules(intrsctn_rules *r) {
    // east can't go straight
    r->east[STRGHT] = 0;
    // south can't go left
    r->soth[LEFT] = 0;
    // west can't go straight or left
    r->west[STRGHT] = 0;
    r->west[LEFT] = 0;
}

void east_to_west_rules(intrsctn_rules *r) {
    // south can't go straight
    r->soth[STRGHT] = 0;
    // west can't go left
    r->west[LEFT] = 0;
    // north can't go straight or left
    r->nrth[STRGHT] = 0;
    r->nrth[LEFT] = 0;
}

void soth_to_nrth_rules(intrsctn_rules *r) {
    // west can't go straight
    r->west[STRGHT] = 0;
    // north can't go left
    r->nrth[LEFT] = 0;
    // east can't go straight or left
    r->east[STRGHT] = 0;
    r->east[LEFT] = 0;
}

void west_to_east_rules(intrsctn_rules *r) {
    // north can't go straight
    r->nrth[STRGHT] = 0;
    // east can't go left
    r->east[LEFT] = 0;
    // south can't go straight or left
    r->soth[STRGHT] = 0;
    r->soth[LEFT] = 0;
}

// For debugging. This might not be helpful, idk. Will delete later.
void print_intrsctn(intrsctn intersection) {
    printf("-------\n");
    if (intersection.nrth) printf("  n   \n");
    else printf("      \n");

    if (intersection.west && intersection.east) printf("w    e\n");
    else if (intersection.west) printf("w     \n");
    else if (intersection.east) printf("     e\n");
    else printf("      \n");

    if (intersection.soth) printf("  s   \n");
    else printf("      \n");
}

// If car has reached it's end point, return 1. Otherwise, return 0.
// If row is even, EW/WE. If row is odd, SN/NS. It doesn't matter which side of street the car is on, but index needs to be adjusted accordingly.
// input should be the even glbl_row_idx and even glbl_col_idx from update_intersections. Assume e_idx if is pointing west or north.
int reached_dest(unsigned long glbl_row_idx, unsigned long glbl_col_idx, unsigned short idx, car *c, unsigned short wn) {
    if ((glbl_row_idx == c->e_row || glbl_row_idx == c->e_row+1) && (glbl_col_idx == c->e_col)
            && ((wn && idx == c->e_idx) || (!wn && idx == ROAD_CAP-1-c->e_idx))) {
        return 1;
    }
    return 0;
}

// make sure to switch which intersections are passed between now and nxt.
// EW streets have even rows numbers and odd columns numbers
// NS streets have odd rows numbers and even column numbers
void update_intersections(unsigned int rpr, unsigned long glbl_row_idx, intrsctn *intrsctn_now, intrsctn *intrsctn_nxt, intrsctn_rules *r) {
    // update nxt based on now. Determine who can move first, corresponding to traffic rules.
    unsigned int block_end = ROAD_CAP-1;
    int east = 0, west = 0, south = 0, left_turn = 0;
    for (size_t i = 0; i < rpr*SIDE_LENGTH; ++i) {
        unsigned long row_idx = glbl_row_idx + 2*(i/SIDE_LENGTH)/*+ 2*(i/(SIDE_LENGTH-glbl_col_idx))*/;
        unsigned long col_idx = 2*(i%SIDE_LENGTH);
        reset_intrsctn(intrsctn_nxt[i], r);
        // Check where car is heading.
        // If coming from north, should not be needing to go top-left or top-right.
        // prioritize going straight, then right, then left, depending on direction of destination
        // Check if car can come from north, and if a car is coming from north.
        if (intrsctn_now[i].nrth && intrsctn_now[i].nrth->go_wn[block_end]) {

            car *c = intrsctn_now[i].nrth->go_wn[block_end];
            // printf("north!!!!!!!!!!!!!!!\n");
            // printf("%lu %lu %d\n", c->e_row, c->e_col, c->e_idx);
            // coming from north, should they turn straight, right, or left?
            // check | c->e_row == row_idx+1? row_idx+1 is NS
            if (c->e_row == row_idx || c->e_row == row_idx+1) {
                // left or right? or arrived? arrived should be checked at end of this function though
                if (c->e_col <= col_idx && intrsctn_nxt[i].west) // go west/right?.   Is my problem that streets aren't correctly set for nxt?
                    west = move_west(&intrsctn_nxt[i], c, r);
                else if (intrsctn_nxt[i].east) { // go east/left?
                    east = move_east(&intrsctn_nxt[i], c, r); // should I be yielding the right of way?
                    // basically no one can go anywhere... I should probably show this another way.
                    // east can't go west, north, or south
                    left_turn = 1;
                }

            } else if (intrsctn_nxt[i].soth) {
                south = move_soth(&intrsctn_nxt[i], c, r); // go straight?
                nrth_to_soth_rules(r);
            }
            move_from_wn(intrsctn_nxt[i].nrth);
            // if (reached_dest(row_idx, col_idx, ROAD_CAP-1, c, 0)) // how should a DES deal with this?
            //     printf("reached destination! %lu %lu %d\n", row_idx, col_idx, ROAD_CAP-1);

        }
        if (intrsctn_now[i].east && intrsctn_now[i].east->go_es[block_end] && !left_turn) {
            car *c = intrsctn_now[i].east->go_es[block_end];
            // printf("east!!!!!!!!!!!!!!!\n");
            // printf("%lu %lu %d\n", c->e_row, c->e_col, c->e_idx);
            if (c->e_col == col_idx) {
                // north or south?
                if (intrsctn_nxt[i].nrth && c->e_row <= row_idx && r->east[RIGHT]) // north/right?
                    move_nrth(&intrsctn_nxt[i], c, r);
                else if (intrsctn_nxt[i].soth && r->east[LEFT]) {
                    south = move_soth(&intrsctn_nxt[i], c, r);
                    left_turn = 1;
                }

            } else if (intrsctn_nxt[i].west && r->east[STRGHT]) {
                west = move_west(&intrsctn_nxt[i], c, r);
                east_to_west_rules(r);
            }
            
            if (!east) move_from_es(intrsctn_nxt[i].east);
        }
        if (intrsctn_now[i].soth && intrsctn_now[i].soth->go_es[block_end] && !left_turn) {

            car *c = intrsctn_now[i].soth->go_es[block_end];
            // printf("south!!!!!!!!!!!!!!!\n");
            // printf("%lu %lu %d\n", c->e_row, c->e_col, c->e_idx);
            if (c->e_row == row_idx || c->e_row == row_idx+1) {
                // east or west?
                if (intrsctn_nxt[i].east && c->e_col >= col_idx && r->soth[RIGHT]) // go east/right?
                    east = move_east(&intrsctn_nxt[i], c, r);
                else if (intrsctn_nxt[i].west && r->soth[LEFT]) {
                    west = move_west(&intrsctn_nxt[i], c, r);
                    left_turn = 1;
                }

            } else if (intrsctn_nxt[i].nrth && r->soth[STRGHT]) {
                move_nrth(&intrsctn_nxt[i], c, r);
                soth_to_nrth_rules(r);
            }
            
            if (!south) move_from_es(intrsctn_nxt[i].soth);

        }
        if (intrsctn_now[i].west && intrsctn_now[i].west->go_wn[block_end] && !left_turn) { // coming from west
            // I can't go
            car *c = intrsctn_now[i].west->go_wn[block_end];
            // printf("west!!!!!!!!!!!!!!!\n");
            // printf("%lu %lu %d\n", c->e_row, c->e_col, c->e_idx);
            if (c->e_col == col_idx) {
                // up or down?
                if (intrsctn_nxt[i].soth && c->e_row >= row_idx && r->west[RIGHT]) // south?
                    south = move_soth(&intrsctn_nxt[i], c, r);
                else if (intrsctn_nxt[i].nrth && r->west[LEFT]) // up?
                    move_nrth(&intrsctn_nxt[i], c, r);

            } else if (intrsctn_nxt[i].east && r->east[STRGHT]) { // try to go straight
                east = move_east(&intrsctn_nxt[i], c, r);
                west_to_east_rules(r);
            }

            if (!west) move_from_wn(intrsctn_nxt[i].west); // ? if another car wasn't already just moved there

        }
        // printf("east: %d\nsouth: %d\nwest: %d\nleft turn: %d\n", east, south, west, left_turn);
        // printf("before\n");
        // print_intrsctn(intrsctn_now[i]);
        // printf("after\n");
        // print_intrsctn(intrsctn_nxt[i]);
    }

}

unsigned long total_grid_dist_to_travel(unsigned long glbl_row_idx, street* sts, unsigned int n){
    unsigned long sum = 0;
    unsigned long r, c;
    unsigned long glbl_col_idx = !(glbl_row_idx%2);
    for(size_t i = 0; i < n; i++)
    {
        r = glbl_row_idx + 2*(i/(SIDE_LENGTH-glbl_col_idx));
        c = glbl_col_idx + 2*(i%(SIDE_LENGTH-glbl_col_idx));
        // printf("%lu %lu\n", r, c);
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            sum += sts[i].go_es[j] ? dist_eloc(sts[i].go_es[j], c, r) : 0;
            sum += sts[i].go_wn[j] ? dist_eloc(sts[i].go_wn[j], c, r) : 0;
            // sum += sts[i].go_es[j] ? 1 : 0;
            // sum += sts[i].go_wn[j] ? 1 : 0;
        }
        
    }
    return sum;
}

unsigned long total_grid_dist_to_travel_ghost(unsigned long glbl_row_idx, street* sts, unsigned int n, int n_or_s){
    // n_or_s 0 -> only count north 1-> only count south
    unsigned long sum = 0;
    unsigned long r, c;
    unsigned long glbl_col_idx = !(glbl_row_idx%2);
    for(size_t i = 0; i < n; i++)
    {
        r = glbl_row_idx + 2*(i/SIDE_LENGTH);
        c = glbl_col_idx + 2*(i%SIDE_LENGTH);
        for(size_t j = 0; j < ROAD_CAP; j++)
        {
            if (n_or_s) {sum += sts[i].go_es[j] ? dist_eloc(sts[i].go_es[j], c, r) : 0;}
            else {sum += sts[i].go_wn[j] ? dist_eloc(sts[i].go_wn[j], c, r) : 0;}
            // if (n_or_s) {sum += sts[i].go_es[j] ? 1 : 0;}
            // else {sum += sts[i].go_wn[j] ? 1 : 0;}

        }
        
    }
    return sum;
}