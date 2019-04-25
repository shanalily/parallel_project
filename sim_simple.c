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

#ifdef BGQ
#define SIZE 32768
#else
#define SIZE 4
#endif

#define NUM_TICKS 3
#define NPI 3
#define FROM_S 1
#define FROM_N 0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

/***************************************************************************/
/* Functions ***********************************************************/
/***************************************************************************/

void mk_cars(int n, float proportion, int* intrsctns, unsigned int g_i, int i_strt){
    int row;
    for(size_t i = i_strt; i < i_strt + n*NPI; i+=NPI)
    {
        row = (i-i_strt)/(NPI*SIZE) + g_i;
        intrsctns[i] = GenVal(row) < proportion;
        if(intrsctns[i]) {
            intrsctns[i+1] = (int) (GenVal(row)*SIZE);
            intrsctns[i+2] = (int) (GenVal(row)*SIZE);;
        }
    }
    
}

int reachdest(int num, int* is, uint g_i, int i_strt) {
    int r, c;
    int sum = 0;
    for(size_t i =i_strt; i < i_strt + num*NPI; i+= NPI) {
        r = (i-i_strt)/(NPI*SIZE) + g_i;
        c = ((i-i_strt)%(NPI*SIZE))/NPI;
        if(is[i] != 0 && is[i+1] == r && is[i+2] == c) {
            sum+=1;
            printf("%d %d has Arrived!\n", r, c);
            is[i] = 0;
        }
    }
    return sum;
}

void exchange_rows(int mpi_myrank, int mpi_commsize, int* is, int rpr){
    MPI_Request north;
    MPI_Request south;
    MPI_Request request_n;
    MPI_Request request_s;

    int* is_south_recv = &is[(rpr+(mpi_myrank!=0))*SIZE*NPI];
    int* is_south_send = &is[(rpr-1+(mpi_myrank!=0))*SIZE*NPI];
    int* is_north_recv = &is[0];
    int* is_north_send = &is[SIZE*NPI];


    if (mpi_myrank != 0) {
        MPI_Irecv(is_north_recv, SIZE*NPI, MPI_INT, mpi_myrank-1, 0, MPI_COMM_WORLD, &north);
    }
    if (mpi_myrank != mpi_commsize - 1) {
        MPI_Irecv(is_south_recv, SIZE*NPI, MPI_INT, mpi_myrank+1, 1, MPI_COMM_WORLD, &south);
    }
    if (mpi_myrank != 0) {
        MPI_Isend(is_north_send, SIZE*NPI, MPI_INT, mpi_myrank-1, 1, MPI_COMM_WORLD, &request_n);
    }
    if (mpi_myrank != mpi_commsize - 1) {
        MPI_Isend(is_south_send, SIZE*NPI, MPI_INT, mpi_myrank+1, 0, MPI_COMM_WORLD, &request_s);
    }

    // printf("a\n");
    if(mpi_myrank != 0) MPI_Wait(&request_n, MPI_STATUS_IGNORE);
    if(mpi_myrank != mpi_commsize -1) MPI_Wait(&request_s, MPI_STATUS_IGNORE);
    // printf("b\n");
    if(mpi_myrank != 0) MPI_Wait(&north, MPI_STATUS_IGNORE);
    // printf("c\n");

    if(mpi_myrank != mpi_commsize -1) MPI_Wait(&south, MPI_STATUS_IGNORE);
    // printf("d\n");

}

void update_i(int i, int n, int s, int w, int e, int* i_now, int* i_nxt){
    // indexes
    
    int n_i = i-NPI*SIZE;
    int s_i = i+NPI*SIZE;
    int w_i = i-NPI;
    int e_i = i+NPI;

    // Poll rank below it to see if move is valid

    // to go must be nothing there rn or in the future
    if(n && i_now[n_i] == 0 && i_nxt[n_i] == 0) {
        i_nxt[n_i] = i_now[i];
        i_nxt[n_i+1] = i_now[i+1];
        i_nxt[n_i+2] = i_now[i+2];
        
    }
        else if(s && i_now[s_i] == 0 && i_nxt[s_i] == 0) {
        i_nxt[s_i] = i_now[i];
        i_nxt[s_i+1] = i_now[i+1];
        i_nxt[s_i+2] = i_now[i+2];
    }
    else if(w && i_now[w_i] == 0 && i_nxt[w_i] == 0) {
        i_nxt[w_i] = i_now[i];
        i_nxt[w_i+1] = i_now[i+1];
        i_nxt[w_i+2] = i_now[i+2];
    }
    else if(e && i_now[e_i] == 0 && i_nxt[e_i] == 0) {
        i_nxt[e_i] = i_now[i];
        i_nxt[e_i+1] = i_now[i+1];
        i_nxt[e_i+2] = i_now[i+2];
    }
    else {
        i_nxt[i] = i_now[i];
        i_nxt[i+1] = i_now[i+1];
        i_nxt[i+2] = i_now[i+2];
    }
}

void answer_rqsts(MPI_Request* rqsts, int i_strt, int* i_now, int* i_nxt, int rank, int* answered){
    MPI_Request request;
    int idx;
    for(size_t j = 0; j < SIZE; j++)
    {
        if(!answered[j]) {
            int flag;
            MPI_Test(&rqsts[j], &flag, MPI_STATUS_IGNORE);
            idx = i_strt+j*NPI;
            if(flag) {
                // printf("Answering %d\n", j);
                // printf("Sending to receive on tag %d\n", SIZE+j);
                answered[j] = 1;
                if(i_now[idx] == 0 && i_nxt[idx] == 0){
                    i_nxt[idx] = 1;
                    MPI_Isend(&i_nxt[idx], 1, MPI_INT, rank, SIZE+j, MPI_COMM_WORLD, &request);
                }
                else {
                    int tmp = 0;
                    // printf("Sending to receive on tag %d\n", j);
                    MPI_Isend(&tmp, 1, MPI_INT, rank, SIZE+j, MPI_COMM_WORLD, &request);
                }
            }
        }
    }
}

int recv_rqsts(MPI_Request* rqsts, int i_strt, int* i_now, int* i_nxt, int* answers, int n_s, int g_i, int* received){
    int n_recv = 0;
    int idx;
    int r, c, n, s, w, e;
    for(size_t j = 0; j < SIZE; j++)
    {
        if(!received[j]){
            int flag;
            MPI_Test(&rqsts[j], &flag, MPI_STATUS_IGNORE);
            if(flag) {
                // printf("Receiving %d\n", j);
                received[j] = 1;
                n_recv+=1;
                idx = i_strt+j*NPI;
                r = (idx-i_strt)/(NPI*SIZE) + g_i;
                c = ((idx-i_strt)%(NPI*SIZE))/NPI;
                n = r!=0 && i_now[idx+1] < r && ((n_s == 0 && answers[j]) || n_s == 1);
                s = r!=SIZE-1 && i_now[idx+1] > r && ((n_s == 1 && answers[j]) || n_s == 0);
                w = c!=0 && i_now[idx+2] < c;
                e = c!=SIZE-1 && i_now[idx+2] > c;
                update_i(idx, n, s, w, e, i_now, i_nxt);
            }
        }
    }
    return n_recv;
}


void update_intersections(int num, int* i_now, int* i_nxt, unsigned int g_i, int i_strt, int rpr) {
    int r, c;
    int n, s, w, e;
    // MPI receives
    MPI_Request recv_rqsts_n[SIZE];
    MPI_Request recv_rqsts_s[SIZE];
    MPI_Request send_rqsts_n[SIZE];
    MPI_Request send_rqsts_s[SIZE];
    MPI_Request done_n;
    int done_i_n;
    MPI_Request done_s;
    int done_i_s;
    MPI_Request request;
    int* answers_n = calloc(SIZE, sizeof(int)); 
    int* answers_s = calloc(SIZE, sizeof(int));
    int* recv_n = calloc(SIZE, sizeof(int)); 
    int* recv_s = calloc(SIZE, sizeof(int)); 
    int* answered_n = calloc(SIZE, sizeof(int)); 
    int* answered_s = calloc(SIZE, sizeof(int));
    int* received_n = calloc(SIZE, sizeof(int)); 
    int* received_s = calloc(SIZE, sizeof(int));
    int mpi_myrank;
    int mpi_commsize;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    int num_recv = 0;


    if(mpi_myrank != 0){
        // printf("Rank %d Receiving on tag %d\n", mpi_myrank, 2*SIZE+1);
        MPI_Irecv(&done_i_n, 1, MPI_INT, mpi_myrank-1, 2*SIZE+1, MPI_COMM_WORLD, &done_n);
        for(size_t i = 0; i < SIZE; i++)
        {
            MPI_Irecv(&recv_n[i], 1, MPI_INT, mpi_myrank-1, i, MPI_COMM_WORLD, &send_rqsts_n[i]);
            MPI_Irecv(&answers_n[i], 1, MPI_INT, mpi_myrank-1, SIZE+i, MPI_COMM_WORLD, &recv_rqsts_n[i]);
        }
    }
    if(mpi_myrank != mpi_commsize-1){
        // printf("Rank %d Receiving on tag %d\n", mpi_myrank, 2*SIZE+2);
        MPI_Irecv(&done_i_s, 1, MPI_INT, mpi_myrank+1, 2*SIZE+2, MPI_COMM_WORLD, &done_s);
        for(size_t i = 0; i < SIZE; i++)
        {
            MPI_Irecv(&recv_s[i], 1, MPI_INT, mpi_myrank+1, i, MPI_COMM_WORLD, &send_rqsts_s[i]);
            MPI_Irecv(&answers_s[i], 1, MPI_INT, mpi_myrank+1, SIZE+i, MPI_COMM_WORLD, &recv_rqsts_s[i]);
        }
    }


    

    for(size_t i =i_strt; i < i_strt + num*NPI; i+= NPI) {

        // answer requests FOR BOTH TOP and BOTTM
        //   IF SPOT is open tell them and mark it used

        if(mpi_myrank != 0){
            answer_rqsts(send_rqsts_n, i_strt, i_now, i_nxt, mpi_myrank-1, answered_n);
        }
        if(mpi_myrank != mpi_commsize-1){
            answer_rqsts(send_rqsts_s, i_strt+rpr*SIZE*NPI, i_now, i_nxt, mpi_myrank+1, answered_s);
        }
        // check for answered requests
        //   IF SPOT is open use it like normal
        //   IF SPOT is closed set n/s as not available and make normal update
        if(mpi_myrank != 0){
            num_recv -= recv_rqsts(recv_rqsts_n, i_strt, i_now, i_nxt, answers_n, 0, g_i, received_n);
        }
        if(mpi_myrank != mpi_commsize-1){
            num_recv -= recv_rqsts(recv_rqsts_s, i_strt+rpr*SIZE*NPI, i_now, i_nxt, answers_s, 0, g_i+rpr*SIZE*NPI, received_s);
        }

        r = (i-i_strt)/(NPI*SIZE) + g_i;
        c = ((i-i_strt)%(NPI*SIZE))/NPI;
        if(i_now[i] != 0) {
            // wants to go
            n = r!=0 && i_now[i+1] < r;
            s = r!=SIZE-1 && i_now[i+1] > r;
            w = c!=0 && i_now[i+2] < c;
            e = c!=SIZE-1 && i_now[i+2] > c;

            if(n && (r - g_i) == 0) {
                // save north for interprocess
                int tmp = (int) i;
                // printf("Sending to answer on tag %d\n", c);
                MPI_Isend(&tmp, 1, MPI_INT, mpi_myrank-1, c, MPI_COMM_WORLD, &request);
                num_recv += 1;
            }
            else if (s && (r - g_i) == (rpr - 1)){
                // save south for interprocess
                int tmp = (int) i;
                // printf("Sending to answer on tag %d\n", c);
                MPI_Isend(&tmp, 1, MPI_INT, mpi_myrank+1, c, MPI_COMM_WORLD, &request);
                num_recv += 1;
            }

            update_i(i, n, s, w, e, i_now, i_nxt);

            
        }


    }
    // printf("Rank %d: num_recv %d\n", mpi_myrank, num_recv);

    int flag_n = 0;
    int flag_s = 0;
    int done_recv = 0;
    while(!flag_n || !flag_s || !done_recv) {
        if(mpi_myrank != 0){
            answer_rqsts(send_rqsts_n, i_strt, i_now, i_nxt, mpi_myrank-1, answered_n);
        }
        if(mpi_myrank != mpi_commsize-1){
            answer_rqsts(send_rqsts_s, i_strt+rpr*SIZE*NPI, i_now, i_nxt, mpi_myrank+1, answered_s);
        }
        // check for answered requests
        //   IF SPOT is open use it like normal
        //   IF SPOT is closed set n/s as not available and make normal update
        if(!done_recv && mpi_myrank != 0){
            num_recv -= recv_rqsts(recv_rqsts_n, i_strt, i_now, i_nxt, answers_n, 0, g_i, received_n);
        }
        if(!done_recv && mpi_myrank != mpi_commsize-1){
            num_recv -= recv_rqsts(recv_rqsts_s, i_strt+rpr*SIZE*NPI, i_now, i_nxt, answers_s, 0, g_i+rpr*SIZE*NPI, received_s);
        }

        if(num_recv == 0 && done_recv == 0){
            done_recv = 1;
            if(mpi_myrank != 0){
                // printf("Rank %d sending on tag %d\n", mpi_myrank, 2*SIZE+2);
                MPI_Isend(&done_i_s, 1, MPI_INT, mpi_myrank-1, 2*SIZE+2, MPI_COMM_WORLD, &request);
            }
            if(mpi_myrank != mpi_commsize-1){
                // printf("Rank %d sending on tag %d\n", mpi_myrank, 2*SIZE+1);
                MPI_Isend(&done_i_n, 1, MPI_INT, mpi_myrank+1, 2*SIZE+1, MPI_COMM_WORLD, &request);
            }
            for(size_t i = 0; i < SIZE; i++)
            {
                if(mpi_myrank != 0){
                    if(received_n[i] == 0) {
                        // printf("%d\n", i);
                        // MPI_Request_free(&recv_rqsts_n[i]);
                        int junk = 0;
                        MPI_Isend(&junk, 1, MPI_INT, mpi_myrank-1, i, MPI_COMM_WORLD, &request);
                    }
                }
                if(mpi_myrank != mpi_commsize-1){
                    if(received_s[i] == 0) {
                        // printf("%d\n", i);
                        // MPI_Request_free(&recv_rqsts_s[i]);
                        int junk = 0;
                        MPI_Isend(&junk, 1, MPI_INT, mpi_myrank+1, i, MPI_COMM_WORLD, &request);
                    }
                }
            }
        }
        if(mpi_myrank != 0 && !flag_n){
            MPI_Test(&done_n, &flag_n, MPI_STATUS_IGNORE);
            flag_s = 1;
            // printf("Rank %d testing done %d\n", mpi_myrank, flag_n);
        }
        if(mpi_myrank != mpi_commsize-1 && !flag_s){
            MPI_Test(&done_s, &flag_s, MPI_STATUS_IGNORE);
            flag_n = 1;
            // printf("Rank %d testing done %d\n", mpi_myrank, flag_s);

        }
    }
    // printf("Rank %d Free from loop\n", mpi_myrank);

    for(size_t i = 0; i < SIZE; i++)
    {
        if(mpi_myrank != 0){
            if(answered_n[i] == 0) {
                // MPI_Request_free(&send_rqsts_n[i]);
                int junk = 0;
                MPI_Isend(&junk, 1, MPI_INT, mpi_myrank-1, i+SIZE, MPI_COMM_WORLD, &request);
            }
        }
        if(mpi_myrank != mpi_commsize-1){
            if(answered_s[i] == 0) {
                // MPI_Request_free(&send_rqsts_s[i]);
                int junk = 0;
                MPI_Isend(&junk, 1, MPI_INT, mpi_myrank+1, i+SIZE, MPI_COMM_WORLD, &request);
            }
        }
    }

    
    free(answers_n); 
    free(answers_s);
    free(recv_n); 
    free(recv_s); 
    free(answered_n); 
    free(answered_s);
    free(received_n); 
    free(received_s);
}

void prnt_ints(unsigned int rpr, unsigned int glbl_index, int* intrsctns, int g_b) {
    for(size_t i = 0; i < rpr; i++)
    {
        printf("Row %lu: %d %d %d %d %d %d %d %d %d %d %d %d\n", glbl_index+i-g_b, 
            intrsctns[i*SIZE*NPI+0], intrsctns[i*SIZE*NPI+1], intrsctns[i*SIZE*NPI+2],
            intrsctns[i*SIZE*NPI+3], intrsctns[i*SIZE*NPI+4], intrsctns[i*SIZE*NPI+5],
            intrsctns[i*SIZE*NPI+6], intrsctns[i*SIZE*NPI+7], intrsctns[i*SIZE*NPI+8],
            intrsctns[i*SIZE*NPI+9], intrsctns[i*SIZE*NPI+10], intrsctns[i*SIZE*NPI+11]);
    }
}

void clear_is(int n, int* is){
    for(size_t i =0; i < n*NPI; i+= 1) {
        is[i] = 0;
    }
}

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
    unsigned int rpr = SIZE/mpi_commsize;
    float proportion = 0.5;
    unsigned int glbl_index = rpr*mpi_myrank;
    int i_strt = (mpi_myrank != 0)*NPI*SIZE;
    InitDefault();

    int ghosts = (mpi_myrank != 0) + (mpi_myrank != mpi_commsize - 1);
    // allocate one array with interstions
    // each intersection gets 3 values. Empty/Not Empty. Final Dest m, n
    int* intrsctns_now = calloc(NPI*(rpr + ghosts)*SIZE, sizeof(int));
    int* intrsctns_nxt = calloc(NPI*(rpr + ghosts)*SIZE, sizeof(int));

    mk_cars(rpr*SIZE, proportion, intrsctns_now, glbl_index, i_strt);

    // prnt_ints(rpr + ghosts, glbl_index, intrsctns_now, mpi_myrank != 0);
    

    MPI_Barrier( MPI_COMM_WORLD );
    if (mpi_myrank == 0) {
        g_start_cycles = GetTimeBase();
    }
    MPI_Barrier( MPI_COMM_WORLD );


    for (size_t i = 1; i < NUM_TICKS; ++i) {

        // reach destination
        reachdest(rpr*SIZE, intrsctns_now, glbl_index, i_strt);
        // do the exchange 
        exchange_rows(mpi_myrank, mpi_commsize, intrsctns_now, rpr);
        update_intersections(rpr*SIZE, intrsctns_now, intrsctns_nxt, glbl_index, i_strt ,rpr);

        // clear rows
        int* tmp = intrsctns_now;
        intrsctns_now = intrsctns_nxt;
        intrsctns_nxt = tmp;

        MPI_Barrier( MPI_COMM_WORLD );

        clear_is(SIZE*(rpr + ghosts), intrsctns_nxt);
        if(i==2 || i==1) {
            prnt_ints(rpr + ghosts, glbl_index, intrsctns_now, mpi_myrank != 0);
        }

    }


    MPI_Barrier( MPI_COMM_WORLD );
    if (mpi_myrank == 0) {
        g_end_cycles = GetTimeBase();
        g_time_in_secs = ((double)(g_end_cycles - g_start_cycles))/g_processor_frequency;
        printf("Sim time: %f\n", g_time_in_secs);
    }
    MPI_Barrier( MPI_COMM_WORLD );

    free( intrsctns_now );
    free( intrsctns_nxt );

    MPI_Finalize();
    return 0;
}