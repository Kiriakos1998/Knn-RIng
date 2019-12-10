#include "knnring.h"
// #include <mpi.h>
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include "time.h"
#include <stdio.h>
#include "functions.c"
#define changePointer(T, x, y) \
  do { \
    T *tmp = *x; \
    *x = *y; \
    *y = tmp; \
  } while (0);
knnresult distrAllkNN(double * X, int n, int d, int k)
{double comp_time=0;
  double comm_time=0;
// value hols the result that each process is going to return
    knnresult value;
int i; // just a counter

    int proccesses, currentRank, toSend, toReceive;
    MPI_Request requests[2];
    MPI_Status statuses[2];


manageMpiStuff(&proccesses,&currentRank,&toSend,&toReceive);
comp_time-=clock();
int fixingValue = currentRank==0?(proccesses-1)*n : (currentRank-1)*n;
// array responsible for proccessing and sending data to the othrer proccesses
comp_time+=clock();
    double *sender = (double *) malloc(d*n*sizeof(double));
//array responsible for receiving data while proccess is still in progress
    double *receiver = (double *) malloc(d*n*sizeof(double));
// send data to the next process
comm_time-=clock();
    MPI_Isend(X, d*n, MPI_DOUBLE, toSend, 1, MPI_COMM_WORLD, &requests[0]);
// receive from the previous
    MPI_Irecv(sender, d*n, MPI_DOUBLE, toReceive, 1, MPI_COMM_WORLD, &requests[1]);
comm_time+=clock();
comp_time-=clock();
    value = kNN(X, X, n, n, d, k);    //calculate the first result
manageIndexes(&value,fixingValue); // manage its indexes according to the fixingValue
comp_time+=clock();
comm_time-=clock();
    MPI_Waitall(2, requests, statuses); // Wait until all proccesses are done
comm_time+=clock();
// loop proccesess-1 times so every corpus has been compared with every query
    for(i=0; i<proccesses; ++i)
    {
// send data to the next procces
comm_time-=clock();
        MPI_Isend(sender, d*n, MPI_DOUBLE, toSend, 1, MPI_COMM_WORLD, &requests[0]);
// receive data from the previous procces
        MPI_Irecv(receiver, d*n, MPI_DOUBLE, toReceive, 1, MPI_COMM_WORLD, &requests[1]);
        comm_time+=clock();
// creare temporary result
knnresult temporary;
// calculate the fixingValue for the indexes
comp_time-=clock();
        fixingValue = calculateFixingValue(fixingValue,proccesses,n); //(idOffset - 1) modulo numtasks*n
// calculate the temporary result
        temporary = kNN(sender, X, n, n, d, k);
// settle the indexes
        manageIndexes(&temporary,fixingValue);
// compare the final result that is to be returned
        compareAndReplace(&value, &temporary);
// free the allocated memory for the ndist array and nidx array
        free(temporary.ndist);
        free(temporary.nidx);
        comp_time+=clock();
// make sure all handlers have finished threir tasks
comm_time-=clock();
for (int j=0;j<2;j++)
MPI_Wait(&requests[j], &statuses[j]);
// change the sender pointer with the receiver
        changePointer(double,&sender, &receiver);
        comm_time+=clock();
    }
// free allocated memory
if(currentRank== 0) {
    printf("TOTAL:%lf\n",(double) (comp_time + comm_time)/CLOCKS_PER_SEC);
    printf("Computations:%lf\n",(double)comp_time/CLOCKS_PER_SEC);
    printf("Communications:%lf\n",(double)comm_time/CLOCKS_PER_SEC);
}
    free(sender);
    free(receiver);
    return value;
}
