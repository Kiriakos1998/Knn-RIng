#include "knnring.h"
#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "functions.c"
#include "time.h"
#include <stdio.h>
// define function to swap two pointers of any type
#define changePointer(T, x, y) \
  do { \
    T *tmp = *x; \
    *x = *y; \
    *y = tmp; \
  } while (0);











knnresult distrAllkNN(double * X, int n, int d, int k)
{double comp_time=0;
  double comm_time=0;
    knnresult value;// holds the final results and it is updated through the loop
    // tha simulates the ring proccess
  // declare variable such as the number of proccess , the rank of each proccess
  // the rank of the proccess that has to send the query it has(toSend)
  // and the rank of proccess that expects to receive data from(toReceive)
    int procceses, currentRank, toSend, toReceive;
manageMpiStuff(&procceses,&currentRank,&toSend,&toReceive);// set the
// correct values to these variables
comp_time-=clock();

     int  fixingValue = (currentRank-1)*n;// caltulate the fixingValue
  // check if the currentRank is 0
      if(currentRank == 0) fixingValue = (procceses-1)*n;// special case when
      // rank is zero





    value = kNN(X, X, n, n, d, k); // calculate the first result
manageIndexes(&value,fixingValue); // manage its indexes)
comp_time+=clock();
double *rec=NULL; // array that will receive data in the odd proccess because
// these procceses will first receive and then send data
    double *sender = (double *) malloc(d * n * sizeof(double));
    // array that will be used for proccesing data and sending

 rec= (double* )malloc (d*n*sizeof(double));// allocate memory for rec
comm_time-=clock();
if((currentRank%2)==0) // check if process is even
 {
// first send data and then receive
   MPI_Send(X,d*n,MPI_DOUBLE,toSend,1,MPI_COMM_WORLD);
 MPI_Recv(sender,d*n,MPI_DOUBLE,toReceive,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
 }
 else {
   // if it is odd first receive data
   MPI_Recv(rec,d*n,MPI_DOUBLE,toReceive,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
// then send
   MPI_Send(X,d*n,MPI_DOUBLE,toSend,1,MPI_COMM_WORLD);
// and change pointers of sender and rec in order sender to point in new data
changePointer(double,&sender,&rec);
}
comm_time+=clock();

    for(int i=1;i<procceses; i++)
    {
comp_time-=clock();
// calculate the fixingValue in order to manage the indexes
fixingValue = calculateFixingValue(fixingValue,procceses,n);
        //calculate the knnresult for this query
        knnresult newResult = kNN(sender, X, n, n, d, k);
        // fix the indexing of the distances
        manageIndexes(&newResult,fixingValue);
        // compare the two result so only the k smallest distances are kept
        compareAndReplace(&value, &newResult);
// send and receive query . Send to the next proccess. receive from the previous
comp_time+=clock();
comm_time-=clock();
if((currentRank%2)==0)// check if process is even
 {
// first send data
MPI_Send(sender,d*n,MPI_DOUBLE,toSend,1,MPI_COMM_WORLD);
// then receive
 MPI_Recv(sender,d*n,MPI_DOUBLE,toReceive,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
 }
 else {
// first receive
   MPI_Recv(rec,d*n,MPI_DOUBLE,toReceive,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
// then send
   MPI_Send(sender,d*n,MPI_DOUBLE,toSend,1,MPI_COMM_WORLD);
// change pointers
changePointer(double,&sender,&rec);
}
comm_time+=clock();

    }

        if(currentRank== 0) {
            printf("TOTAL:%lf\n",(double) (comp_time + comm_time)/CLOCKS_PER_SEC);
            printf("Computations:%lf\n",(double)comp_time/CLOCKS_PER_SEC);
            printf("Communications:%lf\n",(double)comm_time/CLOCKS_PER_SEC);
        }

    free(sender);
  free(rec);
    return value;
  }
