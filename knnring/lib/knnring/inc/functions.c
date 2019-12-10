#ifndef FUNCTIONS_C
#define FUNCTIONS_C
// define macro function that will copy an array of any type to another array
#define copy(T,copyDest,  copyOrig, size){memcpy(copyDest, copyOrig, size *sizeof(T));}
// function that will manage indexes because by default they are 0 based
void manageIndexes(knnresult* result, int fixingValue)
{
	for (int i=0; i<result->m*result->k; i++) // loop through the k*m array
		result->nidx[i] += fixingValue;// add the fixingValue
    // converts the 0 based indexes according to the number of proccess
}
// function that will evaluate given the number of procceses and the size of corpus
// the fixingValue that needs to be added to the indexes

int calculateFixingValue(int fixingValue,int procceses , int n)
{
  return ((fixingValue - n) < 0) ? (fixingValue - n)%(procceses*n) + procceses*n : (fixingValue - n)%(procceses*n);
}
// function to compare the results of the received query and update the final
// result of each proccess
void compareAndReplace(knnresult* finalResult, knnresult* newResult){
  // we will use the distances array and Ids array so when we finish
  // mergining the ndist array will still be sorted
  int finalResultCounter, newResultCounter; // counter to help with
  // indexing in the comparison of the two ndis  arrays
      int m = newResult->m; // store somewhere the m value
      int k = newResult->k; // store the k value also
      double distances[k*m]; // declare double array that is going to be used
      // to copy ndist array from the finalResult
      int Ids[k*m]; // declare int array that is going to be used to copy
      // nid array from finalResult

      ;   //copies the data of store
    copy (double,distances,finalResult->ndist,k*m);// copy ndist from finalResult
    // to distances
    copy(int,Ids,finalResult->nidx,k*m); // copy nidx from finalResult to
    // distances

  // loop through all points of query
      for(int i=0; i<m; i++)
      {
          finalResultCounter= newResultCounter = 0;
  // we start by pointing to zero index
          //we only need to loop through the k neighbors
          for(int j=0; j<k; j++)
          {
              // check which neighbor is closer
              if(distances[i + finalResultCounter*m] < newResult->ndist[i + newResultCounter*m])
              // if it's the one from distances array
              {
                  finalResult->ndist[i + j*m] = distances[i + finalResultCounter*m];
                  // set  this distance
                  finalResult->nidx[i + j*m] = Ids[i + finalResultCounter*m];
                  // and update which index this distance has
                  finalResultCounter++; // increase counter which
                  // indicates that we can move to the next neighbor of
                  // distances array
              }
              else // if newResult has smaller distance
              {
                  finalResult->ndist[i + j*m] = newResult->ndist[i + newResultCounter*m];
                  // set this value to the final result
                  finalResult->nidx[i + j*m] = newResult->nidx[i + newResultCounter*m];
                  // and update which index this distance has
                  newResultCounter++;// increase counter which
                  // indicates that we can move to the next neighbor of
                  // newResult ndist array

              }
          }
      }
  }
	void manageMpiStuff (int * proccesses , int * currentRank , int * toSend , int *toReceive){

	    MPI_Comm_size(MPI_COMM_WORLD, proccesses);// take the num of proccesses
	    MPI_Comm_rank(MPI_COMM_WORLD, currentRank); // take the rank of the proccess



	// the n proccess expects to receive from the n-1 proccess except the 0 proccess
	// that expects to receive from n-1 proccess
	*toReceive= *currentRank==0 ? *proccesses-1:*currentRank-1;
	// in the same way proccess n sends to proccess n+1 except of n-1 proccess that
	// sends to 0 proccess
	*toSend=(*currentRank+1)%*proccesses;
	}
#endif
