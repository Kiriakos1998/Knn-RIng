

#include <cblas.h>
#include "knnring.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
// define function that takes as input the type of two number int float double
// and the two numbers and swaps them
#define SWAP(T, x, y) \
  do { \
    T tmp = x; \
    x = y; \
    y = tmp; \
  } while (0);
// function that will sort the distances for one point of corpus X from all
// the points of the query Y
void sort(double* distances, int* indexes, int begin, int terminal, int
  step);
// index  is the pointer of an integer array and size its size
void initIndexes(int *indexes,int size){
int i; // counter
  for(i=0;i<size;i++)// loop
  indexes[i]=i; // set its index value to be the index
}
double *calculateDistances(double *X, double *Y, int n, int m, int d)

{
//allocate memory for the squares of each point and the final array of distances
double *squareX = malloc(n * sizeof(double)); // allocate memory for X^2 array
double *squareY = malloc(m * sizeof(double)); // allocate memory for Y^2 array
  double *finalDistances = malloc(n*m * sizeof(double)); // allocate memory for
    // the distances
    double roott;// variable to store the distance^2 value
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, d,-2.0, X,
      n, Y, m,0.0, finalDistances, m);
// calculate -2*X*Y array


    for (int i = 0; i < n; i++)
        {
          squareX[i] = cblas_ddot(d, &X[i], n,&X [i], n);
          // calculate the X^2 value for each point in corpus
}
    for (int i = 0; i < m; i++)
      {
         squareY[i] = cblas_ddot(d, &Y[i], m, &Y[i], m);
         // calculate the Y^2 value for each point in query
}
    for (int i = 0; i < n; i++)
{
  for (int j = 0; j < m; j++)
          {
            roott=sqrt(fabs(finalDistances[i*m + j] + squareX[i] + squareY[j]));
            // calculate the distance square X^2+Y^2-2*X*Y
         finalDistances[i*m + j] =roott;// calculate the final distance

            if(isnan(roott)) // check if the value is nan, same points will not
            finalDistances[i*m+j]=0; // lead to 0 distance but due to underflow
            // will have nan value
}
}
    free(squareX); // free  memory  allocated for  X^2
    free(squareY); // free memory allocated for Y^2

    return finalDistances; // return the distances
}

////////////////////////////////////////////////////////////////////////////


knnresult kNN(double *X, double *Y, int n, int m, int d, int k)
{
  double *finalDistances = calculateDistances(X, Y, n, m, d);

  knnresult value; // create a knnresult instance
  value.k = k; // set value k value to be k
  value.m = m; // set value m value to be m
// initialize nidx allocating the neccessery memory
  value.nidx = (int *) malloc(k * m * sizeof(int));

  int *indexes = (int *) malloc(n * sizeof(int));
  //array to hold information about
  //indexing, is going to be copied to value.nidx
  for(int i=0; i<m; i++)
  {
    initIndexes(indexes,n);// initialize the indexes

      sort(finalDistances + i, indexes, 0, n - 1, m);   // sort all the columns
      // from the smaller to the bigger distances so indexes from 0 to k-1
      // hold the k nearest neighbors
      // indexes holds the real indexes of these neighbors

     for(int j=0; j<k; j++)
         value.nidx[m*j + i] = indexes[j];    //update indexes to the value
         // that is to be returned
  }
  free(indexes); // free memory array indexes no longer needed


  // realloc only the memory that is needed
  value.ndist = (double *) realloc(finalDistances, k * m * sizeof(double));

  return value;

}
// distances is the array that holds the finalDistances
// indexes the array that will store the indexes
// begin shows which is the first element of the array to quickSort
// terminal shows which is the last element
// step is parametre that helps us navigate through distances array to
// take the correct data for each points, step takes the m value
void sort(double* distances, int* indexes, int begin, int terminal, int step)
{

int i; // declare counter
    if(begin < terminal) // check if the begin is less than terminal
    {
        int temporary = begin; // declare  temporary variable to help with the
        //partion of the array in the recruision
        double pivot = distances[terminal*step];// set the pivot value
        for ( i = begin; i <= terminal; i++) // loop through the array to be
        {//sorted
            if (distances[i*step] <= pivot) // check if pivot is greter or equal

            {
                SWAP(double,distances[ i*step], distances[temporary*step]);
                // swap the numbers
                SWAP(int,indexes[i], indexes[temporary]);
                // swap theri indexes
                temporary++; // increase temporary
            }
          }
        temporary--; // decrease it
// in the end temporary shows how many values of the array are smaller than
//pivot, by recruision each time the pivot will have on its right all the larger
// number and on its left all the smaller untill the array is sorted
        sort(distances, indexes, begin, temporary - 1, step);
        // sort the left part of the array
        sort(distances, indexes, temporary+ 1, terminal, step);
// sort the right part of the array

}
}
