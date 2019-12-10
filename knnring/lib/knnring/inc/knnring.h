#ifndef KNNRING_H
#define KNNRING_H

typedef struct knnresult{
int *nidx; // Indices (0-based) of nearest neighbors
double *ndist; //  Distance of nearest neighbors
int m;//Number of query points
int k; // Number of nearest neighbors
  } knnresult;


knnresult kNN(double * X, double * Y, int n, int m, int d, int k);
knnresult distrAllkNN(double* X, int n, int d, int k);



#endif
