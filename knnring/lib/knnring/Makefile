
SHELL := /bin/bash
CC = gcc-7
MPICC = mpicc
CFLAGS = -Wall -lopenblas -lblas -O3
INCLUDES = -I ./inc





lib: knnring_sequential.o knnring_synchronous.o knnring_asynchronous.o
	ar rcs lib/knnring_sequential.a lib/knnring_sequential.o

#	COMMENT OUT THE LINE BELLOW TO RUN SYNCHRONOYS MPI 
#	rm -f lib/knnring_mpi.a; ar rcs lib/knnring_mpi.a lib/knnring_synchronous.o lib/knnring_sequential.o

#	COMMENT OUT THE LINE BELLOW TO RUN ASYNCHRONOYS MPI
	rm -f lib/knnring_mpi.a ;ar rcs lib/knnring_mpi.a lib/knnring_asynchronous.o lib/knnring_sequential.o

	rm lib/*.o

knnring_synchronous.o:
	$(MPICC) $(CFLAGS) $(INCLUDES) -c src/knnring_synchronous.c -o lib/knnring_synchronous.o

knnring_asynchronous.o:
	$(MPICC) $(CFLAGS) $(INCLUDES) -c src/knnring_asynchronous.c -o lib/knnring_asynchronous.o

knnring_sequential.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c src/knnring_sequential.c -o lib/knnring_sequential.o
