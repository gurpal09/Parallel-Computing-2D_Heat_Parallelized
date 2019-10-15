/*
 * Purpose: This code shows that the MPI also computes the correct solution as
 * the serial version. The results are printed onto the screen as well as saved
 * Computes for a 10x10 square plate with 1000 time steps
 *
 * Author: Gurpal Singh
 * Date: 04/28/2017
 *
 * to compile: mpicc -std=c00 MPI_plot.c -O3 -mpi_plot.exe
 * to execute: mpirun -np 4 ./mpi_plot.exe 1000
 *
 */

//#include "timer.h"
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#define NX 10   /* includes boundary points on both end */
#define NY NX   /* includes boundary points on both end */
#define LX 1.0f /* length of the domain in x-direction  */
#define LY 1.0f /* length of the domain in x-direction  */
#define dx (REAL)(LX / ((REAL)(NX)))
#define dy (REAL)(LY / ((REAL)(NY)))
#define dt (REAL) 0.5*dx*dx
#define ISTART 0         /* mesh zone to print values to the screen */
#define IEND ISTART + 10 /* mesh zone to print values to the screen */
#define JSTART 0         /* mesh zone to print values to the screen */
#define JEND JSTART + 10 /* mesh zone to print values to the screen */

#define IC (i + j * NX)        /* (i,j)   */
#define IM1 (i + j * NX - 1)   /* (i-1,j) */
#define IP1 (i + j * NX + 1)   /* (i+1,j) */
#define JM1 (i + (j - 1) * NX) /* (i,j-1) */
#define JP1 (i + (j + 1) * NX) /* (i,j+1) */

#ifndef RESTRICT
#define RESTRICT
#endif

#ifndef SINGLE
typedef double REAL;
typedef int    INT;
#define PI 3.14159265358979323846
#else
typedef float REAL;
typedef int   INT;
#define PI 3.1415927f
#endif

// Function to solve equation using Finite Difference Scheme
void solveWave(REAL *RESTRICT unew, const REAL *RESTRICT u, INT ny)
{
    INT i, j;

    for (j = 1; j < ny; j++) {
        for (i = 1; i < NX; i++) {
            unew[IC] = (((u[IP1] - 2.0f * u[IC] + u[IM1]) / (dx * dx))
                        + ((u[JP1] - 2.0f * u[IC] + u[JM1]) / (dy * dy)))
                       * dt
                       + u[IC];
        }
    }
}

// Function to Initialize Boundary Conditions
void initWave(REAL *RESTRICT u, REAL *RESTRICT unew, REAL *RESTRICT x, REAL *RESTRICT y, INT ny)
{
    INT i, j;
    for (j = 0; j < ny; j++) {
        for (i = 0; i < NX; i++) {
            if (j == 0) {
                u[IC]    = 100;
                unew[IC] = 100;
            } else {
                u[IC]    = 0;
                unew[IC] = 0;
            }
        }
    }
}

// Function to create Mesh
void meshGrid(REAL *RESTRICT x, REAL *RESTRICT y)
{
    INT  i, j;
    REAL a;

    for (j = 0; j < NY; j++) {
        a = dx * ((REAL) j);
        for (i = 0; i < NX; i++) {
            x[IC] = dx * ((REAL) i);
            y[IC] = a;
        }
    }
}

// Function to write Output to file
void writeOutput(REAL *RESTRICT phi)
{
    INT   i, j;
    FILE *output;
    output = fopen("MPI_plot_result.txt", "a+");
    fprintf(output,"***************MPI RESULT***************\n");

    for (j = JSTART; j < JEND+1; j++) {
        for (i = ISTART; IEND+1 < NX; i++) {
            fprintf(output, "%.2f\t", phi[IC]);
        }
	fprintf(output,"\n\n");
    }
    fclose(output);
}

// Function for printing to display
void print2Display(REAL *phi, INT ny)
{
    INT i, j;

    for (j = JSTART; j < JEND+1; j++) {
        for (i = ISTART; i < IEND+1; i++) {
            printf("%.2f\t", phi[IC]);
        }
        printf("\n\n");
    }
}

// Function for Decomposing via 1D slicing
void decomposeMesh_1D(const int N, const int nProcs, const int myRank, int *start, int *end,
                      const int nGhostLayers)
{
    int remainder = N % nProcs;
    if (remainder == 0) {
        *start = 0;
        *end   = (N / nProcs) + nGhostLayers - 1;
    } else {
        *start               = 0;
        int pointsPerProcess = (N - remainder) / nProcs + 1;
        if (myRank == (nProcs - 1))
            *end = (N - pointsPerProcess * (nProcs - 1)) + nGhostLayers - 1;
        else
            *end = pointsPerProcess + nGhostLayers - 1;
    }
}

// Function for exchanging ghost cells
void exchange_SendRecv(double *phi, const int start, const int end, int src, int dest,
                       const int myRank, const int nProcs, MPI_Comm comm1D)
{
    int tag0 = 0;
    int tag1 = 1;
    int e    = (end - 1) * NX;
    int s    = 0;

    MPI_Sendrecv(&phi[e], NX, MPI_DOUBLE, dest, tag0, &phi[s], NX, MPI_DOUBLE, src, tag0, comm1D,
                 MPI_STATUS_IGNORE);

    s = (start + 1) * NX;
    e = end * NX;

    // NOTE WE ARE NOW SWAPPING THE SRC & DEST
    int tmp = dest;
    dest    = src;
    src     = tmp;
    MPI_Sendrecv(&phi[s], NX, MPI_DOUBLE, dest, tag1, &phi[e], NX, MPI_DOUBLE, src, tag1, comm1D,
                 MPI_STATUS_IGNORE);
}

int main(INT argc, char *argv[])
{
    int nProcs; /* number of processes */
    int myRank; /* process rank */
    int src;    /* handles for communication, source process id */
    int dest;   /* handles for communication, destination process id */
    int start;  /* start index for each partial domain */
    int end;    /* end index for each partial domain */

    MPI_Init(&argc, &argv);                 /* initialize MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); /* get the number of processes */

    int nDims = 1; // dimension of Cartesian decomposition 1 => slices
    int dimension[nDims];
    int isPeriodic[nDims];
    int reorder = 1; // allow system to optimize(reorder) the mapping of processes to physical cores

    dimension[0]  = nProcs;
    isPeriodic[0] = 0; // periodicty of each dimension

    MPI_Comm comm1D; // define a communicator that would be assigned a new topology
    MPI_Cart_create(MPI_COMM_WORLD, nDims, dimension, isPeriodic, reorder, &comm1D);
    MPI_Comm_rank(comm1D, &myRank); /* get the rank of a process after REORDERING! */
    MPI_Cart_shift(comm1D, 0, 1, &src,
                   &dest); /* Let MPI find out the rank of processes for source and destination */

    int nGhostLayers;
    if (myRank == 0 || myRank == nProcs - 1) {
        nGhostLayers = 1;
    } else {
        nGhostLayers = 2;
    }

    decomposeMesh_1D(NY, nProcs, myRank, &start, &end, nGhostLayers);

    int ny = (end - start) + 1;

    double *u    = calloc(NX * ny, sizeof(*u));
    double *unew = calloc(NX * ny, sizeof(*unew));

    MPI_Barrier(comm1D); // make sure all processes initialized their portion of the problem

    if (argc != 2) {
        perror("Command-line usage: executableName <# time steps>");
        exit(1);
    }

    INT nTimeSteps = atoi(argv[1]);

    REAL *x = calloc(NX * NY, sizeof *x);
    REAL *y = calloc(NX * NY, sizeof *y);
    REAL *tmp;

    meshGrid(x, y);

    //Initialize the BC only on the top slice	
    if (myRank == 0) {
        initWave(u, unew, x, y, ny);
    }

    
    //Loop for Calculation
    for (INT n = 1; n <= nTimeSteps; n++) {
        solveWave(unew, u, ny);

        // Pointer Swap
        exchange_SendRecv(unew, start, end, src, dest, myRank, nProcs, comm1D);
        MPI_Barrier(comm1D);
        tmp  = u;
        u    = unew;
        unew = tmp;
    }

    MPI_Barrier(comm1D);

    int *new = calloc(4, sizeof(*new));
    
    //Gathering ny into an array needed for size of arrays
    MPI_Gather(&ny, 1, MPI_INT, new, 1, MPI_INT, 0, comm1D);

    //Allocate memory for result	
    double *result = calloc(NX * NY + 6 * NX, sizeof(*result));

    int num_no_ghost;
    if (myRank == 0 || myRank == nProcs - 1) {
        num_no_ghost = NX * ny - NX;
    } else {
        num_no_ghost = NX * ny - 2 * NX;
    }

    int *new_no_ghost = calloc(nProcs, sizeof(*new_no_ghost));
    MPI_Gather(&num_no_ghost, 1, MPI_INT, new_no_ghost, 1, MPI_INT, 0, comm1D);

    int *displs_no_ghost = calloc(nProcs, sizeof(*displs_no_ghost));
    int  sum             = 0;
    for (int i = 0; i < nProcs; i++) {
        displs_no_ghost[i] = sum;
        sum                = sum + new_no_ghost[i];
    }

    MPI_Gatherv(u, num_no_ghost, MPI_DOUBLE, result, new_no_ghost, displs_no_ghost, MPI_DOUBLE, 0,
                comm1D);

    if (myRank == 0) {
        printf("Result:\n\n");
        print2Display(result, NY);
        writeOutput(result);
    }

    free(unew);
    free(u);
    free(x);
    free(y);
    free(result);
    free(new);
    free(displs_no_ghost);
    free(new_no_ghost);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

