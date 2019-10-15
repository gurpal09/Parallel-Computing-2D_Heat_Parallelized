/*
 * Serial code to time the finite difference method. Input #timesteps
 *
 *
 *
 * precise time measurements are enabled with the GET_TIME macros that are defined
 * in timer.h
 *
 * Author: Gurpal Singh
 * Date: 04/28/2017
 *
 * to compile:
 * gcc -O3 -std=c99 -lm -DRESTRICT=restrict SERIAL_FD.c -o SERIAL_FD.exe
 * to execute: ./SERIAL_FD.exe <number of time steps>
 *
 */

#include "timer.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#define NX 20000 /* includes boundary points on both end */

#define NY NX   /* includes boundary points on both end */
#define LX 1.0f /* length of the domain in x-direction  */
#define LY 1.0f /* length of the domain in x-direction  */
#define dx (REAL)(LX / ((REAL)(NX)))
#define dy (REAL)(LY / ((REAL)(NY)))
#define dt (REAL)(0.5 * dx * dx)

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

// Function to Solve wave using Finite Difference Scheme
void solveWave(REAL *RESTRICT unew, const REAL *RESTRICT u)
{
    INT i, j;

    for (j = 1; j < NY - 1; j++) {
        for (i = 1; i < NX - 1; i++) {
            unew[IC] = (((u[IP1] - 2.0f * u[IC] + u[IM1]) / (dx * dx))
                        + ((u[JP1] - 2.0f * u[IC] + u[JM1]) / (dy * dy)))
                       * dt
                       + u[IC];
        }
    }
}

// Function to Initialize BC
void initWave(REAL *RESTRICT u, REAL *RESTRICT unew, REAL *RESTRICT x, REAL *RESTRICT y)
{
    // Applying the Boundary Conditions
    INT i, j;
    for (j = 0; j < NY; j++) {
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

// Function to create 2D Mesh
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

// Function for Printing to display
void print2Display(REAL *phi)
{
    INT i, j;

    for (j = JSTART; j < JEND; j++) {
        for (i = ISTART; i < IEND; i++) {
            printf("%.2f\t", phi[IC]);
        }
        printf("\n\n");
    }
}

int main(INT argc, char *argv[])
{
    if (argc != 2) {
        perror("Command-line usage: executableName <# time steps>");
        exit(1);
    }

    INT nTimeSteps = atoi(argv[1]);

    // Dynamic Memory Allocation
    REAL *x = calloc(NX * NY, sizeof *x);
    REAL *y = calloc(NX * NY, sizeof *y);

    REAL *unew  = calloc(NX * NY, sizeof *unew);
    REAL *u     = calloc(NX * NY, sizeof *u);
    REAL *uold  = calloc(NX * NY, sizeof *uold);
    REAL *exact = calloc(NX * NY, sizeof *exact);
    REAL *tmp;

    // Create Mesh
    meshGrid(x, y);
    // Set BC
    initWave(u, unew, x, y);

    // Start Timing
    StartTimer();

    // TimeMarch
    for (INT n = 1; n <= nTimeSteps; n++) {
        solveWave(unew, u);

        // Pointer Swap
        tmp  = u;
        u    = unew;
        unew = tmp;
    }

    // Stop Timing
    double elapsedTime = GetTimer();

    // Print The Numerical Result
    /*
    printf("|||||||----NUMERICAL SOLUTION----|||||||||\n");
    print2Display(u);
    printf("||||||||||||||||||||||||||||||\n");
    */

    // Print Elapsed Time
    printf("Time elapsed = %f s\n", elapsedTime);

    // Save Results to a File
    FILE *file;
    if (fopen("Serial_Timing.txt", "r") == NULL) {
        file = fopen("Serial_Timing.txt", "a+");
        fprintf(file, "Serial Code Timing Results:\n");	
        fprintf(file, "# of TimeSteps: %d\ni",nTimeSteps);
        fprintf(file, "Size:\t Time:\n");
        fprintf(file, "%d\t%.4f\n", NX, elapsedTime);
        fclose(file);
    }

    else {
        file = fopen("Serial_Timing.txt", "a");
        fprintf(file, "%d\t%.4f\n", NX, elapsedTime);
        fclose(file);
    }

    free(unew);
    free(u);
    free(exact);
    free(x);
    free(y);

    return EXIT_SUCCESS;
}

