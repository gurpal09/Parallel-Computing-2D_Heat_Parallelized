/*
 * Numerical and analytical solution of the 2D Heat equation with conductivity = 1
 * Only compare the solutions for small increments (i.e 10x10) due to printing to screen
 *
 *
 * Author: Gurpal Singh
 * Date: 04/28/2017
 *
 *
 * to compile: gcc -O3 -std=c99 -lm -DRESTRICT=restrict NUM_vs_AN.c -o NUM_vs_AN.exe
 * to execute: ./NUM_vs_AN.exe <number of time steps>
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

#define NX 10   /* includes boundary points on both end */
#define NY 10   /* includes boundary points on both end */
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

//Function for Finite Difference Scheme
void solveWave(REAL *RESTRICT unew, const REAL *RESTRICT u)
{
    INT i, j;
    for (j = 1; j < NY; j++) {
        for (i = 1; i < NX; i++) {
            unew[IC] = (((u[IP1] - 2.0f * u[IC] + u[IM1]) / (dx * dx))
                        + ((u[JP1] - 2.0f * u[IC] + u[JM1]) / (dy * dy)))
                       * dt
                       + u[IC];
        }
    }
}

//Function for Setting Boundary Conditions
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

//Function for Analytical Solution 
void analyticalSoln(INT nTimeSteps, REAL *RESTRICT uAnalytical, REAL *RESTRICT XX,
                    REAL *RESTRICT YY, INT timesteps)
{
    INT i, j, n, limit;
    limit     = 101;
    double T0 = 100;
    double L  = LX;
    REAL   x, y, sum;

    for (j = JSTART; j < JEND + 1; j++) {
        for (i = ISTART; i < IEND + 1; i++) {
            x = i * dx;
            y = j * dy;
            sum = 0.0;

            for (n = 1; n < limit; n = n + 2) {
                sum = sum + ((4.0 * T0 / (n * PI)) * sin((n * PI * x) / L) * sinh((n * PI * y) / L)
                             / sinh(n * PI));
            }
            uAnalytical[IC] = sum;
        }
    }

    //Printing the Analytical Soln in reverse
    for (j = JEND + 1; j-- > 0;) {
        for (i = IEND + 1; i-- > 0;) {
            printf("%.2f\t", uAnalytical[IC]);
        }
        printf("\n\n");
    }
    
    //Writing it to A File	
    FILE *output;
    output = fopen("Numerical_vs_Analytical.txt", "a+");

    fprintf(output, "2D Conduction Result Comparison\n");
    fprintf(output, "Square Plate Length: %.2f\n", LX);
    fprintf(output, "Number of Timesteps: %d\n", timesteps);
    fprintf(output, "\n\n");
    fprintf(output, "***************Analytical Solution***************\n");
    for (j = JEND + 1; j-- > 0;) {
        for (i = IEND + 1; i-- > 0;) {
            fprintf(output, "%.2f\t", uAnalytical[IC]);
        }
        fprintf(output, "\n\n");
    }
    fprintf(output, "\n\n");
    fclose(output);
}

//Function for mesh
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

//Function to Write Output
void writeOutput(REAL *RESTRICT phi)
{
    INT   i, j;
    FILE *output;
    output = fopen("Numerical_vs_Analytical.txt", "a+");

    fprintf(output, "***************Numerical Solution***************\n");

    for (j = JSTART; j < JEND + 1; j++) {
        for (i = ISTART; i < IEND + 1; i++) {
            fprintf(output, "%.2f\t", phi[IC]);
        }
        fprintf(output, "\n\n");
    }
    fclose(output);
}

//Function for printing to display
void print2Display(REAL *phi)
{
    INT i, j;
    for (j = JSTART; j < JEND + 1; j++) {
        for (i = ISTART; i < IEND + 1; i++) {
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

    //Dynamic Memory Allocation	
    REAL *x = calloc(NX * NY, sizeof *x);
    REAL *y = calloc(NX * NY, sizeof *y);

    REAL *unew  = calloc(NX * NY + NX + NY, sizeof *unew);
    REAL *u     = calloc(NX * NY + NX + NY, sizeof *u);

    REAL *exact = calloc(NX * NY + NX + NY, sizeof *exact);
    REAL *tmp;

    //Create mesh	
    meshGrid(x, y);

    //Set BC
    initWave(u, unew, x, y);
    writeOutput(u);
    
    //Printing Initial Temp 
    printf("|||||||----Initialized----|||||||||\n");
    print2Display(u);
    printf("||||||||||||||||||||||||||||||\n");

    //Numerical Calculation
    for (INT n = 1; n <= nTimeSteps; n++) {
        solveWave(unew, u);

        // Pointer Swap
        tmp  = u;
        u    = unew;
        unew = tmp;
    }

    // Analytical Solution
    printf("|||||||----ANALYTICAL SOLUTION----|||||||||\n");
    analyticalSoln(nTimeSteps, exact, x, y, nTimeSteps);
    printf("||||||||||||||||||||||||||||||\n");
    
    //Print Numerical Solution
    printf("|||||||----NUMERICAL SOLUTION----|||||||||\n");
    print2Display(u);
    printf("||||||||||||||||||||||||||||||\n");

    //Clean up
    writeOutput(u);
    free(unew);
    free(u);
    free(exact);
    free(x);
    free(y);

    return EXIT_SUCCESS;
}

