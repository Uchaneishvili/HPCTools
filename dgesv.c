#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int _NUM = 8;

// Solves the system of linear equations Ax = B using Gaussian elimination with partial pivoting
// Parameters:
//   n: Size of the system
//   a: Coefficient matrix A (input), overwritten with the upper triangular form (output)
//   b: Right-hand side vector B (input), overwritten with the solution vector x (output)

void my_dgesv(int n, double *a, double *b) {
    int i, j, k;
    double scalingFactor;



    // Array to store the pivot order    
    int *rowOrder = (int *)malloc(n * sizeof(int));

    // Initialize pivot order to the identity permutation
    for (i = 0; i < n; i++) {
        rowOrder[i] = i;
    }

    // Gaussian elimination with partial pivoting
    for (i = 0; i < n; i++) {
        scalingFactor = fabs(a[rowOrder[i] * n + i]);
        int pivotRow = i;

        // Search for a larger pivot element in the remaining rows
        for (j = i + 1; j < n; j++) {
            double abs_val = fabs(a[rowOrder[j] * n + i]);
            if (abs_val > scalingFactor) {
                scalingFactor = abs_val;
                pivotRow = j;
            }
        }

        // Swap rows in the A matrix using the pivot order
        int temp = rowOrder[i];
        rowOrder[i] = rowOrder[pivotRow];
        rowOrder[pivotRow] = temp;

        // Swap corresponding rows in the B matrix

        #pragma vector always
        #pragma loop count (n)
        for (j = 0; j < n; j++) {
            double temp_a = a[rowOrder[i] * n + j];
            a[rowOrder[i] * n + j] = a[rowOrder[pivotRow] * n + j];
            a[rowOrder[pivotRow] * n + j] = temp_a;
        }

        // Swap rows in the B matrix using the same pivot order

        #pragma vector always
        #pragma loop count (n)
        for (j = 0; j < n; j++) {
            double temp_b = a[rowOrder[i] * n + j];
            b[rowOrder[i] * n + j] = b[rowOrder[pivotRow] * n + j];
            b[rowOrder[pivotRow] * n + j] = temp_b;
        }

        scalingFactor = a[rowOrder[i] * n + i];


        #pragma vector always
        #pragma loop count (n)
        for (j = i; j < n; j++) {
            a[rowOrder[i] * n + j] /= scalingFactor;
        }


        #pragma vector always
        #pragma loop count (n)
        for (j = 0; j < n; j++) {
            b[rowOrder[i] * n + j] /= scalingFactor;
        }

        // Eliminate other rows using the pivot row
        for (k = 0; k < n; k++) {
            if (k != i) {
                scalingFactor = a[rowOrder[k] * n + i];


                #pragma vector always
                #pragma loop count (n)
                for (j = i; j < n; j++) {
                    a[rowOrder[k] * n + j] -= scalingFactor * a[rowOrder[i] * n + j];
                }


                #pragma vector always
                #pragma loop count (n)
                for (j = 0; j < n; j++) {
                    b[rowOrder[k] * n + j] -= scalingFactor * b[rowOrder[i] * n + j];
                }
            }
        }
    }

    // Free the allocated memory for pivot order
    free(rowOrder);
}
