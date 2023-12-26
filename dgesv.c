#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int _NUM = 4;

// Solves the system of linear equations Ax = B using Gaussian elimination with partial pivoting
// Parameters:
//   n: Size of the system
//   a: Coefficient matrix A (input), overwritten with the upper triangular form (output)
//   b: Right-hand side vector B (input), overwritten with the solution vector x (output)

void my_dgesv(int n, double *a, double *b) {
    int i, j, k;
    double scalingFactor;

    double* _a = (double*) __builtin_assume_aligned(a, _NUM);
    double* _b = (double*) __builtin_assume_aligned(b, _NUM);

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

        #pragma omp simd aligned(_a:N)
        for (j = 0; j < n; j++) {
            double temp_a = a[rowOrder[i] * n + j];
            _a[rowOrder[i] * n + j] = _a[rowOrder[pivotRow] * n + j];
            _a[rowOrder[pivotRow] * n + j] = temp_a;
        }

        // Swap rows in the B matrix using the same pivot order

        #pragma omp simd aligned(_b:N)

        for (j = 0; j < n; j++) {
            double temp_b = a[rowOrder[i] * n + j];
            _b[rowOrder[i] * n + j] = _b[rowOrder[pivotRow] * n + j];
            _b[rowOrder[pivotRow] * n + j] = temp_b;
        }

        scalingFactor = a[rowOrder[i] * n + i];


        #pragma omp simd aligned(_a:N)
        for (j = i; j < n; j++) {
            _a[rowOrder[i] * n + j] /= scalingFactor;
        }


        #pragma omp simd aligned(_b:N)
        for (j = 0; j < n; j++) {
            _b[rowOrder[i] * n + j] /= scalingFactor;
        }

        // Eliminate other rows using the pivot row
        for (k = 0; k < n; k++) {
            if (k != i) {
                scalingFactor = a[rowOrder[k] * n + i];


                #pragma omp simd aligned(_a:N)
                for (j = i; j < n; j++) {
                    _a[rowOrder[k] * n + j] -= scalingFactor * _a[rowOrder[i] * n + j];
                }


                #pragma omp simd aligned(_a:N)

                for (j = 0; j < n; j++) {
                    b[rowOrder[k] * n + j] -= scalingFactor * _b[rowOrder[i] * n + j];
                }
            }
        }
    }

    // Free the allocated memory for pivot order
    free(rowOrder);
}
