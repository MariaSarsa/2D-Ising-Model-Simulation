#include <stdio.h>
#include <stdlib.h>

// Function to construct data blocks used for system thermalization checks and error estimation (Cv, X)
// x: input data array (energy, magnetization, etc.)
// xb: output array storing the average values of each block
// n_b: number of blocks
// b_s: size of each block (number of elements per block)

void Build_Blocks(double *x, double *xb, int n_b, int b_s)
{
    int i, j;
    for(i = 0; i < n_b; i++)
    {
        xb[i] = 0;
        for(j = 0; j < b_s; j++)
        {
            xb[i] += x[i * b_s + j];
        }
        xb[i] /= b_s; // Average of the current block
    }
}