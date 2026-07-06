/** Function that clears and zeroes out the measurement and statistics variables **/

#include <stdio.h>
#include <stdlib.h>

int Nmed;

// Note: Structural statistical values are passed by pointer references to reset variables globally
void Initialize(double *e, double *M, double *e_avg, double *M_avg, double *des_e, double *des_M, double *des_Cc, double *des_X, double *e_sq, double *M_sq)
{
    int i;

    for (i = 0; i < Nmed; i++)
    {
        e[i] = 0;
        M[i] = 0;
    }
    *e_avg = *M_avg = *des_e = *des_M = *des_Cc = *des_X = *e_sq = *M_sq = 0;
}