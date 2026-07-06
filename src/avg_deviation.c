/** Function that calculates: mean, error (standard deviation of the mean), and mean of squares for a given array (energies or magnetizations) **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int L;

void avg_deviation (double *s, double *average, double *sigma, double *sq, int N)
{
    int i;

    for (i = 0; i < N; i++)
    {
        *average += s[i];
        *sq += s[i] * (s[i]);
    }

    *average = *average / N;
    *sq = *sq / N;
    // Standard error calculation using statistical variance divided by sqrt(N-1)
    *sigma = sqrt((*sq - *average * (*average))) / sqrt(N - 1.0);
}