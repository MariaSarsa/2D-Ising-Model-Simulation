/** Function that computes the total absolute magnetization per spin of a configuration 's' **/

#include <stdio.h>
#include <stdlib.h>

extern int L;

double magnetization (int *s)
{
    int i;
    double total_mag;
    total_mag = 0;

    for (i = 0; i < L * L; i++)
        total_mag += s[i];

    return total_mag / (L * L);
}