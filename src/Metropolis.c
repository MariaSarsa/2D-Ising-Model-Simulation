/** Function that performs a Monte-Carlo iteration using the Metropolis algorithm **/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define NormRANu (2.3283063671E-10F)

int L;

// Definition of variables and functions for the Parisi-Rapuano random number generator
float Random(void);
void ini_ran(int SEMILLA);

void Metropolis(int *s, int *xp, int *xm, int *yp, int *ym, double *prob)
{
    int x, y, n, Index;
    n = 0;
    Index = 0;

    for(y = 0; y < L; y++)
    {
        for(x = 0; x < L; x++)
        {
            // Calculate the index based on the spin interaction with its 4 nearest neighbors:
            // xp/xm (x plus/minus, i.e., East/West) and yp/ym (y plus/minus, i.e., North/South)
            Index = s[n] * (s[n + xp[x]] + s[n + yp[y]] + s[n + xm[x]] + s[n + ym[y]]) / 2 + 2;
            
            // Acceptance criteria via Metropolis choice
            if(Random() < prob[Index])
            {
                s[n] = -s[n]; // Flip the spin
            }
            n++;
        }
    }
}