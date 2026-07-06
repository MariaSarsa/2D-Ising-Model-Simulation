/** Function that computes the total energy per spin of a given configuration 's' **/

#include <stdio.h>
#include <stdlib.h>

int L;

double energy (int *s, int *xp, int *yp)
{
    int i, j, n;
    double total_energy;

    n = 0;
    total_energy = 0;

    for (j = 0; j < L; j++)
    {
        for(i = 0; i < L; i++)
        {
            // Sum interactions with positive offset steps (periodic boundary conditions managed by xp and yp)
            total_energy += s[n] * (s[(n + xp[i])] + s[(n + yp[j])]);
            n++;
        }
    }
    // Normalization factor: negative sign from H = -J*sum(s_i*s_j), divided by 2*L^2 to avoid double counting bonds
    return (-total_energy / (2 * L * L));
}