/** Function that generates the Boltzmann transition probabilities for each beta **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void probability(double *prob,double beta)
{
    prob[0]=exp(beta*8);
    prob[1]=exp(beta*4);
    prob[2]=1;
    prob[3]=exp(beta*(-4));
    prob[4]=exp(beta*(-8));
}
