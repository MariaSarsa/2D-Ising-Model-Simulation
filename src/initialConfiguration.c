/** Function that generates the initial lattice spin configuration and saves it to array 's' **/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define NormRANu (2.3283063671E-10F)

int L;

// Definition of variables and functions for the Parisi-Rapuano random number generator
float Random(void);
void ini_ran(int SEMILLA);

void initialConfiguration (int flag, int *s)
{
    int i, j, n, ChessSize;
    FILE *file_in;

    n = 0;
    ChessSize = 0;

    /** Initialize the initial configuration vector to 0 **/
    for (i = 0; i < L * L; i++)
        s[i] = 0;

    switch(flag)
    {
    case 0: // Random configuration (High Temperature limit T -> infinity)
    {
        for (i = 0; i < (L * L); i++)
        {
            if ((Random() > 0.5))
                s[i] = 1;
            else
                s[i] = -1;
        }
        break;
    }

    case 1: // Frozen configuration (Ferromagnetic ground state, all spins aligned)
    {
        if(Random() > 0.5)
        {
            for (i = 0; i < (L * L); i++)
                s[i] = 1;
            printf("Spin: 1\n");
        }
        else
        {
            for (i = 0; i < (L * L); i++)
                s[i] = -1;
            printf("Spin: -1\n");
        }
        break;
    }

    case 2: // Read a previously saved initial configuration file
    {
        file_in = fopen("confinicial.dat", "rb");
        if(file_in != NULL) {
            fread(s, sizeof(s[0]), L * L, file_in);
            fclose(file_in);
        }
        break;
    }

    case 3: // Anti-ferromagnetic chessboard pattern
    {
        printf("Enter an even number:\n");
        scanf("%d", &ChessSize);

        for (j = 0; j < ChessSize; j++)
        {
            for(i = 0; i < ChessSize; i++)
            {
                if ((i + j) % 2 == 0) s[n] = 1;
                else s[n] = -1;
                n++;
            }
        }
        break;
    }
    }
}