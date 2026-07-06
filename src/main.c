/** Main simulation file for the 2D Ising Model **/

/** Given an input parameters file ("input.txt") containing the temperature interval parameters (betas),
    number of measurements, and maximum thermalization steps. 
    Note: The lattice linear dimension size L is adjusted manually via line 25.
    
    The program returns two output files:
    - Measures.txt (stores final calculated average energy and magnetization along with their respective error bounds per beta step)
    - Cv_X.txt (stores calculated specific heat capacity and magnetic susceptibility alongside error bars per beta step)

    In debugging/testing mode, it can also produce:
    - check.txt: outputs exact completed thermalization cycles tracking convergence step targets per beta increment.
    - An individual file per beta value capturing transient energy and magnetization traces to inspect raw thermalization pathways.
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Header.h"
#define NormRANu (2.3283063671E-10F)

// Constants for estimating Specific Heat (Cv) and Susceptibility (X) block errors
// Tr_Max depends directly on input size configuration scales
#define Tr_Max 1000
#define Size_B1_Max (Tr_Max/100)

// Data arrays mapped for structural error extraction processing blocks
double e_new[Tr_Max], m_new[Tr_Max], dsigma[Tr_Max];
double d_blo[Tr_Max], d_blo_new[Tr_Max], d_final[Tr_Max];

// Global parameters defining global network geometry
extern int L = 4;
extern int Nmed;

// State management definitions matching the Parisi-Rapuano shift register generator 
unsigned int rueda[256];
unsigned int ir1;
unsigned char ind_rand, ig1, ig2, ig3;
extern float Random(void);
extern void ini_ran(int SEMILLA);

int main()
{
    /** Variable Declarations **/
    int Nter, N_steps, N_betas, N_m, M_ter;
    int i, direction;
    int block_size, n_blocks, n_b; // Structural block handling variables for X and Cv calculations
    int *xp, *xm, *yp, *ym;         // Neighbors tracking offset tables (Periodic Boundaries)

    int *s; // Primary 1D vector capturing full 2D spatial spin layouts
    char filename[256]; // Buffer tracking dynamically rendered output logs per beta

    double *e, e_avg, des_e, *M, M_avg, des_M, Cv, des_Cv, X, des_X, e_sq, M_sq, med, sq; // Primary thermodynamics observables
    double beta, b_initial, b_final, b_delta, prob[5]; // Temperature scanning controls
    double *results; // Master 2D summary matrix wrapped as a flat continuous vector

    // Thermalization block control tracking structures
    double *e_b, *m_b, *e_b_new, *m_b_new; 
    double precision; // Maximum allowed epsilon tolerance identifying successful steady state convergence
    double difference_e, difference_m;
    int num_blocks; // Target configuration counting windows matched for comparison
    int s_blocks;   // Size dimensions per comparison window block
    int processed_blocks; 
    int thermalized; // Boolean flag showing if steady state equilibrium bounds met prior to hitting Nter constraints

    /** File Handler Initialization **/
    FILE *f_input = fopen("input.txt", "r"); 
    FILE *f_check = fopen("check.txt", "w"); 
    fclose(f_check);
    FILE *f_out = fopen("Measures.txt", "w"); 
    fclose(f_out);
    FILE *f_out2 = fopen("Cv_X.txt", "w"); 
    fclose(f_out2);

    /** Parse Input Parametric Control Data **/
    if(f_input != NULL) {
        fscanf(f_input, "%lf %lf %lf %d %d", &b_initial, &b_final, &b_delta, &Nmed, &Nter);
        fclose(f_input);
    } else {
        printf("Error: input.txt file missing.\n");
        return -1;
    }

    /** Setup Pseudo Random Number Generator Seed **/
    ini_ran(time(NULL));

    /** Heap Memory Array Allocation **/
    s = ((int*)malloc((L * L) * sizeof(int)));
    xp = ((int*)malloc(L * sizeof(int)));
    xm = ((int*)malloc(L * sizeof(int)));
    yp = ((int*)malloc(L * sizeof(int)));
    ym = ((int*)malloc(L * sizeof(int)));
    e = ((double*)malloc(Nmed * sizeof(double)));
    M = ((double*)malloc(Nmed * sizeof(double)));

    /** Generate Startup Seed Configurations **/
    initialConfiguration(0, s); 

    /** Initialize Periodic Boundary Condition Mapping Offsets (Neighbor indexing) **/
    for (i = 0; i < (L - 1); i++)
    {
        xp[i] = 1;
        yp[i] = L;
    }
    xp[L - 1] = -(L - 1);
    yp[L - 1] = -L * (L - 1);

    for (i = 1; i < L; i++)
    {
        xm[i] = -1;
        ym[i] = -L;
    }
    xm[0] = (L - 1);
    ym[0] = L * (L - 1);

    /** Set Scanning Steps Count and Initial Temperature Coordinates **/
    N_steps = (b_final - b_initial) / b_delta;
    beta = b_initial;

    results = (double*)malloc(7 * N_steps * sizeof(double)); 

    /** Hysteresis Control Loops (Traversing system temperatures in forward and reverse legs) **/
    for(direction = 0; direction < 2; direction++)
    {
        for(N_betas = 0; N_betas < N_steps; N_betas++)
        {
            probability(prob, beta);

            /// Thermalization Phase: Evaluate initial transient steps to decay historical memory artifacts
            f_check = fopen("check.txt", "a");
            fprintf(f_check, "\n-------------\n");
            fprintf(f_check, "Beta:%1.2f\n", beta);

            /// For granular tracking: Create separate transient tracking data sets per temperature tier
            // sprintf(filename, "%1.2f.txt", beta);
            // fterma = fopen(filename, "w");

            precision = 0.01;   // Threshold constraint defining convergence equilibrium
            num_blocks = 2;    // Comparative block tracking registers
            s_blocks = 50;     // Number of operations stored in single evaluation arrays

            e_b = (double *)malloc(num_blocks * s_blocks * sizeof(double));
            m_b = (double *)malloc(num_blocks * s_blocks * sizeof(double));
            e_b_new = (double *)malloc(num_blocks * sizeof(double));
            m_b_new = (double *)malloc(num_blocks * sizeof(double));

            /** THERMALIZATION EVALUATION ENGINE **/
            thermalized = 0;
            processed_blocks = 0;
            i = -1;

            for(M_ter = 0; M_ter < Nter; M_ter++)
            {
                Metropolis(s, xp, xm, yp, ym, prob);

                i++; // Build sequential blocks tracking progress history
                e_b[i] = energy(s, xp, yp);
                m_b[i] = magnetization(s);

                if (i == (num_blocks * s_blocks - 1)) // Once comparison window buffers fill -> evaluate delta trends
                {
                    processed_blocks += num_blocks;
                    Build_Blocks(e_b, e_b_new, num_blocks, s_blocks);
                    Build_Blocks(m_b, m_b_new, num_blocks, s_blocks);

                    difference_e = e_b_new[0] - e_b_new[1]; 
                    difference_m = e_b_new[0] - e_b_new[1];

                    if (difference_e < 0) difference_e = -1.0 * difference_e; // Extract absolute variations
                    if (difference_m < 0) difference_m = -1.0 * difference_m;

                    // If fluctuations settle beneath target epsilon thresholds, stabilization criteria met
                    if ((difference_e < precision) && (difference_m < precision))
                    {
                        fprintf(f_check, "Thermalization steps:%d\nBlocks processed:%d\n", M_ter, processed_blocks);
                        thermalized++;
                        break;
                    }
                    else i = -1; // Reset window indices to keep searching if stability checks fail
                }
            }

            if(thermalized == 0)
                fprintf(f_check, "Thermalization steps:%d\nBlocks processed:%d\n", Nter, 0); // System failed to converge within hard step caps
            fclose(f_check);

            // Free transient analysis memory allocation bounds
            free(e_b); free(m_b); free(e_b_new); free(m_b_new);

            Initialize(e, M, &e_avg, &M_avg, &des_e, &des_M, &des_Cv, &des_X, &e_sq, &M_sq); 

            /** OBSERVABLE SYSTEM MEASUREMENTS **/
            for(N_m = 0; N_m < Nmed; N_m++)
            {
                Metropolis(s, xp, xm, yp, ym, prob);

                // Sample raw equilibrium values
                e[N_m] = energy(s, xp, yp);
                M[N_m] = magnetization(s);
            }

            /** COMPUTE STATISTICAL MEANS AND VARIANCE ERROR BOUNDS **/
            avg_deviation(e, &e_avg, &des_e, &e_sq, Nmed);

            for(i = 0; i < Nmed; i++)
            {
                M[i] = fabs(M[i]); // Compute absolute magnetization tracking order trends safely
            }
            avg_deviation(M, &M_avg, &des_M, &M_sq, Nmed);

            // Archive evaluation points inside master data logs
            results[7 * N_betas + 0] = beta;
            results[7 * N_betas + 1] = e_avg;
            results[7 * N_betas + 2] = des_e;
            results[7 * N_betas + 3] = M_avg;
            results[7 * N_betas + 4] = des_M;
            results[7 * N_betas + 5] = e_sq;
            results[7 * N_betas + 6] = M_sq;

            // Commit structured outputs to Measures.txt logs
            f_out = fopen("Measures.txt", "a");
            fprintf(f_out, "%lf %lf %lf %lf %lf\n", beta, e_avg, des_e, M_avg, des_M);
            fclose(f_out);

            // Derive secondary critical indicators: Specific Heat capacity (Cv) and Magnetic Susceptibility (X) via fluctuation-dissipation relations
            Cv = 2.0 * L * L * (results[7 * N_betas + 5] - results[7 * N_betas + 1] * results[7 * N_betas + 1]);
            X  = L * L * (results[7 * N_betas + 6] - results[7 * N_betas + 3] * results[7 * N_betas + 3]);

            // Data Transformation mapping higher moments for parsing block standard deviations
            for (i = 0; i < Nmed; i++)
            {
                e_new[i] = e[i] * e[i];
                m_new[i] = M[i] * M[i];
            }

            // Execute Blocking Analysis for Energy Fluctuation variance (Cv Error estimation)
            for(block_size = 2; block_size < Size_B1_Max; block_size++)
            {
                n_blocks = (int)(Nmed / block_size);

                Build_Blocks(e, d_blo, n_blocks, block_size);
                Build_Blocks(e_new, d_blo_new, n_blocks, block_size);

                for(n_b = 0; n_b < n_blocks; n_b++)
                    d_final[n_b] = 2.0 * L * L * (d_blo_new[n_b] - d_blo[n_b] * d_blo[n_b]);

                avg_deviation(d_final, &med, &sq, &des_Cv, n_blocks);
            }

            // Execute Blocking Analysis for Magnetization Fluctuation variance (X Error estimation)
            for(block_size = 2; block_size < Size_B1_Max; block_size++)
            {
                n_blocks = (int)(Nmed / block_size);

                Build_Blocks(M, d_blo, n_blocks, block_size);
                Build_Blocks(m_new, d_blo_new, n_blocks, block_size);

                for(n_b = 0; n_b < n_blocks; n_b++)
                    d_final[n_b] = L * L * (d_blo_new[n_b] - d_blo[n_b] * d_blo[n_b]);

                avg_deviation(d_final, &med, &sq, &des_X, n_blocks);
            }

            // Append specific heat and magnetic susceptibility measurements to file
            f_out2 = fopen("Cv_X.txt", "a");
            fprintf(f_out2, "%lf %lf %lf %lf %lf\n", beta, Cv, des_Cv, X, des_X);
            fclose(f_out2);

            beta += b_delta; // Increment or decrement beta tracking index bounds
        }

        b_delta = -b_delta; // Reverse scan tracking delta step paths to trace structural hysteresis pathways
    }

    /** Clean Up Heap Workspace Array Registers **/
    free(s);
    free(xp);
    free(yp);
    free(xm);
    free(ym);
    free(results);
    free(e);
    free(M);

    printf("Execution Finished Successfully.\n");
    return 0;
}

/** Pseudo Random Shift-Register Generator via the Parisi-Rapuano Method **/
float Random (void)
{
    float r;

    ig1 = ind_rand - 24;
    ig2 = ind_rand - 55;
    ig3 = ind_rand - 61;

    rueda[ind_rand] = rueda[ig1] + rueda[ig2];
    ir1 = (rueda[ind_rand] ^ rueda[ig3]);

    ind_rand++;
    r = ir1 * NormRANu;

    return r;
}

/** Seed loader mapping configuration parameters setup initializes random register wheels **/
void ini_ran (int SEMILLA)
{
    int INI, FACTOR, SUM, i;

    srand(SEMILLA);

    INI = SEMILLA;
    FACTOR = 67397;
    SUM = 7364893;

    for (i = 0; i < 256; i++)
    {
        INI = (INI * FACTOR + SUM);
        rueda[i] = INI;
    }
    ind_rand = ig1 = ig2 = ig3 = 0;
}