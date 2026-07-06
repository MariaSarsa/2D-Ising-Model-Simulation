#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

void initialConfiguration (int flag, int *s);
void probability(double *prob, double beta);
void Initialize(double *e, double *M, double *e_avg, double *M_avg, double *des_e, double *des_M, double *des_Cc, double *des_X, double *e_sq, double *M_sq);
void Metropolis(int *s, int *xp, int *xm, int *yp, int *ym, double *prob);
double energy (int *s, int *xp, int *yp);
double magnetization (int *s);
void avg_deviation (double *s, double *average, double *sigma, double *sq, int N);
void Build_Blocks(double *x, double *xb, int n_b, int b_s);

#endif // HEADER_H_INCLUDED
