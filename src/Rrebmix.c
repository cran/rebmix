#include <stdio.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "rngmixf.h"
#include "rebmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

/* Runs RNGMIX template file stream in R. */

void RRNGMIX(char **file, /* File stream. */
             int  *Error) /* Error code. */
{
    *Error = RunRNGMIXTemplateFile(file[0]);
} /* RRNGMIX */

/* Runs REBMIX template file stream in R. */

void RREBMIX(char **file, /* File stream. */
             int  *Error) /* Error code. */
{
    *Error = RunREBMIXTemplateFile(file[0]);
} /* RREBMIX */

/* Returns k-nearest neighbour empirical densities in R. */

void RdensKNearestNeighbourXY(int    *n,     /* Total number of independent observations. */
                              double *x,     /* Pointer to the input array x. */
                              double *y,     /* Pointer to the input array y. */
                              double *p,     /* Pointer to the output array p. */
                              int    *k,     /* k-nearest neighbours. */
                              double *RMIN,  /* Minimum radius of the hypersphere. */
                              double *hx,    /* Normalizing vector. */
                              double *hy,    /* Normalizing vector. */
                              int    *Error) /* Error code. */
{
    double *Dk = NULL;
    double Dc, R, C;
    int    i, j, K, l, m;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (double*)malloc(K * sizeof(double));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * Pi * (*hx) * (*hy));

    for (i = 0; i < *n; i++) {
        for (j = 0; j < K; j++) Dk[j] = DBL_MAX;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (x[i] - x[j]) / (*hx); Dc = R * R;
            R = (y[i] - y[j]) / (*hy); Dc += R * R;

            for (l = 0; l < K; l++) {
                if (Dc <= Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    Dk[l] = Dc;

                    break;
                }
            }
        }

        if (Dk[K - 1] > DBL_MIN) {
            R = sqrt(Dk[K - 1]); if (R < *RMIN) R = *RMIN; 
        }
        else 
            R = *RMIN;

        p[i] = C / (R * R);
    }

E0: if (Dk) free(Dk);
} /* RdensKNearestNeighbourXY */

/* Returns Parzen window empirical densities in R. */

void RdensParzenWindowXY(int    *n,     /* Total number of independent observations. */
                         double *x,     /* Pointer to the input array x. */
                         double *y,     /* Pointer to the input array y. */
                         double *p,     /* Pointer to the output array p. */
                         double *hx,    /* Side of the hypersquare. */
                         double *hy,    /* Side of the hypersquare. */
                         int    *Error) /* Error code. */
{
    int    i, j;
    double C, rx, ry;

    *Error = *n < 1; if (*Error) return;

    C = 1.0 / (*hx) / (*hy) / (*n); rx = 0.5 * (*hx); ry = 0.5 * (*hy);
    
    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if ((fabs(x[j] - x[i]) <= rx) && (fabs(y[j] - y[i]) <= ry)) {
                p[i] += C; p[j] += C;
            }
        }
    }
} /* RdensParzenWindowXY */

/* Returns histogram empirical densities in R. */

void RdensHistogramXY(int    *k,     /* Total number of bins. */ 
                      int    *n,     /* Total number of independent observations. */
                      double *x,     /* Pointer to the input array x. */
                      double *y,     /* Pointer to the input array y. */
                      double *p,     /* Pointer to the output array p. */
                      double *x0,    /* Origin. */ 
                      double *y0,    /* Origin. */ 
                      double *hx,    /* Side of the hypersquare. */
                      double *hy,    /* Side of the hypersquare. */
                      int    *Error) /* Error code. */
{
    int    i, j, m = *k - 1;
    double C, rx, ry;

    *Error = *n < 1; if (*Error) return;
    
    C = 1.0 / (*hx) / (*hy) / (*n); rx = 0.5 * (*hx); ry = 0.5 * (*hy);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + 0.5); 
        
        if (j < 0) j = 0; else if (j > m) j = m;

        x[*k] = (*x0) + j * (*hx);

        j = (int)floor((y[i] - (*y0)) / (*hy) + 0.5); 
        
        if (j < 0) j = 0; else if (j > m) j = m;

        y[*k] = (*y0) + j * (*hy);

        for (j = 0; j < *k; j++) {
            if ((fabs(x[j] - x[*k]) > rx) || (fabs(y[j] - y[*k]) > ry)) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} /* RdensHistogramXY */

/* Returns k-nearest neighbour empirical densities in R. */

void RdensKNearestNeighbourX(int    *n,     /* Total number of independent observations. */
                             double *x,     /* Pointer to the input array x. */
                             double *p,     /* Pointer to the output array p. */
                             int    *k,     /* k-nearest neighbours. */
                             double *RMIN,  /* Minimum radius of the hypersphere. */
                             double *hx,    /* Normalizing vector. */
                             int    *Error) /* Error code. */
{
    double *Dk = NULL;
    double Dc, R, C;
    int    i, j, K, l, m;

    *Error = *n < 1; if (*Error) return;

    K = *k; if (K > 1) K -= 1; else K = 1; 
    
    Dk = (double*)malloc(K * sizeof(double));

    *Error = NULL == Dk; if (*Error) goto E0;

    C = (*k) / ((*n) * 2.0 * (*hx));

    for (i = 0; i < *n; i++) {
        for (j = 0; j < K; j++) Dk[j] = DBL_MAX;

        for (j = 0; j < *n; j++) if (i != j) {
            R = (x[i] - x[j]) / (*hx); Dc = R * R;

            for (l = 0; l < K; l++) {
                if (Dc <= Dk[l]) {
                    for (m = K - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    Dk[l] = Dc;

                    break;
                }
            }
        }

        if (Dk[K - 1] > DBL_MIN) {
            R = sqrt(Dk[K - 1]); if (R < *RMIN) R = *RMIN; 
        }
        else 
            R = *RMIN;

        p[i] = C / R;
    }

E0: if (Dk) free(Dk);
} /* RdensKNearestNeighbourX */

/* Returns Parzen window empirical densities in R. */

void RdensParzenWindowX(int    *n,     /* Total number of independent observations. */
                        double *x,     /* Pointer to the input array x. */
                        double *p,     /* Pointer to the output array p. */
                        double *hx,    /* Side of the hypersquare. */
                        int    *Error) /* Error code. */
{
    int    i, j;
    double C, rx;

    *Error = *n < 1; if (*Error) return;

    C = 1.0 / (*hx) / (*n); rx = 0.5 * (*hx);
    
    for (i = 0; i < *n; i++) {
        for (j = i; j < *n; j++) {
            if (fabs(x[j] - x[i]) <= rx) {
                p[i] += C; p[j] += C;
            }
        }
    }
} /* RdensParzenWindowX */

/* Returns histogram empirical densities in R. */

void RdensHistogramX(int    *k,     /* Total number of bins. */ 
                     int    *n,     /* Total number of independent observations. */
                     double *x,     /* Pointer to the input array x. */
                     double *p,     /* Pointer to the output array p. */
                     double *x0,    /* Origin. */ 
                     double *hx,    /* Side of the hypersquare. */
                     int    *Error) /* Error code. */
{
    int    i, j, m = *k - 1;
    double C, rx;

    *Error = *n < 1; if (*Error) return;
    
    C = 1.0 / (*hx) / (*n); rx = 0.5 * (*hx);

    *k = 0;

    for (i = 0; i < *n; i++) {
        j = (int)floor((x[i] - (*x0)) / (*hx) + 0.5); 
        
        if (j < 0) j = 0; else if (j > m) j = m;

        x[*k] = (*x0) + j * (*hx);

        for (j = 0; j < *k; j++) {
            if (fabs(x[j] - x[*k]) > rx) goto S0;
                
            p[j] += C; goto S1;
S0:;    }

        p[*k] = C; (*k)++; 
S1:;}
} /* RdensHistogramX */
