#include <math.h>
#include <float.h> 
#include <stdio.h>
#include <ctype.h>
#include <time.h>

#if (_REBMIXDLL)
#include <windows.h>
#endif

#include "rebmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

/* Returns the value log(Gamma(y)) for y > 0. See http://www.nrbook.com/a/bookcpdf/c6-1.pdf */

FLOAT Gammaln(FLOAT y)
{
    FLOAT x, z, Tmp, Ser;
    int   j;

    static FLOAT Cof[6] = {(FLOAT)76.18009172947146, -(FLOAT)86.50532032941677,
                           (FLOAT)24.01409824083091, -(FLOAT)1.231739572450155,
                           (FLOAT)0.1208650973866179E-2, -(FLOAT)0.5395239384953E-5};

    static FLOAT Stp = (FLOAT)2.5066282746310005;

    z = x = y; Tmp = x + (FLOAT)5.5; Tmp -= (x + (FLOAT)0.5) * (FLOAT)log(Tmp);

    Ser = (FLOAT)1.000000000190015;

    for (j = 0; j < 6; j++) Ser += Cof[j] / ++z;

    return (-Tmp + (FLOAT)log(Stp * Ser / x));
} /* Gammaln */

/* Returns the incomplete gamma function P(a, y) evaluated by its series
   representation as GamSer. Also returns log(Gamma(a)) as Gln. */

int GammaSer(FLOAT a,       /* Constant a > 0. */
             FLOAT y,       /* Variable y > 0. */
             FLOAT *GamSer, /* Incomplete gamma function. */
             FLOAT *Gamln)  /* Log(Gamma(a)). */
{
    int   i;
    FLOAT Sum, Del, ap;
    int   Error = 0;

    *Gamln = Gammaln(a);

    if (y <= FLOAT_MIN) {
        if (y < (FLOAT)0.0) Error = 1; if (Error) goto E0;

        *GamSer = (FLOAT)0.0;
    }
    else {
        ap = a; Sum = (FLOAT)1.0 / a; Del = Sum; Error = 1; i = 1;

        while ((i <= ItMax) && Error) {
            ap += (FLOAT)1.0; Del *= y / ap; Sum += Del;

            if ((FLOAT)fabs(Del) < (FLOAT)fabs(Sum) * Eps) Error = 0;

            i++;
        }

        if (Error) goto E0;

        *GamSer = Sum * (FLOAT)exp(-y + a * log(y) - *Gamln);
    }

E0: return (Error);
} /* GammaSer */

/* Returns the incomplete gamma function Q(a, y) evaluated by its continued
   fraction representation as GamCfg. Also returns log(Gamma(a)) as Gln. */

int GammaCfg(FLOAT a,       /* Constant a > 0. */
             FLOAT y,       /* Variable y > 0. */
             FLOAT *GamCfg, /* Incomplete gamma function. */
             FLOAT *Gamln)  /* Log(Gamma(a)). */
{
    int   i;
    FLOAT Gold, G, Fac, a0, a1, b0, b1, aif, aia, ai;
    int   Error = 0;

    *Gamln = Gammaln(a); 

    if (y <= FLOAT_MIN) {
        if (y < (FLOAT)0.0) Error = 1; if (Error) goto E0;

        *GamCfg = (FLOAT)0.0;
    }
    else {
        G = (FLOAT)0.0; Gold = (FLOAT)0.0; Fac = (FLOAT)1.0;

        a0 = (FLOAT)1.0; a1 = y; b0 = (FLOAT)0.0; b1 = (FLOAT)1.0; Error = 1; i = 1;

        while ((i <= ItMax) && Error) {
            ai = (FLOAT)1.0 * i; aia = ai - a; aif = ai * Fac;

            a0 = (a1 + a0 * aia) * Fac;
            b0 = (b1 + b0 * aia) * Fac;

            a1 = y * a0 + aif * a1;
            b1 = y * b0 + aif * b1;

            if (a1 != (FLOAT)0.0) {
                Fac = (FLOAT)1.0 / a1; G = b1 * Fac;

                if ((FLOAT)fabs((G - Gold) / G) < Eps) Error = 0; else Gold = G;
            }

            i++;
        }

        if (Error) goto E0;

        *GamCfg = (FLOAT)exp(-y + a * log(y) - *Gamln) * G;
    }

E0: return (Error);
} /* GammaCfg */

/* Returns the incomplete gamma function P(a, y). See http://www.nrbook.com/a/bookcpdf/c6-2.pdf */

int GammaP(FLOAT a,     /* Constant a > 0. */
           FLOAT y,     /* Variable y > 0. */
           FLOAT *GamP) /* Incomplete gamma function. */
{
    FLOAT GamSer, GamCfg, Gamln;
    int   Error = 0;

    if ((y < (FLOAT)0.0) || (a <= FLOAT_MIN)) {
        *GamP = (FLOAT)0.0; Error = 1; if (Error) goto E0;
    }
    else
    if (y < a + (FLOAT)1.0) {
        Error = GammaSer(a, y, &GamSer, &Gamln); 

        if (Error) goto E0;
        
        *GamP = GamSer;
    }
    else {
        Error = GammaCfg(a, y, &GamCfg, &Gamln); 

        if (Error) goto E0;
        
        *GamP = (FLOAT)1.0 - GamCfg;
    }

E0: return (Error);
} /* GammaP */

/* Returns the error function erf(y). */

int ErrorF(FLOAT y,     /* Variable y. */
           FLOAT *ErF)  /* Error function. */
{
    FLOAT GamP;
    int   Error = 0;

    Error = GammaP((FLOAT)0.5, y * y, &GamP);

    if (Error) goto E0;

    if (y < (FLOAT)0.0)
        *ErF = -GamP;
    else
        *ErF = +GamP;

E0: return (Error);
} /* ErrorF */

/* Returns component p.d.f or c.d.f. */ 

int ComponentDist(int                      d,            /* Number of independent random variables. */
                  FLOAT                    *Y,           /* Pointer to the input point [y0,...,yd-1]. */
                  MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                  FLOAT                    *CmpDist,     /* Component distribution. */
                  int                      Cumulative)   /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT y, ypb, ErF, Sum, p;
    int   i, j, k, n;
    int   Error = 0;

    *CmpDist = (FLOAT)1.0;

    if (Cumulative) {
        for (i = 0; i < d; i++) {
            switch (MrgDistType[i].ParametricFamily) {
            case pfNormal:
                y = (Y[i] - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

                Error = ErrorF(y, &ErF);

                if (Error) goto E0;

                *CmpDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);

                break;
            case pfLognormal:
                if (Y[i] > FLOAT_MIN) {
                    y = ((FLOAT)log(Y[i]) - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

                    Error = ErrorF(y, &ErF);

                    if (Error) goto E0;

                    *CmpDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);
                }
                else {
                    Error = 1; if (Error) goto E0;
                }
                
                break;
            case pfWeibull:
                if (Y[i] > FLOAT_MIN) {
                    ypb = (FLOAT)exp(MrgDistType[i].Parameter1 * log(Y[i] / MrgDistType[i].Parameter0));

                    *CmpDist *= ((FLOAT)1.0 - (FLOAT)exp(-ypb));
                }
                else {
                    Error = 1; if (Error) goto E0; 
                }

                break;
            case pfBinomial:
                k = (int)Y[i]; n = (int)MrgDistType[i].Parameter0; p = MrgDistType[i].Parameter1;

                if (k < 0)
                    *CmpDist *= (FLOAT)0.0;
                else
                if (k == 0)
                    *CmpDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
                else
                if (k >= n)
                    *CmpDist *= (FLOAT)1.0;
                else 
                if (p == (FLOAT)0.0)
                    *CmpDist *= (FLOAT)1.0;
                else
                if (p == (FLOAT)1.0)
                    *CmpDist *= (FLOAT)0.0;
                else {
                    Sum = (FLOAT)0.0;

                    for (j = 0; j <= k; j++) {
                        Sum += (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(j + (FLOAT)1.0) - Gammaln(n - j + (FLOAT)1.0) +
                               j * (FLOAT)log(p) + (n - j) * (FLOAT)log((FLOAT)1.0 - p));
                    }

                    *CmpDist *= Sum;
                }
            }
        }
    }
    else {
        for (i = 0; i < d; i++) {
            switch (MrgDistType[i].ParametricFamily) {
            case pfNormal:
                y = (Y[i] - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

                *CmpDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Parameter1);

                break;
            case pfLognormal:
                if (Y[i] > FLOAT_MIN) {
                    y = ((FLOAT)log(Y[i]) - MrgDistType[i].Parameter0)/(Sqrt2 * MrgDistType[i].Parameter1);

                    *CmpDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Parameter1) / Y[i];
                }
                else {
                    Error = 1; if (Error) goto E0;
                }

                break;
            case pfWeibull:
                if (Y[i] > FLOAT_MIN) {
                    ypb = (FLOAT)exp(MrgDistType[i].Parameter1 * log(Y[i] /  MrgDistType[i].Parameter0));

                    *CmpDist *= MrgDistType[i].Parameter1 * ypb * (FLOAT)exp(-ypb) / Y[i];
                }
                else {
                    Error = 1; if (Error) goto E0;
                }

                break;
            case pfBinomial:
                k = (int)Y[i]; n = (int)MrgDistType[i].Parameter0; p = MrgDistType[i].Parameter1;

                if (k < 0)
                    *CmpDist *= (FLOAT)0.0;
                else
                if (k == 0)
                    *CmpDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
                else
                if (k == n)
                    *CmpDist *= (FLOAT)pow(p, n);
                else
                if (k > n)
                    *CmpDist *= (FLOAT)0.0;
                else 
                if ((p == (FLOAT)0.0) || (p == (FLOAT)1.0))
                    *CmpDist *= (FLOAT)0.0;
                else
                    *CmpDist *= (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0) +
                                k * (FLOAT)log(p) + (n - k) * (FLOAT)log((FLOAT)1.0 - p));
            }
        }
    }

E0: return (Error);
} /* ComponentDist */

/* Returns component marginal p.d.f or c.d.f. */ 

int ComponentMarginalDist(int                      i,            /* Index of variable y. */
                          FLOAT                    *Y,           /* Pointer to the input point [y0,...,yd-1]. */
                          MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                          FLOAT                    *CmpMrgDist,  /* Component marginal distribution. */
                          int                      Cumulative)   /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT y, ypb, Sum, p;
    int   j, k, n;
    FLOAT ErF;
    int   Error = 0;
    
    *CmpMrgDist = (FLOAT)1.0;

    if (Cumulative) {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            y = (Y[i] - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

            Error = ErrorF(y, &ErF);

            if (Error) goto E0;

            *CmpMrgDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

                Error = ErrorF(y, &ErF);

                if (Error) goto E0;

                *CmpMrgDist *= (FLOAT)0.5 * ((FLOAT)1.0 + ErF);
            }
            else {
                Error = 1; if (Error) goto E0;
            }
                
            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                ypb = (FLOAT)exp(MrgDistType[i].Parameter1 * log(Y[i] / MrgDistType[i].Parameter0));

                *CmpMrgDist *= ((FLOAT)1.0 - (FLOAT)exp(-ypb));
            }
            else {
                Error = 1; if (Error) goto E0; 
            }

            break;
        case pfBinomial:
            k = (int)Y[i]; n = (int)MrgDistType[i].Parameter0; p = MrgDistType[i].Parameter1;

            if (k < 0)
                *CmpMrgDist *= (FLOAT)0.0;
            else
            if (k == 0)
                *CmpMrgDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
            else
            if (k >= n)
                *CmpMrgDist *= (FLOAT)1.0;
            else
            if (p == (FLOAT)0.0)
                *CmpMrgDist *= (FLOAT)1.0;
            else
            if (p == (FLOAT)1.0)
                *CmpMrgDist *= (FLOAT)0.0;
            else {
                Sum = (FLOAT)0.0;

                for (j = 0; j <= k; j++) {
                    Sum += (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(j + (FLOAT)1.0) - Gammaln(n - j + (FLOAT)1.0) +
                           j * (FLOAT)log(p) + (n - j) * (FLOAT)log((FLOAT)1.0 - p));
                }

                *CmpMrgDist *= Sum;
            }
        }
    }
    else {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            y = (Y[i] - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

            *CmpMrgDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Parameter1);

            break;
        case pfLognormal:
            if (Y[i] > FLOAT_MIN) {
                y = ((FLOAT)log(Y[i]) - MrgDistType[i].Parameter0) / (Sqrt2 * MrgDistType[i].Parameter1);

                *CmpMrgDist *= (FLOAT)exp(-(y * y)) / (Sqrt2Pi * MrgDistType[i].Parameter1) / Y[i];
            }
            else {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfWeibull:
            if (Y[i] > FLOAT_MIN) {
                ypb = (FLOAT)exp(MrgDistType[i].Parameter1 * log(Y[i] /  MrgDistType[i].Parameter0));

                *CmpMrgDist *= MrgDistType[i].Parameter1 * ypb * (FLOAT)exp(-ypb) / Y[i];
            }
            else {
                Error = 1; if (Error) goto E0;
            }

            break;
        case pfBinomial:
            k = (int)Y[i]; n = (int)MrgDistType[i].Parameter0; p = MrgDistType[i].Parameter1;

            if (k < 0)
                *CmpMrgDist *= (FLOAT)0.0;
            else
            if (k == 0)
                *CmpMrgDist *= (FLOAT)pow((FLOAT)1.0 - p, n);
            else
            if (k == n)
                *CmpMrgDist *= (FLOAT)pow(p, n);
            else
            if (k > n)
                *CmpMrgDist *= (FLOAT)0.0;
            else
            if ((p == (FLOAT)0.0) || (p == (FLOAT)1.0))
                *CmpMrgDist *= (FLOAT)0.0;
            else
                *CmpMrgDist *= (FLOAT)exp(Gammaln(n + (FLOAT)1.0) - Gammaln(k + (FLOAT)1.0) - Gammaln(n - k + (FLOAT)1.0) +
                               k * (FLOAT)log(p) + (n - k) * (FLOAT)log((FLOAT)1.0 - p));
        }
    }

E0: return (Error);
} /* ComponentMarginalDist */

/* Returns mixture p.d.f or c.d.f. */ 

int MixtureDist(int                      d,             /* Number of independent random variables. */
                FLOAT                    *Y,            /* Pointer to the input point [y0,...,yd-1]. */
                int                      c,             /* Number of components. */ 
                FLOAT                    *W,            /* Component weights. */
                MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                FLOAT                    *MixDist,      /* Mixture distribution. */
                int                      Cumulative)    /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT CmpDist;
    int   i;
    int   Error = 0;

    *MixDist = (FLOAT)0.0;

    for (i = 0; i < c; i++) {
        Error = ComponentDist(d, Y, MrgDistType[i], &CmpDist, Cumulative);

        if (Error) goto E0;

        *MixDist += W[i] * CmpDist;
    }

E0: return (Error);
} /* MixtureDist */

/* Returns mixture marginal p.d.f or c.d.f. */ 

int MixtureMarginalDist(int                      i,             /* Index of variable y. */  
                        FLOAT                    *Y,            /* Pointer to the input point [y0,...,yd-1]. */
                        int                      c,             /* Number of components. */ 
                        FLOAT                    *W,            /* Component weights. */
                        MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                        FLOAT                    *MixMrgDist,   /* Mixture marginal distribution. */
                        int                      Cumulative)    /* Set 1 if c.d.f. or 0 if p.d.f. */
{
    FLOAT CmpMrgDist;
    int   j;
    int   Error = 0;

    *MixMrgDist = (FLOAT)0.0;

    for (j = 0; j < c; j++) {
        Error = ComponentMarginalDist(i, Y, MrgDistType[j], &CmpMrgDist, Cumulative);

        if (Error) goto E0;

        *MixMrgDist += W[j] * CmpMrgDist;
    }

E0: return (Error);
} /* MixtureMarginalDist */

/* Returns information criterion for k-nearest neighbour. */ 

int InformationCriterionKNN(InformationCriterionType_e ICType,        /* Information criterion type. */
                            int                        n,             /* Total number of independent observations. */
                            int                        d,             /* Number of independent random variables. */ 
                            FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                            int                        c,             /* Number of components. */ 
                            FLOAT                      *W,            /* Component weights. */
                            MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                            FLOAT                      *IC,           /* Information criterion. */
                            FLOAT                      *logL)         /* log-likelihood. */
{
    int   i, j, M;
    FLOAT EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    M = 2 * c * d + c - 1;

    *IC = *logL = EN = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist < FLOAT_MIN)
            *logL += (FLOAT)log(FLOAT_MIN);
        else
            *logL += (FLOAT)log(MixDist);

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 
                
                if (tau > FLOAT_MIN) {
                    EN -= tau * (FLOAT)log(tau); PC += tau * tau; 
                }
            }

            break;
        default:
            EN = PC = 0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M;

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * M;

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * M;

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * ((FLOAT)1.0 + (M + 1) / (n - M - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + M * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + M * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * M * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + (M - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + M * (FLOAT)log((FLOAT)n);
    }

E0: return (Error);
} /* InformationCriterionKNN */

/* Returns information criterion for Parzen window. */ 

int InformationCriterionPW(InformationCriterionType_e ICType,        /* Information criterion type. */
                           int                        n,             /* Total number of independent observations. */
                           int                        d,             /* Number of independent random variables. */ 
                           FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                           int                        c,             /* Number of components. */ 
                           FLOAT                      *W,            /* Component weights. */
                           MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                           FLOAT                      *IC,           /* Information criterion. */
                           FLOAT                      *logL)         /* log-likelihood. */
{
    int   i, j, M;
    FLOAT EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    M = 2 * c * d + c - 1;

    *IC = *logL = EN = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < n; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist < FLOAT_MIN)
            *logL += (FLOAT)log(FLOAT_MIN);
        else
            *logL += (FLOAT)log(MixDist);

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 
                
                if (tau > FLOAT_MIN) {
                    EN -= tau * (FLOAT)log(tau); PC += tau * tau; 
                }
            }

            break;
        default:
            EN = PC = 0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M;

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * M;

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * M;

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * ((FLOAT)1.0 + (M + 1) / (n - M - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + M * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + M * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * M * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + (M - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + M * (FLOAT)log((FLOAT)n);
    }

E0: return (Error);
} /* InformationCriterionPW */

/* Returns information criterion for histogram. */ 

int InformationCriterionH(InformationCriterionType_e ICType,        /* Information criterion type. */
                          int                        k,             /* Total number of bins. */
                          int                        n,             /* Total number of independent observations. */
                          int                        d,             /* Number of independent random variables. */ 
                          FLOAT                      **Y,           /* Pointer to the input points [y0,...,yd-1,kl]. */
                          int                        c,             /* Number of components. */ 
                          FLOAT                      *W,            /* Component weights. */
                          MarginalDistributionType   **MrgDistType, /* Marginal distribution type. */
                          FLOAT                      *IC,           /* Information criterion. */
                          FLOAT                      *logL)         /* log-likelihood. */
{
    int   i, j, M;
    FLOAT EN, PW, K, PC, CmpDist, MixDist, tau;
    int   Error = 0;

    M = 2 * c * d + c - 1;

    *IC = *logL = EN = PW = K = PC = (FLOAT)0.0;

    for (i = 0; i < k; i++) {
        Error = MixtureDist(d, Y[i], c, W, MrgDistType, &MixDist, 0);

        if (Error) goto E0;

        if (MixDist < FLOAT_MIN)
            *logL += Y[i][d] * (FLOAT)log(FLOAT_MIN);
        else
            *logL += Y[i][d] * (FLOAT)log(MixDist);

        switch (ICType) {
        case icAWE: case icCLC: case icICL: case icICLBIC:
            for (j = 0; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                tau = W[j] * CmpDist / MixDist; 

                if (tau > FLOAT_MIN) {
                    EN -= Y[i][d] * tau * (FLOAT)log(tau); PC += Y[i][d] * tau * tau; 
                }
            }

            break;
        default:
            EN = PC = 0.0;
        }
    }

    switch (ICType) {
    case icAIC: /* AIC - Akaike information criterion Akaike (1973). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M;

        break;
    case icAIC3: /* AIC3 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)3.0 * M;

        break;
    case icAIC4: /* AIC4 - Modified Akaike information criterion Smith & Spiegelhalter (1980). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)4.0 * M;

        break;
    case icAICc: /* AICc - Akaike second-order corrected information criterion for small sample sizes Hurvich & Tsai (1989). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * ((FLOAT)1.0 + (M + 1) / (n - M - (FLOAT)1.0));

        break;
    case icBIC: /* BIC - Bayesian information criterion Schwarz (1978). */
        *IC = -(FLOAT)2.0 * (*logL) + M * (FLOAT)log((FLOAT)n);

        break;
    case icCAIC: /* CAIC - Consistent Akaike information criterion Bozdogan (1987). */
        *IC = -(FLOAT)2.0 * (*logL) + M * ((FLOAT)log((FLOAT)n) + (FLOAT)1.0);

        break;
    case icHQC: /* HQC - Hannan-Quinn information criterion Hannan & Quinn (1979). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log(log((FLOAT)n));

        break;
    case icMDL2: /* MDL2 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icMDL5: /* MDL5 - Minimum description length Liang et al. (1992). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)5.0 * M * (FLOAT)log((FLOAT)n);

        break;
    case icAWE: /* AWE - Approximate weight of evidence criterion Banfield & Raftery (1993). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * M * ((FLOAT)3.0 / (FLOAT)2.0 + (FLOAT)log((FLOAT)n));

        break;
    case icCLC: /* CLC - Classification likelihood criterion Biernacki & Govaert (1997). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN;

        break;
    case icICL: /* ICL - Integrated classification likelihood Biernacki et al. (1998). */
        for (j = 0; j < c; j++) {
            PW += W[j] * (FLOAT)log(W[j]); K += (FLOAT)Gammaln(W[j] * n + (FLOAT)0.5);
        }

        K += Gammaln(c * (FLOAT)0.5) - c * Gammaln((FLOAT)0.5) - Gammaln(n + c * (FLOAT)0.5);

        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + (FLOAT)2.0 * n * PW - (FLOAT)2.0 * K + (M - c + (FLOAT)1.0) * (FLOAT)log((FLOAT)n);

        break;
    case icPC: /* PC - Partition coeficient Bezdek (1981). */
        *IC = PC; 

        break;
    case icICLBIC: /* ICL-BIC - Integrated classification likelihood criterion Biernacki et al. (1998). */
        *IC = -(FLOAT)2.0 * (*logL) + (FLOAT)2.0 * EN + M * (FLOAT)log((FLOAT)n);
    }

E0: return (Error);
} /* InformationCriterionH */

/* Preprocessing of observations for k-nearest neighbour. */

int PreprocessingKNN(int   k,    /* k-nearest neighbours. */
                     FLOAT RMIN, /* Minimum radius of the hypersphere. */
                     FLOAT *h,   /* Normalizing vector. */
                     int   n,    /* Total number of independent observations. */
                     int   d,    /* Number of independent random variables. */ 
                     FLOAT **Y)  /* Pointer to the input array [y0,...,yd-1,kl,V,R]. */
{
    FLOAT *Dk = NULL;
    FLOAT Dc, R, V, Vn;
    int   i, j, l, m;
    int   Error = 0;

    if (k > 1) k -= 1; else k = 1; 
    
    Dk = (FLOAT*)malloc(k * sizeof(FLOAT));

    Error = NULL == Dk; if (Error) goto E0;

    Vn = (FLOAT)exp(d * LogPi / (FLOAT)2.0 - Gammaln((FLOAT)1.0 + d / (FLOAT)2.0));

    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) Dk[j] = FLOAT_MAX;

        for (j = 0; j < n; j++) if (i != j) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                R = (Y[i][l] - Y[j][l]) / h[l]; Dc += R * R;
            }

            for (l = 0; l < k; l++) {
                if (Dc <= Dk[l]) {
                    for (m = k - 1; m > l; m--) Dk[m] = Dk[m - 1];

                    Dk[l] = Dc;

                    break;
                }
            }
        }

        if (Dk[k - 1] > FLOAT_MIN) {
            R = (FLOAT)sqrt(Dk[k - 1]); if (R < RMIN) R = RMIN; 
        }
        else 
            R = RMIN;

        V = Vn * (FLOAT)exp(d * log(R));

        for (j = 0; j < d; j++) V *= h[j];

        Y[i][d] = (FLOAT)1.0; Y[i][d + 1] = V; Y[i][d + 2] = R;
    }

E0: if (Dk) free(Dk);

    return (Error);
} /* PreprocessingKNN */

/* Preprocessing of observations for Parzen window. */

int PreprocessingPW(FLOAT *h,   /* Sides of the hypersquare. */
                    int   n,    /* Total number of independent observations. */
                    int   d,    /* Number of independent random variables. */ 
                    FLOAT **Y)  /* Pointer to the input array [y0,...,yd-1,kl,k]. */
{
    int i, j, k;
    int Error = n < 1;

    for (i = 0; i < n; i++) {
        Y[i][d] = (FLOAT)1.0; Y[i][d + 1] = (FLOAT)0.0;

        for (j = i; j < n; j++) {
            for (k = 0; k < d; k++) if ((FLOAT)fabs(Y[i][k] - Y[j][k]) > (FLOAT)0.5 * h[k]) goto S0;

            Y[i][d + 1] += (FLOAT)1.0; Y[j][d + 1] += (FLOAT)1.0;
S0:;    }
    }

    return (Error);
} /* PreprocessingPW */


/* Preprocessing of observations for histogram. */

int PreprocessingH(FLOAT                  *h,          /* Sides of the hypersquare. */
                   FLOAT                  *y0,         /* Origin. */
                   ParametricFamilyType_e *ParFamType, /* Parametric family types. */
                   int                    *k,          /* Total number of bins. */
                   int                    n,           /* Total number of independent observations. */
                   int                    d,           /* Number of independent random variables. */ 
                   FLOAT                  **X,         /* Pointer to the input points [x0,...,xd-1]. */
                   FLOAT                  **Y)         /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int i, j, l, m = *k - 1;
    int Error = n < 1;

    *k = 0;

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            l = (int)floor((X[i][j] - y0[j]) / h[j] + (FLOAT)0.5);

            switch (ParFamType[j]) {
            case pfNormal: case pfLognormal: case pfWeibull:
                if (l < 0) l = 0; else if (l > m) l = m; 

                break;
            case pfBinomial:
                break;
            }

            Y[*k][j] = y0[j] + l * h[j];
        }

        for (j = 0; j < *k; j++) {
            for (l = 0; l < d; l++) if (fabs(Y[j][l] - Y[*k][l]) > (FLOAT)0.5 * h[l]) goto S0;

            Y[j][d] += (FLOAT)1.0; goto S1;
S0:;    }

        Y[*k][d] = (FLOAT)1.0; (*k)++;
S1:;}

    return (Error);
} /* PreprocessingH */

/* Global mode detection for k-nearest neighbour. */

int GlobalModeKNN(FLOAT *h,    /* Normalizing vector. */
                  FLOAT *ymin, /* minimum y. */
                  int   *m,    /* Global mode. */
                  int   n,     /* Total number of independent observations. */
                  int   d,     /* Number of independent random variables. */ 
                  FLOAT **Y)   /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int   i, j, l;
    FLOAT Dc, Ri, Rj;
    int   Error = 0;

    j = 0; 

    for (i = 1; i < n; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] / Y[i][d + 1] > Y[j][d] / Y[j][d + 1]) {
            j = i;
        }
        else
        if (Y[i][d] / Y[i][d + 1] == Y[j][d] / Y[j][d + 1]) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Ri = (Y[i][l] - ymin[l]) / h[l]; Dc += Ri * Ri;
            }

            Ri = Dc; Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Rj = (Y[j][l] - ymin[l]) / h[l]; Dc += Rj * Rj;
            }

            Rj = Dc;

            if (Ri < Rj) j = i;
        } 
    }

    *m = j;

    return (Error);
} /* GlobalModeKNN */

/* Global mode detection for Parzen window. */

int GlobalModePW(FLOAT *h,    /* Sides of the hypersquare. */
                 FLOAT *ymin, /* minimum y. */
                 int   *m,    /* Global mode. */
                 int   n,     /* Total number of independent observations. */
                 int   d,     /* Number of independent random variables. */ 
                 FLOAT **Y)   /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int   i, j, l;
    FLOAT Dc, Ri, Rj;
    int   Error = 0;

    j = 0; 

    for (i = 1; i < n; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] * Y[i][d + 1] > Y[j][d] * Y[j][d + 1]) {
            j = i;
        }
        else
        if (Y[i][d] * Y[i][d + 1] == Y[j][d] * Y[j][d + 1]) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Ri = (Y[i][l] - ymin[l]) / h[l]; Dc += Ri * Ri;
            }

            Ri = Dc; Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Rj = (Y[j][l] - ymin[l]) / h[l]; Dc += Rj * Rj;
            }

            Rj = Dc;

            if (Ri < Rj) j = i;
        } 
    }

    *m = j;

    return (Error);
} /* GlobalModePW */

/* Global mode detection for histogram. */

int GlobalModeH(FLOAT *h,    /* Sides of the hypersquare. */
                FLOAT *ymin, /* Minimum y. */
                int   *m,    /* Global mode. */
                int   k,     /* Total number of bins. */
                int   d,     /* Number of independent random variables. */ 
                FLOAT **Y)   /* Pointer to the input array [y0,...,yd-1,kl]. */
{
    int   i, j, l;
    FLOAT Dc, Ri, Rj;
    int   Error = 0;

    j = 0; 

    for (i = 1; i < k; i++) if (Y[i][d] > FLOAT_MIN) {
        if (Y[i][d] > Y[j][d]) {
            j = i;
        }
        else
        if (Y[i][d] == Y[j][d]) {
            Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Ri = (Y[i][l] - ymin[l]) / h[l]; Dc += Ri * Ri;
            }

            Ri = Dc; Dc = (FLOAT)0.0;

            for (l = 0; l < d; l++) {
                Rj = (Y[j][l] - ymin[l]) / h[l]; Dc += Rj * Rj;
            }

            Rj = Dc;

            if (Ri < Rj) j = i;
        } 
    }

    *m = j;

    return (Error);
} /* GlobalModeH */

/* Returns rough normal parameters. */

int RoughNormalParameters(FLOAT ym,   
                          FLOAT fm,
                          FLOAT *Mean,
                          FLOAT *Stdev)
{
    int Error = 0;

    *Mean = ym; *Stdev = (FLOAT)1.0 / (Sqrt2Pi * fm);

    return (Error);
} /* RoughNormalParameters */

/* Returns rough lognormal parameters. */

int RoughLognormalParameters(FLOAT ym,   
                             FLOAT fm,
                             FLOAT *Mean,
                             FLOAT *Stdev)
{
    FLOAT dStdev, Prd, Tmp;
    int   i;
    int   Error = 0;

    #if (_REBMIXDLL)
    __try {
    #endif
    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    *Stdev = (FLOAT)1.0;

    i = 1; Prd = Sqrt2Pi * ym * fm; Error = 1;
    while ((i <= ItMax) && Error) {
        Tmp = (FLOAT)exp(-(FLOAT)0.5 * (*Stdev) * (*Stdev));

        dStdev = (Prd * (*Stdev) - Tmp) / (Prd + (*Stdev) * Tmp);

        (*Stdev) -= dStdev;

        #if (_REBMIXEXE || _REBMIXR)
        if (IsNan(dStdev) || IsInf(dStdev) || (*Stdev <= (FLOAT)0.0)) {
            Error = 1; goto E0;
        }
        #endif

        if ((FLOAT)fabs(dStdev / (*Stdev)) < Eps) Error = 0;

        i++;
    }

    if (Error) goto E0;

    *Mean = (FLOAT)log(ym) + (*Stdev) * (*Stdev);

    #if (_REBMIXDLL)
    }
    __except(EXCEPTION_EXECUTE_HANDLER) {
        _clearfp(); Error = 1; goto E0;
    }
    #endif

E0: return (Error);
} /* RoughLognormalParameters */

/* Returns rough Weibull parameters. */

int RoughWeibullParameters(FLOAT ym,   
                           FLOAT fm,
                           FLOAT *Theta,
                           FLOAT *Beta)
{
    FLOAT Cor, dBeta, Prd, Tmp;
    int   i;
    int   Error = 0;

    #if (_REBMIXDLL)
    __try {
    #endif
    Error = ym <= FLOAT_MIN; if (Error) goto E0;

    *Beta = (FLOAT)1.0;

    i = 1; Prd = (FLOAT)exp((FLOAT)1.0) * ym * fm; Error = 1;
    while ((i <= ItMax) && Error) {
        Tmp = (FLOAT)exp((FLOAT)1.0 / (*Beta)); Cor = *Beta - (FLOAT)1.0;

        dBeta = (Cor * Tmp - Prd) / (Tmp * ((FLOAT)1.0 - Cor / (*Beta) / (*Beta)));

        *Beta -= dBeta;

        #if (_REBMIXEXE || _REBMIXR)
        if (IsNan(dBeta) || IsInf(dBeta) || (*Beta <= (FLOAT)0.0)) {
            Error = 1; goto E0;
        }
        #endif

        if ((FLOAT)fabs(dBeta / (*Beta)) < Eps) Error = 0; 

        i++;
    }

    if (Error) goto E0;

    *Theta = ym * (FLOAT)exp(log((*Beta) / ((*Beta) - (FLOAT)1.0)) / (*Beta));
    #if (_REBMIXDLL)
    }
    __except(EXCEPTION_EXECUTE_HANDLER) {
        _clearfp(); Error = 1; goto E0;
    }
    #endif

E0: return (Error);
} /* RoughWeibullParameters */

/* Returns rough binomial parameters. */

int RoughBinomialParameters(FLOAT ym,
                            FLOAT fm,
                            FLOAT n,
                            FLOAT *p)
{
    int Error = 0;

    *p = ym / n;

    if (*p == (FLOAT)0.0) {
        *p = (FLOAT)1.0 - (FLOAT)pow(fm, (FLOAT)1.0 / n);
    }
    else
    if (*p == (FLOAT)1.0) {
        *p = (FLOAT)pow(fm, (FLOAT)1.0 / n);
    }

    return (Error);
} /* RoughBinomialParameters */

/* Rough component parameter estimation for k-nearest neighbours. */

int RoughEstimationKNN(int                      n,             /* Total number of independent observations. */
                       int                      d,             /* Number of independent random variables. */ 
                       FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                       int                      k,             /* k-nearest neighbours. */
                       FLOAT                    *h,            /* Normalizing vector. */
                       FLOAT                    nl,            /* Total number of observations in class l. */
                       int                      m,             /* Mode index. */
                       MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                       PestraintsType_e         ResType)       /* Restraints type. */  
{
    const int nf = 15;
    const int ny = 5;

    int                      i, j, l, o, p;
    RoughParameterType       *Mode = NULL;
    FLOAT                    CmpMrgDist, Dc, epsilon, emax, f_lm, dP, Tmp, R;
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    CP[15], CM[15];
    int                      Error = 0;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    TmpDistType = (MarginalDistributionType*)malloc(d * sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    /* Rigid restraints. */

    f_lm = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        if (d > 1) {
            Mode[i].k_lm = (FLOAT)0.0;

            for (j = 0; j < n; j++) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < d; l++) if (i != l) {
                    R = (Y[j][l] - Y[m][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][d + 2]) goto S0;

                Mode[i].k_lm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].k_lm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].f_lm = Y[m][d] * k / (Mode[i].k_lm * (FLOAT)2.0 * Y[m][d + 2] * h[i]); f_lm *= Mode[i].f_lm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] * k / (nl * Y[m][d + 1] * f_lm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].f_lm *= epsilon;

        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = Y[m][d + 2] * h[i];

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = Y[m][d + 2] * h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = Y[m][d + 2] * h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfBinomial:
            Error = 1; goto E0;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    if (nl > Y[m][d] * k) {
        memcpy(TmpDistType, MrgDistType, d * sizeof(MarginalDistributionType));

        for (i = 0; i < d; i++) {
            /* Mode position optimisation. */

            Mode[i].ymin = Mode[i].ym - Mode[i].dy; Mode[i].ymax = Mode[i].ym + Mode[i].dy;

            Mode[i].y = Mode[i].ymin; Mode[i].dy = (FLOAT)2.0 * Mode[i].dy / (ny - (FLOAT)1.0);

            j = 0; Tmp = FLOAT_MAX;

            for (l = 0; l < ny; l++) {
                CP[l] = (FLOAT)0.0;

                switch (TmpDistType[i].ParametricFamily) {
                case pfNormal:
                    Error = RoughNormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                    break;
                case pfLognormal:
                    Error = RoughLognormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1);

                    break;
                case pfWeibull:
                    Error = RoughWeibullParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                    break;
                case pfBinomial:
                    Error = 1; goto E0;
                }

                for (o = 0; o < n; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                    Dc = (FLOAT)0.0;

                    for (p = 0; p < d; p++) if (i != p) {
                        R = (Y[o][p] - Y[m][p]) / h[p]; Dc += R * R;
                    }

                    R = (FLOAT)sqrt(Dc);

                    if (R > Y[m][d + 2]) goto S1;

                    Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                    if (Error) goto E0;

                    dP = (FLOAT)1.0 - Mode[i].k_lm * (FLOAT)2.0 * Y[o][d + 2] * h[i] * CmpMrgDist / Y[o][d] / k;

                    if (dP > CP[l]) CP[l] = dP;
S1:;            } 

                if (CP[l] < Tmp) {
                    j = l; Tmp = CP[l];
                }

                Mode[i].y += Mode[i].dy;
            }

            if ((j == 0) || (j == ny - 1))
                Mode[i].y = Mode[i].ym;
            else
            if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                Mode[i].y = (FLOAT)0.5 * Mode[i].dy * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]) + Mode[i].ymin + Mode[i].dy * j;

            /* Mode empirical density optimisation. */

            Mode[i].y_lmin = FLOAT_MAX; Mode[i].y_lmax = -FLOAT_MAX;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                Dc = (FLOAT)0.0;

                for (l = 0; l < d; l++) if (i != l) {
                    R = (Y[j][l] - Y[m][l]) / h[l]; Dc += R * R;
                }

                R = (FLOAT)sqrt(Dc);

                if (R > Y[m][d + 2]) goto S2;

                if (Y[j][i] < Mode[i].y_lmin) Mode[i].y_lmin = Y[j][i];
                if (Y[j][i] > Mode[i].y_lmax) Mode[i].y_lmax = Y[j][i];
S2:;        } 

            if (Mode[i].y_lmax == Mode[i].y_lmin) goto E0;

            Mode[i].f_lmin = (FLOAT)1.0 / (Mode[i].y_lmax - Mode[i].y_lmin); Mode[i].f_lmax = Mode[i].f_lm;

            if (Mode[i].f_lmax < Mode[i].f_lmin) goto E0;

            Mode[i].f = Mode[i].f_lmax; Mode[i].df = (Mode[i].f_lmax - Mode[i].f_lmin) / (nf - (FLOAT)1.0); emax = (FLOAT)0.0;

            j = 0; Tmp = FLOAT_MAX;

            for (l = 0; l < nf; l++) {
                CP[l] = (FLOAT)0.0; CM[l] = (FLOAT)0.0;

                switch (TmpDistType[i].ParametricFamily) {
                case pfNormal:
                    Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                    break;
                case pfLognormal:
                    Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - 3.09023230616779 * TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + 3.09023230616779 * TmpDistType[i].Parameter1);

                    break;
                case pfWeibull:
                    Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                    break;
                case pfBinomial:
                    Error = 1; goto E0;
                }

                for (o = 0; o < n; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                    Dc = (FLOAT)0.0;

                    for (p = 0; p < d; p++) if (i != p) {
                        R = (Y[o][p] - Y[m][p]) / h[p]; Dc += R * R;
                    }

                    R = (FLOAT)sqrt(Dc);

                    if (R > Y[m][d + 2]) goto S3;

                    Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                    if (Error) goto E0;

                    dP = (FLOAT)1.0 - Mode[i].k_lm * (FLOAT)2.0 * Y[o][d + 2] * h[i] * CmpMrgDist / Y[o][d] / k;

                    if (dP > CP[l])
                        CP[l] = dP;
                    else
                    if (dP < CM[l])
                        CM[l] = dP;
S3:;            } 

                CM[l] = -CM[l];

                if (CP[l] > emax) emax = CP[l];

                if (CP[l] < Tmp) {
                    j = l; Tmp = CP[l];
                }

                if (CM[l] >= emax / ((FLOAT)1.0 - emax)) {
                    if (j == l) j = 0; break;
                }

                Mode[i].f -= Mode[i].df;
            }

            if ((j == 0) || (j == nf - 1))
                Mode[i].f = Mode[i].f_lmax;
            else
            if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                Mode[i].f = Mode[i].f_lmax - Mode[i].df * j - (FLOAT)0.5 * Mode[i].df * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]);

            switch (MrgDistType[i].ParametricFamily) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = 1; goto E0;
            }
        }
    }

E0: if(TmpDistType) free(TmpDistType);

    if(Mode) free(Mode);

    return (Error);
} /* RoughEstimationKNN */

/* Rough component parameter estimation for Parzen window. */

int RoughEstimationPW(int                      n,             /* Total number of independent observations. */
                      int                      d,             /* Number of independent random variables. */ 
                      FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                      FLOAT                    *h,            /* Sides of the hypersquare. */
                      FLOAT                    nl,            /* Total number of observations in class l. */
                      int                      m,             /* Mode index. */
                      MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                      PestraintsType_e         ResType)       /* Restraints type. */  
{
    const int nf = 15;
    const int ny = 5;

    int                      i, j, l, o, p;
    RoughParameterType       *Mode = NULL;
    FLOAT                    CmpMrgDist, epsilon, emax, f_lm, dP, Tmp, V;
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    CP[15], CM[15];
    int                      Error = 0;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    TmpDistType = (MarginalDistributionType*)malloc(d * sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    /* Rigid restraints. */

    f_lm = (FLOAT)1.0; V = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        V *= h[i];

        if (d > 1) {
            Mode[i].k_lm = (FLOAT)0.0;

            for (j = 0; j < n; j++) {
                for (l = 0; l < d; l++) if ((i != l) && ((FLOAT)fabs(Y[j][l] - Y[m][l]) > (FLOAT)0.5 * h[l])) goto S0;

                Mode[i].k_lm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].k_lm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].f_lm = Y[m][d] * Y[m][d + 1] / (Mode[i].k_lm * h[i]); f_lm *= Mode[i].f_lm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] * Y[m][d + 1] / (nl * V * f_lm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].f_lm *= epsilon;

        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = (FLOAT)0.5 * h[i];

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = (FLOAT)0.5 * h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = (FLOAT)0.5 * h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].f_lm, MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = (FLOAT)0.0;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    if (nl > Y[m][d] * Y[m][d + 1]) {
        memcpy(TmpDistType, MrgDistType, d * sizeof(MarginalDistributionType));

        for (i = 0; i < d; i++) {
            /* Mode position optimisation. */

            Mode[i].ymin = Mode[i].ym - Mode[i].dy; Mode[i].ymax = Mode[i].ym + Mode[i].dy;

            Mode[i].y = Mode[i].ymin; Mode[i].dy = (FLOAT)2.0 * Mode[i].dy / (ny - (FLOAT)1.0);

            if (Mode[i].dy != (FLOAT)0.0) {
                j = 0; Tmp = FLOAT_MAX;

                for (l = 0; l < ny; l++) {
                    CP[l] = (FLOAT)0.0;

                    switch (TmpDistType[i].ParametricFamily) {
                    case pfNormal:
                        Error = RoughNormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                        Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                        break;
                    case pfLognormal:
                        Error = RoughLognormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - 3.09023230616779 * TmpDistType[i].Parameter1);
                        Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + 3.09023230616779 * TmpDistType[i].Parameter1);

                        break;
                    case pfWeibull:
                        Error = RoughWeibullParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                        Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                        break;
                    case pfBinomial:
                        Error = RoughBinomialParameters(Mode[i].y, Mode[i].f_lm, TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);
                                                       
                        if (Error) goto E0;

                        Mode[i].b_lmin = (FLOAT)0.0;
                        Mode[i].b_lmax = TmpDistType[i].Parameter0;
                    }

                    for (o = 0; o < n; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                        for (p = 0; p < d; p++) if ((i != p) && ((FLOAT)fabs(Y[o][p] - Y[m][p]) > (FLOAT)0.5 * h[p])) goto S1;

                        Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                        if (Error) goto E0;

                        dP = (FLOAT)1.0 - Mode[i].k_lm * h[i] * CmpMrgDist / Y[o][d] / Y[o][d + 1];

                        if (dP > CP[l]) CP[l] = dP;
S1:;                }

                    if (CP[l] < Tmp) {
                        j = l; Tmp = CP[l];
                    }

                    Mode[i].y += Mode[i].dy;
                }

                if ((j == 0) || (j == ny - 1))
                    Mode[i].y = Mode[i].ym;
                else
                if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                    Mode[i].y = (FLOAT)0.5 * Mode[i].dy * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]) + Mode[i].ymin + Mode[i].dy * j;
            }

            /* Mode empirical density optimisation. */

            Mode[i].y_lmin = FLOAT_MAX; Mode[i].y_lmax = -FLOAT_MAX;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                for (l = 0; l < d; l++) if ((i != l) && ((FLOAT)fabs(Y[j][l] - Y[m][l]) > (FLOAT)0.5 * h[l])) goto S2;

                if (Y[j][i] < Mode[i].y_lmin) Mode[i].y_lmin = Y[j][i];
                if (Y[j][i] > Mode[i].y_lmax) Mode[i].y_lmax = Y[j][i];
S2:;        }

            if (Mode[i].y_lmax == Mode[i].y_lmin) goto E0;

            Mode[i].f_lmin = (FLOAT)1.0 / (Mode[i].y_lmax - Mode[i].y_lmin); Mode[i].f_lmax = Mode[i].f_lm;

            if (Mode[i].f_lmax < Mode[i].f_lmin) goto E0;

            Mode[i].f = Mode[i].f_lmax; Mode[i].df = (Mode[i].f_lmax - Mode[i].f_lmin) / (nf - (FLOAT)1.0); emax = (FLOAT)0.0;

            j = 0; Tmp = FLOAT_MAX;

            for (l = 0; l < nf; l++) {
                CP[l] = (FLOAT)0.0; CM[l] = (FLOAT)0.0;

                switch (TmpDistType[i].ParametricFamily) {
                case pfNormal:
                    Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                    break;
                case pfLognormal:
                    Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - 3.09023230616779 * TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + 3.09023230616779 * TmpDistType[i].Parameter1);

                    break;
                case pfWeibull:
                    Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                    break;
                case pfBinomial:
                    Error = RoughBinomialParameters(Mode[i].y, Mode[i].f, TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);
                                                       
                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)0.0;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0;
                }

                for (o = 0; o < n; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                    for (p = 0; p < d; p++) if ((i != p) && ((FLOAT)fabs(Y[o][p] - Y[m][p]) > (FLOAT)0.5 * h[p])) goto S3;

                    Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                    if (Error) goto E0;

                    dP = (FLOAT)1.0 - Mode[i].k_lm * h[i] * CmpMrgDist / Y[o][d] / Y[o][d + 1];

                    if (dP > CP[l])
                        CP[l] = dP;
                    else
                    if (dP < CM[l])
                        CM[l] = dP;
S3:;            }

                CM[l] = -CM[l];

                if (CP[l] > emax) emax = CP[l];

                if (CP[l] < Tmp) {
                    j = l; Tmp = CP[l];
                }

                if (CM[l] >= emax / ((FLOAT)1.0 - emax)) {
                    if (j == l) j = 0; break;
                }

                Mode[i].f -= Mode[i].df;
            }

            if ((j == 0) || (j == nf - 1))
                Mode[i].f = Mode[i].f_lmax;
            else
            if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                Mode[i].f = Mode[i].f_lmax - Mode[i].df * j - (FLOAT)0.5 * Mode[i].df * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]);

            switch (MrgDistType[i].ParametricFamily) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].y, Mode[i].f, MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);
                                                       
                if (Error) goto E0;
            }
        }
    }

E0: if(TmpDistType) free(TmpDistType);

    if(Mode) free(Mode);

    return (Error);
} /* RoughEstimationPW */

/* Rough component parameter estimation for histogram. */

int RoughEstimationH(int                      k,             /* Total number of bins. */
                     int                      d,             /* Number of independent random variables. */ 
                     FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl]. */
                     FLOAT                    *h,            /* Sides of the hypersquare. */
                     FLOAT                    nl,            /* Total number of observations in class l. */
                     int                      m,             /* Mode index. */
                     MarginalDistributionType *MrgDistType,  /* Marginal distribution type. */
                     PestraintsType_e         ResType)       /* Restraints type. */  
{
    const int nf = 15;
    const int ny = 5;

    int                      i, j, l, o, p;
    RoughParameterType       *Mode = NULL;
    FLOAT                    CmpMrgDist, epsilon, emax, f_lm, dP, Tmp, V;
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    CP[15], CM[15];
    int                      Error = 0;

    Mode = (RoughParameterType*)malloc(d * sizeof(RoughParameterType));

    Error = NULL == Mode; if (Error) goto E0;

    TmpDistType = (MarginalDistributionType*)malloc(d * sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    /* Rigid restraints. */

    f_lm = (FLOAT)1.0; V = (FLOAT)1.0;

    for (i = 0; i < d; i++) {
        V *= h[i];

        if (d > 1) {
            Mode[i].k_lm = (FLOAT)0.0;

            for (j = 0; j < k; j++) {
                for (l = 0; l < d; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S0;

                Mode[i].k_lm += Y[j][d];
S0:;        }
        }
        else
            Mode[i].k_lm = nl;

        Mode[i].ym = Y[m][i]; Mode[i].f_lm = Y[m][d] / (Mode[i].k_lm * h[i]); f_lm *= Mode[i].f_lm;
    }

    epsilon = (FLOAT)exp(log(Y[m][d] / (nl * V * f_lm)) / d);

    for (i = 0; i < d; i++) {
        if (epsilon < (FLOAT)1.0) Mode[i].f_lm *= epsilon;

        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            Error = RoughNormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = h[i];

            break;
        case pfLognormal:
            Error = RoughLognormalParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfWeibull:
            Error = RoughWeibullParameters(Mode[i].ym, Mode[i].f_lm, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = h[i]; Tmp = Mode[i].ym - Mode[i].dy;

            if (Tmp <= (FLOAT)0.0) Mode[i].dy = Mode[i].ym * ((FLOAT)1.0 - Eps);

            break;
        case pfBinomial:
            Error = RoughBinomialParameters(Mode[i].ym, Mode[i].f_lm, MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

            if (Error) goto E0;

            Mode[i].dy = (FLOAT)0.0;
        }
    }

    if (ResType == rtRigid) goto E0;

    /* Loose restraints. */

    if (nl > Y[m][d]) {
        memcpy(TmpDistType, MrgDistType, d * sizeof(MarginalDistributionType));

        for (i = 0; i < d; i++) {
            /* Mode position optimisation. */

            Mode[i].ymin = Mode[i].ym - Mode[i].dy; Mode[i].ymax = Mode[i].ym + Mode[i].dy;

            Mode[i].y = Mode[i].ymin; Mode[i].dy = (FLOAT)2.0 * Mode[i].dy / (ny - (FLOAT)1.0);

            if (Mode[i].dy != (FLOAT)0.0) {
                j = 0; Tmp = FLOAT_MAX;

                for (l = 0; l < ny; l++) {
                    CP[l] = (FLOAT)0.0;

                    switch (TmpDistType[i].ParametricFamily) {
                    case pfNormal:
                        Error = RoughNormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                        Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                        break;
                    case pfLognormal:
                        Error = RoughLognormalParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - 3.09023230616779 * TmpDistType[i].Parameter1);
                        Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + 3.09023230616779 * TmpDistType[i].Parameter1);

                        break;
                    case pfWeibull:
                        Error = RoughWeibullParameters(Mode[i].y, Mode[i].f_lm, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                        if (Error) goto E0;

                        Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                        Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                        break;
                    case pfBinomial:
                        Error = RoughBinomialParameters(Mode[i].y, Mode[i].f_lm, TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);
                                                       
                        if (Error) goto E0;

                        Mode[i].b_lmin = (FLOAT)0.0;
                        Mode[i].b_lmax = TmpDistType[i].Parameter0;
                    }

                    for (o = 0; o < k; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                        for (p = 0; p < d; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S1;

                        Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                        if (Error) goto E0;

                        dP = (FLOAT)1.0 - Mode[i].k_lm * h[i] * CmpMrgDist / Y[o][d];

                        if (dP > CP[l]) CP[l] = dP;
S1:;                }

                    if (CP[l] < Tmp) {
                        j = l; Tmp = CP[l];
                    }

                    Mode[i].y += Mode[i].dy;
                }

                if ((j == 0) || (j == ny - 1))
                    Mode[i].y = Mode[i].ym;
                else
                if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                    Mode[i].y = (FLOAT)0.5 * Mode[i].dy * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]) + Mode[i].ymin + Mode[i].dy * j;
            }

            /* Mode empirical density optimisation. */

            Mode[i].y_lmin = FLOAT_MAX; Mode[i].y_lmax = -FLOAT_MAX;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                for (l = 0; l < d; l++) if ((i != l) && (Y[j][l] != Y[m][l])) goto S2;

                if (Y[j][i] < Mode[i].y_lmin) Mode[i].y_lmin = Y[j][i];
                if (Y[j][i] > Mode[i].y_lmax) Mode[i].y_lmax = Y[j][i];
S2:;        }

            if (Mode[i].y_lmax == Mode[i].y_lmin) goto E0;

            Mode[i].f_lmin = (FLOAT)1.0 / (Mode[i].y_lmax - Mode[i].y_lmin); Mode[i].f_lmax = Mode[i].f_lm;

            if (Mode[i].f_lmax < Mode[i].f_lmin) goto E0;

            Mode[i].f = Mode[i].f_lmax; Mode[i].df = (Mode[i].f_lmax - Mode[i].f_lmin) / (nf - (FLOAT)1.0); emax = (FLOAT)0.0;

            j = 0; Tmp = FLOAT_MAX;

            for (l = 0; l < nf; l++) {
                CP[l] = (FLOAT)0.0; CM[l] = (FLOAT)0.0;

                switch (TmpDistType[i].ParametricFamily) {
                case pfNormal:
                    Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 - (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 + (FLOAT)3.09023230616779 * TmpDistType[i].Parameter1;

                    break;
                case pfLognormal:
                    Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)exp(TmpDistType[i].Parameter0 - 3.09023230616779 * TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = (FLOAT)exp(TmpDistType[i].Parameter0 + 3.09023230616779 * TmpDistType[i].Parameter1);

                    break;
                case pfWeibull:
                    Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);

                    if (Error) goto E0;

                    Mode[i].b_lmin = TmpDistType[i].Parameter0 * (FLOAT)exp(log(0.00100050033358) / TmpDistType[i].Parameter1);
                    Mode[i].b_lmax = TmpDistType[i].Parameter0 * (FLOAT)exp(log(6.90775527898214) / TmpDistType[i].Parameter1);

                    break;
                case pfBinomial:
                    Error = RoughBinomialParameters(Mode[i].y, Mode[i].f, TmpDistType[i].Parameter0, &TmpDistType[i].Parameter1);
                                                       
                    if (Error) goto E0;

                    Mode[i].b_lmin = (FLOAT)0.0;
                    Mode[i].b_lmax = TmpDistType[i].Parameter0;
                }

                for (o = 0; o < k; o++) if ((Y[o][d] > FLOAT_MIN) && (Y[o][i] > Mode[i].b_lmin) && (Y[o][i] < Mode[i].b_lmax)) {
                    for (p = 0; p < d; p++) if ((i != p) && (Y[o][p] != Y[m][p])) goto S3;

                    Error = ComponentMarginalDist(i, Y[o], TmpDistType, &CmpMrgDist, 0);

                    if (Error) goto E0;

                    dP = (FLOAT)1.0 - Mode[i].k_lm * h[i] * CmpMrgDist / Y[o][d];

                    if (dP > CP[l])
                        CP[l] = dP;
                    else
                    if (dP < CM[l])
                        CM[l] = dP;
S3:;            }

                CM[l] = -CM[l];

                if (CP[l] > emax) emax = CP[l];

                if (CP[l] < Tmp) {
                    j = l; Tmp = CP[l];
                }

                if (CM[l] >= emax / ((FLOAT)1.0 - emax)) {
                    if (j == l) j = 0; break;
                }

                Mode[i].f -= Mode[i].df;
            }

            if ((j == 0) || (j == nf - 1))
                Mode[i].f = Mode[i].f_lmax;
            else
            if ((CP[j + 1] > CP[j]) && (CP[j - 1] > CP[j]))
                Mode[i].f = Mode[i].f_lmax - Mode[i].df * j - (FLOAT)0.5 * Mode[i].df * (CP[j - 1] - CP[j + 1]) / (CP[j + 1] - (FLOAT)2.0 * CP[j] + CP[j - 1]);

            switch (MrgDistType[i].ParametricFamily) {
            case pfNormal:
                Error = RoughNormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfLognormal:
                Error = RoughLognormalParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfWeibull:
                Error = RoughWeibullParameters(Mode[i].y, Mode[i].f, &MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);

                if (Error) goto E0;

                break;
            case pfBinomial:
                Error = RoughBinomialParameters(Mode[i].y, Mode[i].f, MrgDistType[i].Parameter0, &MrgDistType[i].Parameter1);
                                                       
                if (Error) goto E0;
            }
        }
    }

E0: if(TmpDistType) free(TmpDistType);

    if(Mode) free(Mode);

    return (Error);
} /* RoughEstimationH */

/* Enhanced component parameter estimation for k-nearest neighbours. */

int EnhancedEstimationKNN(int                      n,             /* Total number of independent observations. */
                          int                      d,             /* Number of independent random variables. */ 
                          FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,V,R]. */
                          FLOAT                    RMIN,          /* Minimum radius of the hypersphere. */
                          FLOAT                    nl,            /* Total number of observations in class l. */
                          MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP, DMIN;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    DMIN = (FLOAT)2.0 * RMIN;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            TmpDistType[i].ParametricFamily = pfNormal;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Parameter0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Parameter0 /= nl;

            for(j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Parameter0; 

                TmpDistType[i].Parameter1 += Y[j][d] * T[0] * T[0];
            }

            T[0] = TmpDistType[i].Parameter1 / nl;

            if (T[0] < DMIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(T[0]);

            break;
        case pfLognormal:
            TmpDistType[i].ParametricFamily = pfLognormal;

            for (j = 0; j < n; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Parameter0 += T[0]; 
                    TmpDistType[i].Parameter1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Parameter0 /= nl; 
            TmpDistType[i].Parameter1 = TmpDistType[i].Parameter1 / nl - TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Parameter0 + TmpDistType[i].Parameter1) * ((FLOAT)exp(TmpDistType[i].Parameter1) - (FLOAT)1.0);

            if (T[0] < DMIN) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(TmpDistType[i].Parameter1);

            break;
        case pfWeibull:
            #if (_REBMIXDLL)
            __try {
            #endif
            TmpDistType[i].ParametricFamily = pfWeibull;

            TmpDistType[i].Parameter1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Parameter1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Parameter1 * TmpDistType[i].Parameter1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Parameter1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Parameter1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Parameter1 <= (FLOAT)0.0)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Parameter1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Parameter0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Parameter1);

            T[1] = TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = T[1] * ((FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Parameter1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Parameter1)));

            if (T[0] < DMIN) {
                Error = 1; if (Error) goto E0;
            }
            #if (_REBMIXDLL)
            }
            __except(EXCEPTION_EXECUTE_HANDLER) {
                _clearfp(); Error = 1; goto E0;
            }
            #endif

            break;
        case pfBinomial:
            Error = 1; goto E0;
        }
    }

    memcpy(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationKNN */

/* Enhanced component parameter estimation for Parzen window. */

int EnhancedEstimationPW(int                      n,             /* Total number of independent observations. */
                         int                      d,             /* Number of independent random variables. */ 
                         FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                         FLOAT                    *h,            /* Sides of the hypersquare. */
                         FLOAT                    nl,            /* Total number of observations in class l. */
                         MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            TmpDistType[i].ParametricFamily = pfNormal;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Parameter0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Parameter0 /= nl;

            for(j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Parameter0; 

                TmpDistType[i].Parameter1 += Y[j][d] * T[0] * T[0];
            }

            T[0] = TmpDistType[i].Parameter1 / nl;

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(T[0]);

            break;
        case pfLognormal:
            TmpDistType[i].ParametricFamily = pfLognormal;

            for (j = 0; j < n; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Parameter0 += T[0]; 
                    TmpDistType[i].Parameter1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Parameter0 /= nl; 
            TmpDistType[i].Parameter1 = TmpDistType[i].Parameter1 / nl - TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Parameter0 + TmpDistType[i].Parameter1) * ((FLOAT)exp(TmpDistType[i].Parameter1) - (FLOAT)1.0);

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(TmpDistType[i].Parameter1);

            break;
        case pfWeibull:
            #if (_REBMIXDLL)
            __try {
            #endif
            TmpDistType[i].ParametricFamily = pfWeibull;

            TmpDistType[i].Parameter1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < n; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Parameter1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Parameter1 * TmpDistType[i].Parameter1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Parameter1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Parameter1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Parameter1 <= (FLOAT)0.0)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Parameter1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Parameter0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Parameter1);

            T[1] = TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = T[1] * ((FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Parameter1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Parameter1)));

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }
            #if (_REBMIXDLL)
            }
            __except(EXCEPTION_EXECUTE_HANDLER) {
                _clearfp(); Error = 1; goto E0;
            }
            #endif

            break;
        case pfBinomial:
            TmpDistType[i].ParametricFamily = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < n; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Parameter1 = T[0] / TmpDistType[i].Parameter0 / nl;
        }
    }

    memcpy(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationPW */

/* Enhanced component parameter estimation for histogram. */

int EnhancedEstimationH(int                      k,             /* Total number of bins. */
                        int                      d,             /* Number of independent random variables. */ 
                        FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1,kl,k]. */
                        FLOAT                    *h,            /* Sides of the hypersquare. */
                        FLOAT                    nl,            /* Total number of observations in class l. */
                        MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    MarginalDistributionType *TmpDistType = NULL;
    FLOAT                    A[4], T[2];
    int                      i, j, l;
    FLOAT                    dP;
    int                      Error = 0;

    TmpDistType = (MarginalDistributionType*)calloc(d, sizeof(MarginalDistributionType));

    Error = NULL == TmpDistType; if (Error) goto E0;

    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            TmpDistType[i].ParametricFamily = pfNormal;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                TmpDistType[i].Parameter0 += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Parameter0 /= nl;

            for(j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] = Y[j][i] - TmpDistType[i].Parameter0; 

                TmpDistType[i].Parameter1 += Y[j][d] * T[0] * T[0];
            }

            T[0] = TmpDistType[i].Parameter1 / nl;

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(T[0]);

            break;
        case pfLognormal:
            TmpDistType[i].ParametricFamily = pfLognormal;

            for (j = 0; j < k; j++) {
                if ((Y[j][d] > FLOAT_MIN) && (Y[j][i] > FLOAT_MIN)) {
                    T[0] = Y[j][d] * (FLOAT)log(Y[j][i]);

                    TmpDistType[i].Parameter0 += T[0]; 
                    TmpDistType[i].Parameter1 += T[0] * (FLOAT)log(Y[j][i]);
                }
            }

            TmpDistType[i].Parameter0 /= nl; 
            TmpDistType[i].Parameter1 = TmpDistType[i].Parameter1 / nl - TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = (FLOAT)exp((FLOAT)2.0 * TmpDistType[i].Parameter0 + TmpDistType[i].Parameter1) * ((FLOAT)exp(TmpDistType[i].Parameter1) - (FLOAT)1.0);

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }

            TmpDistType[i].Parameter1 = (FLOAT)sqrt(TmpDistType[i].Parameter1);

            break;
        case pfWeibull:
            #if (_REBMIXDLL)
            __try {
            #endif
            TmpDistType[i].ParametricFamily = pfWeibull;

            TmpDistType[i].Parameter1 = (FLOAT)1.0;

            j = 1; Error = 1;
            while ((j <= ItMax) && Error) {
                memset(&A, 0, 4 * sizeof(FLOAT));

                for (l = 0; l < k; l++) {
                    if ((Y[l][d] > FLOAT_MIN) && (Y[l][i] > FLOAT_MIN)) {
                        T[0] = (FLOAT)log(Y[l][i]);
                        T[1] = (FLOAT)exp(T[0] * TmpDistType[i].Parameter1);

                        A[0] += Y[l][d] * T[0];
                        A[1] += Y[l][d] * T[1] * T[0];
                        A[2] += Y[l][d] * T[1];
                        A[3] += Y[l][d] * T[1] * T[0] * T[0];
                    }
                }

                A[0] /= nl; T[0] = A[1] / A[2]; T[0] *= T[0]; T[1] = TmpDistType[i].Parameter1 * TmpDistType[i].Parameter1;

                dP = ((FLOAT)1.0 / TmpDistType[i].Parameter1 + A[0] - A[1] / A[2]) / (T[0] - A[3] / A[2] - (FLOAT)1.0 / T[1]);

                TmpDistType[i].Parameter1 -= dP;

                #if (_REBMIXEXE || _REBMIXR)
                if (IsNan(dP) || IsInf(dP) || (TmpDistType[i].Parameter1 <= (FLOAT)0.0)) {
                    Error = 1; goto E0;
                }
                #endif

                if ((FLOAT)fabs(dP / TmpDistType[i].Parameter1) < Eps) Error = 0;

                j++;
            }

            if (Error) goto E0;

            A[2] /= nl;

            TmpDistType[i].Parameter0 = (FLOAT)exp(log(A[2]) / TmpDistType[i].Parameter1);

            T[1] = TmpDistType[i].Parameter0 * TmpDistType[i].Parameter0;

            T[0] = T[1] * ((FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / TmpDistType[i].Parameter1)) - (FLOAT)exp((FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / TmpDistType[i].Parameter1)));

            if (T[0] < h[i] * h[i]) {
                Error = 1; if (Error) goto E0;
            }
            #if (_REBMIXDLL)
            }
            __except(EXCEPTION_EXECUTE_HANDLER) {
                _clearfp(); Error = 1; goto E0;
            }
            #endif

            break;
        case pfBinomial:
            TmpDistType[i].ParametricFamily = pfBinomial;

            T[0] = (FLOAT)0.0;

            for (j = 0; j < k; j++) if (Y[j][d] > FLOAT_MIN) {
                T[0] += Y[j][d] * Y[j][i];
            }

            TmpDistType[i].Parameter1 = T[0] / TmpDistType[i].Parameter0 / nl;
       }
    }

    memcpy(MrgDistType, TmpDistType, d * sizeof(MarginalDistributionType));

E0: if (TmpDistType) free(TmpDistType);

    return (Error);
} /* EnhancedEstimationH */

/* Component mean and variance calculation. */

int MeanVarianceCalculation(int                      d,            /* Number of independent random variables. */ 
                            MarginalDistributionType *MrgDistType, /* Marginal distribution type. */
                            FLOAT                    *Mean,        /* Mean. */
                            FLOAT                    *Variance)    /* Variance. */
{
    int i;
    int Error = 0;
    
    for (i = 0; i < d; i++) {
        switch (MrgDistType[i].ParametricFamily) {
        case pfNormal:
            Mean[i] = MrgDistType[i].Parameter0;

            Variance[i] = MrgDistType[i].Parameter1 * MrgDistType[i].Parameter1 + MrgDistType[i].Parameter0 * MrgDistType[i].Parameter0;

            break;
        case pfLognormal:
            Mean[i] = (FLOAT)exp(MrgDistType[i].Parameter0 + (FLOAT)0.5 * MrgDistType[i].Parameter1 * MrgDistType[i].Parameter1); 
            
            Variance[i] = (FLOAT)exp((FLOAT)2.0 * (MrgDistType[i].Parameter0 + MrgDistType[i].Parameter1 * MrgDistType[i].Parameter1));

            break;
        case pfWeibull:
            Mean[i] = MrgDistType[i].Parameter0 * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType[i].Parameter1));

            Variance[i] = MrgDistType[i].Parameter0 * MrgDistType[i].Parameter0 * (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)2.0 / MrgDistType[i].Parameter1));

            break;
        case pfBinomial:
            Mean[i] = MrgDistType[i].Parameter0 * MrgDistType[i].Parameter1;

            Variance[i] = (FLOAT)0.0;
        }
    }

    return(Error);
} /* MeanVarianceCalculation */

/* Returns Bayes Weibull parameters. */

int BayesWeibullParameters(FLOAT                    Mean,          /* Mean. */
                           FLOAT                    Variance,      /* Variance. */
                           MarginalDistributionType *MrgDistType)  /* Marginal distribution type. */
{
    FLOAT A;
    FLOAT xl, xm, xh, fl, fm, fh, dx;
    FLOAT i;
    int   Error = 0;
    
    A = (FLOAT)log(Variance / Mean / Mean); xl = 0.001; xh = (FLOAT)10.0;

    fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
    fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        if (fl * fh < (FLOAT)0.0)
            Error = 0;
        else
        if ((FLOAT)fabs(fl) < (FLOAT)fabs(fh)) {
            xl += (FLOAT)1.6 * (xl - xh);
            fl = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xl) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xl);
        }
        else {
            xh += (FLOAT)1.6 * (xh - xl);
            fh = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xh) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xh);
        }

        i++;
    }

    if (Error) goto E0;

    /* Root must be bracketed for bisection. */

    if (fl < (FLOAT)0.0) {
        MrgDistType->Parameter1 = xl; dx = xh - xl;
    }
    else {
        MrgDistType->Parameter1 = xh; dx = xl - xh;
    }

    i = 1; Error = 1;
    while ((i <= ItMax) && Error) {
        dx = (FLOAT)0.5 * dx; xm = MrgDistType->Parameter1 + dx;

        fm = A - Gammaln((FLOAT)1.0 + (FLOAT)2.0 / xm) + (FLOAT)2.0 * Gammaln((FLOAT)1.0 + (FLOAT)1.0 / xm);

        if (fm <= FLOAT_MIN) MrgDistType->Parameter1 = xm;

        if (((FLOAT)fabs(dx) < Eps) || ((FLOAT)fabs(fm) < Eps)) Error = 0;

        i++;
    }

    if (Error) goto E0;

    MrgDistType->Parameter0 = Mean / (FLOAT)exp(Gammaln((FLOAT)1.0 + (FLOAT)1.0 / MrgDistType->Parameter1));

E0: return (Error);
} /* BayesWeibullParameters */

/* Bayes classification of the remaining observations for k-nearest neighbour. */

int BayesClassificationKNN(int                      n,             /* Total number of independent observations. */
                           int                      d,             /* Number of independent random variables. */
                           FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                           int                      c,             /* Number of components. */ 
                           FLOAT                    *W,            /* Component weights. */
                           MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                           FLOAT                    **Mean,        /* Means. */
                           FLOAT                    **Variance)    /* Variances. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < n; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                Mean[l][j] += dW * (Y[i][j] - Mean[l][j]) / W[l];

                Variance[l][j] += dW * (Y[i][j] * Y[i][j] - Variance[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParametricFamily) {
        case pfNormal: 
            MrgDistType[i][j].Parameter0 = Mean[i][j]; 
            
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(Variance[i][j] - MrgDistType[i][j].Parameter0 * MrgDistType[i][j].Parameter0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Parameter0 = (FLOAT)2.0 * (FLOAT)log(Mean[i][j]) - (FLOAT)0.5 * (FLOAT)log(Variance[i][j]);
            
           
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(log(Variance[i][j]) - (FLOAT)2.0 * log(Mean[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(Mean[i][j], Variance[i][j], &MrgDistType[i][j]);

            break;
        case pfBinomial:
            Error = 1; goto E0;
        }
    }

E0: return (Error);
} /* BayesClassificationKNN */

/* Bayes classification of the remaining observations for Parzen window. */

int BayesClassificationPW(int                      n,             /* Total number of independent observations. */
                          int                      d,             /* Number of independent random variables. */
                          FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                          int                      c,             /* Number of components. */ 
                          FLOAT                    *W,            /* Component weights. */
                          MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                          FLOAT                    **Mean,        /* Means. */
                          FLOAT                    **Variance)    /* Variances. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < n; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                Mean[l][j] += dW * (Y[i][j] - Mean[l][j]) / W[l];

                Variance[l][j] += dW * (Y[i][j] * Y[i][j] - Variance[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParametricFamily) {
        case pfNormal: 
            MrgDistType[i][j].Parameter0 = Mean[i][j]; 
            
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(Variance[i][j] - MrgDistType[i][j].Parameter0 * MrgDistType[i][j].Parameter0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Parameter0 = (FLOAT)2.0 * (FLOAT)log(Mean[i][j]) - (FLOAT)0.5 * (FLOAT)log(Variance[i][j]);
            
           
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(log(Variance[i][j]) - (FLOAT)2.0 * log(Mean[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(Mean[i][j], Variance[i][j], &MrgDistType[i][j]);

            break;
        case pfBinomial:
            MrgDistType[i][j].Parameter1 = Mean[i][j] / MrgDistType[i][j].Parameter0;
        }
    }

E0: return (Error);
} /* BayesClassificationPW */

/* Bayes classification of the remaining observations for histogram. */

int BayesClassificationH(int                      k,             /* Total number of bins. */
                         int                      n,             /* Total number of independent observations. */
                         int                      d,             /* Number of independent random variables. */
                         FLOAT                    **Y,           /* Pointer to the input points [y0,...,yd-1]. */
                         int                      c,             /* Number of components. */ 
                         FLOAT                    *W,            /* Component weights. */
                         MarginalDistributionType **MrgDistType, /* Marginal distribution type. */
                         FLOAT                    **Mean,        /* Means. */
                         FLOAT                    **Variance)    /* Variances. */
{
    int   i, j, l;
    FLOAT CmpDist, Max, Tmp, dW;
    int   Error = 0;

    for (i = 0; i < k; i++) {
        if (Y[i][d] > FLOAT_MIN) {
            l = 0;

            Error = ComponentDist(d, Y[i], MrgDistType[l], &CmpDist, 0);

            if (Error) goto E0;

            Max = W[l] * CmpDist;

            for (j = 1; j < c; j++) {
                Error = ComponentDist(d, Y[i], MrgDistType[j], &CmpDist, 0);

                if (Error) goto E0;

                Tmp = W[j] * CmpDist;

                if (Tmp > Max) { 
                    l = j; Max = Tmp; 
                }
            }

            dW = Y[i][d] / n; W[l] += dW;

            for (j = 0; j < d; j++) {
                Mean[l][j] += dW * (Y[i][j] - Mean[l][j]) / W[l];

                Variance[l][j] += dW * (Y[i][j] * Y[i][j] - Variance[l][j]) / W[l];
            }
        }
    }

    for (i = 0; i < c; i++) for (j = 0; j < d; j++) { 
        switch (MrgDistType[i][j].ParametricFamily) {
        case pfNormal: 
            MrgDistType[i][j].Parameter0 = Mean[i][j]; 
            
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(Variance[i][j] - MrgDistType[i][j].Parameter0 * MrgDistType[i][j].Parameter0);

            break;
        case pfLognormal: 
            MrgDistType[i][j].Parameter0 = (FLOAT)2.0 * (FLOAT)log(Mean[i][j]) - (FLOAT)0.5 * (FLOAT)log(Variance[i][j]);
            
           
            MrgDistType[i][j].Parameter1 = (FLOAT)sqrt(log(Variance[i][j]) - (FLOAT)2.0 * log(Mean[i][j]));

            break;
        case pfWeibull:
            BayesWeibullParameters(Mean[i][j], Variance[i][j], &MrgDistType[i][j]);

            break;
        case pfBinomial:
            MrgDistType[i][j].Parameter1 = Mean[i][j] / MrgDistType[i][j].Parameter0;
        }
    }

E0: return (Error);
} /* BayesClassificationH */

/* REBMIX algorithm for k-nearest neighbours. */

int REBMIXKNN(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
              OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    FLOAT                    **Y = NULL;
    FLOAT                    *h = NULL, *ymin = NULL, *ymax = NULL;
    FLOAT                    *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                    *W = NULL;
    MarginalDistributionType **Theta = NULL; 
    FLOAT                    **Mean = NULL, **Variance = NULL;
    int                      c = 0, i, I, j, J, l, m;
    FLOAT                    Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL;
    int                      Error = 0, Stop = 0;
    #if(_DEBUG)
    int                      o;
    FILE                     *fp0 = NULL;
    #endif

    #if(_DEBUG)
    if ((fp0 = fopen("debug.txt", "w")) == NULL) {
        Error = 1; goto E0;
    }
    #endif

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType.d + 3) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < InpParType.d; j++) Y[i][j] = OutParType->X[i][j];
    }

    /* Allocation and normalizing vector h calculation. */

    h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    ymin = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    ymax = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        ymin[i] = ymax[i] = OutParType->X[0][i];

        for (j = 1; j < OutParType->n; j++) {
            if (OutParType->X[j][i] < ymin[i]) ymin[i] = OutParType->X[j][i];
            if (OutParType->X[j][i] > ymax[i]) ymax[i] = OutParType->X[j][i];
        }

        h[i] = ymax[i] - ymin[i];
    }

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].ParametricFamily = InpParType.ParFamType[j];
        }

        if (InpParType.Par0 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter0 = InpParType.Par0[j];
        }

        if (InpParType.Par1 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter1 = InpParType.Par1[j];
        }
    }

    Mean = (FLOAT**)malloc(InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Mean; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Mean[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Mean[i]; if (Error) goto E0;
    }

    Variance = (FLOAT**)realloc(Variance, InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Variance; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Variance[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Variance[i]; if (Error) goto E0;
    }

    OutParType->IC = FLOAT_MAX; OutParType->W = NULL; OutParType->Theta = NULL;

    OutParType->h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    for (i = 0; i < InpParType.kmax; i++) {
        /* Preprocessing of observations. */

        Error = PreprocessingKNN(InpParType.k[i], InpParType.RMIN, h, OutParType->n, InpParType.d, Y);

        if (Error) goto E0;

        #if(_DEBUG)
        fprintf(fp0, "%s\t%d\n%s", "k", InpParType.k[i], "h");

        for (j = 0; j < InpParType.d; j++) fprintf(fp0, "\t%E", h[j]);

        fprintf(fp0, "\n");

        for (j = 0; j < OutParType->n; j++) {
            fprintf(fp0, "%s%d%s", "Y[", j, "]");

            for (o = 0; o < InpParType.d + 3; o++) fprintf(fp0, "\t%E", Y[j][o]);

            fprintf(fp0, "\n");
        }
        #endif

        Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * (l + 1) * Dmin) {
                /* Global mode detection. */

                Error = GlobalModeKNN(h, ymin, &m, OutParType->n, InpParType.d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, OutParType->n * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) {
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationKNN(OutParType->n, InpParType.d, Y, InpParType.k[i], h, nl, m, Theta[l], InpParType.ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < OutParType->n; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType.d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType.d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType.d] - nl * fl * Y[j][InpParType.d + 1] / InpParType.k[i];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType.d]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax *= ((FLOAT)1.0 - InpParType.ar);

                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationKNN(OutParType->n, InpParType.d, Y, InpParType.RMIN, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < OutParType->n; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < OutParType->n; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                } 

                /* Component mean and variance calculation. */

                Error = MeanVarianceCalculation(InpParType.d, Theta[l], Mean[l], Variance[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < OutParType->n; j++) Y[j][InpParType.d] = R[j];

                Stop = (c >= InpParType.k[i]) || (c >= InpParType.cmax) || (c * Dmin < InpParType.D);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */
            
            Error = BayesClassificationKNN(OutParType->n, InpParType.d, Y, c, W, Theta, Mean, Variance);

            if (Error) goto E0;

            for (j = 0; j < OutParType->n; j++) Y[j][InpParType.d] = (FLOAT)1.0;

            Error = InformationCriterionKNN(InpParType.ICType, OutParType->n, InpParType.d, Y, c, W, Theta, &IC, &logL);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                OutParType->k = InpParType.k[i];

                memcpy(OutParType->h, h, InpParType.d * sizeof(FLOAT));  
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->c = c; 

                memcpy(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memcpy(OutParType->Theta[j], Theta[j], InpParType.d * sizeof(MarginalDistributionType));  
                }
            }

            #if(_DEBUG)
            fprintf(fp0, "%s\t%d\t%s\t%E\t%s\t%E\t%s\t%E\n", "c", c, "IC", IC, "Dmin", Dmin, "logL", logL);
            #endif 

            Dmin = Dmin * c / (c + 1); J++;
        }
        while (!Stop);
    }

E0: OutParType->W = (FLOAT*)realloc(OutParType->W, OutParType->c * sizeof(FLOAT));

    if (OutParType->Theta) {
        for (i = OutParType->c; i < InpParType.cmax; i++) {
            if (OutParType->Theta[i]) free(OutParType->Theta[i]);
        }

        OutParType->Theta = (MarginalDistributionType**)realloc(OutParType->Theta, OutParType->c * sizeof(MarginalDistributionType*));
    }
    
    if (Variance) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Variance[i]) free(Variance[i]);
        }
         
        free(Variance);
    }

    if (Mean) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Mean[i]) free(Mean[i]);
        }
         
        free(Mean);
    }

    if (Theta) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if(ymax) free(ymax);

    if(ymin) free(ymin);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    #if(_DEBUG)
    if (fp0) fclose(fp0);
    #endif

    return (Error);
} /* REBMIXKNN */

/* REBMIX algorithm for Parzen window. */

int REBMIXPW(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
             OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    FLOAT                    **Y = NULL;
    FLOAT                    *h = NULL, *ymin = NULL, *ymax = NULL;
    FLOAT                    *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                    *W = NULL;
    MarginalDistributionType **Theta = NULL; 
    FLOAT                    **Mean = NULL, **Variance = NULL;
    int                      c = 0, i, I, j, J, l, m;
    FLOAT                    V, Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL;
    int                      Error = 0, Stop = 0;
    #if(_DEBUG)
    int                      o; 
    FILE                     *fp0 = NULL;
    #endif

    #if(_DEBUG)
    if ((fp0 = fopen("debug.txt", "w")) == NULL) {
        Error = 1; goto E0;
    }
    #endif

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType.d + 2) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;

        for (j = 0; j < InpParType.d; j++) Y[i][j] = OutParType->X[i][j];
    }

    /* Allocation and normalizing vector h calculation. */

    h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    ymin = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    ymax = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        ymin[i] = ymax[i] = OutParType->X[0][i];

        for (j = 1; j < OutParType->n; j++) {
            if (OutParType->X[j][i] < ymin[i]) ymin[i] = OutParType->X[j][i];
            if (OutParType->X[j][i] > ymax[i]) ymax[i] = OutParType->X[j][i];
        }
    }

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].ParametricFamily = InpParType.ParFamType[j];
        }

        if (InpParType.Par0 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter0 = InpParType.Par0[j];
        }

        if (InpParType.Par1 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter1 = InpParType.Par1[j];
        }
    }

    Mean = (FLOAT**)malloc(InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Mean; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Mean[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Mean[i]; if (Error) goto E0;
    }

    Variance = (FLOAT**)realloc(Variance, InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Variance; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Variance[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Variance[i]; if (Error) goto E0;
    }

    OutParType->IC = FLOAT_MAX; OutParType->W = NULL; OutParType->Theta = NULL;

    OutParType->h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    for (i = 0; i < InpParType.kmax; i++) {
        /* Preprocessing of observations. */

        V = (FLOAT)1.0;
        
        for (j = 0; j < InpParType.d; j++) {
            switch (InpParType.ParFamType[j]) {
            case pfNormal: case pfLognormal: case pfWeibull:
                h[j] = (ymax[j] - ymin[j]) / InpParType.k[i]; V *= h[j]; 

                break;
            case pfBinomial:
                h[j] = (FLOAT)1.0; V *= h[j];
            }
        }

        Error = PreprocessingPW(h, OutParType->n, InpParType.d, Y);

        if (Error) { 
            Error = 0; goto E1;
        }

        #if(_DEBUG)
        fprintf(fp0, "%s\t%d\n%s", "k", InpParType.k[i], "h");

        for (j = 0; j < InpParType.d; j++) fprintf(fp0, "\t%E", h[j]);

        fprintf(fp0, "\n");

        for (j = 0; j < OutParType->n; j++) {
            fprintf(fp0, "%s%d%s", "Y[", j, "]");

            for (o = 0; o < InpParType.d + 2; o++) fprintf(fp0, "\t%E", Y[j][o]);

            fprintf(fp0, "\n");
        }
        #endif

        Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * (l + 1) * Dmin) {
                /* Global mode detection. */

                Error = GlobalModePW(h, ymin, &m, OutParType->n, InpParType.d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, OutParType->n * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) {
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationPW(OutParType->n, InpParType.d, Y, h, nl, m, Theta[l], InpParType.ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < OutParType->n; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType.d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType.d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType.d] - nl * fl * V / Y[j][InpParType.d + 1];

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType.d]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax *= ((FLOAT)1.0 - InpParType.ar);

                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationPW(OutParType->n, InpParType.d, Y, h, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < OutParType->n; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < OutParType->n; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                }

                /* Component mean and variance calculation. */

                Error = MeanVarianceCalculation(InpParType.d, Theta[l], Mean[l], Variance[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < OutParType->n; j++) Y[j][InpParType.d] = R[j];

                Stop = (c >= InpParType.k[i]) || (c >= InpParType.cmax) || (c * Dmin < InpParType.D);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */

            Error = BayesClassificationPW(OutParType->n, InpParType.d, Y, c, W, Theta, Mean, Variance);

            if (Error) goto E0;

            for (j = 0; j < OutParType->n; j++) Y[j][InpParType.d] = (FLOAT)1.0;

            Error = InformationCriterionPW(InpParType.ICType, OutParType->n, InpParType.d, Y, c, W, Theta, &IC, &logL);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                OutParType->k = InpParType.k[i];

                memcpy(OutParType->h, h, InpParType.d * sizeof(FLOAT));  
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->c = c; 

                memcpy(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memcpy(OutParType->Theta[j], Theta[j], InpParType.d * sizeof(MarginalDistributionType));  
                }
            }

            #if(_DEBUG)
            fprintf(fp0, "%s\t%d\t%s\t%E\t%s\t%E\t%s\t%E\n", "c", c, "IC", IC, "Dmin", Dmin, "logL", logL);
            #endif  

            Dmin = Dmin * c / (c + 1); J++;
        }
        while (!Stop);
E1:; 
    }

E0: OutParType->W = (FLOAT*)realloc(OutParType->W, OutParType->c * sizeof(FLOAT));

    if (OutParType->Theta) {
        for (i = OutParType->c; i < InpParType.cmax; i++) {
            if (OutParType->Theta[i]) free(OutParType->Theta[i]);
        }

        OutParType->Theta = (MarginalDistributionType**)realloc(OutParType->Theta, OutParType->c * sizeof(MarginalDistributionType*));
    }
    
    if (Variance) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Variance[i]) free(Variance[i]);
        }
         
        free(Variance);
    }

    if (Mean) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Mean[i]) free(Mean[i]);
        }
         
        free(Mean);
    }

    if (Theta) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if(ymax) free(ymax);

    if(ymin) free(ymin);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    #if(_DEBUG)
    if (fp0) fclose(fp0);
    #endif

    return (Error);
} /* REBMIXPW */

/* REBMIX algorithm for histogram. */

int REBMIXH(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
            OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    FLOAT                    **Y = NULL;
    FLOAT                    *h = NULL, *y0 = NULL, *ymin = NULL, *ymax = NULL;
    FLOAT                    *R = NULL, *E = NULL, *Epsilon = NULL;
    FLOAT                    *K = NULL;
    FLOAT                    *W = NULL;
    MarginalDistributionType **Theta = NULL; 
    FLOAT                    **Mean = NULL, **Variance = NULL;
    int                      c = 0, i, I, j, J, k, l, m;
    FLOAT                    V, Dmin, r, nl, elp, eln, epsilonlmax, fl, Dl, f, IC, logL;
    int                      Error = 0, Stop = 0;
    #if(_DEBUG)
    int                      o;
    FILE                     *fp0 = NULL;
    #endif

    #if(_DEBUG)
    if ((fp0 = fopen("debug.txt", "w")) == NULL) {
        Error = 1; goto E0;
    }
    #endif

    /* Allocation and initialisation. */

    Y = (FLOAT**)malloc(OutParType->n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < OutParType->n; i++) {
        Y[i] = (FLOAT*)malloc((InpParType.d + 1) * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;
    }

    /* Allocation and normalizing vector h calculation. */

    h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == h; if (Error) goto E0;

    y0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == y0; if (Error) goto E0;

    ymin = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymin; if (Error) goto E0;

    ymax = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == ymax; if (Error) goto E0;

    for (i = 0; i < InpParType.d; i++) {
        ymin[i] = ymax[i] = OutParType->X[0][i];

        for (j = 1; j < OutParType->n; j++) {
            if (OutParType->X[j][i] < ymin[i]) ymin[i] = OutParType->X[j][i];
            if (OutParType->X[j][i] > ymax[i]) ymax[i] = OutParType->X[j][i];
        }
    }

    R = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == R; if (Error) goto E0;

    E = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == E; if (Error) goto E0;

    Epsilon = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == Epsilon; if (Error) goto E0;

    K = (FLOAT*)malloc(OutParType->n * sizeof(FLOAT));

    Error = NULL == K; if (Error) goto E0;

    W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == W; if (Error) goto E0;

    Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == Theta[i]; if (Error) goto E0;

        for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].ParametricFamily = InpParType.ParFamType[j];
        }

        if (InpParType.Par0 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter0 = InpParType.Par0[j];
        }

        if (InpParType.Par1 != NULL) for (j = 0; j < InpParType.d; j++) {
            Theta[i][j].Parameter1 = InpParType.Par1[j];
        }
    }

    Mean = (FLOAT**)malloc(InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Mean; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Mean[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Mean[i]; if (Error) goto E0;
    }

    Variance = (FLOAT**)realloc(Variance, InpParType.cmax * sizeof(FLOAT*));

    Error = NULL == Variance; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        Variance[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Variance[i]; if (Error) goto E0;
    }

    OutParType->IC = FLOAT_MAX; OutParType->W = NULL; OutParType->Theta = NULL;

    OutParType->h = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->h; if (Error) goto E0;

    OutParType->y0 = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

    Error = NULL == OutParType->y0; if (Error) goto E0;

    OutParType->W = (FLOAT*)malloc(InpParType.cmax * sizeof(FLOAT));

    Error = NULL == OutParType->W; if (Error) goto E0;

    OutParType->Theta = (MarginalDistributionType**)malloc(InpParType.cmax * sizeof(MarginalDistributionType*));

    Error = NULL == OutParType->Theta; if (Error) goto E0;

    for (i = 0; i < InpParType.cmax; i++) {
        OutParType->Theta[i] = (MarginalDistributionType*)malloc(InpParType.d * sizeof(MarginalDistributionType));

        Error = NULL == OutParType->Theta[i]; if (Error) goto E0;
    }

    for (i = 0; i < InpParType.kmax; i++) {
        /* Preprocessing of observations. */

        k = InpParType.k[i]; V = (FLOAT)1.0; 
        
        for (j = 0; j < InpParType.d; j++) {
            switch (InpParType.ParFamType[j]) {
            case pfNormal: case pfLognormal: case pfWeibull:
                h[j] = (ymax[j] - ymin[j]) / InpParType.k[i]; y0[j] = ymin[j] + (FLOAT)0.5 * h[j]; V *= h[j]; 

                break;
            case pfBinomial:
                h[j] = (FLOAT)1.0; y0[j] = ymin[j]; V *= h[j];
            }
        }

        Error = PreprocessingH(h, y0, InpParType.ParFamType, &InpParType.k[i], OutParType->n, InpParType.d, OutParType->X, Y);

        if (Error) { 
            Error = 0; goto E1;
        }

        #if(_DEBUG)
        fprintf(fp0, "%s\t%d\n%s", "k", k, "h");

        for (j = 0; j < InpParType.d; j++) fprintf(fp0, "\t%E", h[j]);

        fprintf(fp0, "\n");

        for (j = 0; j < InpParType.k[i]; j++) {
            fprintf(fp0, "%s%d%s", "Y[", j, "]");

            for (o = 0; o < InpParType.d + 1; o++) fprintf(fp0, "\t%E", Y[j][o]);

            fprintf(fp0, "\n");
        }
        #endif

        for (j = 0; j < InpParType.k[i]; j++) K[j] = Y[j][InpParType.d];

        Dmin = (FLOAT)0.25; J = 1;

        /* Outer loop. */

        do {
            l = 0; r = (FLOAT)OutParType->n; nl = (FLOAT)OutParType->n;

            /* Middle loop. */

            while (nl / OutParType->n > (FLOAT)2.0 * (l + 1) * Dmin) {
                /* Global mode detection. */

                Error = GlobalModeH(h, ymin, &m, InpParType.k[i], InpParType.d, Y);

                if (Error) goto E0;

                I = 1; W[l] = nl / OutParType->n; memset(R, 0, InpParType.k[i] * sizeof(FLOAT));

                /* Inner loop. */

                while (I <= ItMax) { 
                    /* Rough component parameter estimation. */

                    Error = RoughEstimationH(InpParType.k[i], InpParType.d, Y, h, nl, m, Theta[l], InpParType.ResType);

                    if (Error) goto E0;

                    elp = eln = epsilonlmax = (FLOAT)0.0;

                    for (j = 0; j < InpParType.k[i]; j++) {
                        E[j] = Epsilon[j] = (FLOAT)0.0;

                        if ((Y[j][InpParType.d] > FLOAT_MIN) || (R[j] > FLOAT_MIN)) {
                            Error = ComponentDist(InpParType.d, Y[j], Theta[l], &fl, 0);

                            if (Error) goto E0;

                            E[j] = Y[j][InpParType.d] - nl * fl * V;

                            if (E[j] > (FLOAT)0.0) {
                                Epsilon[j] = E[j] / Y[j][InpParType.d]; 
                                
                                if (Epsilon[j] > epsilonlmax) epsilonlmax = Epsilon[j]; 
                                
                                elp += E[j];
                            }
                            else {
                                if (E[j] < -R[j]) E[j] = -R[j]; eln -= E[j];
                            }
                        }
                    }

                    Dl = elp / nl; epsilonlmax *= ((FLOAT)1.0 - InpParType.ar);

                    if (Dl <= Dmin / W[l]) {
                        /* Enhanced component parameter estimation. */

                        EnhancedEstimationH(InpParType.k[i], InpParType.d, Y, h, nl, Theta[l]);

                        break;
                    }
                    else {
                        for (j = 0; j < InpParType.k[i]; j++) if (Epsilon[j] > epsilonlmax) {
                            Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                        }

                        if (eln > FLOAT_MIN) {
                            elp = elp / Dl - nl; if (eln > elp) f = elp / eln; else f = (FLOAT)1.0;

                            for (j = 0; j < InpParType.k[i]; j++) if (E[j] < (FLOAT)0.0) {
                                E[j] *= f; Y[j][InpParType.d] -= E[j]; R[j] += E[j]; nl -= E[j];
                            }
                        }

                        W[l] = nl / OutParType->n;
                    }

                    I++;
                }

                /* Component mean and variance calculation. */

                Error = MeanVarianceCalculation(InpParType.d, Theta[l], Mean[l], Variance[l]);

                if (Error) goto E0;

                c = ++l;

                r -= nl; nl = r; for (j = 0; j < InpParType.k[i]; j++) Y[j][InpParType.d] = R[j];

                Stop = (c >= InpParType.k[i]) || (c >= InpParType.cmax) || (c * Dmin < InpParType.D);

                if (Stop) break;
            }

            /* Bayes classification of the remaining observations. */

            Error = BayesClassificationH(InpParType.k[i], OutParType->n, InpParType.d, Y, c, W, Theta, Mean, Variance);

            if (Error) goto E0;

            for (j = 0; j < InpParType.k[i]; j++) Y[j][InpParType.d] = K[j];

            Error = InformationCriterionH(InpParType.ICType, InpParType.k[i], OutParType->n, InpParType.d, Y, c, W, Theta, &IC, &logL);
            
            if (Error) goto E0;

            if (IC < OutParType->IC) {
                OutParType->k = k;

                memcpy(OutParType->h, h, InpParType.d * sizeof(FLOAT));

                memcpy(OutParType->y0, y0, InpParType.d * sizeof(FLOAT)); 
                
                OutParType->IC = IC; OutParType->logL = logL; OutParType->c = c; 

                memcpy(OutParType->W, W, c * sizeof(FLOAT));  

                for (j = 0; j < c; j++) {
                    memcpy(OutParType->Theta[j], Theta[j], InpParType.d * sizeof(MarginalDistributionType));  
                }
            }

            #if(_DEBUG)
            fprintf(fp0, "%s\t%d\t%s\t%E\t%s\t%E\t%s\t%E\n", "c", c, "IC", IC, "Dmin", Dmin, "logL", logL);
            #endif  

            Dmin = Dmin * c / (c + 1); J++;
        }
        while (!Stop);

E1:     InpParType.k[i] = k;
    }

E0: OutParType->W = (FLOAT*)realloc(OutParType->W, OutParType->c * sizeof(FLOAT));

    if (OutParType->Theta) {
        for (i = OutParType->c; i < InpParType.cmax; i++) {
            if (OutParType->Theta[i]) free(OutParType->Theta[i]);
        }

        OutParType->Theta = (MarginalDistributionType**)realloc(OutParType->Theta, OutParType->c * sizeof(MarginalDistributionType*));
    }
    
    if (Variance) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Variance[i]) free(Variance[i]);
        }
         
        free(Variance);
    }

    if (Mean) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Mean[i]) free(Mean[i]);
        }
         
        free(Mean);
    }

    if (Theta) {
        for (i = 0; i < InpParType.cmax; i++) {
            if (Theta[i]) free(Theta[i]);
        }
         
        free(Theta);
    }

    if (W) free(W);

    if (K) free(K);

    if (Epsilon) free(Epsilon);

    if (E) free(E);

    if (R) free(R);

    if(ymax) free(ymax);

    if(ymin) free(ymin);

    if(y0) free(y0);

    if (h) free(h);

    if (Y) {
        for (i = 0; i < OutParType->n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    #if(_DEBUG)
    if (fp0) fclose(fp0);
    #endif

    return (Error);
} /* REBMIXH */

/* Reads input data from the file stream. */

int ReadREBMIXDataFile(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
                       OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    char line[65536];
    char *pchar = NULL;
    int  i, j, BufSize = 0;
    FILE *fp = NULL;
    int  Error = 0;

    memset(OutParType, 0, sizeof(OutputREBMIXParameterType));

    if ((fp = fopen(InpParType.curr, "r")) == NULL) {
        Error = 1; goto E0;
    }

    OutParType->X = (FLOAT**)malloc((BufSize + BufInc) * sizeof(FLOAT*));

    for (i = BufSize; i < BufSize + BufInc; i++) {
        OutParType->X[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));
    }

    BufSize += BufInc;

    OutParType->n = 0;

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n");

        if (!pchar) goto S0;

        j = 0;
        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] == ',') {
                line[j] = '.'; j++;
            }
            else
            if (pchar[i] != ' ') {
                line[j] = pchar[i]; j++;
            }
        }

        line[j] = '\0';

        if (!j) goto S0;

        pchar = strtok(pchar, "\t");

        if (OutParType->n == BufSize) {
            OutParType->X = (FLOAT**)realloc(OutParType->X, (BufSize + BufInc) * sizeof(FLOAT*));

            for (i = BufSize; i < BufSize + BufInc; i++) {
                OutParType->X[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));
            }

            BufSize += BufInc;
        }

        i = 0;
        while (pchar) {
            OutParType->X[OutParType->n][i] = (FLOAT)atof(pchar); pchar = strtok(NULL, "\t"); i++;
        }

        OutParType->n++;
    }

    for (i = OutParType->n; i < BufSize; i++) {
        if (OutParType->X[i]) free(OutParType->X[i]); 
    }

    OutParType->X = (FLOAT**)realloc(OutParType->X, OutParType->n * sizeof(FLOAT*));

E0: if (fp) fclose(fp);

    return (Error);
} /* ReadREBMIXDataFile */

/* Writes input and output parameters into the file stream. */

int WriteREBMIXParameterFile(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
                             OutputREBMIXParameterType OutParType)  /* Output parameters. */
{
    int  i, j;
    char line[65536];
    char mode[2];
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp0 = NULL, *fp1 = NULL;
    int  Error = 0;

    if (InpParType.curr == InpParType.open[0])
        strcpy(mode, "w");
    else
        strcpy(mode, "a");

    strcpy(path, InpParType.save); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_1", ext);

    if ((fp0 = fopen(path, mode)) == NULL) {
        Error = 1; goto E0;
    }

    strcpy(path, InpParType.save); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_2", ext);

    if ((fp1 = fopen(path, mode)) == NULL) {
        Error = 1; goto E0;
    }

    if (!strcmp(mode, "w")) {
        fprintf(fp0, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "Dataset",
                                                       "Preprocessing",
                                                       "D",
                                                       "cmax",
                                                       "InformationCriterion",
                                                       "ar",
                                                       "Restraints", 
                                                       "c");

        switch (InpParType.PreType) {
        case poHistogram:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < InpParType.d; i++) {
                if (InpParType.d == 1)
                    sprintf(line, "%s", "y0");
                else
                    sprintf(line, "%s%d", "y0", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            for (i = 0; i < InpParType.d; i++) {
                if (InpParType.d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poParzenWindow:
            fprintf(fp0, "\t%s", "k");

            for (i = 0; i < InpParType.d; i++) {
                if (InpParType.d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }

            break;
        case poKNearestNeighbour:
            fprintf(fp0, "\t%s\t%s", "k",
                                     "Rmin");

            for (i = 0; i < InpParType.d; i++) {
                if (InpParType.d == 1)
                    sprintf(line, "%s", "h");
                else
                    sprintf(line, "%s%d", "h", i + 1);
                  
                fprintf(fp0, "\t%s", line);
            }
        }

        fprintf(fp0, "\t%s\t%s\t%s\n", "tc",
                                       "IC",
                                       "logL");

        fprintf(fp1, "%s\t%s", "Dataset",
                               "w");

        for (i = 0; i < InpParType.d; i++) {
            if (InpParType.d == 1)
                fprintf(fp1, "\t%s\t%s\t%s", "pdf", "theta1", "theta2");
            else
                fprintf(fp1, "\t%s%d\t%s%d\t%s%d", "pdf", i + 1, "theta1.", i + 1, "theta2.", i + 1);
        }

        fprintf(fp1, "\n");
    }

    strcpy(path, InpParType.curr); 

    pchar = strrchr(path, '\\');

    if (!pchar) {
        pchar = strrchr(path, '/');
    }

    if (pchar) {
        strcpy(path, pchar + 1);
    }

    pchar = strrchr(path, '.'); 
        
    if (pchar) pchar[0] = '\0';

    fprintf(fp0, "%s", path);

    switch (InpParType.PreType) {
    case poHistogram: 
        strcpy(line, "histogram");

        break;
    case poParzenWindow: 
        strcpy(line, "Parzen window");

        break;
    case poKNearestNeighbour:
        strcpy(line, "k-nearest neighbour");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%E\t%d", InpParType.D,
                             InpParType.cmax);

    switch (InpParType.ICType) {
    case icAIC:
        strcpy(line, "AIC");

        break;
    case icAIC3:
        strcpy(line, "AIC3");

        break;
    case icAIC4:
        strcpy(line, "AIC4");

        break;
    case icAICc:
        strcpy(line, "AICc");

        break;
    case icBIC:
        strcpy(line, "BIC");

        break;
    case icCAIC:
        strcpy(line, "CAIC");

        break;
    case icHQC:
        strcpy(line, "HQC");

        break;
    case icMDL2:
        strcpy(line, "MDL2");

        break;
    case icMDL5:
        strcpy(line, "MDL5");

        break;
    case icAWE:
        strcpy(line, "AWE");

        break;
    case icCLC:
        strcpy(line, "CLC");

        break;
    case icICL:
        strcpy(line, "ICL");

        break;
    case icPC:
        strcpy(line, "PC");

        break;
    case icICLBIC:
        strcpy(line, "ICL-BIC");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%E", InpParType.ar);

    switch (InpParType.ResType) {
    case rtRigid:
        strcpy(line, "rigid");

        break;
    case rtLoose:
        strcpy(line, "loose");
    }

    fprintf(fp0, "\t%s", line);

    fprintf(fp0, "\t%d", OutParType.c);

    switch (InpParType.PreType) {
    case poHistogram:
        fprintf(fp0, "\t%d", OutParType.k);

        for (i = 0; i < InpParType.d; i++) {
            fprintf(fp0, "\t%E", OutParType.y0[i]);
        }

        for (i = 0; i < InpParType.d; i++) {
            fprintf(fp0, "\t%E", OutParType.h[i]);
        }

        break;
    case poParzenWindow:
        fprintf(fp0, "\t%d", OutParType.k);

        for (i = 0; i < InpParType.d; i++) {
            fprintf(fp0, "\t%E", OutParType.h[i]);
        }

        break;
    case poKNearestNeighbour:
        fprintf(fp0, "\t%d\t%E", OutParType.k,
                                 InpParType.RMIN);

        for (i = 0; i < InpParType.d; i++) {
            fprintf(fp0, "\t%E", OutParType.h[i]);
        }
    }

    fprintf(fp0, "\t%d\t%E\t%E\n", (int)OutParType.ClcTime,
                                   OutParType.IC,
                                   OutParType.logL);

    for (i = 0; i < OutParType.c; i++) {
        fprintf(fp1, "%s\t%E", path,
                               OutParType.W[i]);

        for (j = 0; j < InpParType.d; j++) switch (OutParType.Theta[i][j].ParametricFamily) {
            case pfNormal:
                fprintf(fp1, "\t%s\t%E\t%E", "normal", 
                                             OutParType.Theta[i][j].Parameter0,
                                             OutParType.Theta[i][j].Parameter1);

                break;
            case pfLognormal:
                fprintf(fp1, "\t%s\t%E\t%E", "lognormal", 
                                             OutParType.Theta[i][j].Parameter0,
                                             OutParType.Theta[i][j].Parameter1);
                break;
            case pfWeibull:
                fprintf(fp1, "\t%s\t%E\t%E", "Weibull", 
                                             OutParType.Theta[i][j].Parameter0,
                                             OutParType.Theta[i][j].Parameter1);
                break;
            case pfBinomial:
                fprintf(fp1, "\t%s\t%E\t%E", "binomial", 
                                             OutParType.Theta[i][j].Parameter0,
                                             OutParType.Theta[i][j].Parameter1);
        }

        fprintf(fp1, "\n");
    }

E0: if (fp0) fclose(fp0);
    if (fp1) fclose(fp1);

    if (OutParType.X) {
        for (i = 0; i < OutParType.n; i++) {
            if (OutParType.X[i]) free(OutParType.X[i]);
        }
         
        free(OutParType.X);
    }

    if (OutParType.Theta) {
        for (i = 0; i < OutParType.c; i++) {
            if (OutParType.Theta[i]) free(OutParType.Theta[i]);
        }
         
        free(OutParType.Theta);
    }

    if (OutParType.W) free(OutParType.W);

    if (OutParType.y0) free(OutParType.y0);

    if (OutParType.h) free(OutParType.h);

    return (Error);
} /* WriteREBMIXParameterFile */

/* REBMIX algorithm. */

int REBMIX(InputREBMIXParameterType  InpParType,   /* Input parameters. */ 
           OutputREBMIXParameterType *OutParType)  /* Output parameters. */
{
    long StartTime, EndTime;
    int  Error = 0;

    StartTime = clock();

    switch (InpParType.PreType) {
    case poHistogram:
        Error = REBMIXH(InpParType, OutParType);

        if (Error) goto E0;

        break;
    case poParzenWindow:
        Error = REBMIXPW(InpParType, OutParType);

        if (Error) goto E0;

        break;
    case poKNearestNeighbour:
        Error = REBMIXKNN(InpParType, OutParType);

        if (Error) goto E0;
    }

    EndTime = clock(); OutParType->ClcTime = 1000 * (EndTime - StartTime) / CLOCKS_PER_SEC;

E0: return (Error);
} /* REBMIX */

/* Runs REBMIX template file stream. */

int RunREBMIXTemplateFile(char *file)
{
    char                      line[65536], ident[65536], list[65536];
    char                      *pchar = NULL, *rchar = NULL;
    int                       i, imin, imax, iinc, j, k, isI;
    FLOAT                     isF;
    InputREBMIXParameterType  InpParType;
    OutputREBMIXParameterType OutParType;
    FILE                      *fp = NULL;
    int                       Error = 0;

    memset(&InpParType, 0, sizeof(InputREBMIXParameterType));

    /* Recommended values. */

    InpParType.ar = (FLOAT)0.1;
    InpParType.RMIN = (FLOAT)0.001;
    InpParType.ResType = rtLoose;

    memset(&OutParType, 0, sizeof(OutputREBMIXParameterType));

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

    #if (_REBMIXEXE)
    printf("REBMIX Version 2.2.1\n");
    #elif (_REBMIXR)
    Rprintf("REBMIX Version 2.2.1\n");
    R_FlushConsole();
    #endif

S0: while (fgets(line, 2048, fp) != NULL) {
        pchar = strtok(line, "\n"); 
        
        pchar = strtok(pchar, "=");

        if (!pchar) goto S0;

        j = 0;

        for (i = 0; i < (int)strlen(pchar); i++) {
            if (pchar[i] != ' ') {
                ident[j] = (char)toupper(pchar[i]); j++;
            }
        }

        ident[j] = '\0';

        j = 0; list[j] = '\0'; imin = 0; imax = 0;

        while((pchar = strtok(NULL, ",")) != NULL) {
            if (!strcmp(ident, "DATASET") || !strcmp(ident, "SAVE")) {
                for (i = 0; pchar[i] != '\0'; i++) {
                    list[j] = pchar[i]; j++;

                    if (pchar[i] != ' ') {
                        imax = i; if (!imin) imin = i;
                    }
                }

                j = imax + 1 - imin;

                for (i = 0; i < j; i++) {
                    list[i] = list[imin + i];
                }
            }
            else {
                for (i = 0; pchar[i] != '\0'; i++) if (pchar[i] != '[' && pchar[i] != ']') {
                    if (pchar[i] == ' ') {
                    }
                    else {
                        list[j] = (char)toupper(pchar[i]); j++;
                    }
                }
            }

            list[j] = '\t'; j++;
        }

        if (!j) goto S0; else list[j - 1] = '\0';

        pchar = strtok(list, "\t");

        if (!strcmp(ident, "RUN")) {
            for (k = 0; k < InpParType.o; k++) {
                InpParType.curr = InpParType.open[k];

                Error = ReadREBMIXDataFile(InpParType, &OutParType);

                if (Error) goto E0;

                #if (_REBMIXEXE)
                printf("Dataset = %s\n", InpParType.curr);
                #elif (_REBMIXR)
                Rprintf("Dataset = %s\n", InpParType.curr);
                R_FlushConsole();
                #endif

                Error = REBMIX(InpParType, &OutParType);

                if (Error) goto E0;

                Error = WriteREBMIXParameterFile(InpParType, OutParType);

                if (Error) goto E0;
            }
        }
        else
        if (!strcmp(ident, "DATASET")) {
            InpParType.open = (char**)realloc(InpParType.open, (InpParType.o + 1) * sizeof(char*));

            Error = NULL == InpParType.open; if (Error) goto E0;

            InpParType.open[InpParType.o] = (char*)malloc((strlen(pchar) + 1) * sizeof(char));

            Error = NULL == InpParType.open[InpParType.o]; if (Error) goto E0;

            strcpy(InpParType.open[InpParType.o], pchar); InpParType.o++;
        }
        else
        if (!strcmp(ident, "PREPROCESSING")) {
            if (!strcmp(pchar, "HISTOGRAM"))
                InpParType.PreType = poHistogram;
            else
            if (!strcmp(pchar, "PARZENWINDOW"))
                InpParType.PreType = poParzenWindow;
            else
            if (!strcmp(pchar, "K-NEARESTNEIGHBOUR"))
                InpParType.PreType = poKNearestNeighbour;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "D")) {
            InpParType.D = isF = (FLOAT)atof(pchar);

            Error = (isF <= (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;
        } else
        if (!strcmp(ident, "CMAX")) {
            InpParType.cmax = isI = (int)atol(pchar);

            Error = isI <= 0; if (Error) goto E0;
        } else
        if (!strcmp(ident, "INFORMATIONCRITERION")) {
            if (!strcmp(pchar, "AIC"))
                InpParType.ICType = icAIC;
            else
            if (!strcmp(pchar, "AIC3"))
                InpParType.ICType = icAIC3;
            else
            if (!strcmp(pchar, "AIC4"))
                InpParType.ICType = icAIC4;
            else
            if (!strcmp(pchar, "AICC"))
                InpParType.ICType = icAICc;
            else
            if (!strcmp(pchar, "BIC"))
                InpParType.ICType = icBIC;
            else
            if (!strcmp(pchar, "CAIC"))
                InpParType.ICType = icCAIC;
            else
            if (!strcmp(pchar, "HQC"))
                InpParType.ICType = icHQC;
            else
            if (!strcmp(pchar, "MDL2"))
                InpParType.ICType = icMDL2;
            else
            if (!strcmp(pchar, "MDL5"))
                InpParType.ICType = icMDL5;
            else
            if (!strcmp(pchar, "AWE"))
                InpParType.ICType = icAWE;
            else
            if (!strcmp(pchar, "CLC"))
                InpParType.ICType = icCLC;
            else
            if (!strcmp(pchar, "ICL"))
                InpParType.ICType = icICL;
            else
            if (!strcmp(pchar, "PC"))
                InpParType.ICType = icPC;
            else
            if (!strcmp(pchar, "ICL-BIC"))
                InpParType.ICType = icICLBIC;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "PDF")) {
            i = 0;

            while (pchar) {
                InpParType.ParFamType = (ParametricFamilyType_e*)realloc(InpParType.ParFamType, (i + 1) * sizeof(ParametricFamilyType_e));

                Error = NULL == InpParType.ParFamType; if (Error) goto E0;

                if (!strcmp(pchar, "NORMAL"))
                    InpParType.ParFamType[i] = pfNormal; 
                else
                if (!strcmp(pchar, "LOGNORMAL"))
                    InpParType.ParFamType[i] = pfLognormal;
                else
                if (!strcmp(pchar, "WEIBULL"))
                    InpParType.ParFamType[i] = pfWeibull;
                else
                if (!strcmp(pchar, "BINOMIAL"))
                    InpParType.ParFamType[i] = pfBinomial;
                else {
                    Error = 1; goto E0;
                }
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "THETA1")) {
            i = 0;

            while (pchar) {
                InpParType.Par0 = (FLOAT*)realloc(InpParType.Par0, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.Par0; if (Error) goto E0;

                InpParType.Par0[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "THETA2")) {
            i = 0;

            while (pchar) {
                InpParType.Par1 = (FLOAT*)realloc(InpParType.Par1, (i + 1) * sizeof(FLOAT));

                Error = NULL == InpParType.Par1; if (Error) goto E0;

                InpParType.Par1[i] = (FLOAT)atof(pchar);
                
                pchar = strtok(NULL, "\t"); ++i;
            }

            if ((InpParType.d > 0) && (InpParType.d != i)) {
                Error = 1; goto E0;
            }
            else {
                InpParType.d = i;
            }
        } else
        if (!strcmp(ident, "K")) {
            i = 0;

            while (pchar != NULL) {
                if ((rchar = strrchr(pchar, '-')) != NULL) {
                    imin = (int)atol(pchar); imax = (int)atol(rchar + 1);

                    if (imin > imax) {
                        j = imin; imin = imax; imax = j;
                    }

                    if ((rchar = strrchr(pchar, '+')) != NULL)
                        iinc = (int)atol(rchar + 1);
                    else
                        iinc = 1;

                    InpParType.k = (int*)realloc(InpParType.k, (i + (imax - imin) / iinc + 1) * sizeof(int));

                    Error = NULL == InpParType.k; if (Error) goto E0;

                    for (j = imin; j <= imax; j += iinc) {
                        InpParType.k[i] = isI = j; 

                        Error = isI <= 0; if (Error) goto E0;
                        
                        i++;
                    }
                
                    InpParType.kmax = i;
                }
                else {
                    InpParType.k = (int*)realloc(InpParType.k, (i + 1) * sizeof(int));

                    Error = NULL == InpParType.k; if (Error) goto E0;

                    InpParType.k[i] = isI = (int)atol(pchar);

                    Error = isI <= 0; if (Error) goto E0;
                
                    InpParType.kmax = ++i;
                }

                pchar = strtok(NULL, "\t"); 
            }
        } else
        if (!strcmp(ident, "RMIN")) {
            InpParType.RMIN = isF = (FLOAT)atof(pchar);

            Error = (isF <= (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;
        } else
        if (!strcmp(ident, "AR")) {
            InpParType.ar = isF = (FLOAT)atof(pchar);

            Error = (isF <= (FLOAT)0.0) || (isF > (FLOAT)1.0); if (Error) goto E0;

        } else
        if (!strcmp(ident, "RESTRAINTS")) {
            if (!strcmp(pchar, "RIGID"))
                InpParType.ResType = rtRigid;
            else
            if (!strcmp(pchar, "LOOSE"))
                InpParType.ResType = rtLoose;
            else {
                Error = 1; goto E0;
            }
        } else
        if (!strcmp(ident, "SAVE")) {
            InpParType.save = (char*)realloc(InpParType.save, (strlen(pchar) + 1) * sizeof(char));

            Error = NULL == InpParType.save; if (Error) goto E0;

            strcpy(InpParType.save, pchar);
        }
    }

E0: if (fp) fclose(fp);

    if (InpParType.save) free(InpParType.save);

    if (InpParType.k) free(InpParType.k);

    if (InpParType.Par1) free(InpParType.Par1);

    if (InpParType.Par0) free(InpParType.Par0);

    if (InpParType.ParFamType) free(InpParType.ParFamType);
    
    if (InpParType.open) {
        for (i = 0; i < InpParType.o; i++) {
            if (InpParType.open[i]) free(InpParType.open[i]);
        }
         
        free(InpParType.open);
    }

    return (Error);
} /* RunREBMIXTemplateFile */
