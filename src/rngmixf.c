#include <math.h>
#include <stdio.h>
#include <ctype.h>

#include "rngmixf.h"

#if (_REBMIXR)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#endif

int   NDevISet = 0;
FLOAT NDevVSet = (FLOAT)0.0;
int   LDevISet = 0;
FLOAT LDevVSet = (FLOAT)0.0;
FLOAT Bn = (FLOAT)(-1.0), Bp = (FLOAT)(-1.0), Be, Bg, Bplog, Bpc, Bpclog;
long  IY = 0;
long  IV[NTAB];

/* Minimal random number generator of Park and Miller with Bays-Durham shuffle and added
   safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
   values). Call with IDum a negative integer to initialize; thereafter do not alter IDum
   between successive deviates in a sequence. RNMX should approximate the largest floating 
   value that is less than 1. See http://www.nrbook.com/a/bookcpdf/c7-1.pdf */

FLOAT Ran1(int *IDum)
{
    int   j, k;
    FLOAT Tmp;

    if (*IDum <= 0 || !IY) {
        *IDum = (-(*IDum) < 1) ? 1 : -(*IDum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = *IDum / IQ;

            *IDum = IA * (*IDum - k * IQ) - IR * k;

            if (*IDum < 0) *IDum += IM; 
            
            if (j < NTAB) IV[j] = *IDum;
        }

        IY = IV[0];
    }

    k = *IDum / IQ; *IDum = IA * (*IDum - k * IQ) - IR * k;

    if (*IDum < 0) *IDum += IM;

    j = IY / NDIV; IY = IV[j]; IV[j] = *IDum;

    if ((Tmp = AM * IY) > RNMX) return (RNMX); else return (Tmp);
} /* Ran1 */

/* Writes input parameters into the file stream. */

int WriteRNGMIXParameterFile(InputRNGMIXParameterType InpParType)  /* Input parameters. */ 
{
    char path[FILENAME_MAX];
    char ext[FILENAME_MAX];
    char *pchar = NULL;
    FILE *fp = NULL;
    int  Error = 0;

    strcpy(path, InpParType.save); 
        
    pchar = strrchr(path, '.'); 
        
    if (pchar) {
        strcpy(ext, pchar); pchar[0] = '\0';
    }
    else {
        strcpy(ext, "");
    }
        
    sprintf(path, "%s%s%s", path, "_1", ext);

    if ((fp = fopen(path, "w")) == NULL) {
        Error = 1; goto E0;
    }

    fprintf(fp, "%s\n", "rseed");

    fprintf(fp, "%d\n", InpParType.IDum);

E0: if (fp) fclose(fp);

    return (Error);
} /* WriteRNGMIXParameterFile */

/* Returns random sample of independent observations. */

int RNGMIX(InputRNGMIXParameterType InpParType) /* Input parameters. */
{
    FLOAT **Y = NULL;
    FLOAT C[8];
    FLOAT y, p;
    int   i, j, k, l, m, n = 0;
    FILE  *fp = NULL;
    int   Error = 0;

    if ((fp = fopen(InpParType.curr, "w")) == NULL) {
        Error = 1; goto E0;
    }

    n = 0; for (i = 0; i < InpParType.c; i++) n += InpParType.n[i];
    
    Y = (FLOAT**)malloc(n * sizeof(FLOAT*));

    Error = NULL == Y; if (Error) goto E0;

    for (i = 0; i < n; i++) {
        Y[i] = (FLOAT*)malloc(InpParType.d * sizeof(FLOAT));

        Error = NULL == Y[i]; if (Error) goto E0;
    }

    l = 0;
    for (i = 0; i < InpParType.c; i++) {
        for (j = 0; j < InpParType.n[i]; j++) {
            for (k = 0; k < InpParType.d; k++) {
                switch (InpParType.Theta[i][k].ParametricFamily) {
                case pfNormal:
                    if (NDevISet == 0) {
                        do {
                            C[0] = (FLOAT)2.0 * Ran1(&InpParType.IDum) - (FLOAT)1.0;
                            C[1] = (FLOAT)2.0 * Ran1(&InpParType.IDum) - (FLOAT)1.0;

                            C[2] = C[0] * C[0] + C[1] * C[1];
                        }
                        while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                        C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * log(C[2]) / C[2]);

                        y = C[3] * C[0];

                        NDevISet = 1; NDevVSet = C[3] * C[1];
                    }
                    else {
                        y = NDevVSet; NDevISet = 0;
                    }

                    Y[l][k] = InpParType.Theta[i][k].Parameter1 * y + 
                              InpParType.Theta[i][k].Parameter0;

                    break;
                case pfLognormal:
                    if (LDevISet == 0) {
                        do {
                            C[0] = (FLOAT)2.0 * Ran1(&InpParType.IDum) - (FLOAT)1.0;
                            C[1] = (FLOAT)2.0 * Ran1(&InpParType.IDum) - (FLOAT)1.0;

                            C[2] = C[0] * C[0] + C[1] * C[1];
                        }
                        while ((C[2] >= (FLOAT)1.0) || (C[2] == (FLOAT)0.0));

                        C[3] = (FLOAT)sqrt(-(FLOAT)2.0 * log(C[2]) / C[2]);

                        y = C[3] * C[0];

                        LDevISet = 1; LDevVSet = C[3] * C[1];
                    }
                    else {
                        y = LDevVSet; LDevISet = 0;
                    }

                    Y[l][k] = (FLOAT)exp(InpParType.Theta[i][k].Parameter1 * y + 
                                         InpParType.Theta[i][k].Parameter0);

                    break;
                case pfWeibull:
                    Y[l][k] = InpParType.Theta[i][k].Parameter0 * 
                              (FLOAT)exp(log(log((FLOAT)1.0 / Ran1(&InpParType.IDum))) / 
                              InpParType.Theta[i][k].Parameter1);

                    break;
                case pfBinomial:
                    if (InpParType.Theta[i][k].Parameter1 < (FLOAT)0.5) {
                        p = InpParType.Theta[i][k].Parameter1;
                    }
                    else {
                        p = (FLOAT)1.0 - InpParType.Theta[i][k].Parameter1;
                    }

                    C[0] = InpParType.Theta[i][k].Parameter0 * p;
                    if ((int)InpParType.Theta[i][k].Parameter0 < 25) {
                        Y[l][k] = (FLOAT)0.0;
                        
                        for (m = 0; m < (int)InpParType.Theta[i][k].Parameter0; m++) {
                            if (Ran1(&InpParType.IDum) < p) ++Y[l][k];
                        }
                    }
                    else 
                    if (C[0] < (FLOAT)1.0) {
                        C[1] = (FLOAT)exp(-C[0]); C[2] = (FLOAT)1.0;

                        for (m = 0; m < (int)InpParType.Theta[i][k].Parameter0; m++) {
                            C[2] *= Ran1(&InpParType.IDum); if (C[2] < C[1]) break;
                        }

                        if (m > (int)InpParType.Theta[i][k].Parameter0) {
                            Y[l][k] = InpParType.Theta[i][k].Parameter0;
                        }
                        else {
                            Y[l][k] = m;
                        }
                    }
                    else {
                        if (InpParType.Theta[i][k].Parameter0 != Bn) {
                            Be = InpParType.Theta[i][k].Parameter0;
                            Bg = Gammaln(Be + (FLOAT)1.0);
                            Bn = InpParType.Theta[i][k].Parameter0;
                        }

                        if (p != Bp) {
                            Bpc = (FLOAT)1.0 - p; 
                            Bplog = (FLOAT)log(p); 
                            Bpclog = (FLOAT)log(Bpc);
                            Bp = p;
                        }

                        C[3] = (FLOAT)sqrt((FLOAT)2.0 * C[0] * Bpc);

                        do {
                            do {
                                C[4] = Pi * Ran1(&InpParType.IDum);

                                C[5] = (FLOAT)tan(C[4]);

                                C[6] = C[3] * C[5] + C[0];
                            } while ((C[6] < (FLOAT)0.0) || (C[6] >= Be + (FLOAT)1.0));


                            C[6] = (FLOAT)floor(C[6]);

                            C[7] = (FLOAT)1.2 * C[3] * ((FLOAT)1.0 + C[5] * C[5]) *
                                   (FLOAT)exp(Bg - Gammaln(C[6] + (FLOAT)1.0) -
                                   Gammaln(Be - C[6] + (FLOAT)1.0) + 
                                   C[6] * Bplog + (Be - C[6]) * Bpclog);

                        } while (Ran1(&InpParType.IDum) > C[7]);

                        Y[l][k] = C[6];
                    }

                    if (p != InpParType.Theta[i][k].Parameter1) {
                        Y[l][k] = InpParType.Theta[i][k].Parameter0 - Y[l][k];
                    }
                }
            }

            fprintf(fp, "%E", Y[l][0]); 
        
            for (k = 1; k < InpParType.d; k++) fprintf(fp, "\t%E", Y[l][k]); 
        
            fprintf(fp, "\n");

            l++;
        }
    }

E0: if (fp) fclose(fp);
    
    if (Y) {
        for (i = 0; i < n; i++) {
            if (Y[i]) free(Y[i]);
        }
         
        free(Y);
    }

    return (Error);
} /* RNGMIX */

/* Runs RNGMIX template file stream. */

int RunRNGMIXTemplateFile(char *file)  /* File stream. */
{
    InputRNGMIXParameterType InpParType;
    int                      i, imin, imax, j, k, isI;
    FLOAT                    isF;
    char                     line[65536], ident[65536], list[65536];
    char                     *pchar = NULL;
    FILE                     *fp = NULL;
    int                      Error = 0;

    memset(&InpParType, 0, sizeof(InputRNGMIXParameterType));

    if ((fp = fopen(file, "r")) == NULL) {
        Error = 1; goto E0;
    }

    #if (_REBMIXEXE)
    printf("RNGMIX Version 2.2.1\n");
    #elif (_REBMIXR)
    Rprintf("RNGMIX Version 2.2.1\n");
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
            Error = WriteRNGMIXParameterFile(InpParType);

            if (Error) goto E0;

            for (k = 0; k < InpParType.o; k++) {
                InpParType.curr = InpParType.open[k]; 

                #if (_REBMIXEXE)
                printf("Dataset = %s\n", InpParType.curr);
                #elif (_REBMIXR)
                Rprintf("Dataset = %s\n", InpParType.curr);
                R_FlushConsole();
                #endif

                Error = RNGMIX(InpParType);

                if (Error) goto E0;

                InpParType.IDum--;
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
        if (!strcmp(ident, "RSEED")) {
            InpParType.IDum = isI = (int)atol(pchar);

            Error = isI >= 0; if (Error) goto E0;
        } else
        if (!strcmp(ident, "NTHETA")) {
            InpParType.n = (int*)realloc(InpParType.n, (InpParType.c + 1) * sizeof(int));

            Error = NULL == InpParType.n; if (Error) goto E0;

            InpParType.n[InpParType.c] = isI = (int)atol(pchar);

            Error = isI < 1; if (Error) goto E0;

            pchar = strtok(NULL, "\t"); 

            while (pchar) {
                InpParType.Theta = (MarginalDistributionType**)realloc(InpParType.Theta, (InpParType.c + 1) * sizeof(MarginalDistributionType*));

                Error = NULL == InpParType.Theta; if (Error) goto E0;

                InpParType.d = 0; InpParType.Theta[InpParType.c] = NULL;

                while (pchar) {
                    InpParType.Theta[InpParType.c] = (MarginalDistributionType*)realloc(InpParType.Theta[InpParType.c], (InpParType.d + 1) * sizeof(MarginalDistributionType));

                    Error = NULL == InpParType.Theta[InpParType.c]; if (Error) goto E0;

                    if (!strcmp(pchar, "NORMAL")) {
                        InpParType.Theta[InpParType.c][InpParType.d].ParametricFamily = pfNormal; 

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter0 = (FLOAT)atof(pchar);

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter1 = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "LOGNORMAL")) {
                        InpParType.Theta[InpParType.c][InpParType.d].ParametricFamily = pfLognormal;

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter0 = (FLOAT)atof(pchar);

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter1 = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "WEIBULL")) {
                        InpParType.Theta[InpParType.c][InpParType.d].ParametricFamily = pfWeibull; 

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter0 = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter1 = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;
                    }
                    else
                    if (!strcmp(pchar, "BINOMIAL")) {
                        InpParType.Theta[InpParType.c][InpParType.d].ParametricFamily = pfBinomial; 

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter0 = isF = (FLOAT)atof(pchar);

                        Error = isF <= (FLOAT)0.0; if (Error) goto E0;

                        pchar = strtok(NULL, "\t"); 

                        InpParType.Theta[InpParType.c][InpParType.d].Parameter1 = isF = (FLOAT)atof(pchar);

                        Error = (isF < (FLOAT)0.0) || (isF > (FLOAT)1.0) ; if (Error) goto E0;
                    }
                    else {
                        Error = 1; goto E0;
                    }

                    InpParType.d++;

                    pchar = strtok(NULL, "\t");
                }

                InpParType.c++;
            }
        } else
        if (!strcmp(ident, "SAVE")) {
            InpParType.save = (char*)realloc(InpParType.save, (strlen(pchar) + 1) * sizeof(char));

            Error = NULL == InpParType.save; if (Error) goto E0;

            strcpy(InpParType.save, pchar);
        }
    }

E0: if (fp) fclose(fp);

    if (InpParType.open) {
        for (i = 0; i < InpParType.o; i++) {
            if (InpParType.open[i]) free(InpParType.open[i]);
        }
         
        free(InpParType.open);
    }

    if (InpParType.n) free(InpParType.n);

    if (InpParType.Theta) {
        for (i = 0; i < InpParType.c; i++) {
            if (InpParType.Theta[i]) free(InpParType.Theta[i]);
        }

        free(InpParType.Theta);
    }

    if (InpParType.save) free(InpParType.save);

    return (Error);
} /* RunRNGMIXTemplateFile */
