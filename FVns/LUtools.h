//
// Created by tskoepli on 4/9/2024.
//

#ifndef FVNS_LUTOOLS_H
#define FVNS_LUTOOLS_H

void LUPDecompose(double **A, int N, double Tol, int *P);
void LUPSolve(double **A, int *P, double *b, int N, double *x);

#endif //FVNS_LUTOOLS_H
