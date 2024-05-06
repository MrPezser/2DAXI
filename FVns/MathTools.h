//
// Created by Tsail on 5/5/2024.
//

#ifndef FVNS_MATHTOOLS_H
#define FVNS_MATHTOOLS_H

#include "Indexing.h"

double *vecadd(  const double*a , const double*b, size_t n, double* vout );
double *vecaddself(    double*a , const double*b, size_t n);
double *vecsub(  const double*a , const double*b, size_t n, double* vout );
double *vecscale(const double *a, const double k, size_t n, double *out);
double *veccopy(double *a, const double *b, size_t n);
double vecinner( const double *a, const double *b, size_t n);
void vecinitialize( double* a, size_t n);
double norm22(   const double *vec, size_t n);
double norm2(   const double *vec, size_t n);
double normL2( const double *vec, const double*dx, size_t n);
double *matvec(  const double *A, const double *x, size_t n, double *out);
double vecmean( const double *vec, size_t n);

#endif //FVNS_MATHTOOLS_H
