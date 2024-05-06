//
// Created by Tsail on 5/5/2024.
//

#include "MathTools.h"

double *vecadd(  const double*a , const double*b, size_t n, double* vout ) {
    //be careful
    for (unsigned i=0; i < n; i++) {
        vout[i] = a[i] + b[i];
    }
    return vout;
}
double *vecaddself(    double*a , const double*b, size_t n) {
    //be careful
    for (unsigned i=0; i < n; i++) {
        a[i] += b[i];
    }
    return a;
}
double *vecsub(  const double*a , const double*b, size_t n, double* vout ) {
    //be careful
    for (unsigned i=0; i < n; i++) {
        vout[i] = a[i] - b[i];
    }
    return vout;
}
double *vecscale(const double *a, const double k, size_t n, double *out) {
    for (unsigned i=0;i<n;i++) {
        out[i] = k*a[i];
    }
    return out;
}
double *veccopy(double *a, const double *b, size_t n) {
    for (int i=0; i<n ; i++) {
        a[i] = b[i];
    }
    return a;
}
double vecinner( const double *a, const double *b, size_t n) {
    double out = 0.0;
    for (int i=0; i < n; i++) {
        out += a[i] * b[i];
    }
    return out;
}
void vecinitialize( double* a, size_t n){
    for (int i=0; i<n; i++){
        a[i] = 0.0;
    }
}
double norm22(   const double *vec, size_t n){
    double sum = 0;
    for (unsigned i=0; i<n; ++i ) {
        sum += vec[i]*vec[i];
    }
    return sum;
}
double norm2(   const double *vec, size_t n){
    double sum = 0;
    for (unsigned i=0; i<n; ++i ) {
        sum += vec[i]*vec[i];
    }
    return sqrt(sum);
}
double normL2( const double *vec, const double*dx, size_t n) {
    double sum = 0.0;
    for (int i=0;i<n;i++) {
        sum += ((vec[i]*vec[i]) * dx[i]);
    }
    double out = sqrt(sum);
    return out;
}
double *matvec(  const double *A, const double *x, size_t n, double *out) {
    // Performs a matrix vector multiplication
    for (int i=0;i<n;i++) {
        out[i] = 0;
        for (int j=0;j<n;j++) {
            out[i] += A[IJ(i,j,n)] * x[j];
        }
    }
    return out;
}
double vecmean( const double *vec, size_t n) {
    double sum = 0.0;
    for (int i=0;i<n;i++) {
        sum += vec[i];
    }
    return (sum / n);
}