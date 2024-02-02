//
// Created by tskoepli on 1/26/2024.
//

#ifndef FVEULER_EULERFLUX_H
#define FVEULER_EULERFLUX_H

void getPrimatives(double gam, const double *unkel, double *rho, double *u, double *v, double *p, double *c, double *M);
void LeerFlux(double gam, double normx, double normy, double* uLeft, double* uRight, double* fout);




#endif //FVEULER_EULERFLUX_H
