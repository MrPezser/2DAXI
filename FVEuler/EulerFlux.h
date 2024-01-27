//
// Created by tskoepli on 1/26/2024.
//

#ifndef FVEULER_EULERFLUX_H
#define FVEULER_EULERFLUX_H

void getPrimatives(double gam, double *unkel, double *rho, double *u, double *v, double *p, double *c, double *M);
void LeerFlux(const double gam, const double vx, const double vy, const double nx, const double ny,
              const double Mn, const double rho, const double c, double *fout, int isplus);
#endif //FVEULER_EULERFLUX_H
