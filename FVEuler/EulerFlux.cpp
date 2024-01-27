//
// Created by tskoepli on 1/26/2024.
//

#include <cmath>
#include "EulerFlux.h"

void getPrimatives(const double gam, const double *unkel, double *rho, double *u, double *v, double *p, double *c, double *M) {
    // Gets the primative variables at a specific node
    rho[0]       = unkel[0];
    double rhou  = unkel[1];
    double rhov  = unkel[2];
    double rhoe  = unkel[3];

    //Density Limiter
    rho[0] = fmax(1e-8, rho[0]);

    //break down into primatives
    u[0] = rhou / rho[0];
    v[0] = rhov / rho[0];
    double v2 = (u[0]*u[0] + v[0]*v[0]);

    p[0] = (gam - 1) * (rhoe - (0.5 * rho[0] * v2));

    //Pressure Limit
    p[0] = fmax(1e-8, p[0]);

    c[0] = sqrt(gam * p[0] / rho[0]);
    M[0] = sqrt(v2) / c[0];
}

double F1pm(const double M, const double rho, const double c, const int isplus){
    double fout = NAN;
    if (isplus == 1){
        if (M<=-1.0) {
            //F1+
            fout = 0.0;
            return fout;
        }
        if (M>=1.0) {
            fout = rho * M * c;
            return fout;
        }
        fout =  0.25*rho*c*(M+1)*(M+1);
        return fout;
    }
    if (isplus == 0) {
        if (M <= -1.0) {
            //F1-
            fout = -rho * M * c;
            return fout;
        }
        if (M >= 1.0) {
            fout = 0.0;
            return fout;
        }
        fout = -0.25*rho*c*(M-1)*(M-1);
        return fout;
    }
    return fout;
}
void LeerFluxPart(const double gam, const double vx, const double vy, const double nx, const double ny,
              const double Mn, const double rho, const double c, double *fout, int isplus){
    double F1;
    F1 = F1pm(Mn, rho, c, isplus);
    //Calculate velocity information
    double vn = Mn * c;
    double vmags = vx*vx + vy*vy;
    if (isplus == 1) {
        if (Mn >= 1){
            fout[0] = rho*vn;
            double p = fmax(1e-8,c*c*rho/(gam));
            fout[1] = rho*vn*vx + p*nx;
            fout[2] = rho*vn*vy + p*ny;
            fout[3] = vn*( (p/(gam-1)) + 0.5*rho*vmags + p);
            return;
        }
        if (Mn <= -1) {
            fout[0] = 0.0;
            fout[1] = 0.0;
            fout[2] = 0.0;
            fout[3] = 0.0;
            return;
        }
        fout[0] = F1;
        fout[1] = F1 * (vx + nx*(-vn + 2*c)/gam);
        fout[2] = F1 * (vy + ny*(-vn + 2*c)/gam);
        double A = ((gam-1) * vn) + 2*c;
        fout[3] = F1 * ( 0.5*(vmags - vn*vn) + A*A*0.5/(gam*gam - 1.0)) ;
        return;
    }
    if (isplus == 0) {
        if (Mn <= -1){
            fout[0] = rho*vn;
            double p = fmax(1e-8,c*c*rho/(gam));
            fout[1] = rho*vn*vx + p*nx;
            fout[2] = rho*vn*vy + p*ny;
            fout[3] = vn*( (p/(gam-1)) + 0.5*rho*vmags + p);
            return;
        }
        if (Mn >= 1) {
            fout[0] = 0.0;
            fout[1] = 0.0;
            fout[2] = 0.0;
            fout[3] = 0.0;
            return;
        }
        fout[0] = F1;
        fout[1] = F1 * (vx + nx*(-vn - 2*c)/gam);
        fout[2] = F1 * (vy + ny*(-vn - 2*c)/gam);
        double A = ((gam-1) * vn) - 2*c;
        fout[3] = F1 * ( 0.5*(vmags - vn*vn) + A*A * 0.5 / (gam*gam - 1.0)) ;
        return;
    }
}

void LeerFlux(const double gam, const double vx, const double vy, const double nx, const double ny,
                  const double Mn, const double rho, const double c, double *fout, int isplus) {


}