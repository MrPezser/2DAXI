//
// Created by tskoepli on 1/26/2024.
//

#include <cmath>
#include <cstdio>
#include "EulerFlux.h"
#include "StateVariables.h"

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

void LeerFlux(const double gam, double normx, double normy, double* uLeft, State varL, double* uRight, State varR,
              double* fout) {
    double fPlus[4],fMnus[4],MnL,MnR;

    //Calculate normal mach number
    MnL = (uLeft[1]*normx + uLeft[2]*normy)/varL.a;
    MnR = (uRight[1]*normx + uLeft[2]*normy)/varR.a;

    //Calculate positive and negative fluxes
    LeerFluxPart(gam, uLeft[1], uLeft[2], normx, normy, MnL, uLeft[0], varL.a, &(fPlus[0]), 1);
    LeerFluxPart(gam, uRight[1],uRight[2],normx, normy, MnR, uRight[0],varR.a, &(fMnus[0]), 0);

    //Find the combined face flux
    fout[0] = fPlus[0] + fMnus[0];
    fout[1] = fPlus[1] + fMnus[1];
    fout[2] = fPlus[2] + fMnus[2];
    fout[3] = fPlus[3] + fMnus[3];

    ASSERT(!_isnan(fout[0]), "nan flux[0]")
    if (_isnan(fout[0]) or _isnan(fout[1]) or _isnan(fout[2]) or _isnan(fout[3])) {
        printf("leerflux isnan\n");
        exit(0);
    }
}


void LDFSS(double normx, double normy, double len, double yface, double* uLeft, State varL, double* uRight, State varR,
           double* flux, double* parr) {

//--------------------------------------------------------------------
//----- inviscid flux contribution (LDFSS)
//
//    rho - density
//    p - pressure
//    u - velocity
//    ho - stagnation enthalpy
//    ys - mass fractions
//    a - sound speed
//     area - interface area
//     dx   - mesh spacing for cell
//     res  - residual vector
//     ev - interface flux vector
// --------------------------------------------------------------------
    //unk = [rho, u, v, T]

    double ahalf = 0.5 * (varL.a + varR.a);

    // Flux Calculation
    double xml = (uLeft[1] *normx +  uLeft[2]*normy)/ahalf;
    double xmr = (uRight[1]*normx + uRight[2]*normy)/ahalf;

    double all = 0.5*(1.0 + sign(xml));
    double alr = 0.5*(1.0 - sign(xmr));

    double btl = -fmax(0.0,1.0-double(int(fabs(xml))));
    double btr = -fmax(0.0,1.0-double(int(fabs(xmr))));

    double xmml =  0.25*(xml+1.0)*(xml+1.0);
    double xmmr = -0.25*(xmr-1.0)*(xmr-1.0);

    double xmhalf = sqrt(0.5*(xml*xml + xmr*xmr));
    double xmc = 0.25*btl*btr*(xmhalf - 1.0)*(xmhalf - 1.0);

    double delp = varL.p - varR.p;
    double psum = varL.p + varR.p;

    double xmcp = xmc * fmax(0.0,(1.0 - (delp/psum + 2.0*fabs(delp)/varL.p)));
    double xmcm = xmc * fmax(0.0,(1.0 + (delp/psum - 2.0*fabs(delp)/varR.p)));
    double cvlp = all*(1.0+btl)*xml - btl*xmml;
    double cvlm = alr*(1.0+btr)*xmr - btr*xmmr;
    double cep = cvlp - xmcp;
    double cem = cvlm + xmcm;

    double fml = len*uLeft[0]*ahalf*cep;
    double fmr = len*uRight[0]*ahalf*cem;

    double ppl = 0.25*(xml+1.0)*(xml+1.0)*(2.0-xml);
    double ppr = 0.25*(xmr-1.0)*(xmr-1.0)*(2.0+xmr);

    double pnet = (all*(1.0+btl) - btl*ppl)*varL.p
                + (alr*(1.0+btr) - btr*ppr)*varR.p;

    parr[0] = len * pnet * normy;

    if (not IAXI){yface = 1.0;}

    flux[0] = yface*(fml + fmr);                                        //continuity
    flux[1] = yface*(fml*uLeft[1] + fmr*uRight[1] + len*pnet*normx);        //x momentum
    flux[2] = yface*(fml*uLeft[2] + fmr*uRight[2]);                     //y momentum w/o  "+pnet*normy"
    flux[3] = yface*(fml*varL.h0  + fmr*(varR.h0));                     //total energy
}