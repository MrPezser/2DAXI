//
// Functions to aid in the addition of DG for higher order extension
// Created by Tsail on 5/29/2024.
//

#include "DGP1Tools.h"

void get_u_val(const double* unk, const double* ux, const double* uy, double xsi, double eta, double* uout){
    if (ACCUR == 1) {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k] + (xsi * ux[k]) + (eta * uy[k]);
        }
    } else {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k];
        }
    }
}


void DGP1_volume_integral(int nx, int ny, double vol, double* xfa, double* yfa, double* geoel, double* unk, State* ElemVar,
                          double* duxdt, double* duydt){
    // NOTE: Includes axisymetric flux modification
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++) {
            int iu = IJK(i,j,0,nx-1,NVAR);
            int iel = IJ(i,j,nx-1);

            //Calculate cell centered flux
            double Fx[NVAR], Fy[NVAR];
            double* unkij = &(unk[iu]);
            State var = ElemVar[iel];

            //Find Flux Vector
            double ycc = geoel[IJK(i,j,2,nx-1, 3)];
            double rho = unkij[0];
            Fx[0] = ycc *  rho * unkij[1];
            Fx[1] = ycc * (rho * unkij[1] * unkij[1] + var.p);
            Fx[2] = ycc *  rho * unkij[1] * unkij[2];
            Fx[3] = ycc *  rho * unkij[1] * var.h0;//unkij[1] * (rho*var.e + var.p);

            Fy[0] = ycc *  rho * unkij[2];
            Fy[1] = ycc *  rho * unkij[1] * unkij[2];
            Fy[2] = ycc * (rho * unkij[2] * unkij[2] + var.p);
            Fy[3] = ycc *  rho * unkij[2] * var.h0;//unkij[2] * (rho*var.e + var.p);

            // Derivatives of Basis Functions
            double xL,xR,xD,xU,yL,yR,yD,yU;
            double db1dx, db1dy, db2dx,db2dy;
            xL = xfa[IJK(i,j,1,nx,2)];
            xR = xfa[IJK(i+1,j,1,nx,2)];
            xD = xfa[IJK(i,j,0,nx,2)];
            xU = xfa[IJK(i,j+1,0,nx,2)];

            yL = yfa[IJK(i,j,1,nx,2)];
            yR = yfa[IJK(i+1,j,1,nx,2)];
            yD = yfa[IJK(i,j,0,nx,2)];
            yU = yfa[IJK(i,j+1,0,nx,2)];

            //d(xsi)/dx
            //d(xsi)/dy
                //xsi plane normal
            double nxsi[3];
            nxsi[0] = -0.5*(yU-yD);
            nxsi[1] =  0.5*(xU-xD);
            nxsi[2] =  -0.25 * ((xU-xD) * (yL-yR) - (yU-yD) * (xL-xR));
            ASSERT(!_isnan(nxsi[2]), "Flat basis function dummy.")
            db1dx = nxsi[0] / nxsi[2];
            db1dy = nxsi[1] / nxsi[2];

            //d(eta)/dx
            //d(eta)/dy
                //eta plane normal
            double neta[3];
            neta[0] = -0.5*(yL-yR);
            neta[1] =  0.5*(xL-xR);
            neta[2] =  -0.25 * ((xL-xR) * (yD-yU) - (yL-yR) * (xD-xU));
            ASSERT(!_isnan(neta[2]), "Flat basis function dummy.")
            db2dx = neta[0] / neta[2];
            db2dy = neta[1] / neta[2];

            //Assemble volume contribution
            for (int kvar=0; kvar<NVAR; kvar++){
                duxdt[iu+kvar] += 3.0 * (Fx[kvar]*db1dx + Fy[kvar]*db1dy);
                duydt[iu+kvar] += 3.0 * (Fx[kvar]*db2dx + Fy[kvar]*db2dy);
            }

        }
    }
}