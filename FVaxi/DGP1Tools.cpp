//
// Functions to aid in the addition of DG for higher order extension
// Created by Tsail on 5/29/2024.
//

#include "DGP1Tools.h"
#include "EulerFlux.h"

void get_u_val(const double* unk, State var, Thermo air,const double* ux, const double* uy, double xsi, double eta, double* uout){
    if (ACCUR == 1) {
        for (int k = 0; k < NVAR; k++) {
            if (k==0) {

                //find pressure and derivatives
                double p = var.p;
                double dpdxsi, dpdeta, pout;
                dpdxsi = ux[0] * air.Rs[0] * unk[3] +
                         unk[0] * air.Rs[0] * ux[3];  // = (dp/drho)(drho/dxsi) + (dp/dT)(dT/dxsi)
                dpdeta = uy[0] * air.Rs[0] * unk[3] +
                         unk[0] * air.Rs[0] * uy[3];  // = (dp/drho)(drho/deta) + (dp/dT)(dT/deta)

                //extrapolate with pressure instead of density
                pout = p + (xsi*dpdxsi) + (eta*dpdeta);

                uout[k] = pout / (air.Rs[0]*(unk[3] + (xsi * ux[3]) + (eta * uy[3])));

            } else{
                uout[k] = unk[k] + (xsi * ux[k]) + (eta * uy[k]);
            }
        }
    } else {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k];
        }
    }
}

void get_u_val_standardrecon(const double* unk,const double* ux, const double* uy, double xsi, double eta, double* uout){
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
            db1dx = nxsi[0] / nxsi[2];
            db1dy = nxsi[1] / nxsi[2];

            //d(eta)/dx
            //d(eta)/dy
                //eta plane normal
            double neta[3];
            neta[0] = -0.5*(yL-yR);
            neta[1] =  0.5*(xL-xR);
            neta[2] =  -0.25 * ((xL-xR) * (yD-yU) - (yL-yR) * (xD-xU));
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


void DGP1_xsi_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                        const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                        double* rhsel, double* rhselx, double* rhsely){
    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uLFace[NVAR], uRFace[NVAR], parr;
    State varL = State();
    State varR = State();
    varL.Initialize(uLFace);
    varR.Initialize(uRFace);

    //Cell centers used for adding in pressure term
    double ycL = yCenter[0];
    double ycR = yCenter[1];

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsiL, etaL, xsiR, etaR;

    //1st point
    xsiL = 1.0;
    etaL = 1.0/sqrt(3.0);
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0;
    etaR = 1.0/sqrt(3.0);
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    xsiL = 1.0;
    etaL = -1.0/sqrt(3.0);
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0;
    etaR = -1.0/sqrt(3.0);
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

}

void DGP1_eta_face_integral(int ieL, int ieR, int iuL, int iuR,double* unk, State* ElemVar, double* ux, double* uy,
                            const double* yCenter, Thermo air, double rFace, double* fNormal, double len,
                            double* rhsel, double* rhselx, double* rhsely){
    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uLFace[NVAR], uRFace[NVAR], parr;
    State varL = State();
    State varR = State();
    varL.Initialize(uLFace);
    varR.Initialize(uRFace);

    //Cell centers used for adding in pressure term
    double ycL = yCenter[0];
    double ycR = yCenter[1];

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsiL, etaL, xsiR, etaR;

    //1st point
    xsiL = 1.0/sqrt(3.0);
    etaL = -1.0;
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = 1.0/sqrt(3.0);
    etaR = 1.0;
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    xsiL = -1.0/sqrt(3.0);
    etaL = -1.0;
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0/sqrt(3.0);
    etaR = 1.0;
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

}

void DGP1_boundary_face_integral(int ieIn, int ieEx, int iuIn, int iuEx,double* unk, State* ElemVar, double* ux, double* uy,
                            int iFaceType, double* unkExt, State* EVExt,
                            double yCenter, Thermo air, double rFace, double* fNormal, double len,
                            double* rhsel, double* rhselx, double* rhsely){
    //Input left and right variable/state information
    //Output addition of flux contribution to respective elemets
    double fflux[NVAR], uInFace[NVAR], uExFace[NVAR], parr;
    State varIn = State();
    State varEx = EVExt[ieEx];
    varIn.Initialize(uInFace);
    varEx.Initialize(uExFace);

    //Cell centers used for adding in pressure term
    double ycIn = yCenter;

    //DG extension (xsi flux on vertical face)
    double weight = 0.5;
    double xsi, eta;

    //1st point
    switch (iFaceType) {
        case 1 : //horizontal face - bottom boundary
            xsi = 1.0 / sqrt(3.0);
            eta = -1.0;
        case 2 : //vertical face - right boundary
            xsi = 1.0;
            eta = 1.0 / sqrt(3.0);
        case 3 : //horizontal face - top boundary
            xsi = 1.0 / sqrt(3.0);
            eta = 1.0;
        case 4 : //vertical face - left boundary
            xsi = -1.0;
            eta = 1.0 / sqrt(3.0);
        default:
            printf("No Boundary Face Type Specified\n");
    }


    get_u_val(&(unk[iuIn]), ElemVar[ieIn], air, &(ux[iuIn]), &(uy[iuIn]), xsi, eta, uInFace);
    varIn.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uInFace, varIn,
          uExFace, varEx, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

    //2nd point
    xsiL = -1.0/sqrt(3.0);
    etaL = -1.0;
    get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLFace);
    varL.UpdateState(air);
    xsiR = -1.0/sqrt(3.0);
    etaR = 1.0;
    get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRFace);
    varR.UpdateState(air);

    //Find interface flux
    LDFSS(fNormal[0], fNormal[1], len, rFace, uLFace, varL, uRFace, varR, fflux, &parr);

    //Add flux contribution to elements
    for (int kvar=0; kvar<NVAR; kvar++){
        if (kvar == NSP+1){ //Axisymmetric pressure correction
            rhsel[iuL + kvar]  -= weight * (fflux[kvar] + (ycL*parr));
            rhsel[iuR + kvar]  += weight * (fflux[kvar] + (ycR*parr));
            rhselx[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar] + (ycL*parr)) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar] + (ycR*parr)) * etaR;
        } else {
            rhsel[iuL + kvar]  -= weight * fflux[kvar];
            rhsel[iuR + kvar]  += weight * fflux[kvar];
            rhselx[iuL + kvar] -= weight * (fflux[kvar]) * xsiL;
            rhselx[iuR + kvar] += weight * (fflux[kvar]) * xsiR;
            rhsely[iuL + kvar] -= weight * (fflux[kvar]) * etaL;
            rhsely[iuR + kvar] += weight * (fflux[kvar]) * etaR;
        }
    }
    ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")

}