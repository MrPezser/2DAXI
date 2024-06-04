//
// Created by Tsail on 1/28/2024.
//

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "SpatialDiscretization.h"
#include "Indexing.h"
#include "BoundaryConditions.h"
#include "EulerFlux.h"
#include "StateVariables.h"
#include "DGP1Tools.h"

void generate_ghost_cells(int nx, int ny, double* unk, double* ux, double* uy, State* ElemVar, Thermo air, int* ibound,
                            double* geofa, double* uFS, double* uGBot, double* uGTop, double* uGLeft, double* uGRight,
                            State* BotVar, State* TopVar, State* LeftVar, State* RightVar){
    double unkelij[NVAR];
    State varij = State();

    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[i];
        int iint = IJK(i,0,0, nx-1, 4);
        int iel = IJ(i, 0, nx-1);

        //DG extension
        double xsi, eta;
        xsi = 0.0;
        eta = -1.0;
        get_u_val(&(unk[iint]), ElemVar[iel], air, &(ux[iint]), &(uy[iint]), xsi, eta, unkelij);
        varij.Initialize(unkelij);
        varij.UpdateState(air);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];
        //==========Ghost State
        BotVar[i].Initialize(&(uGBot[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,unkelij, varij,
                       &(uGBot[IJ(0,i,NVAR)]));
        BotVar[i].UpdateState(air);
    }
    //right side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[j+ nx-1];
        int iint = IJK(nx-2,j,0, nx-1, NVAR);
        int iel = IJ(nx-2, j, nx-1);

        //DG extension
        double xsi, eta;
        xsi = 1.0;
        eta = 0.0;
        get_u_val(&(unk[iint]), ElemVar[iel], air, &(ux[iint]), &(uy[iint]), xsi, eta, unkelij);
        varij.Initialize(unkelij);
        varij.UpdateState(air);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(0, j, 4,nx,6)];
        normy = geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        RightVar[j].Initialize(&(uGRight[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,unkelij, varij,
                       &(uGRight[IJ(0,j,NVAR)]));
        RightVar[j].UpdateState(air);
    }

    //top side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        int ib = nx-2-i;
        btype = ibound[ib+nx+ny-2];
        int iint = IJK(i,ny-2,0, nx-1, NVAR);
        int iel = IJ(i, ny-2, nx-1);

        //DG extension
        double xsi, eta;
        xsi = 0.0;
        eta = 1.0;
        get_u_val(&(unk[iint]), ElemVar[iel], air, &(ux[iint]), &(uy[iint]), xsi, eta, unkelij);
        varij.Initialize(unkelij);
        varij.UpdateState(air);

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];
        //==========Ghost State
        TopVar[i].Initialize(&(uGTop[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, unkelij, varij,
                       &(uGTop[IJ(0,i,NVAR)]));
        TopVar[i].UpdateState(air);
    }

    //left side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        int jb = (ny-2)-j;
        btype = ibound[jb+(2*nx)+ny-3];
        int iint = IJK(0,j,0, nx-1, 4);
        int iel = IJ(0, j, nx-1);

        //DG extension
        double xsi, eta;
        xsi = -1.0;
        eta =  0.0;
        get_u_val(&(unk[iint]), ElemVar[iel], air, &(ux[iint]), &(uy[iint]), xsi, eta, unkelij);
        varij.Initialize(unkelij);
        varij.UpdateState(air);

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        LeftVar[j].Initialize(&(uGLeft[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, unkelij, varij,
                       &(uGLeft[IJ(0,j,NVAR)]));
        LeftVar[j].UpdateState(air);
    }
}

void viscous(int nx, double normy, double normx, double* uLeft, State& varL, double* uRight, State varR, double* dc, double* visc_contrib){

    if (IVISC) {
        // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
        ///Need to make axisymmatric modification and add extra termsnn/JE
        printf("Need to update viscous fluxes for axisymmetric and DP higher order\n");
        exit(0);
        double tau[4], Sij[4], st[4], dc2i, trace, vflux[2];
        double uL, vL, uR, vR;

        double mu = 0.5 * (varL.mu + varR.mu);

        dc2i = 1.0 / (dc[0] * dc[0] + dc[1] * dc[1]);

        uL = uLeft[1];
        vL = uLeft[2];
        uR = uRight[1];
        vR = uRight[2];
        //stress tensor
        st[0] = (uR - uL) * dc[0] * dc2i;//  du/dx
        st[1] = (uR - uL) * dc[1] * dc2i;//  du/dy
        st[2] = (vR - vL) * dc[0] * dc2i;//  dv/dx
        st[3] = (vR - vL) * dc[1] * dc2i;//  dv/dy
        trace = (1.0 / 3.0) * (st[0] + st[3]);

        // Viscous strain rate
        Sij[0] = st[0] - trace;                //S_1,1
        Sij[1] = 0.5 * (st[1] + st[2]);         //S_1,2
        Sij[2] = Sij[1];                        //S_2,1
        Sij[3] = st[1] - trace;                //S_2,2

        //Viscous Stress  (Stoke's Law for Monoatoms)
        tau[0] = 2 * mu * Sij[0];
        tau[1] = 2 * mu * Sij[1];
        tau[2] = 2 * mu * Sij[2];
        tau[3] = 2 * mu * Sij[3];

        // Viscous Contributions to Flux
        vflux[0] = tau[0] * normx + tau[2] * normy;
        vflux[1] = tau[1] * normx + tau[3] * normy;

        visc_contrib[0] = vflux[0];
        visc_contrib[1] = vflux[1];
        visc_contrib[2] = ((uL * tau[0] + vL * tau[2]) * normx + (uL * tau[1] + vL * tau[3]) * normy);

        visc_contrib[3] = vflux[0];
        visc_contrib[4] = vflux[1];
        visc_contrib[5] = ((uR * tau[0] + vR * tau[2]) * normx + (uR * tau[1] + vR * tau[3]) * normy);

        for (int iv = 0; iv < 6; iv++) {
            if (_isnan(visc_contrib[iv])) {
                printf("asdfaed\n");
            }
            ASSERT(!_isnan(visc_contrib[iv]), "Visc Flux NaN");
        }
    } else {

        visc_contrib[0] = 0.0;
        visc_contrib[1] = 0.0;
        visc_contrib[2] = 0.0;
        visc_contrib[3] = 0.0;
        visc_contrib[4] = 0.0;
        visc_contrib[5] = 0.0;
    }
}

void calc_dudt(int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* yfa, double* xfa, double* unk, double* ux, double* uy, double* dudt, double* duxdt, double* duydt) {
    int nelem = (nx-1)*(ny-1);
    double *rhsel, *rhselx, *rhsely, parr;
    rhsel  = (double*)malloc(NVAR*nelem*sizeof(double));
    rhselx = (double*)malloc(NVAR*nelem*sizeof(double));
    rhsely = (double*)malloc(NVAR*nelem*sizeof(double));
    for(int i=0; i<NVAR*nelem; i++) {
        rhsel[i]  = 0.0;
        rhselx[i] = 0.0;
        rhsely[i] = 0.0;
    }


    //Calculate boundary cell state (ghost state)
    double *uGBot, *uGRight, *uGTop, *uGLeft;
    State *BotVar, *TopVar, *RightVar, *LeftVar;
    if (ACCUR == 0) {
        uGBot   = (double*)malloc((NVAR * (nx - 1))*sizeof(double));
        uGRight = (double*)malloc((NVAR * (ny - 1))*sizeof(double));
        uGTop   = (double*)malloc((NVAR * (nx - 1))*sizeof(double));
        uGLeft  = (double*)malloc((NVAR * (ny - 1))*sizeof(double));
        BotVar   = (State*)malloc((NVAR * (nx - 1))*sizeof(State));
        RightVar = (State*)malloc((NVAR * (ny - 1))*sizeof(State));
        TopVar   = (State*)malloc((NVAR * (nx - 1))*sizeof(State));
        LeftVar  = (State*)malloc((NVAR * (ny - 1))*sizeof(State));
        generate_ghost_cells(nx, ny, unk, ux, uy, ElemVar, air, ibound, geofa, uFS, uGBot, uGTop, uGLeft, uGRight,
                             BotVar, TopVar, LeftVar, RightVar);
    } else if (ACCUR==1){
        uGBot   = (double*)malloc(2*(NVAR * (nx - 1))*sizeof(double));
        uGRight = (double*)malloc(2*(NVAR * (ny - 1))*sizeof(double));
        uGTop   = (double*)malloc(2*(NVAR * (nx - 1))*sizeof(double));
        uGLeft  = (double*)malloc(2*(NVAR * (ny - 1))*sizeof(double));
        BotVar   = (State*)malloc(2*((nx - 1))*sizeof(State));
        RightVar = (State*)malloc(2*((ny - 1))*sizeof(State));
        TopVar   = (State*)malloc(2*((nx - 1))*sizeof(State));
        LeftVar  = (State*)malloc(2*((ny - 1))*sizeof(State));
        DGP1_ghost_cell_generator(nx, ny, unk, ux, uy, ElemVar, air, ibound, geofa, uFS, uGBot, uGTop, uGLeft, uGRight,
                BotVar, TopVar, LeftVar, RightVar);
    } else {
        printf("ACCUR must be 1 or 0.");
        exit(0);
    }

    //====================Evaluate Flux Contributions====================
    //dudt = sum(flux_in * face_length) / volume
    //==========Fully Interior Faces


    // I|xi Fluxes
    for (int i=1; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            // ~~~~~~~~~~ Inviscid fluxes ~~~~~~~~~~
            //face above point i,j
            double len, fNormal[2], fflux[NVAR], rFace, uLeft[NVAR], uRight[NVAR], yCenter[2];
            State varL = State();
            State varR = State();
            len = geofa[IJK(i,j,3,nx,6)];
            fNormal[0] = geofa[IJK(i,j,4,nx,6)];
            fNormal[1] = geofa[IJK(i,j,5,nx,6)];
            rFace = yfa[IJK(i,j,1,nx,2)];

            int iuL = IJK(i-1,j,0,nx-1,NVAR);
            int iuR = IJK(i  ,j,0,nx-1,NVAR);
            int ieL = IJ(i-1, j, nx-1);
            int ieR = IJ(i,   j, nx-1);

            yCenter[0] = geoel[IJK(i-1, j, 2, nx-1, 3)];
            yCenter[1] = geoel[IJK(i,   j, 2, nx-1, 3)];

            DGP1_xsi_face_integral(ieL, ieR, iuL, iuR, unk, ElemVar, ux, uy, yCenter, air,
                                   rFace, fNormal, len, rhsel, rhselx, rhsely);

            /*
            //DG extension (xsi flux on vertical face)
            double xsiL, etaL, xsiR, etaR;
            xsiL = 1.0;
            etaL = 0.0;
            get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLeft);
            varL.Initialize(uLeft);
            varL.UpdateState(air);
            xsiR = -1.0;
            etaR = 0.0;
            get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRight);
            varR.Initialize(uRight);
            varR.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, len, yface, uLeft, varL, uRight, varR, fflux, &parr);

            //Add pressure term in
            double ycL = geoel[IJK(i-1, j, 2, nx-1, 3)];
            double ycR = geoel[IJK(i,   j, 2, nx-1, 3)];

            //Add flux contribution to elements
            rhsel[iuL  ] -= fflux[0];
            rhsel[iuL+1] -= fflux[1];
            rhsel[iuL+2] -=(fflux[2] + (ycL*parr));
            rhsel[iuL+3] -= fflux[3];

            rhsel[iuR  ] += fflux[0];
            rhsel[iuR+1] += fflux[1];
            rhsel[iuR+2] +=(fflux[2] + (ycR*parr));
            rhsel[iuR+3] += fflux[3];



            if (ACCUR == 1){
                rhselx[iuL  ] -= (fflux[0]) * xsiL;
                rhselx[iuL+1] -= (fflux[1]) * xsiL;
                rhselx[iuL+2] -= (fflux[2] + (ycL*parr)) * xsiL;
                rhselx[iuL+3] -= (fflux[3]) * xsiL;

                rhselx[iuR  ] += (fflux[0]) * xsiR;
                rhselx[iuR+1] += (fflux[1]) * xsiR;
                rhselx[iuR+2] += (fflux[2] + (ycR*parr)) * xsiR;
                rhselx[iuR+3] += (fflux[3]) * xsiR;

                //rhsely[iuL  ] -= (fflux[0]) * etaL;
                //rhsely[iuL+1] -= (fflux[1]) * etaL;
                //rhsely[iuL+2] -= (fflux[2] + (ycL*parr)) * etaL;
                //rhsely[iuL+3] -= (fflux[3]) * etaL;

                //rhsely[iuR  ] += (fflux[0]) * etaR;
                //rhsely[iuR+1] += (fflux[1]) * etaR;
                //rhsely[iuR+2] += (fflux[2] + (ycR*parr)) * etaR;
                //rhsely[iuR+3] += (fflux[3]) * etaR;
            }

            ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")
            */
            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double vflux[6], dc[2];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i,j,1,nx-1,3)] - geoel[IJK(i-1,j,1,nx-1,3)];
            dc[1] = geoel[IJK(i,j,2,nx-1,3)] - geoel[IJK(i-1,j,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL+1] += len * vflux[0];
            rhsel[iuL+2] += len * vflux[1];
            rhsel[iuL+3] += len * vflux[2];

            rhsel[iuR+1] -= len * vflux[3];
            rhsel[iuR+2] -= len * vflux[4];
            rhsel[iuR+3] -= len * vflux[5];
             */
        }
    }

    // J|eta fluxes
    for (int i=0; i<nx-1; i++){
        for (int j=1; j<ny-1; j++){
            //faces to the left/above point i,j
            double len, fNormal[2], fflux[NVAR], rFace, uLeft[NVAR], uRight[NVAR], yCenter[2];
            State varL = State();
            State varR = State();
            len = geofa[IJK(i,j,0,nx,6)];
            fNormal[0] = geofa[IJK(i,j,1,nx,6)];
            fNormal[1] = geofa[IJK(i,j,2,nx,6)];
            rFace = yfa[IJK(i,j,0,nx,2)];

            int iuL = IJK(i,j  ,0,nx-1,NVAR);
            int iuR = IJK(i,j-1,0,nx-1,NVAR);
            int ieL = IJ(i, j  , nx-1);
            int ieR = IJ(i, j-1, nx-1);

            yCenter[0] = geoel[IJK(i, j,   2, nx-1, 3)];
            yCenter[1] = geoel[IJK(i, j-1, 2, nx-1, 3)];

            DGP1_eta_face_integral(ieL, ieR, iuL, iuR, unk, ElemVar, ux, uy, yCenter, air,
                                   rFace, fNormal, len, rhsel, rhselx, rhsely);

            /*
            //DG extension (eta flux on 'horizontal' face)
            double xsiL, etaL, xsiR, etaR;
            xsiL =  0.0;
            etaL = -1.0;
            get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLeft);
            varL.Initialize(uLeft);
            varL.UpdateState(air);
            xsiR = 0.0;
            etaR = 1.0;
            get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRight);
            varR.Initialize(uRight);
            varR.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            double normx = fNormal[0];
            double normy = fNormal[1];
            double yface = rFace;
            LDFSS(normx, normy, len, yface, uLeft, varL, uRight, varR, fflux, &parr);

            //Add pressure term in
            double ycL = geoel[IJK(i, j,   2, nx-1, 3)];
            double ycR = geoel[IJK(i, j-1, 2, nx-1, 3)];

            //Add flux contribution to elements
            rhsel[iuL  ] -= fflux[0];
            rhsel[iuL+1] -= fflux[1];
            rhsel[iuL+2] -=(fflux[2] + (ycL*parr));
            rhsel[iuL+3] -= fflux[3];

            rhsel[iuR  ] += fflux[0];
            rhsel[iuR+1] += fflux[1];
            rhsel[iuR+2] +=(fflux[2] + (ycR*parr));
            rhsel[iuR+3] += fflux[3];

            if (ACCUR == 1){
                //rhselx[iuL  ] -= (fflux[0]) * xsiL;
                //rhselx[iuL+1] -= (fflux[1]) * xsiL;
                //rhselx[iuL+2] -= (fflux[2] + (ycL*parr)) * xsiL;
                //rhselx[iuL+3] -= (fflux[3]) * xsiL;

                //rhselx[iuR  ] += (fflux[0]) * xsiR;
                //rhselx[iuR+1] += (fflux[1]) * xsiR;
                //rhselx[iuR+2] += (fflux[2] + (ycR*parr)) * xsiR;
                //rhselx[iuR+3] += (fflux[3]) * xsiR;

                rhsely[iuL  ] -= (fflux[0]) * etaL;
                rhsely[iuL+1] -= (fflux[1]) * etaL;
                rhsely[iuL+2] -= (fflux[2] + (ycL*parr)) * etaL;
                rhsely[iuL+3] -= (fflux[3]) * etaL;

                rhsely[iuR  ] += (fflux[0]) * etaR;
                rhsely[iuR+1] += (fflux[1]) * etaR;
                rhsely[iuR+2] += (fflux[2] + (ycR*parr)) * etaR;
                rhsely[iuR+3] += (fflux[3]) * etaR;
            }

            ASSERT(!_isnan(fflux[0]+fflux[1]+fflux[2]+fflux[3]),"NAN in flux splitting")
            */

            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double vflux[6], dc[2];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i,j-1,1,nx-1,3)] - geoel[IJK(i,j,1,nx-1,3)];
            dc[1] = geoel[IJK(i,j-1,2,nx-1,3)] - geoel[IJK(i,j,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL+1] += len * vflux[0];
            rhsel[iuL+2] += len * vflux[1];
            rhsel[iuL+3] += len * vflux[2];

            rhsel[iuR+1] -= len * vflux[3];
            rhsel[iuR+2] -= len * vflux[4];
            rhsel[iuR+3] -= len * vflux[5];
             */
        }
    }



    for (int iu=0; iu<nelem*NVAR; iu++){
        if (_isnan(rhsel[iu])){
            printf("NAN RHSEL\n");
        }
    }

    //==========Boundary faces
    double len, normx, normy, fflux[NVAR], yface, uRight[NVAR], uLeft[NVAR];
    State varL = State();
    State varR = State();

    // I|xi fluxes
    for (int j=0; j<ny-1; j++){
        //--------------------  LEFT BOUNDARY
        len = geofa[IJK(0, j, 3, nx, 6)];
        normx = geofa[IJK(0, j, 4, nx, 6)];
        normy = geofa[IJK(0, j, 5, nx, 6)];
        yface = yfa[IJK(0, j, 1, nx, 2)];
        int iuL = IJ(0, j, NVAR);
        int iuR = IJK(0, j, 0, nx - 1, NVAR);
        int ieL = j;
        int ieR = IJ(0, j, nx - 1);
        //double ycL = 2.0 * geoel[IJK(0, j, 2, nx-1, 3)] - geoel[IJK(0, j+1, 2, nx-1, 3)]; //extrapolate to get gost y coord
        double ycR = geoel[IJK(0, j, 2, nx - 1, 3)];

        if (ACCUR==1){
            int iFaceType = 4; //number CCW starting from bot
            int ieEx = IJ(0,j,2);
            int iuEx = IJK(0,j,0,2,NVAR);
            double yCenter = ycR;
            double fNormal[2] = {normx, normy};
            DGP1_boundary_face_integral(ieR, ieEx, iuR, iuEx, unk, ElemVar, ux, uy, iFaceType, uGLeft, LeftVar,
                    yCenter, air, yface, fNormal, len, rhsel, rhselx, rhsely);
        } else {

            //DG extension (xsi flux on vertical face)
            double xsiR, etaR;
            xsiR = -1.0;
            etaR = 0.0;
            get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRight);
            varR.Initialize(uRight);
            varR.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, len, yface, &(uGLeft[iuL]), LeftVar[j], uRight, varR, fflux, &parr);

            //Add flux contribution to elements
            rhsel[iuR] += fflux[0];
            rhsel[iuR + 1] += fflux[1];
            rhsel[iuR + 2] += (fflux[2] + (ycR * parr));
            rhsel[iuR + 3] += fflux[3];
        }

            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6]{};
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(1,j,1,nx-1,3)] - geoel[IJK(0,j,1,nx-1,3)];
            dc[1] = geoel[IJK(1,j,2,nx-1,3)] - geoel[IJK(0,j,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(uGLeft[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR+1] -= len * vflux[3];
            rhsel[iuR+2] -= len * vflux[4];
            rhsel[iuR+3] -= len * vflux[5];
             */



            //--------------------  RIGHT BOUNDARY
            len = geofa[IJK(nx - 1, j, 3, nx, 6)];
            normx = geofa[IJK(nx - 1, j, 4, nx, 6)];
            normy = geofa[IJK(nx - 1, j, 5, nx, 6)];
            yface = yfa[IJK(nx - 1, j, 1, nx, 2)];
            iuL = IJK(nx - 2, j, 0, nx - 1, NVAR);
            iuR = IJ(0, j, NVAR);
            ieL = IJ(nx - 2, j, nx - 1);
            ieR = j;
            double ycL = geoel[IJK(nx - 2, j, 2, nx - 1, 3)];

        if (ACCUR==1){
            int iFaceType = 2; //number CCW starting from bot = 1
            int ieEx = IJ(0,j,2);
            int iuEx = IJK(0,j,0,2,NVAR);
            double yCenter = ycL;
            double fNormal[2] = {normx, normy};
            DGP1_boundary_face_integral(ieL, ieEx, iuL, iuEx, unk, ElemVar, ux, uy, iFaceType, uGRight, RightVar,
                                        yCenter, air, yface, fNormal, len, rhsel, rhselx, rhsely);
        } else {

            //DG extension (xsi flux on vertical face)
            double xsiL, etaL;
            xsiL = 1.0;
            etaL = 0.0;
            get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLeft);
            varL.Initialize(uLeft);
            varL.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, len, yface, uLeft, varL, &(uGRight[iuR]), RightVar[j], fflux, &parr);



            //Add flux contribution to elements
            rhsel[iuL] -= fflux[0];
            rhsel[iuL + 1] -= fflux[1];
            rhsel[iuL + 2] -= (fflux[2] + (ycL * parr));
            rhsel[iuL + 3] -= fflux[3];
        }
            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(nx-2,j,1,nx-1,3)] - geoel[IJK(nx-3,j,1,nx-1,3)];
            dc[1] = geoel[IJK(nx-2,j,2,nx-1,3)] - geoel[IJK(nx-3,j,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(uGRight[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL+1] += len * vflux[0];
            rhsel[iuL+2] += len * vflux[1];
            rhsel[iuL+3] += len * vflux[2];
             */

    }


    // J|eta fluxes
    for (int i=0; i<nx-1; i++) {
        //BOTTOM BOUNDARY
        len   = geofa[IJK(i,0,0,nx,6)];
        normx = geofa[IJK(i,0,1,nx,6)];
        normy = geofa[IJK(i,0,2,nx,6)];
        yface = yfa[IJK(i,0,0,nx,2)];

        int iuL = IJK(i, 0, 0, nx - 1, NVAR);
        int iuR = IJ(0,i,NVAR);
        int ieL = IJ(i,0,nx-1);
        int ieR = i;
        double ycL = geoel[IJK(i, 0, 2, nx - 1, 3)];

        if (ACCUR==1){
            int iFaceType = 1; //number CCW starting from bot = 1
            int ieEx = IJ(0,i,2);
            int iuEx = IJK(0,i,0,2,NVAR);
            double yCenter = ycL;
            double fNormal[2] = {normx, normy};
            DGP1_boundary_face_integral(ieL, ieEx, iuL, iuEx, unk, ElemVar, ux, uy, iFaceType, uGBot, BotVar,
                                        yCenter, air, yface, fNormal, len, rhsel, rhselx, rhsely);
        } else {

            //DG extension (eta flux on 'horizontal' face)
            double xsiL, etaL;
            xsiL = 0.0;
            etaL = -1.0;
            get_u_val(&(unk[iuL]), ElemVar[ieL], air, &(ux[iuL]), &(uy[iuL]), xsiL, etaL, uLeft);
            varL.Initialize(uLeft);
            varL.UpdateState(air);
            //xsi = 0.0;
            //eta = 1.0;
            //get_u_val(&(unk[iuR]), &(ux[iuR]), &(uy[iuR]), xsi, eta, uRight);
            //varR.Initialize(uRight);
            //varR.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, len, yface, uLeft, varL, &(uGBot[iuR]), BotVar[i], fflux, &parr);

            //Add flux contribution to elements
            rhsel[iuL] -= fflux[0];
            rhsel[iuL + 1] -= fflux[1];
            rhsel[iuL + 2] -= (fflux[2] + ycL * parr);
            rhsel[iuL + 3] -= fflux[3];
        }
            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i,0,1,nx-1,3)] - geoel[IJK(i,1,1,nx-1,3)];
            dc[1] = geoel[IJK(i,0,2,nx-1,3)] - geoel[IJK(i,1,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(uGBot[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL+1] += len * vflux[0];
            rhsel[iuL+2] += len * vflux[1];
            rhsel[iuL+3] += len * vflux[2];
             */

        //TOP BOUNDARY
        len = geofa[IJK(i, ny - 1, 0, nx, 6)];
        normx = geofa[IJK(i, ny - 1, 1, nx, 6)];
        normy = geofa[IJK(i, ny - 1, 2, nx, 6)];
        yface = yfa[IJK(i, ny - 1, 0, nx, 2)];

        iuL = IJ(0, i, NVAR);
        iuR = IJK(i, ny - 2, 0, nx - 1, NVAR);
        ieL = i;
        ieR = IJ(i, ny - 2, nx - 1);

        double ycR = geoel[IJK(i, ny - 2, 2, nx - 1, 3)];

        if (ACCUR==1){
            int iFaceType = 3; //number CCW starting from bot = 1
            int ieEx = IJ(0,i,2);
            int iuEx = IJK(0,i,0,2,NVAR);
            double yCenter = ycR;
            double fNormal[2] = {normx, normy};
            DGP1_boundary_face_integral(ieR, ieEx, iuR, iuEx, unk, ElemVar, ux, uy, iFaceType, uGTop, TopVar,
                                        yCenter, air, yface, fNormal, len, rhsel, rhselx, rhsely);
        } else {

            //DG extension (eta flux on 'horizontal' face)
            double xsiR, etaR;
            xsiR = 0.0;
            etaR = 1.0;
            get_u_val(&(unk[iuR]), ElemVar[ieR], air, &(ux[iuR]), &(uy[iuR]), xsiR, etaR, uRight);
            varR.Initialize(uRight);
            varR.UpdateState(air);

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, len, yface, &(uGTop[iuL]), TopVar[i], uRight, varR, fflux, &parr);


            //Add flux contribution to elements
            rhsel[iuR] += fflux[0];
            rhsel[iuR + 1] += fflux[1];
            rhsel[iuR + 2] += (fflux[2] + (ycR * parr));
            rhsel[iuR + 3] += fflux[3];

        }
            /*
            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i,ny-3,1,nx-1,3)] - geoel[IJK(i,ny-2,1,nx-1,3)];
            dc[1] = geoel[IJK(i,ny-3,2,nx-1,3)] - geoel[IJK(i,ny-2,2,nx-1,3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(uGTop[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR+1] -= len * vflux[3];
            rhsel[iuR+2] -= len * vflux[4];
            rhsel[iuR+3] -= len * vflux[5];
             */

    }

    //====================Combine to Find du/dt====================
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){

            int iu = IJK(i,j,0,nx-1,NVAR);
            double vol = geoel[IJK(i,j,0,nx-1,3)];
            dudt[iu  ] = rhsel[iu  ] / vol;
            dudt[iu+1] = rhsel[iu+1] / vol;
            dudt[iu+2] = rhsel[iu+2] / vol;
            dudt[iu+3] = rhsel[iu+3] / vol;

            if (ACCUR==1) {
                duxdt[iu]     = 3.0 * rhselx[iu]     / vol;
                duxdt[iu + 1] = 3.0 * rhselx[iu + 1] / vol;
                duxdt[iu + 2] = 3.0 * rhselx[iu + 2] / vol;
                duxdt[iu + 3] = 3.0 * rhselx[iu + 3] / vol;

                duydt[iu]     = 3.0 * rhsely[iu]     / vol;
                duydt[iu + 1] = 3.0 * rhsely[iu + 1] / vol;
                duydt[iu + 2] = 3.0 * rhsely[iu + 2] / vol;
                duydt[iu + 3] = 3.0 * rhsely[iu + 3] / vol;
            }

            if (_isnan(dudt[iu]) or _isnan(dudt[iu+1]) or _isnan(dudt[iu+2]) or _isnan(dudt[iu+3])) {
                printf("dug\n");
            }

            ASSERT(!_isnan(dudt[iu]), "NaN Encounterd in Flux formulation")

        }
    }
    /// SINCE M IS DIAGONAL MATRIX, THE BELOW ALREADY INCLUDES THE MULTIPLE 3/VOL FROM ITS INVERSION
    DGP1_volume_integral(nx, ny, 1.0, xfa, yfa, geoel, unk, ElemVar, duxdt, duydt);

    free(rhsel);
    free(rhselx);
    free(rhsely);

    free(BotVar);
    free(TopVar);
    free(RightVar);
    free(LeftVar);

    free(uGBot);
    free(uGRight);
    free(uGTop);
    free(uGLeft);
}
