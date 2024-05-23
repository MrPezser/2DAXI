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


void viscous(int nx, double normy, double normx, double* uLeft, State& varL, double* uRight, State varR, double* dc, double* visc_contrib){

    // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
    double tau[4], Sij[4], st[4], dc2i, trace, vflux[2];
    double uL, vL, uR, vR;

    double mu = 0.5*(varL.mu + varR.mu);

    dc2i = 1.0/(dc[0]*dc[0] + dc[1]*dc[1]);

    uL = uLeft[1];
    vL = uLeft[2];
    uR = uRight[1];
    vR = uRight[2];
    //stress tensor
    st[0] = (uR - uL) * dc[0] * dc2i ;//  du/dx
    st[1] = (uR - uL) * dc[1] * dc2i ;//  du/dy
    st[2] = (vR - vL) * dc[0] * dc2i ;//  dv/dx
    st[3] = (vR - vL) * dc[1] * dc2i ;//  dv/dy
    trace = (1.0/3.0)*(st[0] + st[3]);

    // Viscous strain rate
    Sij[0] = st[0] - trace;                //S_1,1
    Sij[1] = 0.5*(st[1] + st[2]);         //S_1,2
    Sij[2] = Sij[1];                        //S_2,1
    Sij[3] = st[1] - trace;                //S_2,2

    //Viscous Stress  (Stoke's Law for Monoatoms)
    tau[0] = 2*mu*Sij[0];
    tau[1] = 2*mu*Sij[1];
    tau[2] = 2*mu*Sij[2];
    tau[3] = 2*mu*Sij[3];

    // Viscous Contributions to Flux
    vflux[0] = tau[0]*normx + tau[2]*normy;
    vflux[1] = tau[1]*normx + tau[3]*normy;

    visc_contrib[0] = vflux[0];
    visc_contrib[1] = vflux[1];
    visc_contrib[2] = ((uL*tau[0] + vL*tau[2])*normx + (uL*tau[1] + vL*tau[3])*normy);

    visc_contrib[3] = vflux[0];
    visc_contrib[4] = vflux[1];
    visc_contrib[5] = ((uR*tau[0] + vR*tau[2])*normx + (uR*tau[1] + vR*tau[3])*normy);

    for (int iv=0; iv<6; iv++){
        if (_isnan(visc_contrib[iv])){
            printf("asdfaed\n");
        }
        ASSERT(!_isnan(visc_contrib[iv]),"Visc Flux NaN");
    }
}

void calc_dudt(int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* unk, double* dudt) {
    int nelem = (nx-1)*(ny-1);
    double* rhsel;
    rhsel = (double*)malloc(NVAR*nelem*sizeof(double));
    for(int i=0; i<NVAR*nelem; i++) {
        rhsel[i] = 0.0;
    }

    //Calculate boundary cell state (ghost state)
    double uGBot[NVAR*(nx-1)], uGRight[NVAR*(ny-1)], uGTop[NVAR*(nx-1)], uGLeft[NVAR*(ny-1)];
    State BotVar[NVAR*(nx-1)], TopVar[NVAR*(nx-1)], RightVar[NVAR*(ny-1)], LeftVar[NVAR*(ny-1)];
    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[i];
        int iint = IJK(i,0,0, nx-1, 4);
        int iel = IJ(i, 0, nx-1);

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];
        //==========Ghost State
        BotVar[i].Initialize(&(uGBot[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,&(unk[iint]), ElemVar[iel],
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

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(0, j, 4,nx,6)];
        normy = geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        RightVar[j].Initialize(&(uGRight[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS,&(unk[iint]), ElemVar[iel],
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

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];
        //==========Ghost State
        TopVar[i].Initialize(&(uGTop[IJ(0,i,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, &(unk[iint]), ElemVar[iel],
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

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        LeftVar[j].Initialize(&(uGLeft[IJ(0,j,NVAR)]));
        boundary_state(btype,air,normx,normy,uFS, &(unk[iint]), ElemVar[iel],
                       &(uGLeft[IJ(0,j,NVAR)]));
        LeftVar[j].UpdateState(air);
    }

    //====================Evaluate Flux Contributions====================
    //dudt = sum(flux_in * face_length) / volume
    //==========Fully Interior Faces

    // I|xi Fluxes
    for (int i=1; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            // ~~~~~~~~~~ Inviscid fluxes ~~~~~~~~~~
            //face above point i,j
            double len, normx, normy, fflux[4];
            len = geofa[IJK(i,j,3,nx,6)];
            normx = geofa[IJK(i,j,4,nx,6)];
            normy = geofa[IJK(i,j,5,nx,6)];

            int iuL = IJK(i-1,j,0,nx-1,NVAR);
            int iuR = IJK(i  ,j,0,nx-1,NVAR);
            int ieL = IJ(i-1, j, nx-1);
            int ieR = IJ(i,   j, nx-1);
            State varL = ElemVar[ieL];
            State varR = ElemVar[ieR];

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL  ] -= len * fflux[0];
            rhsel[iuL+1] -= len * fflux[1];
            rhsel[iuL+2] -= len * fflux[2];
            rhsel[iuL+3] -= len * fflux[3];

            rhsel[iuR  ] += len * fflux[0];
            rhsel[iuR+1] += len * fflux[1];
            rhsel[iuR+2] += len * fflux[2];
            rhsel[iuR+3] += len * fflux[3];


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
        }
    }

    // J|eta fluxes
    for (int i=0; i<nx-1; i++){
        for (int j=1; j<ny-1; j++){
            //faces to the left/above point i,j
            double len, normx, normy, fflux[4];
            len = geofa[IJK(i,j,0,nx,6)];
            normx = geofa[IJK(i,j,1,nx,6)];
            normy = geofa[IJK(i,j,2,nx,6)];

            int iuL = IJK(i,j  ,0,nx-1,NVAR);
            int iuR = IJK(i,j-1,0,nx-1,NVAR);
            int ieL = IJ(i, j  , nx-1);
            int ieR = IJ(i, j-1, nx-1);
            State varL = ElemVar[ieL];
            State varR = ElemVar[ieR];
            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL  ] -= len * fflux[0];
            rhsel[iuL+1] -= len * fflux[1];
            rhsel[iuL+2] -= len * fflux[2];
            rhsel[iuL+3] -= len * fflux[3];

            rhsel[iuR  ] += len * fflux[0];
            rhsel[iuR+1] += len * fflux[1];
            rhsel[iuR+2] += len * fflux[2];
            rhsel[iuR+3] += len * fflux[3];

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
        }
    }

    for (int iu=0; iu<nelem*NVAR; iu++){
        if (_isnan(rhsel[iu])){
            printf("NAN RHSEL\n");
        }
    }

    //==========Boundary faces
    double len, normx, normy, fflux[4];
    // I|xi fluxes
    for (int j=0; j<ny-1; j++){
        //--------------------  LEFT BOUNDARY
        len   = geofa[IJK(0,j,3,nx,6)];
        normx = geofa[IJK(0,j,4,nx,6)];
        normy = geofa[IJK(0,j,5,nx,6)];
        int iuL = IJ(0,j,NVAR);
        int iuR =IJK(0,j,0,nx-1,NVAR);
        int ieL = IJ(0, j, nx-1);
        int ieR = IJ(0, j, nx-1);
        State varL = ElemVar[ieL];
        State varR = ElemVar[ieR];

        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(uGLeft[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuR  ] += len * fflux[0];
        rhsel[iuR+1] += len * fflux[1];
        rhsel[iuR+2] += len * fflux[2];
        rhsel[iuR+3] += len * fflux[3];

        // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
        double dc[2], vflux[6];
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



        //--------------------  RIGHT BOUNDARY
        len   = geofa[IJK(nx-1,j,3,nx,6)];
        normx = geofa[IJK(nx-1,j,4,nx,6)];
        normy = geofa[IJK(nx-1,j,5,nx,6)];
        iuL = IJK(nx-2,j,0,nx-1,NVAR);
        iuR = IJ(0,j,NVAR);
        ieL = IJ(nx-2,j, nx-1);
        ieR = j;
        varL = ElemVar[ieL];
        varR = RightVar[ieR];
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(unk[iuL]), varL, &(uGRight[iuR]), varR, &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuL  ] -= len * fflux[0];
        rhsel[iuL+1] -= len * fflux[1];
        rhsel[iuL+2] -= len * fflux[2];
        rhsel[iuL+3] -= len * fflux[3];

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
    }

    // J|eta fluxes
    for (int i=0; i<nx-1; i++) {
        //BOTTOM BOUNDARY
        len   = geofa[IJK(i,0,0,nx,6)];
        normx = geofa[IJK(i,0,1,nx,6)];
        normy = geofa[IJK(i,0,2,nx,6)];

        int iuL = IJK(i, 0, 0, nx - 1, NVAR);
        int iuR = IJ(0,i,NVAR);
        int ieL = IJ(i,0,nx-1);
        int ieR = i;
        State varL = ElemVar[ieL];
        State varR = BotVar[ieR];

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(unk[iuL]), varL, &(uGBot[iuR]), varR, &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuL]     -= len * fflux[0];
        rhsel[iuL + 1] -= len * fflux[1];
        rhsel[iuL + 2] -= len * fflux[2];
        rhsel[iuL + 3] -= len * fflux[3];

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



        //TOP BOUNDARY
        len   = geofa[IJK(i,ny-1,0,nx,6)];
        normx = geofa[IJK(i,ny-1,1,nx,6)];
        normy = geofa[IJK(i,ny-1,2,nx,6)];

        iuL = IJ(0,i,NVAR);
        iuR = IJK(i,ny-2,0,nx-1,NVAR);
        ieL = i;
        ieR = IJ(i, ny-2, nx-1);
        varL = TopVar[ieL];
        varR = ElemVar[ieR];

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(uGTop[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuR  ] += len * fflux[0];
        rhsel[iuR+1] += len * fflux[1];
        rhsel[iuR+2] += len * fflux[2];
        rhsel[iuR+3] += len * fflux[3];

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

            if (_isnan(dudt[iu])) {
                printf("nan dudt\n");
            }

        }
    }
    free(rhsel);
}