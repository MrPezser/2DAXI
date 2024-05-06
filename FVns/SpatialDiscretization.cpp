//
// Created by Tsail on 1/28/2024.
//

#include "SpatialDiscretization.h"


double find_dt(Thermo& air, int nx, int ny, double CFL, const double* uRef, State& var, double* geofa){
    double vmax, dt, mindx;

    //Max info speed = v+c
    vmax = sqrt(var.v2) + var.a;

    //Min grid spacing
    mindx = 1.0;
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            mindx = fmin(mindx, geofa[IJK(i,j,0,nx,6)]);
            mindx = fmin(mindx, geofa[IJK(i,j,3,nx,6)]);
        }
    }

    dt = CFL * mindx / vmax;
    ASSERT(!_isnan(dt), "NAN dt")
    return dt;
}

void calculate_residual(int nx, int ny, double* res, double* ressum){
    double res2[NVAR] = {0.0};

    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            int iu = IJK(i,j,0,nx-1,NVAR);

            for (int k=0; k<NVAR; k++){
                ASSERT(!_isnan(res[iu+k]),"res NaN")
                res2[k] += res[iu+k]*res[iu+k];
            }
        }
    }

    ressum[0] = sqrt(res2[0]);
    ressum[1] = sqrt(res2[1]);
    ressum[2] = sqrt(res2[2]);
    ressum[3] = sqrt(res2[3]);
}

void viscous(int nx, double normy, double normx, double* uLeft, State& varL, double* uRight, State varR, double* dc, double* visc_contrib){
    if (IVISC) {
        // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
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
               double* geofa, double* unk, BC& bound, double* dudt) {
    int nelem = (nx-1)*(ny-1);
    double* rhsel;
    rhsel = (double*)malloc(NVAR*nelem*sizeof(double));
    for(int i=0; i<NVAR*nelem; i++) {
        rhsel[i] = 0.0;
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
        LDFSS(normx, normy, &(bound.uGLeft[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

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
        viscous(nx, normy, normx, &(bound.uGLeft[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

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
        varR = bound.RightVar[ieR];
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(unk[iuL]), varL, &(bound.uGRight[iuR]), varR, &(fflux[0]));

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
        viscous(nx, normy, normx, &(unk[iuL]), varL, &(bound.uGRight[iuR]), varR, &(dc[0]), &(vflux[0]));

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
        State varR = bound.BotVar[ieR];

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(unk[iuL]), varL, &(bound.uGBot[iuR]), varR, &(fflux[0]));

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
        viscous(nx, normy, normx, &(unk[iuL]), varL, &(bound.uGBot[iuR]), varR, &(dc[0]), &(vflux[0]));

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
        varL = bound.TopVar[ieL];
        varR = ElemVar[ieR];

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
        LDFSS(normx, normy, &(bound.uGTop[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

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
        viscous(nx, normy, normx, &(bound.uGTop[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

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

void calc_dudt_element(int iel, int jel, int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* unk, BC& bound, double* dudt) {
    //====================Evaluate Flux Contributions====================
    // !!!! This just evaluates fluxes on a single element and adjacent elements
    //dudt = sum(flux_in * face_length) / volume

    int nelem = (nx-1)*(ny-1);
    double* rhsel;
    rhsel = (double*)malloc(NVAR*nelem*sizeof(double));
    for(int i=0; i<NVAR*nelem; i++) {
        rhsel[i] = 0.0;
    }
    /*
    //center element
    int ielem = IJ(iel,jel,nx-1);
    int iunkel = IJK(iel,jel,0,nx-1,NVAR);

    //element below
    if (jel > 0) {
        int ielemjm = IJ(iel, jel-1, nx-1);
        int iunkeljm = IJK(iel, jel-1, 0,nx-1, NVAR);
    } else {
        int ielmjm = iel;
        int iunkeljm = IJ(0,iel,NVAR);
    }

    //element above
    if (jel < ny-2) {
        int ielemjp  = IJ(iel, jel+1, nx-1);
        int iunkeljp = IJK(iel, jel+1, 0,nx-1, NVAR);
    } else {
        int ielmjp   = iel;
        int iunkeljp = IJ(0,iel,NVAR);
    }

    int ielem = IJ(iel,jel,nx-1);
    int iunkel = IJK(iel,jel,0,nx-1,NVAR);

    int ielem = IJ(iel,jel,nx-1);
    int iunkel = IJK(iel,jel,0,nx-1,NVAR);
    */

    int ilow = iel-2;
    int ihigh= iel+2;
    int jlow = jel-2;
    int jhigh= jel+2;

    //==========Fully Interior Faces
    // I|xi Fluxes
    for (int i = ilow; i <= ihigh; i++) {
        if (i < 1 or i > nx-2) {continue;}

        for (int j = jlow; j <= jhigh; j++) {
            if (j < 0 or j > ny-2) {continue;}
            // ~~~~~~~~~~ Inviscid fluxes ~~~~~~~~~~
            //face above point i,j
            double len, normx, normy, fflux[4];
            len = geofa[IJK(i, j, 3, nx, 6)];
            normx = geofa[IJK(i, j, 4, nx, 6)];
            normy = geofa[IJK(i, j, 5, nx, 6)];

            int iuL = IJK(i - 1, j, 0, nx - 1, NVAR);
            int iuR = IJK(i, j, 0, nx - 1, NVAR);
            int ieL = IJ(i - 1, j, nx - 1);
            int ieR = IJ(i, j, nx - 1);
            State varL = ElemVar[ieL];
            State varR = ElemVar[ieR];

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL] -= len * fflux[0];
            rhsel[iuL + 1] -= len * fflux[1];
            rhsel[iuL + 2] -= len * fflux[2];
            rhsel[iuL + 3] -= len * fflux[3];

            rhsel[iuR] += len * fflux[0];
            rhsel[iuR + 1] += len * fflux[1];
            rhsel[iuR + 2] += len * fflux[2];
            rhsel[iuR + 3] += len * fflux[3];


            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double vflux[6], dc[2];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i, j, 1, nx - 1, 3)] - geoel[IJK(i - 1, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, j, 2, nx - 1, 3)] - geoel[IJK(i - 1, j, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];

            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
        }
    }

    // J|eta fluxes
    for (int i = ilow; i <= ihigh; i++) {
        if (i < 0 or i > nx-2) {continue;}

        for (int j = jlow; j <= jhigh; j++) {
            if (j < 1 or j > ny-2) {continue;}

            //faces to the left/above point i,j
            double len, normx, normy, fflux[4];
            len = geofa[IJK(i, j, 0, nx, 6)];
            normx = geofa[IJK(i, j, 1, nx, 6)];
            normy = geofa[IJK(i, j, 2, nx, 6)];

            int iuL = IJK(i, j, 0, nx - 1, NVAR);
            int iuR = IJK(i, j - 1, 0, nx - 1, NVAR);
            int ieL = IJ(i, j, nx - 1);
            int ieR = IJ(i, j - 1, nx - 1);
            State varL = ElemVar[ieL];
            State varR = ElemVar[ieR];
            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL] -= len * fflux[0];
            rhsel[iuL + 1] -= len * fflux[1];
            rhsel[iuL + 2] -= len * fflux[2];
            rhsel[iuL + 3] -= len * fflux[3];

            rhsel[iuR] += len * fflux[0];
            rhsel[iuR + 1] += len * fflux[1];
            rhsel[iuR + 2] += len * fflux[2];
            rhsel[iuR + 3] += len * fflux[3];

            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double vflux[6], dc[2];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i, j - 1, 1, nx - 1, 3)] - geoel[IJK(i, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, j - 1, 2, nx - 1, 3)] - geoel[IJK(i, j, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];

            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
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
    if (iel == 0 or iel == 1) {
        for (int j = jlow; j <= jhigh; j++) {
            if (j < 0 or j > (ny - 2)) { continue; }
            //--------------------  LEFT BOUNDARY
            len = geofa[IJK(0, j, 3, nx, 6)];
            normx = geofa[IJK(0, j, 4, nx, 6)];
            normy = geofa[IJK(0, j, 5, nx, 6)];
            int iuL = IJ(0, j, NVAR);
            int iuR = IJK(0, j, 0, nx - 1, NVAR);
            int ieL = IJ(0, j, nx - 1);
            int ieR = IJ(0, j, nx - 1);
            State varL = ElemVar[ieL];
            State varR = ElemVar[ieR];

            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(bound.uGLeft[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuR] += len * fflux[0];
            rhsel[iuR + 1] += len * fflux[1];
            rhsel[iuR + 2] += len * fflux[2];
            rhsel[iuR + 3] += len * fflux[3];

            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(1, j, 1, nx - 1, 3)] - geoel[IJK(0, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(1, j, 2, nx - 1, 3)] - geoel[IJK(0, j, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(bound.uGLeft[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
        }
    }

    if (iel == nx-2 or iel == nx-1) {
        for (int j = jlow; j <= jhigh; j++) {
            if (j < 0 or j > (ny - 2)) { continue; }
            //--------------------  RIGHT BOUNDARY
            len = geofa[IJK(nx - 1, j, 3, nx, 6)];
            normx = geofa[IJK(nx - 1, j, 4, nx, 6)];
            normy = geofa[IJK(nx - 1, j, 5, nx, 6)];
            int iuL = IJK(nx - 2, j, 0, nx - 1, NVAR);
            int iuR = IJ(0, j, NVAR);
            int ieL = IJ(nx - 2, j, nx - 1);
            int ieR = j;
            State varL = ElemVar[ieL];
            State varR = bound.RightVar[ieR];
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(bound.uGRight[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL] -= len * fflux[0];
            rhsel[iuL + 1] -= len * fflux[1];
            rhsel[iuL + 2] -= len * fflux[2];
            rhsel[iuL + 3] -= len * fflux[3];

            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(nx - 2, j, 1, nx - 1, 3)] - geoel[IJK(nx - 3, j, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(nx - 2, j, 2, nx - 1, 3)] - geoel[IJK(nx - 3, j, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(bound.uGRight[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];
        }
    }

    // J|eta fluxes
    if (jel == 0 or jel == 1) {
        for (int i = ilow; i <= ihigh; i++) {
            if (i < 0 or i > (nx - 2)) { continue; }
            //BOTTOM BOUNDARY
            len = geofa[IJK(i, 0, 0, nx, 6)];
            normx = geofa[IJK(i, 0, 1, nx, 6)];
            normy = geofa[IJK(i, 0, 2, nx, 6)];

            int iuL = IJK(i, 0, 0, nx - 1, NVAR);
            int iuR = IJ(0, i, NVAR);
            int ieL = IJ(i, 0, nx - 1);
            int ieR = i;
            State varL = ElemVar[ieL];
            State varR = bound.BotVar[ieR];

            if (iuR < 0 or iuL < 0) {
                printf("wut\n");
            }
            if ((iuR > nelem * NVAR) or (iuL > nelem * NVAR)) {
                printf("wut\n");
            }

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(unk[iuL]), varL, &(bound.uGBot[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL] -= len * fflux[0];
            rhsel[iuL + 1] -= len * fflux[1];
            rhsel[iuL + 2] -= len * fflux[2];
            rhsel[iuL + 3] -= len * fflux[3];

            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i, 0, 1, nx - 1, 3)] - geoel[IJK(i, 1, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, 0, 2, nx - 1, 3)] - geoel[IJK(i, 1, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(unk[iuL]), varL, &(bound.uGBot[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuL + 1] += len * vflux[0];
            rhsel[iuL + 2] += len * vflux[1];
            rhsel[iuL + 3] += len * vflux[2];
        }
    }
    if (jel == ny-2 or jel == ny-1) {
        for (int i = ilow; i <= ihigh; i++) {
            if (i < 0 or i > (nx - 2)) { continue; }
            //TOP BOUNDARY
            len = geofa[IJK(i, ny - 1, 0, nx, 6)];
            normx = geofa[IJK(i, ny - 1, 1, nx, 6)];
            normy = geofa[IJK(i, ny - 1, 2, nx, 6)];

            int iuL = IJ(0, i, NVAR);
            int iuR = IJK(i, ny - 2, 0, nx - 1, NVAR);
            int ieL = i;
            int ieR = IJ(i, ny - 2, nx - 1);
            State varL = bound.TopVar[ieL];
            State varR = ElemVar[ieR];

            if (iuR < 0 or iuL < 0) {
                printf("wut\n");
            }
            if ((iuR > nelem * NVAR) or (iuL > nelem * NVAR)) {
                printf("wut\n");
            }

            //Find interface flux
            //ASSERT(varR.a*varL.a > 0.0, "nonpositive wave speed")
            LDFSS(normx, normy, &(bound.uGTop[iuL]), varL, &(unk[iuR]), varR, &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuR] += len * fflux[0];
            rhsel[iuR + 1] += len * fflux[1];
            rhsel[iuR + 2] += len * fflux[2];
            rhsel[iuR + 3] += len * fflux[3];

            // ~~~~~~~~~~ Viscous fluxes ~~~~~~~~~~
            double dc[2], vflux[6];
            //mirror the next interior cell to the boundary for that flux contrib
            dc[0] = geoel[IJK(i, ny - 3, 1, nx - 1, 3)] - geoel[IJK(i, ny - 2, 1, nx - 1, 3)];
            dc[1] = geoel[IJK(i, ny - 3, 2, nx - 1, 3)] - geoel[IJK(i, ny - 2, 2, nx - 1, 3)];
            varL = ElemVar[ieL];
            varR = ElemVar[ieR];
            viscous(nx, normy, normx, &(bound.uGTop[iuL]), varL, &(unk[iuR]), varR, &(dc[0]), &(vflux[0]));

            //Add flux contribution to elements
            rhsel[iuR + 1] -= len * vflux[3];
            rhsel[iuR + 2] -= len * vflux[4];
            rhsel[iuR + 3] -= len * vflux[5];
        }
    }

    //====================Combine to Find du/dt====================
    for (int i=ilow; i<=ihigh; i++) {
        for (int j = jlow; j <= jhigh; j++) {
            if (i < 0 or i > nx - 2) { continue; }
            if (j < 0 or j > ny - 2) { continue; }

            int iu = IJK(i, j, 0, nx - 1, NVAR);
            double vol = geoel[IJK(i, j, 0, nx - 1, 3)];
            dudt[iu] = rhsel[iu] / vol;
            dudt[iu + 1] = rhsel[iu + 1] / vol;
            dudt[iu + 2] = rhsel[iu + 2] / vol;
            dudt[iu + 3] = rhsel[iu + 3] / vol;

            if (_isnan(dudt[iu])) {
                printf("nan dudt\n");
            }

        }
    }
    free(rhsel);
}