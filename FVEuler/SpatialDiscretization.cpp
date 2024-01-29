//
// Created by Tsail on 1/28/2024.
//

#include <cstdio>
#include <cmath>
#include "SpatialDiscretization.h"
#include "Indexing.h"
#include "BoundaryConditions.h"
#include "EulerFlux.h"


void calc_dudt(int nx, int ny, double gam, double *uFS, int* ibound, double* geoel, double* geofa, double* unk, double* dudt) {
    int nelem = (nx-1)*(ny-1);
    double rhsel[NVAR*nelem];
    for(int i=0;i<NVAR*nelem;i++) rhsel[i] = 0.0;

    //Calculate boundary cell state (ghost state)
    double uGBot[NVAR*(nx-1)], uGRight[NVAR*(ny-1)], uGTop[NVAR*(nx-1)], uGLeft[NVAR*(ny-1)];
    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[i];

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(i, 0, 1,nx,6)];
        normy = geofa[IJK(i, 0, 2,nx,6)];
        //==========Ghost State
        boundary_state(btype,gam,normx,normy,uFS,&(unk[IJK(i,0,0, nx-1, 4)]),&(uGBot[IJ(0,i,NVAR)]));
    }

    //right side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        btype = ibound[j+ nx-1];

        //==========Face Normal
        double normx, normy;
        normx = geofa[IJK(0, j, 4,nx,6)];
        normy = geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        boundary_state(btype,gam,normx,normy,uFS,&(unk[IJK(nx-2,j,0, nx-1, NVAR)]),&(uGRight[IJ(0,j,NVAR)]));
    }

    //top side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        int btype;
        int ib = nx-2-i;
        btype = ibound[ib+nx+ny-2];

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(i, ny-1, 1,nx,6)];
        normy = -geofa[IJK(i, ny-1, 2,nx,6)];
        //==========Ghost State
        boundary_state(btype,gam,normx,normy,uFS, &(unk[IJK(i,ny-2,0, nx-1, NVAR)]), &(uGTop[IJ(0,i,NVAR)]));
    }

    //left side of domain
    for (int j=0; j<(ny-1); j++){
        // left state = interior, right state = ghost
        int btype;
        int jb = (ny-2)-j;
        btype = ibound[jb+(2*nx)+ny-3];

        //==========Face Normal
        double normx, normy;
        normx = -geofa[IJK(0, j, 4,nx,6)];
        normy = -geofa[IJK(0, j, 5,nx,6)];
        //==========Ghost State
        boundary_state(btype,gam,normx,normy,uFS, &(unk[IJK(0,j,0, nx-1, 4)]), &(uGLeft[IJ(0,j,NVAR)]));
    }

    //====================Evaluate Flux Contributions====================
    //dudt = sum(flux_in * face_length) / volume
    //==========Fully Interior Faces
    // I|xi fluxes
    for (int i=1; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            //faces to the left/above point i,j
            double len, normx, normy, fflux[4];
            len = geofa[IJK(i,j,3,nx,6)];
            normx = geofa[IJK(i,j,4,nx,6)];
            normy = geofa[IJK(i,j,5,nx,6)];

            int iuL = IJK(i-1,j,0,nx-1,NVAR);
            int iuR = IJK(i  ,j,0,nx-1,NVAR);

            //Find interface flux
            LeerFlux(gam, normx, normy, &(unk[iuL]), &(unk[iuR]), &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL  ] -= len * fflux[0];
            rhsel[iuL+1] -= len * fflux[1];
            rhsel[iuL+2] -= len * fflux[2];
            rhsel[iuL+3] -= len * fflux[3];

            rhsel[iuR  ] += len * fflux[0];
            rhsel[iuR+1] += len * fflux[1];
            rhsel[iuR+2] += len * fflux[2];
            rhsel[iuR+3] += len * fflux[3];
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
            //Find interface flux
            LeerFlux(gam, normx, normy, &(unk[iuL]), &(unk[iuR]), &(fflux[0]));

            //Add flux contribution to elements
            rhsel[iuL  ] -= len * fflux[0];
            rhsel[iuL+1] -= len * fflux[1];
            rhsel[iuL+2] -= len * fflux[2];
            rhsel[iuL+3] -= len * fflux[3];

            rhsel[iuR  ] += len * fflux[0];
            rhsel[iuR+1] += len * fflux[1];
            rhsel[iuR+2] += len * fflux[2];
            rhsel[iuR+3] += len * fflux[3];
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
        //LEFT BOUNDARY
        len   = geofa[IJK(0,j,3,nx,6)];
        normx = geofa[IJK(0,j,4,nx,6)];
        normy = geofa[IJK(0,j,5,nx,6)];

        int iuL = IJ(0,j,NVAR);
        int iuR = IJK(0  ,j,0,nx-1,NVAR);
        LeerFlux(gam, normx, normy, &(uGLeft[iuL]), &(unk[iuR]), &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuR  ] += len * fflux[0];
        rhsel[iuR+1] += len * fflux[1];
        rhsel[iuR+2] += len * fflux[2];
        rhsel[iuR+3] += len * fflux[3];

        //RIGHT BOUNDARY
        len   = geofa[IJK(nx-1,j,3,nx,6)];
        normx = geofa[IJK(nx-1,j,4,nx,6)];
        normy = geofa[IJK(nx-1,j,5,nx,6)];

        iuL = IJK(nx-2,j,0,nx-1,NVAR);
        iuR = IJ(0,j,NVAR);
        LeerFlux(gam, normx, normy, &(unk[iuL]), &(uGRight[iuR]), &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuL  ] -= len * fflux[0];
        rhsel[iuL+1] -= len * fflux[1];
        rhsel[iuL+2] -= len * fflux[2];
        rhsel[iuL+3] -= len * fflux[3];
    }

    // J|eta fluxes
    for (int i=0; i<nx-1; i++) {
        //BOTTOM BOUNDARY
        len   = geofa[IJK(i,0,0,nx,6)];
        normx = geofa[IJK(i,0,1,nx,6)];
        normy = geofa[IJK(i,0,2,nx,6)];

        int iuL = IJK(i, 0, 0, nx - 1, NVAR);
        int iuR = IJ(0,i,NVAR);

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        LeerFlux(gam, normx, normy, &(unk[iuL]), &(uGBot[iuR]), &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuL]     -= len * fflux[0];
        rhsel[iuL + 1] -= len * fflux[1];
        rhsel[iuL + 2] -= len * fflux[2];
        rhsel[iuL + 3] -= len * fflux[3];

        //TOP BOUNDARY
        len   = geofa[IJK(i,ny-1,0,nx,6)];
        normx = geofa[IJK(i,ny-1,1,nx,6)];
        normy = geofa[IJK(i,ny-1,2,nx,6)];

        iuL = IJ(0,i,NVAR);
        iuR = IJK(i,ny-2,0,nx-1,NVAR);

        if (iuR<0 or iuL < 0){
            printf("wut\n");
        }
        if ((iuR > nelem*NVAR) or (iuL > nelem*NVAR)){
            printf("wut\n");
        }

        //Find interface flux
        LeerFlux(gam, normx, normy, &(uGTop[iuL]), &(unk[iuR]), &(fflux[0]));

        //Add flux contribution to elements
        rhsel[iuR  ] += len * fflux[0];
        rhsel[iuR+1] += len * fflux[1];
        rhsel[iuR+2] += len * fflux[2];
        rhsel[iuR+3] += len * fflux[3];
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

}