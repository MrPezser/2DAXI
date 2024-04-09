/* 2D Navier Stokes solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Feb2024
 */
#include <iostream>
#include <cmath>
#include "FileIO.h"
#include "Indexing.h"
#include "SpatialDiscretization.h"
#include "MeshModule.h"
#include "StateVariables.h"


double find_dt(double gam, int nx, int ny, double CFL, const double* uRef, double* geofa){
    double rho, u, v, rhoe, v2, p, c, vmax, dt, mindx;
    rho = uRef[0];
    u = uRef[1]/rho;
    v = uRef[2]/rho;
    rhoe = uRef[3];

    //Find speed of sound
    v2 = (u*u + v*v);
    p = (gam - 1) * (rhoe - (0.5 * rho * v2));
    c = sqrt(gam * p / rho);

    //Max info speed = v+c
    vmax = sqrt(v2) + c;

    //Min grid spacing
    mindx = 1.0;
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            mindx = fmin(mindx, geofa[IJK(i,j,0,nx,6)]);
            mindx = fmin(mindx, geofa[IJK(i,j,3,nx,6)]);
        }
    }

    return CFL * mindx / vmax;
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

void vec_copy(double n, double* a, double* b){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

int main() {
    //Read in setup file
    double gam, mu, mach, tol, CFL;
    int mxiter;
    gam =1.4;
    mu = 1e-5; // ~ 1/Re
    mach = 3.0;
    tol = 1e-6;
    mxiter = 1e6; //maximum number of iteration before stopping
    CFL = 0.3;

    printf("==================== Loading Mesh ====================\n");
    //==================== Load Mesh ====================
    int nx, ny, npoin, nb, nelem, nface;
    double *x, *y;
    int* ibound;
    double* geoel;
    double* geofa;

    //read in mesh file
    printf("Reading Mesh File..... \n");
    read_mesh(&nx, &ny, &ibound, &x, &y);
    npoin = nx*ny;
    nb = 2*nx + 2*ny;
    nelem = (nx-1)*(ny-1);
    nface = 2*nelem + nx + ny;

    //Find elem volume and centroid, face len and normals
    printf("Calculating Grid Metrics..... \n");
    calc_geoel_geofa(nx, ny, x, y, &geoel, &geofa);

    printf("==================== Initializing ====================\n");
    //==================== Setup for Sim ====================
    auto* unk    = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dudt   = (double*)malloc(NVAR*nelem*sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = 1.0;           //0.246319280397921945993707147708312;
    uFS[1] = 1.0;
    uFS[2] = 0.0;
    uFS[3] = 0.5 + 1 / (gam*(gam-1)*mach*mach);

    //Plenum State
    uBP[0] = 1.0;
    uBP[1] = 0.0;
    uBP[2] = 0.0;
    uBP[3] = uFS[3];//13.387171568521152;

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uFS[0];
        unk[NVAR*ielem+1] = uFS[1];
        unk[NVAR*ielem+2] = uFS[2];
        unk[NVAR*ielem+3] = uFS[3];
    }

    printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    print_state("Initial State", nx, ny, gam, x, y, unk, geoel);

    printf("Calculating Timestep..... \n");
    //Find timestep based off of CFL limit for initial condition (dt = CFL dx / c )
    double dt;

    printf("==================== Starting Solver ====================\n");
    State* ElemVar;//
    ElemVar = (State*)malloc(nelem*sizeof(State));// [nelem];
    //Set up structures for calculating/containing non-state variables on each element
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++) {
            int ie = IJ(i,j,nx-1);
            ElemVar[ie].Initialize(&(unk[IJK(i,j,0,nx-1,NVAR)]));
            ElemVar[ie].UpdateState(gam);
        }
    }



    double res0[NVAR]{};
    double res[4], ressum;
    int iter;
    for (iter=0; iter<mxiter; iter++){
        //Explicit Euler Time Integration
        dt = find_dt(gam, nx, ny, CFL, unk, geofa);
        calc_dudt(nx, ny, gam, mu, ElemVar, uFS, uBP, ibound, geoel, geofa, unk, dudt);

        for (int ielem=0; ielem<nelem; ielem++){
            int iu = NVAR*ielem;
            unk[iu  ] += dudt[iu  ]*dt;
            unk[iu+1] += dudt[iu+1]*dt;
            unk[iu+2] += dudt[iu+2]*dt;
            unk[iu+3] += dudt[iu+3]*dt;
            ElemVar[ielem].UpdateState(gam);
        }
        calculate_residual(nx, ny, dudt, res);
        if (iter==0) {
            for (int i=0; i<NVAR; i++){
                res0[i] = res[i];
            }
        }

        ressum = 0.0;
        for (int i=0; i<NVAR; i++){
            ASSERT(res0[i] > 0.0, "Negative or Zero Residual")
            ressum += res[i] / res0[i];
        }

        if (iter%50 == 0) {
            //printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e  \t\t Abs Residuals:  %12.2e%12.2e%12.2e%12.2e\n", \
                    iter, dt,ressum / res0, res[0], res[1], res[2], res[3]);
            printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e\n", \
                    iter, dt,ressum);
        }
        if (iter > 0 and iter%1000 == 0){
            printf("Saving current Solution\n");
            print_state("Final State", nx, ny, gam, x, y, unk, geoel);
        }
        if (ressum < tol) break;

    }
    printf("==================== Solution Found ====================\n");
    printf("Saving Solution File..... \n");

    print_state("Final State", nx, ny, gam, x, y, unk, geoel);


    printf("Complete.");

    free(ElemVar);
}
