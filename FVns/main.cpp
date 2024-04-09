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
#include "LUtools.h"
#include "Jacobian.h"


double find_dt(double gam, int nx, int ny, double CFL, const double* uRef, State& var, double* geofa){
    double rho, u, v, rhoe, v2, p, c, vmax, dt, mindx;
    rho = uRef[0];
    u = uRef[1];
    v = uRef[2];
    rhoe = rho * var.h;

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
    auto* res   = (double*)malloc(NVAR * nelem * sizeof(double));
    auto* dv   = (double*)malloc(NVAR * nelem * sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = 0.5;
    uFS[1] = 500.0;
    uFS[2] = 0.0;
    uFS[3] = 350;

    //Plenum State
    ///back pressure here

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uFS[0];
        unk[NVAR*ielem+1] = uFS[1];
        unk[NVAR*ielem+2] = uFS[2];
        unk[NVAR*ielem+3] = uFS[3];
    }

    printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    print_state("Initial State", nx, ny, gam, x, y, unk, geoel);

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
    //Same memory to be used for each local matrix (chg this if making parallel)
    auto D = (double**)malloc((NSP+3) * sizeof(double*));
    for (int isp = 0; isp < NSP+3; isp++)
        D[isp] = (double*)malloc( (NSP+3) * sizeof(double));


    double res0[NVAR]{};
    double ressum[4], restotal;
    int iter;
    for (iter=0; iter<mxiter; iter++){
        //Explicit Euler Time Integration

        //Find global timestep based off of CFl condition
        dt = find_dt(gam, nx, ny, CFL, unk, ElemVar[0], geofa);

        //calculate the right hand side residual term (change of conserved quantities)
        calc_dudt(nx, ny, gam, mu, ElemVar, uFS, uBP, ibound, geoel, geofa, unk, res);
        calculate_residual(nx, ny, res, ressum);

        //========== Solve linear system on each element (turns chg in conservatives to change in solution variables)
        int flg = 0;
        for (int i=0; i<nx-1; i++) {
            for (int j=0; j<ny-1; j++) {
                double *unkij = &(unk[IJK(i, j, 0, nx-1, NVAR)]);
                double LUtol = 1e-16;
                int iel = IJ(i,j,nx-1);
                int N = NVAR;
                int P[NVAR]{}; //permutation vector for pivoting

                //Evaluate the jacobian / Implicit matrix
                BuildJacobian(gam, dt, unkij, ElemVar[iel], D);

                //get the rhs block needed
                double *b = &(res[IJK(i, j, 0, nx - 1, NVAR)]);
                double *xLU = &(dv[ IJK(i, j, 0, nx - 1, NVAR)]);
                LUPDecompose(D, N, LUtol, P);
                LUPSolve(D, P, b, N, xLU);
            }
        }
        //perform iteration
        for (int ielem=0; ielem<nelem; ielem++){
            int iu = NVAR*ielem;
            unk[iu  ] += dv[iu  ];
            unk[iu+1] += dv[iu + 1];
            unk[iu+2] += dv[iu + 2];
            unk[iu+3] += dv[iu + 3];
            ElemVar[ielem].UpdateState(gam);
        }



        if (iter==0) {
            for (int i=0; i<NVAR; i++){
                res0[i] = ressum[i];
            }
        }
        restotal = 0.0;
        for (int i=0; i<NVAR; i++){
            ASSERT(res0[i] > 0.0, "Nonpositive Residual")
            restotal += ressum[i] / res0[i];
        }
        if (iter%50 == 0) {
            printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e\n", \
                    iter, dt, restotal);
        }
        if (iter > 0 and iter%1000 == 0){
            printf("Saving current Solution\n");
            print_state("Final State", nx, ny, gam, x, y, unk, geoel);
        }
        if (restotal < tol) break;

    }
    printf("==================== Solution Found ====================\n");
    printf("Saving Solution File..... \n");

    print_state("Final State", nx, ny, gam, x, y, unk, geoel);


    printf("Complete.");

    free(ElemVar);
    free(res);
    free(dv);
}
