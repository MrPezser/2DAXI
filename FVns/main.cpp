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
#include "Thermo.h"
#include "InexactNewtonCg.h"
#include "NewtonLU.h"


void vec_copy(double n, double* a, double* b){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

int main() {
    //Read in setup file
    double p0, u0, tol, CFL, T0, v0, rho0;
    int mxiter;
    tol = 1e-6;
    mxiter = 1e6; //maximum number of iteration before stopping
    CFL = 10.0;//0.8;
    u0 = 10.0;
    T0 = 300;
    rho0 = 1.0;
    v0 = 0.0;

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
    auto* unk  = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* res  = (double*)malloc(NVAR*nelem * sizeof(double));
    auto* dv   = (double*)malloc(NVAR*nelem * sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = rho0;
    uFS[1] = u0;
    uFS[2] = v0;
    uFS[3] = T0;

    //Plenum State
    ///back pressure here

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uFS[0];
        unk[NVAR*ielem+1] = uFS[1];
        unk[NVAR*ielem+2] = uFS[2];
        unk[NVAR*ielem+3] = uFS[3];
    }

    Thermo air = Thermo();

    printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    print_state("Initial State", nx, ny, air, x, y, unk, geoel);

    //Find timestep based off of CFL limit for initial condition (dt = CFL dx / c )
    double dt;
    State* ElemVar;//
    ElemVar = (State*)malloc(nelem*sizeof(State));// [nelem];
    //Set up structures for calculating/containing non-state variables on each element
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++) {
            int ie = IJ(i,j,nx-1);
            double* unkel = &(unk[IJK(i,j,0,nx-1,NVAR)]);

            unkel[0] = uFS[0];
            unkel[1] = uFS[1];
            unkel[2] = uFS[2];
            unkel[3] = uFS[3];

            ElemVar[ie].Initialize(unkel);
            ElemVar[ie].UpdateState(air);
        }
    }
    //Same memory to be used for each local matrix (chg this if making parallel) (explicit only)
    //auto D = (double**)malloc((NVAR) * sizeof(double*));
    //for (int k = 0; k < NVAR; k++)
    //    D[k] = (double*)malloc( (NVAR) * sizeof(double));

    printf("==================== Starting Solver ====================\n");
    BC bound = BC(nx,ny);

    //int isolved = INCG(x, y, nx, ny, CFL, air, ElemVar, bound, uFS,
    //                   ibound, geoel, geofa, unk);
    int isolved = NLU(x, y, nx, ny, CFL, air, ElemVar, bound, uFS, ibound, geoel, geofa, unk);

    if (isolved==0){
        printf("Failed to find satifactory solution\n");
    }

    /*
    double res0[NVAR]{};
    double ressum[NVAR], restotal;
    int iter;
    BC bound = BC(nx,ny);
    for (iter=0; iter<mxiter; iter++){
        //Explicit Euler Time Integration

        //Find global timestep based off of CFl condition
        dt = find_dt(air, nx, ny, CFL, unk, ElemVar[0], geofa);

        //Calculate the ghost cell values
        bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, unk);

        //calculate the right hand side residual term (change of conserved quantities)
        calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, unk, bound,res);
        calculate_residual(nx, ny, res, ressum);

        //========== Solve linear system on each element (turns chg in conservatives to change in solution variables)
        int flg = 0;
        for (int i=0; i<nx-1; i++) {
            for (int j=0; j<ny-1; j++) {
                double *unkij = &(unk[IJK(i, j, 0, nx-1, NVAR)]);
                double LUtol = 1e-16;
                int iel = IJ(i,j,nx-1);
                int N = NVAR;
                int P[NVAR+1]{}; //permutation vector for pivoting

                //Evaluate the jacobian / Implicit matrix
                RegularizationTerm(dt, unkij, ElemVar[iel], D);

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
            ElemVar[ielem].UpdateState(air);
        }



        if (iter==0) {
            for (int i=0; i<NVAR; i++){
                res0[i] = ressum[i];
            }
        }
        restotal = 0.0;
        for (int i=0; i<NVAR; i++){
            ASSERT(ressum[i] > 0.0, "Nonpositive Residual")
            if (res0[i] < 1e-10) res0[i] = ressum[i];
            restotal += ressum[i] / res0[i];
        }
        if (iter%100 == 0) {
            printf("Iter:%7d\tdt:%7.4e \t\t RelativeTotalResisual:  %8.5e\n", \
                    iter, dt, restotal);
        }
        if (iter > 0 and iter%1000 == 0){
            printf("Saving current Solution\n");
            print_state("Final State", nx, ny, air, x, y, unk, geoel);
        }
        if (restotal < tol) break;

    }
    printf("==================== Solution Found ====================\n");
    printf("Saving Solution File..... \n");

    print_state("Final State", nx, ny, air, x, y, unk, geoel);


    printf("Complete.");
    */

    free(ElemVar);
    free(res);
    free(dv);
    free(unk);
    free(geoel);
    free(geofa);
    free(ibound);
    free(x);
    free(y);
}
