//
// Created by Tsail on 5/7/2024.
//

#include <cerrno>
#include "NewtonLU.h"


int NLU(double* xgeo, double* ygeo, int nx, int ny,double CFL, Thermo& air, State* ElemVar, BC& bound, double *uFS, int* ibound, double* geoel,
         double* geofa, double* unkel){
    // funciton for carrying out the inexact newton iterations with LU linear solve
    ///THIS METHOD IS VERY INEFFICIENT AND SHOULD ONLY BE USED TO TEST OTHER PARTS OF CODE SUCH AS JACOBIAN CALCULATIONS
    double dt, ressum[NVAR]{}, rnorm, rnormnew, normmax, norm0, CFL0;
    int nelem, iter, mxitarmijo, nu, iarmijo;

    nelem = (nx-1)*(ny-1);
    nu = nelem*NVAR;
    auto dv     = (double*)malloc(nu*sizeof(double));
    auto x     = (double*)malloc(nu*sizeof(double));
    auto RHS    = (double*)malloc(nu*sizeof(double));

    CFL0 = CFL;

    // ==========  Initializing for Newton's method  ==========
    //Calculate the ghost cell values
    bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, unkel);
    //calculate the right hand side residual term (change of conserved quantities)
    calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, unkel, bound,RHS);
    calculate_residual(nx, ny, RHS, ressum);
    rnorm = norm2(ressum, NVAR);
    norm0 = rnorm;

    auto Jac = (double**)malloc((nu) * sizeof(double*));
    for (int k = 0; k < nu; k++)
        Jac[k] = (double*)malloc( (nu) * sizeof(double));

    //save residual history
    FILE* fres = fopen("../Outputs/res.tec", "w");
    if (fres == nullptr) {
        printf("~~~~~~~~~ Failed to save residual file, error:%d\n", errno);}
    else {
        fprintf(fres, "Residual history\n");
        fprintf(fres, "%d,\t%le,\t%le\n",0,1.0,CFL0);
    }

    for(iter=0; iter<MXITER; iter++) {
        //Calculate Pseudo-Timestep to be used for Jacobian
        //dt = find_dt(air, nx, ny, CFL, unkel, ElemVar[0], geofa);

        //Build jacobian Matrix
        BuildJacobian(nx, ny, CFL, air, ElemVar, uFS, ibound, geoel, geofa, unkel, RHS, bound, Jac);

        //Compute LU factorization
//        double *unkij = &(unkel[IJK(i, j, 0, nx-1, NVAR)]);
        double LUtol = 1e-16;
        vecinitialize(dv,nu);
        int N = nu;
        int P[N+1]; //permutation vector for pivoting

        //calculate LU decomposition
        DUNG("Before LU Decomp")
        ILUPDecompose(Jac, N, LUtol, P);

        //Solve LU system
        DUNG("Before LUSOLVE")
        LUPSolve(Jac, P, RHS, N, dv);

        //Update variables
        (void) vecadd(unkel, dv, nu, unkel);
        for (int i=0; i<nx-1; i++){
            for (int j=0; j<ny-1; j++) {
                int ie = IJ(i,j,nx-1);
                ElemVar[ie].UpdateState(air);
            }
        }


        bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, unkel);
        calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, unkel, bound,RHS);
        calculate_residual(nx, ny, RHS, ressum);
        rnorm = norm2(ressum, NVAR);

        if (iter%1 == 0) {
            printf("Iter:%7d\tCFL:%7.4e\tRelativeTotalResisual:  %8.5e\n", \
                    iter, CFL, rnorm/norm0);
        }
        if (iter >= 0 and iter%1 == 0){
            //printf("Saving current Solution\n\n");
            print_state("Final State", nx, ny, air, xgeo, ygeo, unkel, geoel);
        }
        if (rnorm/norm0 < RESTOL) {
            printf("==================== Solution Found ====================\n");
            printf("Saving Solution File..... \n");
            print_state("Final State", nx, ny, air, xgeo, ygeo, unkel, geoel);
            free(RHS);
            free(x);
            free(dv);
            fclose(fres);
            return 1;
        }
        fprintf(fres, "%d,\t%le,\t%le\n",
                iter+1, rnorm/norm0, CFL);
        CFL = fmin(1.0e5, CFL0 * 0.5 * (1.0 + norm0/rnorm));

    }
    free(RHS);
    free(x);
    free(dv);
    fclose(fres);
    return 0;
}