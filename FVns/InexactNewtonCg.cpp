//
// Created by Tsail on 5/5/2024.
//


#include "InexactNewtonCg.h"
#include "LUtools.h"

int CGsolve(int nx, int ny,double dt, Thermo& air, State* ElemVar, BC& bound, double *uFS, int* ibound, double* geoel,
            double* geofa, double* unkel, double *RHS, double *x) {
    //Solves the Ax = b system using the conjugate gradient method
    // A = function to multiply a vector by the Jacobian of the RHS
    // b = -Right Hand Side
    // x = guess at a descent direction
    int n = (nx-1)*(ny-1)*NVAR;
    double tol = 1e-4;
    unsigned kmx = 20;//n;
    double r[n], s[n], omega[n], dummy[n];
    vecinitialize(r,n);

    JacobianVectorMultiply(nx,ny,dt,air,ElemVar,uFS,ibound,geoel,geofa,unkel,RHS,bound,x,r);
    for (int iu=0; iu<n; iu++){
        r[iu] = RHS[iu] - r[iu];
    }
    //vecinitialize(x,n);
    //(void) veccopy(r, RHS, n);
    (void) veccopy(s, r, n);

    double rho0 = vecinner(r, r, n);
    double rhok = rho0;
    double rhokm1 = rho0;

    printf("CG Start: %d \t\t |R0|: %e\t |PHI|: %le\n", 0, sqrt(rhok), vecinner(x,x,n));

    for (int k=1;k<=kmx;k++) {

        ///    matvec(A, s, n, omega);   //omega = As
        JacobianVectorMultiply(nx,ny,dt,air,ElemVar,uFS,ibound,geoel,geofa,unkel,RHS,bound,s,omega);

        double s_omega = vecinner(s, omega, n);
        //printf("s_omega = %e\n",s_omega);
        if( s_omega < 0.0) {
            printf("CG failed,%d\t",k);//, CGITER=%d\n", k);
            //ASSERT(k>1,"CG FAILED ON FIRST ITERATION, NO DESCENT DIRECTION FOUND")
            return 0;
        }

        double alpha = rhok /s_omega;
        for (int iu=0; iu<n; iu++){
            x[iu] += (alpha * s[iu]);
        }

        for (int iu=0; iu<n; iu++){
            r[iu] += -(alpha * omega[iu]);
        }

        rhokm1 = rhok;
        rhok   = vecinner(r,r,n);
        if (k % 1 == 0) {
            printf("CG Iteration: %d \t |Res|: %e \t |PHI|: %le\n", k, sqrt(rhok/rho0) , vecinner(x,x,n));
        }

        if (sqrt(rhok/rho0) < tol) {
            printf("Conjugate Gradient exit in %d iterations. Res: %e\n", k, sqrt(rhok/rho0));
            return 1;
        }

        double beta = rhok / rhokm1;

        if ( beta > 1.0) {
            for (int iu=0; iu<n; iu++){
                x[iu] -= (alpha * s[iu]);
            }
            //printf("CG failed, Beta > 1, it:%d\n",k);
            printf("%d\t",k);
            return 0;
        }

        for (int iu=0; iu<n; iu++){
            s[iu] = r[iu] +  (beta * s[iu]);
        }
    }
    printf("Conjugate Gradient Did not converge. Iter %d\t Res: %e\n", kmx, sqrt(rhok));
    return 0;
}

int INCG(double* xgeo, double* ygeo, int nx, int ny,double CFL, Thermo& air, State* ElemVar, BC& bound, double *uFS, int* ibound, double* geoel,
          double* geofa, double* unkel){
    // funciton for carrying out the inexact newton-CG method with line search
    double dt, eta, alpha, alp0, ressum[NVAR]{}, rnorm, rnormnew, normmax, norm0;
    int nelem, iter, mxitarmijo, nu, iarmijo;

    nelem = (nx-1)*(ny-1);
    nu = nelem*NVAR;
    auto dv     = (double*)malloc(nu*sizeof(double));
    auto x     = (double*)malloc(nu*sizeof(double));
    auto RHS    = (double*)malloc(nu*sizeof(double));
    auto RHSnew = (double*)malloc(nu*sizeof(double));
    auto unknew = (double*)malloc(nu*sizeof(double));
    auto ELVanew= (State*)malloc(nelem*sizeof(State));



    //Parameter defining sufficient decrease based off of the initial jacobian for each iteration
    ///Hijacked eta to be a relaxation of the requirement for early iterations
    eta = 1.0;//0.5;
    //Starting stepsize for armijo line search, lowering for initial transients might help,
    // want to be ==1 in region of local convergence
    alp0 = 1.0;
    //Max number of iterations for armijo linesearch
    mxitarmijo = 100;

    // ==========  Initializing for Newton's method  ==========
    //Calculate the ghost cell values
    bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, unkel);
    //calculate the right hand side residual term (change of conserved quantities)
    calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, unkel, bound,RHS);
    calculate_residual(nx, ny, RHS, ressum);
    rnorm = norm2(ressum, NVAR);
    norm0 = rnorm;

    for(iter=0; iter<MXITER; iter++) {
        //Calculate Pseudo-Timestep to be used for Jacobian
        dt = find_dt(air, nx, ny, CFL, unkel, ElemVar[0], geofa);

        //Find a descent direction used conjugate gradient
        //(void) vecinitialize(x,nu); //start with x=0 guess, try using explicit step for start?
        //========== Solve linear system on each element (turns chg in conservatives to change in solution variables)
        auto D = (double**)malloc((NVAR) * sizeof(double*));
        for (int k = 0; k < NVAR; k++)
            D[k] = (double*)malloc( (NVAR) * sizeof(double));
        for (int i=0; i<nx-1; i++) {
            for (int j=0; j<ny-1; j++) {
                double *unkij = &(unkel[IJK(i, j, 0, nx-1, NVAR)]);
                double LUtol = 1e-16;
                int iel = IJ(i,j,nx-1);
                int N = NVAR;
                int P[NVAR+1]{}; //permutation vector for pivoting

                //Evaluate the jacobian / Implicit matrix
                RegularizationTerm(dt, unkij, ElemVar[iel], D);

                //get the rhs block needed
                double *b = &(RHS[IJK(i, j, 0, nx - 1, NVAR)]);
                double *xLU = &(dv[ IJK(i, j, 0, nx - 1, NVAR)]);
                LUPDecompose(D, N, LUtol, P);
                LUPSolve(D, P, b, N, xLU);
            }
        }
        (void) veccopy(x, dv, nu);
        CGsolve(nx,ny,dt,air,ElemVar,bound,uFS,ibound,geoel,geofa,unkel,RHS,x);
        for (int k=0; k<NVAR; k++){
            free(D[k]);
        }
        free(D);

        //Define sufficient descent
        // for now just accept some descent
        normmax = rnorm;
        // accept norm increase for first number of iterations > globilazation
        if (iter <= ITGLOBAL) normmax = eta*rnorm;

        //First state guess
        vecinitialize(dv,nu);
        alpha = alp0;
        (void) vecscale(x, alpha, nu, dv);
        (void) vecadd(unkel, dv, nu, unknew);
        for (int i=0; i<nx-1; i++){
            for (int j=0; j<ny-1; j++) {
                int ie = IJ(i,j,nx-1);
                int iunk = IJK(i,j,0,nx-1,NVAR);
                ELVanew[ie].Initialize(&(unknew[iunk]));
                ELVanew[ie].UpdateState(air);
            }
        }
        bound.set_boundary_conditions(nx, ny, air, ELVanew, uFS, ibound, geofa, unknew);
        calc_dudt(nx, ny, air, ELVanew, uFS, ibound, geoel, geofa, unknew, bound,RHSnew);
        calculate_residual(nx, ny, RHSnew, ressum);
        rnormnew = norm2(ressum, NVAR);

        //Find a step size so that we have sufficient descent
        for (iarmijo=0; iarmijo<mxitarmijo; iarmijo++){
            if (rnormnew < normmax) break;
            alpha *= 0.5;

            //First state guess
            (void) vecscale(x, alpha, nu, dv);
            (void) vecadd(unkel, dv, nu, unknew);
            for (int i=0; i<nx-1; i++){
                for (int j=0; j<ny-1; j++) {
                    int ie = IJ(i,j,nx-1);
                    ELVanew[ie].UpdateState(air);
                }
            }
            bound.set_boundary_conditions(nx, ny, air, ELVanew, uFS, ibound, geofa, unknew);
            calc_dudt(nx, ny, air, ELVanew, uFS, ibound, geoel, geofa, unknew, bound,RHSnew);
            calculate_residual(nx, ny, RHSnew, ressum);
            rnormnew = norm2(ressum, NVAR);
        }
        //printf("Armijo Line Search Competed, iarm=%d\n", iarmijo);

        if (iter%1 == 0) {
            printf("\nIter:%7d\tdt:%7.4e\tiarm=%5d\tRelativeTotalResisual:  %8.5e\n\n", \
                    iter, dt,iarmijo, rnormnew/norm0);
        }
        if (iter > 0 and iter%1 == 0){
            //printf("Saving current Solution\n\n");
            print_state("Final State", nx, ny, air, xgeo, ygeo, unknew, geoel);
        }
        if (rnormnew/norm0 < RESTOL) {
            printf("==================== Solution Found ====================\n");
            printf("Saving Solution File..... \n");
            print_state("Final State", nx, ny, air, xgeo, ygeo, unknew, geoel);
            free(ELVanew);
            free(unknew);
            free(RHSnew);
            free(RHS);
            free(x);
            free(dv);
            return 1;
        }

        (void) veccopy(unkel,unknew, nu);
        (void) veccopy(RHS, RHSnew, nu);
        for (int i=0; i<nx-1; i++){
            for (int j=0; j<ny-1; j++) {
                int ie = IJ(i,j,nx-1);
                ElemVar[ie].UpdateState(air);
            }
        }
        rnorm = rnormnew;

    }
    free(ELVanew);
    free(unknew);
    free(RHSnew);
    free(RHS);
    free(x);
    free(dv);
    return 0;
}