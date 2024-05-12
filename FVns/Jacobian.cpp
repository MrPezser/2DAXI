//
// Created by tskoepli on 4/9/2024.
//
#include "Jacobian.h"
void RegularizationTerm(double dt, const double* unk, State& var,double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};

    double dti = 1.0/dt;

    for (int i=0; i<NVAR; i++){
        for (int j=0; j<NVAR; j++){
            D[i][j] = 0.0;
        }
    }

    //d rho_s / d rho_s
    for (int i=0; i<NSP; i++){
        //Top left identity block
        D[i][i] += dti * 1.0;
    }

    // momentum derivatives   d rhou,rhov / d ()
    for (int j=0; j<NSP; j++){
        D[NSP][j]   += dti *  unk[NSP];
        D[NSP+1][j] += dti *  unk[NSP+1];
    }
    D[NSP][NSP]     += dti * unk[0];
    D[NSP+1][NSP+1] += dti * unk[0];

    //total energy derivatives
    for (int isp=0; isp<NSP; isp++){
        D[NSP+2][isp] += dti * (var.e + 0.5*var.v2);       //- (var.rhoCv/var.rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]);
    }
    D[NSP+2][NSP]   += dti * unk[0]* unk[NSP];
    D[NSP+2][NSP+1] += dti * unk[0]* unk[NSP+1];
    D[NSP+2][NSP+2] += dti * unk[0]* var.Cv;

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            ASSERT(!__isnan(D[i][j]), "NaN Jacobian")
            ASSERT(!std::isinf(D[i][j]), "INF Jacobian)")
        }
    }

}
void RegularizationTerm(double CFL, double* dx, const double* unk, State& var,double** D) {
    //bluid the matrix transformation/ jacobian for d(conserv) / d(primative)
    //double D[NSP+3][NSP+3]{0.0};

    //Max info speed = v+c
    double vmax = sqrt(var.v2) + var.a;

    //Min grid spacing
    double mindx = fmin(dx[0], dx[1]);
    double dt = CFL * mindx / vmax;

    double dti = 1.0/dt;

    for (int i=0; i<NVAR; i++){
        for (int j=0; j<NVAR; j++){
            D[i][j] = 0.0;
        }
    }

    //d rho_s / d rho_s
    for (int i=0; i<NSP; i++){
        //Top left identity block
        D[i][i] += dti * 1.0;
    }

    // momentum derivatives   d rhou,rhov / d ()
    for (int j=0; j<NSP; j++){
        D[NSP][j]   += dti *  unk[NSP];
        D[NSP+1][j] += dti *  unk[NSP+1];
    }
    D[NSP][NSP]     += dti * unk[0];
    D[NSP+1][NSP+1] += dti * unk[0];

    //total energy derivatives
    for (int isp=0; isp<NSP; isp++){
        D[NSP+2][isp] += dti * (var.e + 0.5*var.v2);       //- (var.rhoCv/var.rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]);
    }
    D[NSP+2][NSP]   += dti * unk[0]* unk[NSP];
    D[NSP+2][NSP+1] += dti * unk[0]* unk[NSP+1];
    D[NSP+2][NSP+2] += dti * unk[0]* var.Cv;

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            ASSERT(!__isnan(D[i][j]), "NaN Jacobian")
            ASSERT(!std::isinf(D[i][j]), "INF Jacobian)")
        }
    }

}


void BuildJacobian(int nx, int ny, double CFL, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
                            double* geofa, double* unkel, const double* RHS, BC& bound, double** Jout){
    /*
     * Function which exacts a Jacobian matrix, vector multiply
     * This is very memory efficient, but computationally inefficient
     *
     * qin = vector to be multiplied with
     * qout = result of matvec multiplication
     * both are size of unk
     */

    //Loop through elements
    //Start out with regularization term (1/dt)(dudv)
    //Then add contributions from neighboring elements

    DUNG("IN:: BuildJac")

    int i,j, nelem, nu, ielm, iunk, jqin, iqout;
    nelem = (nx - 1) * (ny - 1);
    nu = nelem * NVAR;
    auto uPerturb = (double*)malloc(nu*sizeof(double));
    auto RHSPertu = (double*)malloc(nu*sizeof(double));
    double* unk;
    State var;


    auto D = (double**)malloc((NVAR) * sizeof(double*));
    auto Jip = (double**)malloc((NVAR) * sizeof(double*));
    auto Jim = (double**)malloc((NVAR) * sizeof(double*));
    auto Jjp = (double**)malloc((NVAR) * sizeof(double*));
    auto Jjm = (double**)malloc((NVAR) * sizeof(double*));
    for (int k = 0; k < NVAR; k++) {
        D[k] = (double *) malloc((NVAR) * sizeof(double));
        Jip[k] = (double *) malloc((NVAR) * sizeof(double));
        Jim[k] = (double *) malloc((NVAR) * sizeof(double));
        Jjp[k] = (double *) malloc((NVAR) * sizeof(double));
        Jjm[k] = (double *) malloc((NVAR) * sizeof(double));
    }

    (void)veccopy(uPerturb,unkel,nu);

    //initialize jacobian matrix
    for (int imat=0; imat<nu; imat++){
        for (int jmat=0; jmat<nu; jmat++){
            Jout[imat][jmat] = 0.0;
        }
    }

    for(j=0; j<(ny-1); j++){
        for(i=0; i<(nx-1);i++) {
            //printf("i,j = %5d,%5d\n",i,j);


            ielm = IJ(i, j, nx - 1);
            iunk = IJK(i, j, 0, nx - 1, NVAR);

            unk = &(unkel[iunk]);
            var = ElemVar[ielm];

            //Perturb element
            //each column is dF/d(var_c)
            for (int jvar=0; jvar<NVAR; jvar++) {

                for (int imat=0; imat<NVAR; imat++){
                    for (int jmat=0; jmat<NVAR; jmat++){
                        D[imat][jmat] = 0.0;
                        Jip[imat][jmat] = 0.0;
                        Jim[imat][jmat] = 0.0;
                        Jjp[imat][jmat] = 0.0;
                        Jjm[imat][jmat] = 0.0;
                    }
                }
                //regularization term
                double dx[2];
                dx[0] = geofa[IJK(i,j,0,nx,6)];
                dx[1] = geofa[IJK(i,j,3,nx,6)];
                RegularizationTerm(CFL, dx, unk, var, D);



                //Perturb for finite difference solution
                double delvar = 1e-8 * fabs(unk[jvar]);
                if (delvar <= 1e-8) delvar = 1e-8;
                uPerturb[iunk + jvar] += delvar;

                ElemVar[ielm].Initialize(&(uPerturb[iunk]));
                ElemVar[ielm].UpdateState(air);

                //Find the new residual/RHS value
                bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, uPerturb);

                calc_dudt_element(i, j, nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, uPerturb,
                                  bound, RHSPertu);
                //calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, uPerturb,
                //          bound, RHSPertu);

                //Find the derivatives with finite difference and gather to output vector
                for (int ivar = 0; ivar < NVAR; ivar++) {
                    //same element
                    // dFi / dVj
                    D[ivar][jvar] += -(RHSPertu[iunk + ivar] - RHS[iunk + ivar]) / delvar;
                    if (__isnan(D[ivar][jvar]) or std::isinf(D[ivar][jvar])){
                        //DUNG("bread")
                    }

                    //element neighbor to left
                    if (i > 0) {
                        //int ielmN = IJ(i-1, j, nx - 1);
                        int iunkN = IJK(i - 1, j, 0, nx - 1, NVAR);
                        Jim[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor above
                    if (i < nx - 2) {
                        //int ielmN = IJ(i+1, j, nx - 1);
                        int iunkN = IJK(i + 1, j, 0, nx - 1, NVAR);
                        Jip[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor below
                    if (j > 0) {
                        //int ielmN = IJ(i, j-1, nx - 1);
                        int iunkN = IJK(i, j - 1, 0, nx - 1, NVAR);
                        Jjm[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor above
                    if (j < ny - 2) {
                        //int ielmN = IJ(i, j+1, nx - 1);
                        int iunkN = IJK(i, j + 1, 0, nx - 1, NVAR);
                        Jjp[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }

                }

                jqin = iunk + jvar;  // index for the column being multiplied === element of vector to use
                //Gather this component of the matrix multiply on RHS
                for (int ivar = 0; ivar < NVAR; ivar++) {
                    //same element
                    iqout = iunk + ivar;
                    Jout[iqout][jqin] += D[ivar][jvar];

                    //element neighbor to left
                    if (i > 0) {
                        //int ielmN = IJ(i-1, j, nx - 1);
                        int iunkN = IJK(i - 1, j, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        Jout[iqout][jqin] += Jim[ivar][jvar];
                    }
                    //element neighbor above
                    if (i < nx - 2) {
                        //int ielmN = IJ(i+1, j, nx - 1);
                        int iunkN = IJK(i + 1, j, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        Jout[iqout][jqin] += Jip[ivar][jvar];
                    }
                    //element neighbor below
                    if (j > 0) {
                        //int ielmN = IJ(i, j-1, nx - 1);
                        int iunkN = IJK(i, j - 1, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        Jout[iqout][jqin] += Jjm[ivar][jvar];
                    }
                    //element neighbor above
                    if (j < ny - 2) {
                        //int ielmN = IJ(i, j+1, nx - 1);
                        int iunkN = IJK(i, j + 1, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        Jout[iqout][jqin] += Jjp[ivar][jvar];
                    }

                }

                //reset perturbed variable
                uPerturb[iunk + jvar] = unkel[iunk + jvar];
                veccopy(RHSPertu, RHS, nu);
                ElemVar[ielm].Initialize(unk);
                ElemVar[ielm].UpdateState(air);
            }

            //ElemVar[ielm].Initialize(unk);
            //ElemVar[ielm].UpdateState(air);
        }
    }


    free(uPerturb);
    free(RHSPertu);
    for (int k = 0; k < NVAR; k++) {
        free(D[k]);
        free(Jip[k]);
        free(Jim[k]);
        free(Jjp[k]);
        free(Jjm[k]);
    }
    free(D);
    free(Jip);
    free(Jim);
    free(Jjp);
    free(Jjm);

    DUNG("OUT: BuildJac")
}


void JacobianVectorMultiply(int nx, int ny, double dt, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
                            double* geofa, double* unkel, const double* RHS, BC& bound, double* qin, double* qout){
    /*
     * Function which exacts a Jacobian matrix, vector multiply
     * This is very memory efficient, but computationally inefficient
     *
     * qin = vector to be multiplied with
     * qout = result of matvec multiplication
     * both are size of unk
     */

    //Loop through elements
    //Start out with regularization term (1/dt)(dudv)
    //Then add contributions from neighboring elements

    int i,j, nelem, nu, ielm, iunk, iqin, iqout;
    nelem = (nx - 1) * (ny - 1);
    nu = nelem * NVAR;
    auto uPerturb = (double*)malloc(nu*sizeof(double));
    auto RHSPertu = (double*)malloc(nu*sizeof(double));
    double* unk;
    State var;


    auto D = (double**)malloc((NVAR) * sizeof(double*));
    auto Jip = (double**)malloc((NVAR) * sizeof(double*));
    auto Jim = (double**)malloc((NVAR) * sizeof(double*));
    auto Jjp = (double**)malloc((NVAR) * sizeof(double*));
    auto Jjm = (double**)malloc((NVAR) * sizeof(double*));
    for (int k = 0; k < NVAR; k++) {
        D[k] = (double *) malloc((NVAR) * sizeof(double));
        Jip[k] = (double *) malloc((NVAR) * sizeof(double));
        Jim[k] = (double *) malloc((NVAR) * sizeof(double));
        Jjp[k] = (double *) malloc((NVAR) * sizeof(double));
        Jjm[k] = (double *) malloc((NVAR) * sizeof(double));
    }

    (void)veccopy(uPerturb,unkel,nu);
    (void)vecinitialize(qout,nu);  //set to zero

    for(j=0; j<(ny-1); j++){
        for(i=0; i<(nx-1);i++) {

            ielm = IJ(i, j, nx - 1);
            iunk = IJK(i, j, 0, nx - 1, NVAR);

            unk = &(unkel[iunk]);
            var = ElemVar[ielm];

            //Perturb element
            //each column is dF/d(var_c)
            for (int jvar=0; jvar<NVAR; jvar++) {

                for (int imat=0; imat<NVAR; imat++){
                    for (int jmat=0; jmat<NVAR; jmat++){
                        D[imat][jmat] = 0.0;
                        Jip[imat][jmat] = 0.0;
                        Jim[imat][jmat] = 0.0;
                        Jjp[imat][jmat] = 0.0;
                        Jjm[imat][jmat] = 0.0;
                    }
                }
                //regularization term
                RegularizationTerm(dt, unk, var, D);



                //Perturb for finite difference solution
                double delvar = 1e-8 * fabs(unk[jvar]);
                if (delvar <= 1e-8) delvar = 1e-8;
                uPerturb[iunk + jvar] += delvar;

                ElemVar[ielm].Initialize(&(uPerturb[iunk]));
                ElemVar[ielm].UpdateState(air);

                //Find the new residual/RHS value
                bound.set_boundary_conditions(nx, ny, air, ElemVar, uFS, ibound, geofa, uPerturb);

                //calc_dudt_element(i, j, nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, uPerturb,
                //                  bound, RHSPertu);
                calc_dudt(nx, ny, air, ElemVar, uFS, ibound, geoel, geofa, uPerturb,
                                  bound, RHSPertu);

                //Find the derivatives with finite difference and gather to output vector
                for (int ivar = 0; ivar < NVAR; ivar++) {
                    //same element
                    // dFi / dVj
                    D[ivar][jvar] += -(RHSPertu[iunk + ivar] - RHS[iunk + ivar]) / delvar;
                    if (__isnan(D[ivar][jvar]) or std::isinf(D[ivar][jvar])){
                        //DUNG("bread")
                    }

                    //element neighbor to left
                    if (i > 0) {
                        //int ielmN = IJ(i-1, j, nx - 1);
                        int iunkN = IJK(i - 1, j, 0, nx - 1, NVAR);
                        Jim[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor above
                    if (i < nx - 2) {
                        //int ielmN = IJ(i+1, j, nx - 1);
                        int iunkN = IJK(i + 1, j, 0, nx - 1, NVAR);
                        Jip[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor below
                    if (j > 0) {
                        //int ielmN = IJ(i, j-1, nx - 1);
                        int iunkN = IJK(i, j - 1, 0, nx - 1, NVAR);
                        Jjm[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }
                    //element neighbor above
                    if (j < ny - 2) {
                        //int ielmN = IJ(i, j+1, nx - 1);
                        int iunkN = IJK(i, j + 1, 0, nx - 1, NVAR);
                        Jjp[ivar][jvar] += -(RHSPertu[iunkN + ivar] - RHS[iunkN + ivar]) / delvar;
                    }

                }

                iqin = iunk + jvar;  // index for the column being multiplied === element of vector to use
                //Gather this component of the matrix multiply on RHS
                for (int ivar = 0; ivar < NVAR; ivar++) {
                    //same element
                    iqout = iunk + ivar;
                    qout[iqout] += D[ivar][jvar] * qin[iqin];

                    //element neighbor to left
                    if (i > 0) {
                        //int ielmN = IJ(i-1, j, nx - 1);
                        int iunkN = IJK(i - 1, j, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        qout[iqout] += Jim[ivar][jvar] * qin[iqin];
                    }
                    //element neighbor above
                    if (i < nx - 2) {
                        //int ielmN = IJ(i+1, j, nx - 1);
                        int iunkN = IJK(i + 1, j, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        qout[iqout] += Jip[ivar][jvar] * qin[iqin];
                    }
                    //element neighbor below
                    if (j > 0) {
                        //int ielmN = IJ(i, j-1, nx - 1);
                        int iunkN = IJK(i, j - 1, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        qout[iqout] += Jjm[ivar][jvar] * qin[iqin];
                    }
                    //element neighbor above
                    if (j < ny - 2) {
                        //int ielmN = IJ(i, j+1, nx - 1);
                        int iunkN = IJK(i, j + 1, 0, nx - 1, NVAR);
                        iqout = iunkN + ivar;
                        qout[iqout] += Jjp[ivar][jvar] * qin[iqin];
                    }

                }

                //reset perturbed variable
                uPerturb[iunk + jvar] = unkel[iunk + jvar];
                ElemVar[ielm].Initialize(unk);
                ElemVar[ielm].UpdateState(air);
            }

            //ElemVar[ielm].Initialize(unk);
            //ElemVar[ielm].UpdateState(air);
        }
    }


    free(uPerturb);
    free(RHSPertu);
    for (int k = 0; k < NVAR; k++) {
        free(D[k]);
        free(Jip[k]);
        free(Jim[k]);
        free(Jjp[k]);
        free(Jjm[k]);
    }
    free(D);
    free(Jip);
    free(Jim);
    free(Jjp);
    free(Jjm);
}