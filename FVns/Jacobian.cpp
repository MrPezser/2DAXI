//
// Created by tskoepli on 4/9/2024.
//
#include "Jacobian.h"
#include "Indexing.h"

void BuildJacobian(double  dt,const double* unk,State& var,double** D) {
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
        D[i][i] = dti * 1.0;
    }

    // momentum derivatives
    for (int j=0; j<NSP; j++){
        D[NSP][j]   = dti *  unk[NSP];
        D[NSP+1][j] = dti *  unk[NSP+1];
    }
    D[NSP][NSP]     = dti * unk[0];
    D[NSP+1][NSP+1] = dti * unk[0];

    //total energy derivatives
    for (int isp=0; isp<NSP; isp++){
        D[NSP+2][isp] = dti * var.h0;               //- (var.rhoCv/var.rhoR)*(air.Ruv/air.Mw[isp])*unk[NSP+1]);
    }
    D[NSP+2][NSP]   = dti * unk[0]* unk[NSP];
    D[NSP+2][NSP+1] = dti * unk[0]* unk[NSP+1];
    D[NSP+2][NSP+2] = dti * unk[0]* var.Cv;

    for (int i=0; i<NSP+3; i++){
        for (int j=0; j<NSP+3; j++){
            ASSERT(!_isnan(D[i][j]), "NaN Jacobian")
        }
    }

}