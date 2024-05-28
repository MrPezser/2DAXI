//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_BOUNDARYCONDITIONS_H
#define FVEULER_BOUNDARYCONDITIONS_H



#include "StateVariables.h"
#include "Indexing.h"
#include "EulerFlux.h"

void boundary_state(int btype, Thermo& air,double normx, double normy, const double *uFS,
                    const double* uLeft, State varL, double* uRight);

struct BC {

public:
    double *uGBot, *uGRight, *uGTop, *uGLeft;
    State *BotVar, *TopVar, *RightVar, *LeftVar;
    BC(int nx, int ny) {
        uGBot   = (double*) malloc(NVAR * (nx - 1) * sizeof(double));
        uGRight = (double*) malloc(NVAR * (ny - 1) * sizeof(double));
        uGTop   = (double*) malloc(NVAR * (nx - 1) * sizeof(double));
        uGLeft  = (double*) malloc(NVAR * (ny - 1) * sizeof(double));

        BotVar    = (State*) malloc(NVAR * (nx - 1) * sizeof(State));
        TopVar    = (State*) malloc(NVAR * (nx - 1) * sizeof(State));
        RightVar  = (State*) malloc(NVAR * (ny - 1) * sizeof(State));
        LeftVar   = (State*) malloc(NVAR * (ny - 1) * sizeof(State));
    }
    ~BC(){
        free(uGBot);
        free(uGRight);
        free(uGTop);
        free(uGLeft);
        free(BotVar);
        free(TopVar);
        free(RightVar);
        free(LeftVar);
    }

    void set_boundary_conditions(int nx, int ny, Thermo& air, State* ElemVar, double *uFS,const int* ibound, double* geofa,
                                 double* unk);
};


#endif //FVEULER_BOUNDARYCONDITIONS_H
