//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_BOUNDARYCONDITIONS_H
#define FVEULER_BOUNDARYCONDITIONS_H

#include "StateVariables.h"

void boundary_state(int btype, Thermo& air,double normx, double normy, const double *uFS,
                    const double* uLeft, State varL, double* uRight);

#endif //FVEULER_BOUNDARYCONDITIONS_H
