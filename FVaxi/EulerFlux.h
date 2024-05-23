//
// Created by tskoepli on 1/26/2024.
//

#ifndef FVEULER_EULERFLUX_H
#define FVEULER_EULERFLUX_H

#include "StateVariables.h"

void LeerFlux(const double gam, double normx, double normy, double* uLeft, State varL, double* uRight, State varR, double* fout);
void LDFSS(double normx, double normy, double* uLeft, State varL, double* uRight, State varR, double* fout);



#endif //FVEULER_EULERFLUX_H
