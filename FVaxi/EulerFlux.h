//
// Created by tskoepli on 1/26/2024.
//

#ifndef FVEULER_EULERFLUX_H
#define FVEULER_EULERFLUX_H

#include "StateVariables.h"

void LDFSS(double normx, double normy, double len, double yface, double* uLeft, State varL, double* uRight, State varR,
           double* flux, double* parr);
//void LDFSS(double normx, double normy, double* uLeft, State varL, double* uRight, State varR, double* fout);



#endif //FVEULER_EULERFLUX_H
