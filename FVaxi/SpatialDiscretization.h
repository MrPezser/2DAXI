//
// Created by Tsail on 1/28/2024.
//

#ifndef FVEULER_SPATIALDISCRETIZATION_H
#define FVEULER_SPATIALDISCRETIZATION_H

#include "StateVariables.h"

void calc_dudt(int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* yfa, double* xfa, double* unk, double* ux, double* uy, double* dudt, double* duxdt, double* duydt);

#endif //FVEULER_SPATIALDISCRETIZATION_H
