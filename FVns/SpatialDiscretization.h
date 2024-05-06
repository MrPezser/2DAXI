//
// Created by Tsail on 1/28/2024.
//

#ifndef FVEULER_SPATIALDISCRETIZATION_H
#define FVEULER_SPATIALDISCRETIZATION_H

#include "StateVariables.h"
#include "Indexing.h"
#include "EulerFlux.h"
#include "BoundaryConditions.h"

void calc_dudt(int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* unk, BC& bound, double* dudt);
void calc_dudt_element(int iel, int jel, int nx, int ny, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
               double* geofa, double* unk, BC& bound, double* dudt);
double find_dt(Thermo& air, int nx, int ny, double CFL, const double* uRef, State& var, double* geofa);
void calculate_residual(int nx, int ny, double* res, double* ressum);


#endif //FVEULER_SPATIALDISCRETIZATION_H
