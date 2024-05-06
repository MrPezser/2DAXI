//
// Created by Tsail on 5/5/2024.
//

#ifndef FVNS_INEXACTNEWTONCG_H
#define FVNS_INEXACTNEWTONCG_H


#include "MathTools.h"
#include "Indexing.h"
#include "Thermo.h"
#include "StateVariables.h"
#include "BoundaryConditions.h"
#include "Jacobian.h"
#include "FileIO.h"

int INCG(double* xgeo, double* ygeo, int nx, int ny,double CFL, Thermo& air, State* ElemVar, BC& bound, double *uFS, int* ibound, double* geoel,
         double* geofa, double* unkel);

#endif //FVNS_INEXACTNEWTONCG_H
