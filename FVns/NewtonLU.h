//
// Created by Tsail on 5/7/2024.
//

#ifndef FVNS_NEWTONLU_H
#define FVNS_NEWTONLU_H

#include "Indexing.h"
#include "StateVariables.h"
#include "BoundaryConditions.h"
#include "SpatialDiscretization.h"
#include "MathTools.h"
#include "Jacobian.h"
#include "LUtools.h"
#include "FileIO.h"

int NLU(double* xgeo, double* ygeo, int nx, int ny,double CFL, Thermo& air, State* ElemVar, BC& bound, double *uFS, int* ibound, double* geoel,
        double* geofa, double* unkel);

#endif //FVNS_NEWTONLU_H
