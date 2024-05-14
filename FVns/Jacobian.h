//
// Created by tskoepli on 4/9/2024.
//

#ifndef FVNS_JACOBIAN_H
#define FVNS_JACOBIAN_H

#include "StateVariables.h"
#include "Indexing.h"
#include "MathTools.h"
#include "BoundaryConditions.h"
#include "SpatialDiscretization.h"

void RegularizationTerm(double  dt,const double* unk,State& var,double** D);
void RegularizationTerm(double CFL, double* dx, const double* unk, State& var,double** D);

void BuildJacobian(int nx, int ny, double CFL, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
                   double* geofa, double* unkel, const double* RHS, BC& bound, double** Jout);

void JacobianVectorMultiply(int nx, int ny, double CFL, Thermo& air, State* ElemVar, double *uFS, int* ibound, double* geoel,
                            double* geofa, double* unkel,const double* RHS, BC& bound, const double* qin, double* qout);
#endif //FVNS_JACOBIAN_H
