//
// Created by tskoepli on 4/9/2024.
//

#ifndef FVNS_JACOBIAN_H
#define FVNS_JACOBIAN_H

#include "StateVariables.h"
void BuildJacobian(double  dt,const double* unk,State& var,double** D);

#endif //FVNS_JACOBIAN_H
