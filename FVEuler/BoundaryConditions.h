//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_BOUNDARYCONDITIONS_H
#define FVEULER_BOUNDARYCONDITIONS_H

void boundary_state(int btype, double gam,double normx, double normy, double *uFS, double* uBP, double* uLeft, double* uRight);

#endif //FVEULER_BOUNDARYCONDITIONS_H
