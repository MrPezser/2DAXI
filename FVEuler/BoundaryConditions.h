//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_BOUNDARYCONDITIONS_H
#define FVEULER_BOUNDARYCONDITIONS_H

void boundary_state(int nx, int btype, double gam, double vDOTn, double normx, double normy, double *uFS, double* geofa, double* uLeft, double* uRight);

#endif //FVEULER_BOUNDARYCONDITIONS_H
