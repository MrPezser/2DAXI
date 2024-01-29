//
// Created by Tsail on 1/28/2024.
//

#ifndef FVEULER_SPATIALDISCRETIZATION_H
#define FVEULER_SPATIALDISCRETIZATION_H

void calc_dudt(int nx, int ny, double gam, double *uFS, int* ibound, double* geoel, double* geofa, double* unk, double* dudt);

#endif //FVEULER_SPATIALDISCRETIZATION_H
