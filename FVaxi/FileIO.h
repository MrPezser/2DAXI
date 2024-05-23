//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_FILEIO_H
#define FVEULER_FILEIO_H

#include "Thermo.h"

void read_mesh(int* nx, int* ny, int** ibound, double** x, double** y);
void print_elem_stats(const char *title, int nx, int ny, const double* geoel);
void print_state(const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* geoel );

#endif //FVEULER_FILEIO_H
