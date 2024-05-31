//
// Created by Tsail on 5/29/2024.
//

#ifndef FVAXI_DGP1TOOLS_H
#define FVAXI_DGP1TOOLS_H

#include "Indexing.h"
#include "StateVariables.h"

void get_u_val(const double* unk, const double* ux, const double* uy, double xsi, double eta, double* uout);
void DGP1_volume_integral(int nx, int ny, double vol, double* xfa, double* yfa, double* geoel, double* unk, State* ElemVar,
                          double* duxdt, double* duydt);

#endif //FVAXI_DGP1TOOLS_H
