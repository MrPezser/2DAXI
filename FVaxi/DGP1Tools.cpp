//
// Functions to aid in the addition of DG for higher order extension
// Created by Tsail on 5/29/2024.
//

#include "DGP1Tools.h"

void get_u_val(const double* unk, const double* ux, const double* uy, double xsi, double eta, double* uout){
    if (ACCUR == 1) {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k] + xsi * ux[k] + eta * uy[k];
        }
    } else {
        for (int k = 0; k < NVAR; k++) {
            uout[k] = unk[k];
        }
    }
}