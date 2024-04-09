//
// Created by tskoepli on 1/27/2024.
//

#include <cmath>
#include <cstdio>
#include "BoundaryConditions.h"
#include "EulerFlux.h"



void SubsonInflo(Thermo& air, double vxint, double vyint, double cint, const double *unkel0, State var0, double nx, double ny,
                 double *rhogst, double *vxgst, double *vygst, double *Tgst) {
    //Copied from CFD2 project, works but should be cleaned up
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (unkel0[1] * nx) + (unkel0[2] * ny);
    double vinftyTANx = unkel0[1] - vinftyDOTn * nx;
    double vinftyTANy = unkel0[2] - vinftyDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (air.gam - 1)) * (cint - var0.a);
    //c of ghost cell
    double cgst = var0.a + 0.5 * (air.gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = unkel0[0] * pow(cgst * cgst / (var0.a * var0.a), (air.gam - 1));
    //p  of ghost cell
    double pgst = var0.p * pow(rhogst[0] / unkel0[0], air.gam);
    //finish by finding the conserved variables
    vxgst[0] =  (vinftyTANx + vgstDOTn * nx);
    vygst[0] =  (vinftyTANy + vgstDOTn * ny);
    Tgst[0] = pgst / (rhogst[0]*air.Rs[0]);   ///hardcoded but using this as motivation for more indepth overhaul
}
void SubsonOutfl(Thermo& air, double rhoint, double pint, double vxint, double vyint, double cint, const double *unkel0, State var0,
                 double nx, double ny, double *rhogst, double *vxgst, double *vygst, double *Tgst) {
    //Copied from CFD2 project, works but should be cleaned up
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (unkel0[1] * nx) + (unkel0[2] * ny);
    double vintTANx = vxint - vintDOTn * nx;
    double vintTANy = vyint - vintDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (air.gam - 1)) * (cint - var0.a);
    //c of ghost cell
    double cgst = var0.a + 0.5 * (air.gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = rhoint * pow(cgst * cgst / (cint * cint), (air.gam - 1));
    //p  of ghost cell
    double pgst = pint * pow(rhogst[0] / rhoint, air.gam);
    //finish by finding the rest of the conserved variables
    vxgst[0] = (vintTANx + vgstDOTn * nx);
    vygst[0] = (vintTANy + vgstDOTn * ny);
    Tgst[0] = pgst / (rhogst[0]*air.Rs[0]);
}

void boundary_state(int btype, Thermo& air,double normx, double normy, const double *uFS,
                    const double* uLeft, State varL, double* uRight) {
    //==========Apply Boundary Condition
    double rhoL, uL, vL, vDOTn, uBP[NVAR];

    //==========Find Normal Velocity
    rhoL = uLeft[0];
    uL = uLeft[1];
    vL = uRight[2];
    vDOTn = uL*normx + vL*normy;

    //Wall BC
    if (btype == 0 or btype == 4) {
        //Density and temperature(adiabatic) are constant
        uRight[0] = uLeft[0];
        uRight[3] = uLeft[3];

        if (btype==4) {
            //''ghost'' velocity is mirrored (slip)
            // symmetry BC
            double uR = uL - 2 * vDOTn * normx;
            double vR = vL - 2 * vDOTn * normy;
            uRight[1] = uLeft[0] * uR;
            uRight[2] = uLeft[0] * vR;
            return;
        }

        //''ghost'' velocity is opposite (no slip)
        double uR = -uL;
        double vR = -vL;
        uRight[1] = uLeft[0]*uR;
        uRight[2] = uLeft[0]*vR;

        if (_isnan(normx) or _isnan(normy)){
            printf("Undef. Surface Normal!\n");
        }
    }

    //Freestream, Back Pressure, and Outflow BC
    if (btype == 1 or btype == 2 or btype == 3) {
        double uBound[4];

        uBound[0] = uFS[0];
        uBound[1] = uFS[1];
        uBound[2] = uFS[2];
        uBound[3] = uFS[3];

        if (btype==2) {
            uBound[0] = uBP[0];
            uBound[1] = uBP[1];
            uBound[2] = uBP[2];
            uBound[3] = uBP[3];
        } else if (btype==3){
            uRight[0] = uLeft[0];
            uRight[1] = uLeft[1];
            uRight[2] = uLeft[2];
            uRight[3] = uLeft[3];
            return;
        }


        //get all interior primitives

        //Normal Mach number
        double MDOTn = vDOTn / varL.a;

        if (MDOTn <= -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Inflow - fully determined by free stream
            uRight[0] = uBound[0];
            uRight[1] = uBound[1];
            uRight[2] = uBound[2];
            uRight[3] = uBound[3];
            return;
        }
        if (MDOTn <= 0 && MDOTn > -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Inflow
            double rhoR, uR, vR, TR;
            SubsonInflo(air, uL, vL, varL.a, uBound, varL, normx, normy, &rhoR, &uR, &vR, &TR);
            uRight[0] = rhoR;
            uRight[1] = uR;
            uRight[2] = vR;
            uRight[3] = TR;
            return;
        }
        if (MDOTn >= 1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Outflow - fully determined by interior
            uRight[0] = uLeft[0];
            uRight[1] = uLeft[1];
            uRight[2] = uLeft[2];
            uRight[3] = uLeft[3];
            return;
        }
        if (MDOTn > 0 && MDOTn < 1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Outflow
            double rhoR, uR, vR, TR;
            SubsonOutfl(air, rhoL, varL.p, uL, vL, varL.a, uFS, varL, normx,
                        normy, &rhoR, &uR, &vR, &TR);
            uRight[0] = rhoR;
            uRight[1] = uR;
            uRight[2] = vR;
            uRight[3] = TR;
            return;
        }
    }
}
