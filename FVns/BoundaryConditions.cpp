//
// Created by tskoepli on 1/27/2024.
//

#include <cmath>
#include <cstdio>
#include "BoundaryConditions.h"
#include "EulerFlux.h"



void SubsonInflo(double gam, double vxint, double vyint, double cint, const double *unkel0, State var0, double nx, double ny,
                 double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //Copied from CFD2 project, works but should be cleaned up
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (var0.vx * nx) + (var0.vy * ny);
    double vinftyTANx = var0.vx - vinftyDOTn * nx;
    double vinftyTANy = var0.vy - vinftyDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - var0.a);
    //c of ghost cell
    double cgst = var0.a + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = unkel0[0] * pow(cgst * cgst / (var0.a * var0.a), (gam - 1));
    //p  of ghost cell
    double pgst = var0.p * pow(rhogst[0] / unkel0[0], gam);
    //finish by finding the conserved variables
    rhovxgst[0] =  unkel0[0] * (vinftyTANx + vgstDOTn * nx);
    rhovygst[0] =  unkel0[0] * (vinftyTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}
void SubsonOutfl(double gam, double rhoint, double pint, double vxint, double vyint, double cint, const double *unkel0, State var0,
                 double nx, double ny, double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //Copied from CFD2 project, works but should be cleaned up
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (var0.vx * nx) + (var0.vy * ny);
    double vintTANx = vxint - vintDOTn * nx;
    double vintTANy = vyint - vintDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - var0.a);
    //c of ghost cell
    double cgst = var0.a + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = rhoint * pow(cgst * cgst / (cint * cint), (gam - 1));
    //p  of ghost cell
    double pgst = pint * pow(rhogst[0] / rhoint, gam);
    //finish by finding the rest of the conserved variables
    rhovxgst[0] = unkel0[0] * (vintTANx + vgstDOTn * nx);
    rhovygst[0] = unkel0[0] * (vintTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}

void boundary_state(int btype, double gam,double normx, double normy, const double *uFS, const double* uBP,
                    const double* uLeft, State varL, double* uRight) {
    //==========Apply Boundary Condition
    double rhoL, uL, vL, vDOTn;

    //==========Find Normal Velocity
    rhoL = uLeft[0];
    uL = varL.vx;
    vL = varL.vy;
    vDOTn = uL*normx + vL*normy;

    //Wall BC
    if (btype == 0 or btype == 4) {
        //Density and energy are constant
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
            double rhoR, rhouR, rhovR, rhoeR;
            SubsonInflo(gam, uL, vL, varL.a, uBound, varL, normx, normy, &rhoR, &rhouR, &rhovR, &rhoeR);
            uRight[0] = rhoR;
            uRight[1] = rhouR;
            uRight[2] = rhovR;
            uRight[3] = rhoeR;
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
            double rhoR, rhouR, rhovR, rhoeR;
            SubsonOutfl(gam, rhoL, varL.p, uL, vL, varL.a, uFS, varL, normx, normy, &rhoR, &rhouR, &rhovR,\
                            &rhoeR);
            uRight[0] = rhoR;
            uRight[1] = rhouR;
            uRight[2] = rhovR;
            uRight[3] = rhoeR;
            return;
        }
    }
}
