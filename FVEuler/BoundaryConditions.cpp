//
// Created by tskoepli on 1/27/2024.
//

#include <cmath>
#include <cstdio>
#include "BoundaryConditions.h"
#include "EulerFlux.h"



void SubsonInflo(double gam, double vxint, double vyint, double cint, double *unkel0, double nx, double ny,
                 double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //Copied from CFD2 project, works but should be cleaned up
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    double rhoinfty, vxinfty, vyinfty, pinfty, cinfty, Minfty;
    getPrimatives(gam, unkel0, &rhoinfty, &vxinfty, &vyinfty, &pinfty, &cinfty, &Minfty);
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (vxinfty * nx) + (vyinfty * ny);
    double vinftyTANx = vxinfty - vinftyDOTn * nx;
    double vinftyTANy = vyinfty - vinftyDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - cinfty);
    //c of ghost cell
    double cgst = cinfty + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = rhoinfty * pow(cgst * cgst / (cinfty * cinfty), (gam - 1));
    //p  of ghost cell
    double pgst = pinfty * pow(rhogst[0] / rhoinfty, gam);
    //finish by finding the conserved variables
    rhovxgst[0] = rhoinfty * (vinftyTANx + vgstDOTn * nx);
    rhovygst[0] = rhoinfty * (vinftyTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}
void SubsonOutfl(double gam, double rhoint, double pint, double vxint, double vyint, double cint, double *unkel0,
                 double nx, double ny, double *rhogst, double *rhovxgst, double *rhovygst, double *rhoegst) {
    //INPUT = int and fs values
    //OUTPUT = ghost values
    //Get the free stream in primative variables
    double rhoinfty, vxinfty, vyinfty, pinfty, cinfty, Minfty;
    getPrimatives(gam, unkel0, &rhoinfty, &vxinfty, &vyinfty, &pinfty, &cinfty, &Minfty);
    //Calculate velocity components
    double vintDOTn = (vxint * nx) + (vyint * ny);
    double vinftyDOTn = (vxinfty * nx) + (vyinfty * ny);
    double vintTANx = vxint - vintDOTn * nx;
    double vintTANy = vyint - vintDOTn * ny;
    //Vn at ghost cell
    double vgstDOTn = 0.5 * (vintDOTn + vinftyDOTn) + (1 / (gam - 1)) * (cint - cinfty);
    //c of ghost cell
    double cgst = cinfty + 0.5 * (gam - 1) * (vgstDOTn - vinftyDOTn);
    //rho of ghost cell
    rhogst[0] = rhoint * pow(cgst * cgst / (cint * cint), (gam - 1));
    //p  of ghost cell
    double pgst = pint * pow(rhogst[0] / rhoint, gam);
    //finish by finding the rest of the conserved variables
    rhovxgst[0] = rhoinfty * (vintTANx + vgstDOTn * nx);
    rhovygst[0] = rhoinfty * (vintTANy + vgstDOTn * ny);
    rhoegst[0] = (pgst / (gam - 1)) + 0.5 * (rhovxgst[0] * rhovxgst[0] + rhovygst[0] * rhovygst[0]) / rhogst[0];
}

void boundary_state(int btype, double gam,double normx, double normy, double *uFS, double* uBP, double* uLeft, double* uRight) {
    //==========Apply Boundary Condition
    double rhoL, uL, vL, vDOTn;

    //==========Find Normal Velocity
    rhoL = uLeft[0];
    uL = uLeft[1]/rhoL;
    vL = uLeft[2]/rhoL;
    vDOTn = uL*normx + vL*normy;

    //Wall BC
    if (btype == 0) {
        //Density and energy are constant
        uRight[0] = uLeft[0];
        uRight[3] = uLeft[3];

        //''ghost'' velocity is mirrored
        double uR = uL - 2*vDOTn*normx;
        double vR = vL - 2*vDOTn*normy;
        uRight[1] = uLeft[0]*uR;
        uRight[2] = uLeft[0]*vR;

        if (_isnan(normx) or _isnan(normy)){
            printf("Undef. Surface Normal!\n");
        }
    }

    //Freestream BC
    if (btype == 1 or btype == 2 ){
        double uBound[4];

        uBound[0] = uFS[0];
        uBound[1] = uFS[1];
        uBound[2] = uFS[2];
        uBound[3] = uFS[3];

        if (btype==2){
            uBound[0] = uBP[0];
            uBound[1] = uBP[1];
            uBound[2] = uBP[2];
            uBound[3] = uBP[3];
        }


        //get all interior primitives
        double pL, cL, ML;
        getPrimatives(gam, uLeft, &rhoL, &uL, &vL, &pL, &cL, &ML);

        //Normal Mach number
        double MDOTn = vDOTn / cL;

        if (MDOTn <= -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Inflow - fully determined by free stream
            uRight[0] = uBound[0];
            uRight[1] = uBound[1];
            uRight[2] = uBound[2];
            uRight[3] = uBound[3];
        }
        if (MDOTn <= 0 && MDOTn > -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Inflow
            double rhoR, rhouR, rhovR, rhoeR;
            SubsonInflo(gam, uL, vL, cL, uBound, normx, normy, &rhoR, &rhouR, &rhovR, &rhoeR);
            uRight[0] = rhoR;
            uRight[1] = rhouR;
            uRight[2] = rhovR;
            uRight[3] = rhoeR;
        }
        if (MDOTn >= 1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Outflow - fully determined by interior
            uRight[0] = uLeft[0];
            uRight[1] = uLeft[1];
            uRight[2] = uLeft[2];
            uRight[3] = uLeft[3];
        }
        if (MDOTn > 0 && MDOTn < 1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Outflow
            double rhoR, rhouR, rhovR, rhoeR;
            SubsonOutfl(gam, rhoL, pL, uL, vL, cL, uFS, normx, normy, &rhoR, &rhouR, &rhovR,\
                            &rhoeR);
            uRight[0] = rhoR;
            uRight[1] = rhouR;
            uRight[2] = rhovR;
            uRight[3] = rhoeR;
        }
    }
}
