/* 2D Euler Equations solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#include <cmath>
#include "EulerFlux.h"
#include "FileIO.h"
#include "BoundaryConditions.h"
#include "Indexing.h"

#define NVAR 4

void calc_geoel_geofa(const int nx, const int ny, double* x, double* y, \
            double** geoel, double** geofa) {
    int nelem, npoin;
    nelem = (nx-1)*(ny-1);
    npoin = nx*ny;

    (*geoel) = (double*)malloc(3*nelem*sizeof(double));
    (*geofa) = (double*)malloc(3*2*npoin*sizeof(double));

    //Calculate mesh stats
    //element stats
    for(int ie=0;ie<(nx-1); ie++){
        for(int je=0; je<(ny-1); je++){
            int igeo = IJK(ie, je, 0, nx-1, 3);
            int ip1 = IJ(ie, je, nx);
            int ip2 = IJ(ie+1, je, nx);
            int ip3 = IJ(ie+1, je+1, nx);
            int ip4 = IJ(ie, je+1, nx);

            // find volume of element
            double vol = 0.5 *(x[ip1]*y[ip2] + x[ip2]*y[ip3] + x[ip3]*y[ip4] + x[ip4]*y[ip1] \
                - x[ip2]*y[ip1] - x[ip3]*y[ip2] - x[ip4]*y[ip3] - x[ip1]*y[ip4]);
            (*geoel)[igeo] = fabs(vol);

            //find centroid of element
            (*geoel)[igeo+1] = 0.25*(x[ip1] + x[ip2] + x[ip3] + x[ip4]);
            (*geoel)[igeo+2] = 0.25*(y[ip1] + y[ip2] + y[ip3] + y[ip4]);
        }
    }

    //face stats
    //interior points
    for(int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            int ip0 = IJ(i,j,nx);
            int iph = IJ(i+1,j,nx);
            int ipv = IJ(i,j+1,nx);

            //          Horizontal face 0,1,2
            //length
            double len  = sqrt((x[iph]-x[ip0])*(x[iph]-x[ip0]) + (y[iph]-y[ip0])*(y[iph]-y[ip0]) );
            (*geofa)[IJK(i,j,0,nx,6)] = len;
            //normal
            double normx,normy;
            normx =  (y[iph] - y[ip0])/len;
            normy = -(x[iph] - x[ip0])/len;
            (*geofa)[IJK(i,j,1,nx,6)] = normx;
            (*geofa)[IJK(i,j,2,nx,6)] = normy;

            //          Vertical face 3,4,5
            //length
            len  = sqrt((x[ipv]-x[ip0])*(x[ipv]-x[ip0]) + (y[ipv]-y[ip0])*(y[ipv]-y[ip0]) );
            (*geofa)[IJK(i,j,3,nx,6)] = len;
            //normal
            normx =  (y[ipv] - y[ip0])/len;
            normy = -(x[ipv] - x[ip0])/len;
            (*geofa)[IJK(i,j,4,nx,6)] = normx;
            (*geofa)[IJK(i,j,6,nx,6)] = normy;
        }
    }

}

void calc_dudt(int nx, int ny, double gam, double *uFS, int* ibound, double* geoel, double* geofa, double* unk, double* dudt) {
    int nelem = (nx-1)*(ny-1);
    double rhsel[nelem];

    //Calculate boundary cell state (ghost state)
    double uGBot[NVAR*(nx-1)], uGRight[NVAR*(ny-1)], uGTop[NVAR*(nx-1)], uGLeft[NVAR*(ny-1)];
    //bottom side of domain
    for (int i=0; i<(nx-1); i++){
        // left state = interior, right state = ghost
        double *uLeft, *uRight;
        int btype;
        btype = ibound[i];

        uLeft = &(unk[IJK(i,0,0, nx-1, 4)]);
        uRight = (double*)malloc(4*sizeof(double));

        //==========Find Normal Velocity
        double rhoL, uL, vL, normx, normy, vDOTn;
        rhoL = uLeft[0];
        uL = uLeft[1]/rhoL;
        vL = uLeft[2]/rhoL;
        normx = geofa[IJK(i, 0, 1,nx-1,6)];
        normy = geofa[IJK(i, 0, 2,nx-1,6)];
        vDOTn = uL*normx + vL*normy;
        boundary_state(nx,btype,gam,vDOTn,normx,normy,uFS,geofa,uLeft,uRight);


    }

}



int main() {
    //Read in setup file
    double gam, mach, tol;
    int mxiter;
    gam =1.4;
    mach = 0.5;
    tol = 1e-6;
    mxiter = 1e4; //maximum number of iteration before stopping

    //==================== Load Mesh ====================
    int nx, ny, npoin, nb, nelem, nface;
    double *x, *y;
    int* ibound;
    double* geoel;
    double* geofa;

    //read in mesh file
    read_mesh(&nx, &ny, &ibound, &x, &y);
    npoin = nx*ny;
    nb = 2*nx + 2*ny;
    nelem = (nx-1)*(ny-1);
    nface = 2*nelem + nx + ny;

    //Find elem volume and centroid, face len and normals
    calc_geoel_geofa(nx, ny, x, y, &geoel, &geofa);

    //==================== Setup for Sim ====================
    double res;
    int iter;
    auto* u = (double*)malloc(NVAR*nelem*sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4];
    uFS[0] = 1.0;
    uFS[1] = 1.0;
    uFS[2] = 0.0;
    uFS[3] = 0.5 + 1 / (gam*(gam-1)*mach*mach);
    for (int ielem=0; ielem<nelem; ielem++){
        u[NVAR*ielem]   = uFS[0];
        u[NVAR*ielem+1] = uFS[1];
        u[NVAR*ielem+2] = uFS[2];
        u[NVAR*ielem+3] = uFS[3];
    }


    res = 1.0; //Initialized residual
    iter = 0;

    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    print_state("Initial State", nx, ny, x, y, u, geoel);

}
