//
// Created by tskoepli on 2/2/2024.
//

#include <cstdlib>
#include <cmath>
#include "MeshModule.h"
#include "Indexing.h"

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
    for (int ipt=0; ipt<npoin; ipt++) (*geofa)[ipt] = NAN;
    for(int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
            int ip0 = IJ(i,j,nx);
            double len, normx, normy;

            if (i<nx-1) {
                //          Horizontal face 0,1,2
                int iph = IJ(i+1,j,nx);
                //length
                len = sqrt(((x[iph] - x[ip0]) * (x[iph] - x[ip0])) + ((y[iph] - y[ip0]) * (y[iph] - y[ip0])));
                //normal
                normx = (y[iph] - y[ip0]) / len;
                normy = -(x[iph] - x[ip0]) / len;

                (*geofa)[IJK(i,j,0,nx,6)] = len;
                (*geofa)[IJK(i,j,1,nx,6)] = normx;
                (*geofa)[IJK(i,j,2,nx,6)] = normy;
            }

            if (j<ny-1) {
                //          Vertical face 3,4,5
                int ipv = IJ(i,j+1,nx);
                //length
                len = sqrt((x[ipv] - x[ip0]) * (x[ipv] - x[ip0]) + (y[ipv] - y[ip0]) * (y[ipv] - y[ip0]));
                //normal
                normx = (y[ipv] - y[ip0]) / len;
                normy = -(x[ipv] - x[ip0]) / len;

                (*geofa)[IJK(i,j,3,nx,6)] = len;
                (*geofa)[IJK(i,j,4,nx,6)] = normx;
                (*geofa)[IJK(i,j,5,nx,6)] = normy;
            }
        }
    }

}