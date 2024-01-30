/* 2D Euler Equations solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#include <cmath>
#include "FileIO.h"
#include "Indexing.h"
#include "SpatialDiscretization.h"


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

double find_dt(double gam, int nx, int ny, double CFL, double* uRef, double* geofa){
    double rho, u, v, rhoe, v2, p, c, vmax, dt, mindx;
    rho = uRef[0];
    u = uRef[1]/rho;
    v = uRef[2]/rho;
    rhoe = uRef[3];

    //Find speed of sound
    v2 = (u*u + v*v);
    p = (gam - 1) * (rhoe - (0.5 * rho * v2));
    c = sqrt(gam * p / rho);

    //Max info speed = v+c
    vmax = sqrt(v2) + c;

    //Min grid spacing
    mindx = 1.0;
    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            mindx = fmin(mindx, geofa[IJK(i,j,0,nx,6)]);
            mindx = fmin(mindx, geofa[IJK(i,j,3,nx,6)]);
        }
    }

    return CFL * mindx / vmax;
}

void calculate_residual(int nx, int ny, double* unk, double* unknew, double* res){
    double res2[NVAR] = {0.0};

    for (int i=0; i<nx-1; i++){
        for (int j=0; j<ny-1; j++){
            int iu = IJK(i,j,0,nx-1,NVAR);

            for (int k=0; k<NVAR; k++){
                double diff = unknew[iu+k] - unk[iu+k];
                res2[k] += diff*diff;

                if (_isnan(diff)){
                    printf("oeups (rescalc nan)\n");
                }
            }
        }
    }

    res[0] = sqrt(res2[0]);
    res[1] = sqrt(res2[1]);
    res[2] = sqrt(res2[2]);
    res[3] = sqrt(res2[3]);
}

void vec_copy(double n, double* a, double* b){
    for (int i=0; i<n; i++){
        a[i] = b[i];
    }
}

int main() {
    //Read in setup file
    double gam, mach, tol, CFL;
    int mxiter;
    gam =1.4;
    mach = 3.0;
    tol = 1e-6;
    mxiter = 1e4; //maximum number of iteration before stopping
    CFL = 0.75;

    printf("==================== Loading Mesh ====================\n");
    //==================== Load Mesh ====================
    int nx, ny, npoin, nb, nelem, nface;
    double *x, *y;
    int* ibound;
    double* geoel;
    double* geofa;

    //read in mesh file
    printf("Reading Mesh File..... \n");
    read_mesh(&nx, &ny, &ibound, &x, &y);
    npoin = nx*ny;
    nb = 2*nx + 2*ny;
    nelem = (nx-1)*(ny-1);
    nface = 2*nelem + nx + ny;

    //Find elem volume and centroid, face len and normals
    printf("Calculating Grid Metrics..... \n");
    calc_geoel_geofa(nx, ny, x, y, &geoel, &geofa);

    printf("==================== Initializing ====================\n");
    //==================== Setup for Sim ====================
    double res0 = 1.0;
    double res[4], ressum;
    int iter;
    auto* unk    = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* unknew = (double*)malloc(NVAR*nelem*sizeof(double));
    auto* dudt   = (double*)malloc(NVAR*nelem*sizeof(double));

    //initialize solution on mesh (zero aoa)
    double uFS[4], uBP[4];
    uFS[0] = 0.246319280397921945993707147708312;
    uFS[1] = 1.0;
    uFS[2] = 0.0;
    uFS[3] = 0.5 + 1 / (gam*(gam-1)*mach*mach);

    //Plenum Pressure (mach 3 for CD nozzle, A/A* = 1.68749999)
    uBP[0] = 1.0;
    uBP[1] = 1.0;
    uBP[2] = 0.0;
    uBP[3] = 13.387171568521152;

    for (int ielem=0; ielem<nelem; ielem++){
        unk[NVAR*ielem]   = uBP[0];
        unk[NVAR*ielem+1] = uBP[1];
        unk[NVAR*ielem+2] = uBP[2];
        unk[NVAR*ielem+3] = uBP[3];
    }

    printf("===== Generating Mesh and Initial State Tecplot Files ====\n");
    print_elem_stats("MeshVolumeStats", nx, ny, geoel);
    print_state("Initial State", nx, ny, gam, x, y, unk, geoel);

    printf("Calculating Timestep..... \n");
    //Find timestep based off of CFL limit for initial condition (dt = CFL dx / c )
    double dt = find_dt(gam, nx, ny, CFL, uFS, geofa);
    printf("dt = %f\n", dt);

    printf("==================== Starting Solver ====================\n");

    for (iter=0; iter<mxiter; iter++){
        //calculate dudt
        calc_dudt(nx, ny, gam, uFS, uBP, ibound, geoel, geofa, unk, dudt);

        for (int ielem=0; ielem<nelem; ielem++){
            int iu = NVAR*ielem;
            unknew[iu  ] = unk[iu  ] + dudt[iu  ]*dt;
            unknew[iu+1] = unk[iu+1] + dudt[iu+1]*dt;
            unknew[iu+2] = unk[iu+2] + dudt[iu+2]*dt;
            unknew[iu+3] = unk[iu+3] + dudt[iu+3]*dt;
        }


        calculate_residual(nx, ny, unk, unknew, res);
        ressum = res[0]+res[1]+res[2]+res[3];
        if (iter==0) res0 = ressum;

        if (iter%50 == 0) {
            printf("Iter:%10d Rel Tot Res:  %8.5e\t\t Equation Abs Res:\t%10.2e%10.2e%10.2e%10.2e\n", \
                    iter, ressum / res0, res[0], res[1], res[2], res[3]);
        }

        if (ressum/res0 < tol) break;

        vec_copy(nelem*NVAR, unk, unknew);

    }
    printf("==================== Solution Found ====================\n");
    printf("Saving Solution File..... \n");

    print_state("Final State", nx, ny, gam, x, y, unk, geoel);


    printf("Complete.");
}
