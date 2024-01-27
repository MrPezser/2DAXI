/* 2D Euler Equations solver for a structured quadrilateral grid
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#include <cmath>
#include "EulerFlux.h"

#define NVAR 4
#define IJ(i, j, ni)  (((j)*(ni)) + (i))
#define IJK(i, j, k, ni, nk)  ((((j)*(ni)) + (i))*(nk) + k)

void read_mesh(int* nx, int* ny, int** ibound, double** x, double** y){
    int  npoin, nb = {0};
    //==================== Read in mesh file ====================
    FILE* fmsh = fopen("../Case/mesh.dat","r");
    fscanf(fmsh, "%d%d", nx, ny);

    npoin = (*nx) * (*ny);
    nb = 2*(*nx) + 2*(*ny);

    (*x)        = (double*)malloc(npoin*sizeof(double));
    (*y)        = (double*)malloc(npoin*sizeof(double));
    (*ibound)   = (int*)malloc(nb*sizeof(int));

    //read mesh nodes
    for(int ip=0; ip<npoin; ip++){
        fscanf(fmsh,"%lf%lf", &(*x)[ip], &(*y)[ip]);
    }
    //read BC assignments
    fscanf(fmsh,"%*s");
    for(int ib=0; ib<nb; ib++){
        fscanf(fmsh, "%d", &(*ibound)[ib]);
    }
}

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

void print_elem_stats(const char *title, int nx, int ny, const double* geoel) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;

    FILE *fout = fopen("../Outputs/volstat.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"VOL\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx - 1, ny - 1);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < (ny-1); j++) {
        for (int i=0; i< (ny-1); i++) {
            double x = geoel[IJK(i,j,1,nx-1,3)];
            double y = geoel[IJK(i,j,2,nx-1,3)];
            double vol = geoel[IJK(i,j,0,nx-1,3)];
            fprintf(fout, "%lf,\t%lf,\t%lf\n", x, y, vol);
        }
    }
    fclose(fout);
}

void print_state(const char *title, int nx, int ny, double* x, double* y, double* unk, double* geoel ) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;

    FILE *fout = fopen("../Outputs/solution.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"e\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx-1, ny-1);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < (ny-1); j++) {
        for (int i=0; i< (ny-1); i++) {
            double xp, yp, rho, u, v, e;
            xp = geoel[IJK(i,j,1,nx-1,3)];
            yp = geoel[IJK(i,j,2,nx-1,3)];
            rho = unk[IJK(i,j,0,nx-1,4)];
            u = unk[IJK(i,j,1,nx-1,4)];
            v = unk[IJK(i,j,2,nx-1,4)];
            e = unk[IJK(i,j,3,nx-1,4)];
            fprintf(fout, "%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf\n", xp, yp, rho, u, v, e);
        }
    }
    fclose(fout);
}


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

void boundary_state(int nx, int btype, double gam, double vDOTn, double normx, double normy, double *uFS, double* geofa, double* uLeft, double* uRight) {
    //==========Apply Boundary Condition
    double rhoL, uL, vL;

    rhoL = uLeft[0];
    uL = uLeft[1]/rhoL;
    vL = uLeft[2]/rhoL;

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
    }

    //Freestream BC
    if (btype == 4){
        //get all interior primitives
        double pL, cL, ML;
        getPrimatives(gam, uLeft, &rhoL, &uL, &vL, &pL, &cL, &ML);

        //Normal Mach number
        double MDOTn = vDOTn / cL;

        if (MDOTn <= -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Supersonic Inflow - fully determined by free stream
            uRight[0] = uFS[0];
            uRight[1] = uFS[1];
            uRight[2] = uFS[2];
            uRight[3] = uFS[3];
        }
        if (MDOTn <= 0 && MDOTn > -1) {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~Subsonic Inflow
            double rhoR, rhouR, rhovR, rhoeR;
            SubsonInflo(gam, uL, vL, cL, uFS, normx, normy, &rhoR, &rhouR, &rhovR, &rhoeR);
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



void calc_dudt(int nx, int ny, double gam, double *uFS, int* ibound, double* geoel, double* geofa, double* unk, double* dudt) {
    int nelem = (nx-1)*(ny-1);
    double rhsel[nelem];

    //Calculate boundary cell contributions to flux first
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
