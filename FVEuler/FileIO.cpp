//
// Created by tskoepli on 1/27/2024.
//
#include <iostream>
#include <valarray>
#include "FileIO.h"
#include "Indexing.h"

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

void print_state(const char *title, int nx, int ny, double gam, double* x, double* y, double* unk, double* geoel ) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;

    FILE *fout = fopen("../Outputs/solution.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"e\", \"p\", \"c\", \"M\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx-1, ny-1);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < (ny-1); j++) {
        for (int i=0; i< (ny-1); i++) {
            double xp, yp, rho, u, v, e, p, c, M;
            xp = geoel[IJK(i,j,1,nx-1,3)];
            yp = geoel[IJK(i,j,2,nx-1,3)];

            rho = unk[IJK(i,j,0,nx-1,NVAR)];
            u   = unk[IJK(i,j,1,nx-1,NVAR)]/rho;
            v   = unk[IJK(i,j,2,nx-1,NVAR)]/rho;
            e   = unk[IJK(i,j,3,nx-1,NVAR)]/rho;
            p = (gam - 1) * (rho*e- (0.5 * rho * (u*u + v*v)));
            c = sqrt(gam * p / rho);
            M = sqrt((u*u + v*v)) / c;

            fprintf(fout, "%lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf \n", xp, yp, rho, u, v, e, p, c, M);
        }
    }
    fclose(fout);
}