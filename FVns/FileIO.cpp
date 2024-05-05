//
// Created by tskoepli on 1/27/2024.
//
#include <iostream>
#include <valarray>
#include "FileIO.h"
#include "Indexing.h"
#include "StateVariables.h"

void read_mesh(int* nx, int* ny, int** ibound, double** x, double** y){
    int  npoin, nb = {0};
    //==================== Read in mesh file ====================
    FILE* fmsh = fopen("../../MeshGen/Outputs/mesh.dat","r");//fopen("../Case/mesh.dat","r");
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
        for (int i=0; i< (nx-1); i++) {
            double x = geoel[IJK(i,j,1,nx-1,3)];
            double y = geoel[IJK(i,j,2,nx-1,3)];
            double vol = geoel[IJK(i,j,0,nx-1,3)];
            fprintf(fout, "%lf,\t%lf,\t%lf\n", x, y, vol);
        }
    }
    fclose(fout);
}

void print_state(const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* geoel ) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;

    FILE *fout = fopen("../Outputs/solution.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"T\", \"p\", \"c\", \"M\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx-1, ny-1);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < (ny-1); j++) {
        for (int i=0; i< (nx-1); i++) {
            double xp, yp, rho, u, v, T, M;
            xp = geoel[IJK(i,j,1,nx-1,3)];
            yp = geoel[IJK(i,j,2,nx-1,3)];

            State var;
            var.Initialize(&(unk[IJK(i,j,0,nx-1,NVAR)]));
            var.UpdateState(air);

            rho = unk[IJK(i,j,0,nx-1,NVAR)];
            u   = unk[IJK(i,j,1,nx-1,NVAR)];
            v   = unk[IJK(i,j,2,nx-1,NVAR)];
            T   = unk[IJK(i,j,3,nx-1,NVAR)];
            M = sqrt(var.v2) / var.a;

            fprintf(fout, "%lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf \n",
                    xp, yp, rho, u, v, T, var.p, var.a, M);
        }
    }
    fclose(fout);
}