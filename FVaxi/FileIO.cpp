//
// Created by tskoepli on 1/27/2024.
//
#include <iostream>
#include <valarray>
#include "FileIO.h"
#include "Indexing.h"
#include "StateVariables.h"
#include "DGP1Tools.h"

void read_mesh(int* nx, int* ny, int** ibound, double** x, double** y){
    int  npoin, nb = {0};
    //==================== Read in mesh file ====================
    FILE* fmsh = fopen("../Case/mesh.dat","r"); //fopen("../../MeshGen/Outputs/mesh.dat","r");
    fscanf(fmsh, "%d%d", nx, ny);

    npoin = (*nx) * (*ny);
    nb = 2*(*nx) + 2*(*ny);

    (*x)        = (double*)malloc(npoin*sizeof(double));
    (*y)        = (double*)malloc(npoin*sizeof(double));
    (*ibound)   = (int*)malloc(nb*sizeof(int));

    //read mesh nodes
    for(int ip=0; ip<npoin; ip++){
        fscanf(fmsh,"%lf%lf", &(*x)[ip], &(*y)[ip]);
        //(*y)[ip] += 10000.0;
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

            fprintf(fout, "%.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e \n",
                    xp, yp, rho, u, v, T, var.p, var.a, M);
        }
    }
    fclose(fout);
}
void print_state_DGP1(const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* ux, double* uy, double* geoel ) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;

    FILE *fout = fopen("../Outputs/solution.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"T\", \"p\", \"c\", \"M\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx, ny);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < ny; j++) {
        for (int i=0; i< nx; i++) {
            double xp, yp, rho, u, v, T, M;
            xp = x[IJ(i,j,nx)];
            yp = y[IJ(i,j,nx)];

            //Average out contributions from adjacent cells
            double *upp, *ump, *umm, *upm, unkout[NVAR];
            int iupp, iump, iumm, iupm;
            iupp = IJK(i,j,0,nx-1,NVAR);
            iump = IJK(i-1,j,0,nx-1,NVAR);
            iumm = IJK(i-1,j-1,0,nx-1,NVAR);
            iupm = IJK(i,j-1,0,nx-1,NVAR);
            upp = &(unk[iupp]);
            ump = &(unk[iump]);
            umm = &(unk[iumm]);
            upm = &(unk[iupm]);
            if (i==0 and j==0){
                //BL corner
                get_u_val_standardrecon(upp, &(ux[iupp]), &(uy[iupp]), -1.0, -1.0, unkout);
            } else if (i==0 and j==ny-1) {
                //TL corner
                get_u_val_standardrecon(upm, &(ux[iupm]), &(uy[iupm]), -1.0, 1.0, unkout);
            } else if (i==nx-1 and j==0) {
                //BR corner
                get_u_val_standardrecon(ump, &(ux[iump]), &(uy[iump]), 1.0, -1.0, unkout);
            } else if (i==nx-1 and j==ny-1){
                //TR corner
                get_u_val_standardrecon(umm, &(ux[iumm]), &(uy[iumm]), 1.0, 1.0, unkout);
            } else if (j==0){
                //bottom boundary
                double unk1[4], unk2[4];
                get_u_val_standardrecon(upp, &(ux[iupp]), &(uy[iupp]), -1.0, -1.0, unk1);
                get_u_val_standardrecon(ump, &(ux[iump]), &(uy[iump]),  1.0, -1.0, unk2);
                for (int kvar=0;kvar<NVAR;kvar++) unkout[kvar] = 0.5*(unk1[kvar] + unk2[kvar]);
            } else if (j==ny-1){
                //top boundary
                double unk1[4], unk2[4];
                get_u_val_standardrecon(upm, &(ux[iupm]), &(uy[iupm]), -1.0, 1.0, unk1);
                get_u_val_standardrecon(umm, &(ux[iumm]), &(uy[iumm]),  1.0, 1.0, unk2);
                for (int kvar=0;kvar<NVAR;kvar++) unkout[kvar] = 0.5*(unk1[kvar] + unk2[kvar]);
            } else if (i==0){
                //left boundary
                double unk1[4], unk2[4];
                get_u_val_standardrecon(upp, &(ux[iupp]), &(uy[iupp]), -1.0, -1.0, unk1);
                get_u_val_standardrecon(upm, &(ux[iupm]), &(uy[iupm]), -1.0,  1.0, unk2);
                for (int kvar=0;kvar<NVAR;kvar++) unkout[kvar] = 0.5*(unk1[kvar] + unk2[kvar]);
            } else if (i==nx-1){
                //right boundary
                double unk1[4], unk2[4];
                get_u_val_standardrecon(ump, &(ux[iump]), &(uy[iump]), 1.0, -1.0, unk1);
                get_u_val_standardrecon(umm, &(ux[iumm]), &(uy[iumm]), 1.0,  1.0, unk2);
                for (int kvar=0;kvar<NVAR;kvar++) unkout[kvar] = 0.5*(unk1[kvar] + unk2[kvar]);
            } else {
                //internal
                double unk1[4], unk2[4], unk3[4], unk4[4];
                get_u_val_standardrecon(upp, &(ux[iupp]), &(uy[iupp]), -1.0, -1.0, unk1);
                get_u_val_standardrecon(upm, &(ux[iupm]), &(uy[iupm]), -1.0,  1.0, unk2);
                get_u_val_standardrecon(ump, &(ux[iump]), &(uy[iump]),  1.0, -1.0, unk3);
                get_u_val_standardrecon(umm, &(ux[iumm]), &(uy[iumm]),  1.0,  1.0, unk4);
                for (int kvar=0;kvar<NVAR;kvar++) unkout[kvar] = 0.25*(unk1[kvar] + unk2[kvar] + unk3[kvar] + unk4[kvar]);
            }

            State var;
            var.Initialize(unkout);
            var.UpdateState(air);

            rho = unkout[0];//unk[IJK(i,j,0,nx-1,NVAR)];
            u   = unkout[1];//unk[IJK(i,j,1,nx-1,NVAR)];
            v   = unkout[2];//unk[IJK(i,j,2,nx-1,NVAR)];
            T   = unkout[3];//unk[IJK(i,j,3,nx-1,NVAR)];
            M = sqrt(var.v2) / var.a;

            fprintf(fout, "%.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e,\t %.16e \n",
                    xp, yp, rho, u, v, T, var.p, var.a, M);
        }
    }
    fclose(fout);
}

void print_state_axi(const char *title, int nx, int ny, Thermo& air, double* x, double* y, double* unk, double* geoel ) {
    //Makes a tecplot file of the grid and a setup file for the solver
    //int nb = 2*nx + 2*ny;
    int nelem = nx * ny;
    int nrad = 50;

    FILE *fout = fopen("../Outputs/solution3D.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\", \"Z\",\"rho\", \"u\", \"v\", \"w\", \"T\", \"p\", \"c\", \"M\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, K=%d DATAPACKING=POINT\n", nx-1, ny-1, nrad);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j < (ny-1); j++) {
        for (int i=0; i< (nx-1); i++) {
            for (int k=0; k<nrad; k++) {
                double xp, yp, zp, rho, u, v, w, T, M, psi;
                xp = geoel[IJK(i, j, 1, nx - 1, 3)];
                yp = geoel[IJK(i, j, 2, nx - 1, 3)];

                psi = 2.0 * M_PI * k / nrad;

                zp = yp * sin(psi);
                yp = yp * cos(psi);

                State var;
                var.Initialize(&(unk[IJK(i, j, 0, nx - 1, NVAR)]));
                var.UpdateState(air);

                rho = unk[IJK(i, j, 0, nx - 1, NVAR)];
                u = unk[IJK(i, j, 1, nx - 1, NVAR)];
                v = unk[IJK(i, j, 2, nx - 1, NVAR)];
                w = v * sin(psi);
                v = v * cos(psi);
                T = unk[IJK(i, j, 3, nx - 1, NVAR)];
                M = sqrt(var.v2) / var.a;

                fprintf(fout, "%lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf,\t %lf \n",
                        xp, yp, zp, rho, u, v, w, T, var.p, var.a, M);
            }
        }
    }
    fclose(fout);
}