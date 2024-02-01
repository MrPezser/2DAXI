/* 2D Structured Mesh Generation for Various Ramps
 *
 *
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#define IU(i, j, ni)  (((j)*(ni)) + (i))

//test ramp geometry for now
double ramp_surface(double x, double h, double L) {
    //test geometry: straight - ramp - straight
    if(x < 1) {
        return 0;
    } else if( x > 1+L) {
        return h;
    } else {
        return (h/L)*(x-1);
    }
}

double nozzle_surface(double x, double L1, double L2) {
    //test geometry: Mach 3.0 CD nozzle
    double h =  1.5;//0.81481480 - 0.01481480;
    if(x < 0.0) {
        return 0;
    } else if( x > L1) {
        return h * ((L2-(x-L1))/L2);
    } else {
        return (h/L1)*(x);
    }
}

void printgrid(const char *title, int nx, int ny, double *x, double *y, int* ibound) {
    //Makes a tecplot file of the grid and a setup file for the solver
    int nb = 2*(nx-1) + 2*(ny-1);

    FILE* fout = fopen("../Outputs/grid.tec", "w");
    if (fout == nullptr) printf("oeups\n");

    //printf("\nDisplaying Grid Header\n");
    fprintf(fout, "TITLE = \"%s\"\n", title);
    fprintf(fout, "VARIABLES = \"X\", \"Y\"\n");
    fprintf(fout, "ZONE I=%d, J=%d, DATAPACKING=POINT\n", nx, ny);

    //printf("Printing Coordinate Information\n");
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            int ip = IU(i,j,nx);
            fprintf(fout, "%lf,\t%lf\n", x[ip], y[ip]);
        }
    }
    fclose(fout);


    FILE* fout2 = fopen("../Outputs/mesh.dat", "w");
    if (fout2 == nullptr) printf("oeups\n");
    fprintf(fout2, "%d %d\n", nx, ny);
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            int ip = IU(i,j,nx);
            fprintf(fout2, "%lf \t%lf\n", x[ip], y[ip]);
        }
    }
    //printing out boundary conditions
    fprintf(fout2, "BC\n");
    for (int ib=0; ib<nb; ib++){
        fprintf(fout2,"%d\n", ibound[ib]);
    }

    fclose(fout2);
}

int main() {
    // ========== Input Parameters (change to file input) ==========
    double height, length;
    int nx, ny;
    height = 2.0;
    length = 5.0;
    nx = 101;
    ny = 101;

    /*
    int bias; // 0:no bias, 1:bottom bias, 2:top bias
    bias = 0;
    double bias_ratio; // ratio of consecutive element size for bias option
    bias_ratio = 1.0;
    */
    double y_offset;   // Offset for axisymmetric applications
    y_offset = 0.0;


    /*
     * ==================== Geometry Input ====================
     * Need to represent bottom and top surfaces of geometry
     */
    double ramp_height = 0.25;
    double ramp_length = 1.5;

    /*
     * ==================== Mesh Generation ====================
     * ibound - flag for boundary condition (convention BL corner CCW)
     *          0-wall, 1-characteristic
     */
    int npoin = nx*ny;
    int nbound = 2*(nx-1) + 2*(ny-1);
    double dx, dy, ymin, ymax;
    auto* x = (double*)malloc(npoin*sizeof(double));
    auto* y = (double*)malloc(npoin*sizeof(double));
    auto* ibound = (int*)malloc(nbound*sizeof(int));

    dx = length / (nx-1);

    //define coordinates
    for (int i =0; i<nx; i++){
        ymax = height-nozzle_surface(i*dx, ramp_length, length-ramp_length)+y_offset;
                // ramp_surface(i*dx, ramp_height, ramp_length) + y_offset;
        ymin = y_offset;             // flat top
        dy = (ymax - ymin) / ny;

        for (int j=0; j<ny; j++){
            int ip = IU(i,j,nx);
            x[ip] = i*dx;
            y[ip] = j*dy + ymin;
        }
    }

    //========== Boundary cells
    //initialize with freestream
    for(int ib=0; ib<nbound; ib++)
        ibound[ib] = 1;


    //hardcoded nosecone stagnation point finder - sucks
    int istag=0;
    for(int ix=0; ix<nx; ix++){
        if(y[IU(ix,0,nx)] > 1e-10){
            istag = ix-1;
            break;
        }
    }
    //apply wall boundary condition
    for (int ib = istag; ib < nx-1; ib++) {
        ibound[ib] = 0;     //bottom surface
        ibound[ib+nx+ny-2-istag] = 0; //top surface
    }

    //Back Pressure
    for (int ib = nx-1; ib<nx+ny-2; ib++){
        ibound[ib] = 2;
    }




    // ==================== Output Mesh ====================
    printgrid("2DRampMesh", nx, ny, x, y, ibound);

}