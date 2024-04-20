/* 2D Structured Mesh Generation for Various Ramps
 *
 *
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#include <cmath>

#define IU(i, j, ni)  (((j)*(ni)) + (i))

//test ramp geometry for now
double ramp_surface(double x, double h, double L) {
    //test geometry: straight - ramp - straight
    double l0 = 0.5;//1.0;
    if(x < l0) {
        return 0;
    } else if( x > l0+L) {
        return h;
    } else {
        return (h/L)*(x-l0);
    }
}

double nozzle_surface(double x, double L1, double L2) {
    //test geometry: Mach 3.0 CD nozzle
    double h =  0.555;//0.81481480 - 0.01481480;
    double h2 = 0.365;
    if(x < 0.0) {
        return 0;
    } else if( x > L1) {
        return x*(h-h2)/(L1-L2) + (-h*L2 + h2*L1)/(L1-L2);
        // h * ((L2-(x-L1))/L2);
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
    length = 1.5;//4.0;
    nx = 201;
    ny = 201;
    double bias = 1.0;
    double y_offset;   // Offset for axisymmetric applications
    y_offset = 0.0;


    /*
     * ==================== Geometry Input ====================
     * Need to represent bottom and top surfaces of geometry
     */
    double ramp_height = 0.75;//0.3;
    double ramp_length = 1.0;//1.0;

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

        double xi = i*dx;
        ymax = height + y_offset;
        ymin = ramp_surface(xi,ramp_height, ramp_length);
        dy = (ymax - ymin) / ny;

        for (int j=0; j<ny; j++){
            int ip = IU(i,j,nx);
            x[ip] = xi;
            //y[ip] = j*dy + ymin; //equal spacing
            double xi = (double)j/(ny-1);
            double k=6;
            double f = 1.0 / (1.0 + exp(-k*(xi-0.5)));//(cbrt(xi-0.5) + 2 - cbrt(0.5)) / (2.0);
            double f1 = 1.0 / (1.0 + exp(-k*(1.0-0.5)));
            y[ip] = (bias)*(f/f1)*dy*ny + (1.0-bias)*(j*dy) + ymin; //biased spacing
            if (i==0) printf("y, %lf\n", y[ip]);
        }
    }


    //   ==================== Boundary cells ====================
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
    //full top/bot surf
    for (int ib = 0; ib < nx-1; ib++) {
        if (ib*dx > -1.0) { ///////////////
            ibound[ib] = 0; //bot surf
            if (ib*dx > 2.0) {
                ibound[-ib + 2 * nx + ny - 3] = 3; //top surface
            }
        }
    }
    //Back Pressure (2) or outflow (3)
    for (int ib = nx-1; ib<nx+ny-2; ib++){
        ibound[ib] = 3;
        //ibound[ib+nx+ny-2] = 0;
    }


    // ==================== Output Mesh ====================
    printgrid("2DRampMesh", nx, ny, x, y, ibound);

}