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
    /*
    double l0 = 0.5;//1.0;
    if(x < l0) {
        return 0;
    } else if( x > l0+L) {
        return h;
    } else {
        return (h/L)*(x-l0);
    }*/

    //alternate geometry, converging diverging
    double l0 = L;
    if(x < l0) {
        return (h/L)*(x);
    } else if( x > l0+L) {
        return h;
    } else {
        return h - (h/L)*(x-l0);
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
            fprintf(fout, "%20.16lf,%20.16lf\n", x[ip], y[ip]);
        }
    }
    fclose(fout);


    FILE* fout2 = fopen("../Outputs/mesh.dat", "w");
    if (fout2 == nullptr) printf("oeups\n");
    fprintf(fout2, "%d %d\n", nx, ny);
    for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {
            int ip = IU(i,j,nx);
            fprintf(fout2, "%20.16lf %20.16lf\n", x[ip], y[ip]);
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
    int irefine, nx, ny, nyrefine{};
    height = 1.0;
    length = 2.0;//4.0;
    nx = 101;
    ny = 51;
    double bias = 1.0;
    double y_offset;   // Offset for axisymmetric applications
    y_offset = 0.0;
    irefine = 0; //option to include a tanh distribution for boundary layoe on bottom surface

    /*
     * ==================== Geometry Input ====================
     * Need to represent bottom and top surfaces of geometry
     */
    double ramp_height = 0.25;//0.3;
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
    double factor = 2.0;
    for (int i = 0; i < nx; i++) {
        double xi = i * dx;
        ymax = height + y_offset;
        ymin = ramp_surface(xi, ramp_height, ramp_length);
        dy = (ymax - ymin) / (ny - 1);

        if (irefine==1) {
            double delmin, del, etamin;
            int iconv;
            iconv = 0;
            etamin = 0.05;
            delmin = etamin*sqrt(1.789e-5 * 0.2 / 10.0 );//;1e-6; //adjust to get desired Y+ or eta value

            while (iconv == 0 && i == 0) {
                del = 0.0;
                for (int j = 0; j < ny; j++) {
                    //printf("j:%d \t del:%le\n",j,del);
                    int ip = IU(i, j, nx);
                    x[ip] = xi;
                    if (j == 0) {
                        y[ip] = ymin;
                        del = delmin * pow(factor, j);
                        continue;
                    }
                    y[ip] = y[IU(i, j-1, nx)] + del;
                    //if (i == 0) printf("y, %le\n", y[ip]);

                    del = delmin * pow(factor, j);
                    if (y[ip] >= ymax) {
                        factor = 0.99 * factor + 0.01 * 1.0;
                        iconv = -1;
                        printf("y failed max, %le\n", y[ip]);
                        break;
                    }
                }

                if (iconv == 0) printf("y final max, %lf\n", y[IU(i, ny - 1, nx)]);
                printf("Growth Factor: %lf\n", factor);
                iconv += 1;
            }

            if (i >= 1) {
                for (int j = 0; j < ny; j++) {
                    int ip = IU(i, j, nx);
                    y[ip] = y[IU(i-1,j,nx)];
                    x[ip] = xi;

                }
            }
        } else {

            for (int j = 0; j < ny; j++) {
                int ip = IU(i, j, nx);
                x[ip] = xi;
                y[ip] = (j * dy) + ymin;
                if (i == 0) printf("y, %lf\n", y[ip]);
            }
        }
    }


    //   ==================== Boundary cells ====================
    // 0 = viscous wall
    // 1 = freestream
    // 2 = backpressure (kinda not used)
    // 3 = outflow / extrapolation
    //initialize with freestream
    for(int ib=0; ib<nbound; ib++)
        ibound[ib] = 1;

    //top/bot surf
    for (int ib = 0; ib < nx-1; ib++) {
        int itop = -ib + (2 * nx) + ny - 3 - 1;
        int ibot = ib;

        if (ib*dx > 0.5) { ///////////////
            ibound[ibot] = 4;       //bot surf
        }

        ibound[itop] = 3;   //top surface

    }
    //Back Pressure (2) or outflow (3)
    for (int ib = nx-1; ib<nx+ny-2; ib++){
        int iback = ib;
        int ifront = ib+nx+ny-2;

        ibound[iback] = 3;      //      right/exit
        //ibound[ifront] = 0;
    }


    // ==================== Output Mesh ====================
    printgrid("2DRampMesh", nx, ny, x, y, ibound);

}