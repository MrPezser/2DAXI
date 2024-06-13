/* 2D Structured Mesh Generation for Various Ramps
 *
 *
 *
 * T. Sailor Koeplinger - Jan2024
 */
#include <iostream>
#include <cmath>
#include <errno.h>
#include <string.h>

#define IN2M (0.0254)

#define IU(i, j, ni)  (((j)*(ni)) + (i))

double array_max(size_t n, double* a){
    double maxelem = a[0];
    for (int i=0; i<n; i++){
        maxelem = fmax(maxelem, a[i]);
    }
    return maxelem;
}

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
    if (fout == nullptr) printf("Unable to open grid file.\n%s\n", strerror(errno));

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

void get_nozzle(const int nx, const double* x, double* y){
    // Uses a given nozzle contour to find the maximum y values fr the given set of x points

    FILE* fcont = fopen("../axicontour.dat","r");
    if (fcont == nullptr) printf("Couldn't open nozzle contour file\n");
    int ncon;
    fscanf(fcont, " %d",&ncon);

    auto z = (double*)malloc(ncon*sizeof(double));
    auto r = (double*)malloc(ncon*sizeof(double));

    for (int icon=0; icon<ncon; icon++){
        fscanf(fcont,"%lf,%lf",&z[icon],&r[icon]);

        z[icon] -= z[0];
    }


    for (int ipoin=0; ipoin<nx; ipoin++){
        y[ipoin] = 0.0; //initialize
        //printf("Interpolating to point %d\n",ipoin);

        if (x[ipoin] >= z[ncon-1]) {y[ipoin] = IN2M * r[ncon-1];continue;}


        //find closest point on the left
        double mindel = 999.0;
        int ileft = -1;
        for (int j=0; j<ncon; j++) {
            if (x[ipoin]-z[j] < 0.0) continue;

            mindel = fmin(mindel, x[ipoin]-z[j]);
            if (mindel == x[ipoin]-z[j]) ileft = j;
        }

        if (ileft == -1){
            printf("Filed to find closest neightbor\n");
            exit(0);
        }

        //interpolate to point
        double drdz, dz;
        drdz = (r[ileft+1] - r[ileft]) / (z[ileft+1] - z[ileft]);
        dz = x[ipoin] - z[ileft];

        y[ipoin] = IN2M * (r[ileft] + (drdz*dz));

        /*
         * //Use lagrange polynomial to interpolate, be wary of spurious oscilations
         *
        for (int j=0; j<ncon; j++){
            double l = 1.0;

            for (int m=0; m<ncon; m++){
                if (j==m) continue;

                l *= (x[ipoin] - z[m]) / (z[j] - z[m]);

            }

            y[ipoin] += r[j] * l;
        }
         */

    }

    fclose(fcont);

}

int main() {
    // ========== Input Parameters (change to file input) ==========
    double height, length;
    int irefine, nx, ny, nyrefine{};
    height = 0.5;
    length = 7.75*IN2M - 0.1;//;21.07*IN2M;//1.0;//4.0;
    nx = 201;
    ny = 51;
    double bias = 1.0;
    double y_offset;   // Offset for axisymmetric applications
    y_offset = 0.0;//0.001;
    irefine = 0; //option to include a tanh distribution for boundary layoe on bottom surface

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
    auto* z = (double*)malloc(nx*sizeof(double));
    auto* r = (double*)malloc(nx*sizeof(double));
    auto* ibound = (int*)malloc(nbound*sizeof(int));

    dx = length / (nx-1);

    //Read in nozzle geometry
    for (int i=0; i<nx; i++) z[i] = (i*dx + 0.1)/IN2M;
    get_nozzle(nx, z, r);
    y_offset = 0.0 ;//0.05 * array_max(nx, r);
    /*
    //Modify upper surface to maintain the same area
    for (int i=0; i<nx; i++){
        double yoff2 = sqrt(r[i]*r[i] + y_offset*y_offset) - r[i];
        if (yoff2 < 0){
            printf("Upper surface modification negative\n");
            exit(0);
        }

        r[i] += yoff2;

    }
     */

    //define coordinates
    double factor = 2.0;
    for (int i = 0; i < nx; i++) {
        double xi = i * dx;
        ymax = r[i];//height + y_offset;
        ymin = y_offset; //ramp_surface(xi, ramp_height, ramp_length);
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
                if (i == nx-1) printf("y, %lf\n", y[ip]);
            }
        }
    }


    //   ==================== Boundary cells ====================
    // 0 = viscous wall
    // 1 = freestream
    // 2 = backpressure (kinda not used)
    // 3 = outflow / extrapolation
    // 4 = inviscid wall / symmetry
    // 5 = bernouli inflow
    //initialize with freestream
    for(int ib=0; ib<nbound; ib++)
        ibound[ib] = 1;

    //top/bot surf
    for (int ib = 0; ib < nx-1; ib++) {
        int itop = -ib + 2 * nx + ny - 4;
        int ibot = ib;

        //if (ib*dx > 0.2) ibound[ibot] = 4;       //bot surf
        ibound[ibot] = 4;

        //if (ib*dx > 2.0) ibound[itop] = 3;   //top surface
        ibound[itop] = 4;
    }
    //Back Pressure (2) or outflow (3)
    for (int ib = nx-1; ib<nx+ny-2; ib++){
        int iback = ib;
        int ifront = ib+nx+ny-2;

        ibound[iback] = 3;      //      right/exit
        ibound[ifront] = 5;
    }


    // ==================== Output Mesh ====================
    printgrid("2DRampMesh", nx, ny, x, y, ibound);

}