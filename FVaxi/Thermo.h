//
// Created by tskoepli on 4/9/2024.
//

#ifndef FVNS_THERMO_H
#define FVNS_THERMO_H

#include "Indexing.h"

class Thermo {
private:
    double Al[NSP]{}, Bl[NSP]{}, Cl[NSP]{}, Dl[NSP]{}, El[NSP]{}, Fl[NSP]{}, Gl[NSP]{};
    double Ah[NSP]{}, Bh[NSP]{}, Ch[NSP]{}, Dh[NSP]{}, Eh[NSP]{}, Fh[NSP]{}, Gh[NSP]{};

    void LoadCurveFits(){
        int isp = 0;
        //Function for reading in the therm.dat file and loading in thermo curve fits for 5 species air
        FILE* ftherm = fopen("../therm_air.dat","r");
        ASSERT(ftherm != nullptr, "Unable to open thermo file");

        //get molecular weight from first line
        //fscanf(ftherm, "%*18s %*6c %*1c%*f%*1c %*f %*f %*2f %*s %*f %*f %lf %*d", &(Mw[isp]));
        fscanf(ftherm,"%*[^\n]\n"); //skip 1st line
        Mw[isp] = 28.9;
        //get curve fits
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Ah[isp]),&(Bh[isp]),&(Ch[isp]),&(Dh[isp]),&(Eh[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Fh[isp]),&(Gh[isp]),&(Al[isp]),&(Bl[isp]),&(Cl[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%*15lE%*d",&(Dl[isp]),&(El[isp]),&(Fl[isp]),&(Gl[isp]));

        Rs[isp] = Ruv/Mw[isp];

        fclose(ftherm);
    }

public:
    double Ruv = 8314.34;
    double gam = 1.28;//1.4;
    double Mw[NSP]{}, Rs[NSP]{};
    double CalcEnthalpy(double T);
    double CalcCp(double T);

    Thermo() {
        LoadCurveFits();
    }

};


#endif //FVNS_THERMO_H
