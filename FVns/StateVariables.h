//
// Created by tskoepli on 4/8/2024.
//

#ifndef FVNS_STATEVARIABLES_H
#define FVNS_STATEVARIABLES_H


//
// Created by tskoepli on 3/27/2024.
//

#ifndef PROJECT2_STATEVARIABLES_H
#define PROJECT2_STATEVARIABLES_H

#include <cmath>
#include "Indexing.h"
// class for storing calculated properties in an element

class State {

private:
    const double* unk{};

public:
    double p{NAN},a{NAN}, h{NAN}, h0{NAN}, v2{NAN}, Mw[NSP]{}, Rs[NSP]{}, Cp[NSP]{}, Cv{};
    double Al[NSP]{}, Bl[NSP]{}, Cl[NSP]{}, Dl[NSP]{}, El[NSP]{}, Fl[NSP]{}, Gl[NSP]{};
    double Ah[NSP]{}, Bh[NSP]{}, Ch[NSP]{}, Dh[NSP]{}, Eh[NSP]{}, Fh[NSP]{}, Gh[NSP]{};
    double Ruv = 8314.34;

    State(){
        LoadCurveFits();
    };

    void Initialize(const double* u){
        unk = u;
    }

    void UpdateState(double gam) {
        p=0;
        h0=0.0;
        a=0.0;
        int isp = 0;

        v2 = unk[1]*unk[1] + unk[2]*unk[2];
        p = unk[0]*Rs[isp]*unk[3];
        a = sqrt(gam*Rs[isp]*unk[3]);
        h = CalcEnthalpy(unk[3]);
        h0 = h + 0.5*v2;
        Cp[0] = CalcCp(unk[3]);
        Cv = Cp[0] - Rs[0]; //total cv, not species.... not that it matters now

        ASSERT(!_isnan(p*a), "Error in finding pressure or wavespeed.")
    }

    double CalcEnthalpy(double T){
        // Calculate species specific enthalpy according to curve fit alone
        // isp = index of species
        int isp = 0;
        // T = temperature
        // Cp = output Cp at temperature

        double T2,T3,T4,T5;
        T2 = T*T;
        T3 = T*T2;
        T4 = T*T3;
        T5 = T*T4;

        ASSERT(T > 200.0, "Temp below curve fit")
        ASSERT(T < 6000.0, "Temp above curve fit")

        if (T > 1000.0) {
            return (Rs[isp])*(Ah[isp]*T + Bh[isp]*T2*0.5 + Ch[isp]*T3/3.0 + Dh[isp]*T4*0.25 + Eh[isp]*T5*0.2 + Fh[isp]);
        } else {  //T<1000
            return (Rs[isp])*(Al[isp]*T + Bl[isp]*T2*0.5 + Cl[isp]*T3/3.0 + Dl[isp]*T4*0.25 + El[isp]*T5*0.2 + Fl[isp]);
        }
    }

    double CalcCp(double T){
        int isp = 0;
        double T2,T3,T4,T5;
        T2 = T*T;
        T3 = T*T2;
        T4 = T*T3;

        ASSERT(T > 200.0, "Temp below curve fit")
        ASSERT(T < 6000.0, "Temp above curve fit")

        if (T >= 1000.0) {
            return (Rs[isp])*(Ah[isp] + Bh[isp]*T + Ch[isp]*T2 + Dh[isp]*T3 + Eh[isp]*T4);
        } else {  //T<1000
            return (Rs[isp])*(Al[isp] + Bl[isp]*T + Cl[isp]*T2 + Dl[isp]*T3 + El[isp]*T4);
        }
    }

    void LoadCurveFits(){
        int isp = 0;
        //Function for reading in the therm.dat file and loading in thermo curve fits for 5 species air
        FILE* ftherm = fopen("../therm_airmix.dat","r");

        //get molecular weight from first line
        fscanf(ftherm, "%*18s %*6c %*1c%*f%*1c %*f %*f %*2f %*s %*f %*f %lf %*d", &(Mw[isp]));
        //get curve fits
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Ah[isp]),&(Bh[isp]),&(Ch[isp]),&(Dh[isp]),&(Eh[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%15lE%*d",&(Fh[isp]),&(Gh[isp]),&(Al[isp]),&(Bl[isp]),&(Cl[isp]));
        fscanf(ftherm, "%15lE%15lE%15lE%15lE%*15lE%*d",&(Dl[isp]),&(El[isp]),&(Fl[isp]),&(Gl[isp]));

        Rs[isp] = Ruv/Mw[isp];

        fclose(ftherm);
    }

};

#endif //PROJECT2_STATEVARIABLES_H

#endif //FVNS_STATEVARIABLES_H
