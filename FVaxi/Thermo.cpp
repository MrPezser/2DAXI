//
// Created by tskoepli on 4/9/2024.
//

#include "Thermo.h"

double Thermo::CalcEnthalpy(double T){
    // Calculate species specific enthalpy according to curve fit alone
    // isp = index of species
    int isp = 0;
    // T = temperature
    // Cp = output Cp at temperature

    if (IGAM > 0) { return (T * Rs[isp]/(1.0 - 1.0/(IGAM)) + Fh[isp]);}

    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = T*T2;
    T4 = T*T3;
    T5 = T*T4;

    ASSERT(T > 200.0, "Temp below curve fit")
    ASSERT(T < 10000.0, "Temp above curve fit")

    if (T > 1000.0) {
        return (Rs[isp])*(Ah[isp]*T + Bh[isp]*T2*0.5 + Ch[isp]*T3/3.0 + Dh[isp]*T4*0.25 + Eh[isp]*T5*0.2 + Fh[isp]);
    } else {  //T<1000
        return (Rs[isp])*(Al[isp]*T + Bl[isp]*T2*0.5 + Cl[isp]*T3/3.0 + Dl[isp]*T4*0.25 + El[isp]*T5*0.2 + Fl[isp]);
    }
}

double Thermo::CalcCp(double T){
    int isp = 0;
    if (IGAM > 0) { return Rs[isp]/(1.0 - 1.0/(IGAM));}
    double T2,T3,T4,T5;
    T2 = T*T;
    T3 = T*T2;
    T4 = T*T3;

    ASSERT(T > 200.0, "Temp below curve fit")
    ASSERT(T < 10000.0, "Temp above curve fit")

    if (T >= 1000.0) {
        return (Rs[isp])*(Ah[isp] + Bh[isp]*T + Ch[isp]*T2 + Dh[isp]*T3 + Eh[isp]*T4);
    } else {  //T<1000
        return (Rs[isp])*(Al[isp] + Bl[isp]*T + Cl[isp]*T2 + Dl[isp]*T3 + El[isp]*T4);
    }
}