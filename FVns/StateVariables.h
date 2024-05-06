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
#include "Thermo.h"
// class for storing calculated properties in an element

class State {

private:
    double* unk{};


public:
    double p{NAN},a{NAN}, h{NAN}, h0{NAN}, v2{NAN}, Cp[NSP]{}, Cv{}, mu{}, e{};


    State() = default;

    void Initialize(double* u){
        unk = u;
        // vars = [rho, u, v, T]
    }

    void UpdateState(Thermo& air ) {
        int isp = 0;
        unk[3] = fmin(unk[3], 9999.0);
        unk[3] = fmax(unk[3], 201.0);
        double T = unk[3];

        v2 = unk[1]*unk[1] + unk[2]*unk[2];
        p = unk[0]*air.Rs[isp]*T;
        a = sqrt(air.gam*air.Rs[isp]*T);
        h =air.CalcEnthalpy(T);
        h0 = h + 0.5*v2;
        e = h - air.Rs[0]*T;
        Cp[0] = air.CalcCp(T);
        Cv = Cp[0] - air.Rs[isp]; //total cv, not species.... not that it matters now

        //Sutherland's law for viscosity
        double S, C1;
        S = 110.4;
        C1 = 1.458e-6;
        mu = C1 * pow(T, 1.5) / (T + S);

        CHECKD(a > 0.0, "bad wave speed", a)
        CHECKD(p > 0.0, "bad pressure", p)
        ASSERT(!__isnan(p*a*T), "Error in finding pressure or wavespeed or temperature.")
    }

};

#endif //PROJECT2_STATEVARIABLES_H

#endif //FVNS_STATEVARIABLES_H
