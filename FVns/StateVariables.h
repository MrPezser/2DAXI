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
    const double* unk{};


public:
    double p{NAN},a{NAN}, h{NAN}, h0{NAN}, v2{NAN}, Cp[NSP]{}, Cv{};


    State() = default;

    void Initialize(const double* u){
        unk = u;
    }

    void UpdateState(Thermo& air ) {
        int isp = 0;

        v2 = unk[1]*unk[1] + unk[2]*unk[2];
        p = unk[0]*air.Rs[isp]*unk[3];
        a = sqrt(air.gam*air.Rs[isp]*unk[3]);
        h =air.CalcEnthalpy(unk[3]);
        h0 = h + 0.5*v2;
        Cp[0] = air.CalcCp(unk[3]);
        Cv = Cp[0] - air.Rs[isp]; //total cv, not species.... not that it matters now

        ASSERT(!_isnan(p*a), "Error in finding pressure or wavespeed.")
        CHECKD(a > 0.0, "bad wave speed", a)
        CHECKD(p > 0.0, "bad pressure", p)
    }

};

#endif //PROJECT2_STATEVARIABLES_H

#endif //FVNS_STATEVARIABLES_H
