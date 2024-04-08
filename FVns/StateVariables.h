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
    double p{},a{NAN},h0{}, vx{}, vy{}, v2{};

    State() = default;

    void Initialize(const double* u){
        unk = u;
    }

    void UpdateState(double gam) {
        p=0;
        h0=0.0;
        a=0.0;
        vx = 0.0;
        vy = 0.0;

        vx = unk[1] / unk[0];
        vy = unk[2] / unk[0];

        v2 = vx*vx + vy*vy;
        p = (gam - 1) * (unk[3] - (0.5 * unk[0] * v2));
        a = sqrt(gam*p / unk[0]);

        ASSERT(!_isnan(p*a), "Error in finding pressure or wavespeed.")
    }

};

#endif //PROJECT2_STATEVARIABLES_H

#endif //FVNS_STATEVARIABLES_H
