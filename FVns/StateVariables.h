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
    double* unk{};

public:
    double p{},a{},h0{}, u{}, v{}, M{};

    State() = default;

    void Initialize(double* u){
        unk = u;
    }

    void UpdateState(double gam) {
        p=0;
        h0=0.0;
        a=0.0;
        u = 0.0;
        v = 0.0;
        M = 0.0;

        u = unk[1] / unk[0];
        v = unk[2] / unk[0];

        double v2 = u*u + v*v;
        p = fmax((gam - 1) * (unk[3] - (0.5 * unk[0] * v2)), 1e-8);


        a = sqrt(gam*p / unk[0]);
    }

};

#endif //PROJECT2_STATEVARIABLES_H

#endif //FVNS_STATEVARIABLES_H
