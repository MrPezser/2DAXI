//
// Created by tskoepli on 2/2/2024.
//

#ifndef FVNS_MESHMODULE_H
#define FVNS_MESHMODULE_H

void calc_geoel_geofa(const int nx, const int ny, double* x, double* y, \
            double** geoel, double** geofa);

#endif //FVNS_MESHMODULE_H
