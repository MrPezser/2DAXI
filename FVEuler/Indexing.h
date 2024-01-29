//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_INDEXING_H
#define FVEULER_INDEXING_H

#define IJ(i, j, ni)  (((j)*(ni)) + (i))
#define IJK(i, j, k, ni, nk)  ((((j)*(ni)) + (i))*(nk) + (k))
#define NVAR 4

#endif //FVEULER_INDEXING_H
