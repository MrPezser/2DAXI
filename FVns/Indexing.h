//
// Created by tskoepli on 1/27/2024.
//

#ifndef FVEULER_INDEXING_H
#define FVEULER_INDEXING_H
#include <cstdio>
#include <cstdlib>
#include <cmath>

#define IJ(i, j, ni)  (((j)*(ni)) + (i))
#define IJK(i, j, k, ni, nk)  ((((j)*(ni)) + (i))*(nk) + (k))
#define NVAR 4
#define NSP 1
#define IVISC 1
#define MXITER 1e6
#define ITGLOBAL 1e6
#define RESTOL 1e-6
#define IDUNG 0
#define IIMPLI 1

#define sign(x)  ((std::signbit(x) ?  -1 : 1))
#define ASSERT(cond, msg) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n", __FILE__, __LINE__, #cond, msg); exit(0);}
#define CHECKD(cond, msg, val) if(!(cond)){printf("Failed Assert: %s:%u %s\n %s\n Val:%lf\n", __FILE__, __LINE__, #cond, msg, val); exit(0);}

#if IDUNG
#define DUNG(msg) printf("DUNG: %s\n",msg);
#else
#define DUNG(msg) nullptr;
#endif

#endif //FVEULER_INDEXING_H
