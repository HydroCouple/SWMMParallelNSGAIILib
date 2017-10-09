#ifndef SWMMPARALLELNSGAII_H
#define SWMMPARALLELNSGAII_H

#include <list>
#include <string>
#include "swmmparallelnsgaii_global.h"

#ifdef __cplusplus
extern "C" {
#endif

 void SWMMParallelNSGAII_EXPORT SWMMParallelNSGAII(int gen, int nreal, double *xreal, int nbin, double *xbin, int *nbits,
                                                  int **gene, int nobj, double *obj, int ncon, double *constr, const std::list<std::string>& optionalArgs);

#ifdef __cplusplus
}
#endif

#endif // SWMMPARALLELNSGAII_H
