#ifndef SWMMPARALLELNSGAII_GLOBAL_H
#define SWMMPARALLELNSGAII_GLOBAL_H


#  ifdef Q_OS_WIN
#    define SWMMParallelNSGAII_EXPORT     __declspec(dllexport)
#    define SWMMParallelNSGAII_IMPORT     __declspec(dllimport)
#  else
#    define SWMMParallelNSGAII_EXPORT     __attribute__((visibility("default")))
#    define SWMMParallelNSGAII_IMPORT     __attribute__((visibility("default")))
#    define SWMMParallelNSGAII_HIDDEN     __attribute__((visibility("hidden")))
#  endif

#ifdef __cplusplus
extern "C" {
#endif

typedef  void SWMMParallelNSGAII_EXPORT (*problemDef)(int gen, int nreal, double *xreal, int nbin, double *xbin, int *nbits, int **gene, int nobj, double *obj, int ncon, double *constr, const std::list<std::string>& optionalArgs);

#ifdef __cplusplus
}
#endif

#endif // SWMMPARALLELNSGAII_GLOBAL_H
