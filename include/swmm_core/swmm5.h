//-----------------------------------------------------------------------------
//   swmm5.h
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/24/14  (Build 5.1.001)
//   Author:  L. Rossman
//
//   Prototypes for SWMM5 functions exported to swmm5.dll.
//
//-----------------------------------------------------------------------------
#ifndef SWMM5_H
#define SWMM5_H

#include "swmmparallelnsgaii_global.h"

// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
#define WINDOWS
#endif
#ifdef __WIN32__
#define WINDOWS
#endif

// --- use "C" linkage for C++ programs

#ifdef __cplusplus
extern "C" { 
#endif 

typedef struct Project Project;
typedef double DateTime;
typedef struct TNode TNode;
typedef struct TLink TLink;
typedef struct TSubcatch TSubcatch;

int  SWMMParallelNSGAII_EXPORT  swmm_run(char* f1, char* f2, char* f3);
Project* SWMMParallelNSGAII_EXPORT  swmm_open(char* f1, char* f2, char* f3);
int  SWMMParallelNSGAII_EXPORT   swmm_start(Project* project, int saveFlag);
int  SWMMParallelNSGAII_EXPORT   swmm_step(Project* project, double *elapsedTime);
int  SWMMParallelNSGAII_EXPORT   swmm_end(Project* project);
int  SWMMParallelNSGAII_EXPORT   swmm_report(Project* project);
int  SWMMParallelNSGAII_EXPORT   swmm_getMassBalErr(Project* project, float* runoffErr, float* flowErr,float* qualErr);
int  SWMMParallelNSGAII_EXPORT   swmm_close(Project* project);
int  SWMMParallelNSGAII_EXPORT   swmm_getVersion();

//hydrocouple
double SWMMParallelNSGAII_EXPORT getRoutingStep(Project *project);
int SWMMParallelNSGAII_EXPORT performRoutingIteration(Project *project, double *elapsedTime);

//OpenMI 2.0

int  SWMMParallelNSGAII_EXPORT   swmm_getErrorCode(Project* project);
double SWMMParallelNSGAII_EXPORT swmm_getDateTime(Project* project, char* beginorend);
void SWMMParallelNSGAII_EXPORT  datetime_decodeDateTime(DateTime dateTime, int* y, int* m, int* d, int* h, int* mm, int* s);
char * getErrorMsg(int errorCode);
int SWMMParallelNSGAII_EXPORT getObjectTypeCount(Project* project, int type);

//TNode
TNode* SWMMParallelNSGAII_EXPORT getNode(Project* project, int index);
TNode* SWMMParallelNSGAII_EXPORT getNodeById(Project* project, char* id);
void SWMMParallelNSGAII_EXPORT setNode(Project* project, char* nodeId, char* propertyName, double value);

//TLink
TLink*  SWMMParallelNSGAII_EXPORT getLink(Project* project, int index);
TLink* SWMMParallelNSGAII_EXPORT getLinkById(Project* project, char* id);
void SWMMParallelNSGAII_EXPORT setLink(Project* project, char* linkId, char* propertyName, double value);

//TSubcatch
TSubcatch*  SWMMParallelNSGAII_EXPORT getSubcatch(Project* project, int index);
TSubcatch* SWMMParallelNSGAII_EXPORT getSubcatchById(Project* project, char* id);
void SWMMParallelNSGAII_EXPORT setSubcatch(Project* project, char* subCatchId, char* propertyName, double value);

#ifdef __cplusplus 
}   // matches the linkage specification from above */ 
#endif

#endif
