// swmm5_iface.h
//
// Header file for SWMM 5 interfacing functions
//
// #include this file in any C module that references the functions
// contained in swmm5_iface.c.
//

#ifndef SWMM5_IFACE_H
#define SWMM5_IFACE_H

#include "globals.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct IFaceData
{
   int  SWMM_Nperiods;           // number of reporting periods
   int  SWMM_FlowUnits;          // flow units code
   int  SWMM_Nsubcatch;          // number of subcatchments
   int  SWMM_Nnodes;             // number of drainage system nodes
   int  SWMM_Nlinks;             // number of drainage system links
   int  SWMM_Npolluts;           // number of pollutants tracked
   double SWMM_StartDate;          // start date of simulation
   int  SWMM_ReportStep;         // reporting time step (seconds)
   int SubcatchVars;               // number of subcatch reporting variables
   int NodeVars;                   // number of node reporting variables
   int LinkVars;                   // number of link reporting variables
   int SysVars;                    // number of system reporting variables
   int    StartPos;                // file position where results start
   int    BytesPerPeriod;          // bytes used for results in each period
   FILE *Fout;
}IFaceData;


int    RunSwmmExe(char* cmdLine);
int    RunSwmmDll(char* inpFile, char* rptFile, char* outFile);
IFaceData *OpenSwmmOutFile(char* outFile, int *error);
int GetSwmmResult(IFaceData *faceData, int iType, int iIndex, int vIndex, int period, float* value);
void CloseSwmmOutFile(IFaceData *faceData);

#ifdef __cplusplus
}
#endif

#endif // SWMM5_IFACE_H
