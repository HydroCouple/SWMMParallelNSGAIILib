// swmm5_iface.c
//
// Example code for interfacing SWMM 5 with C/C++ programs.
//
// Remember to #include the file swmm5_iface.h in the calling program.

#include <stdio.h>
#include "swmm5_iface.h"
#include "swmm5.h"
#include "globals.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

//#include "globals.h"

//int    faceData->SWMM_Nperiods;                  // number of reporting periods
//int    faceData->SWMM_FlowUnits;                 // flow units code
//int    faceData->SWMM_Nsubcatch;                 // number of subcatchments
//int    faceData->SWMM_Nnodes;                    // number of drainage system nodes
//int    faceData->SWMM_Nlinks;                    // number of drainage system links
//int    faceData->SWMM_Npolluts;                  // number of pollutants tracked
//double faceData->SWMM_StartDate;                 // start date of simulation
//int    faceData->SWMM_ReportStep;                // reporting time step (seconds)


//int    RunSwmmExe(char* cmdLine);
//int    RunSwmmDll(char* inpFile, char* rptFile, char* outFile);
//int    OpenSwmmOutFile(char* outFile);
//int    GetSwmmResult(int iType, int iIndex, int vIndex, int period, float* value);
//void   CloseSwmmOutFile(void);

static const int IFACESUBCATCH = 0;
static const int IFACENODE     = 1;
static const int IFACELINK     = 2;
static const int IFACESYS      = 3;
static const int RECORDSIZE = 4;       // number of bytes per file record


static void   ProcessMessages(void);

//-----------------------------------------------------------------------------
int RunSwmmExe(char* cmdLine)
//-----------------------------------------------------------------------------
{
  //  int exitCode;
  //  STARTUPINFO si;
  //  PROCESS_INFORMATION  pi;

  //  // --- initialize data structures
  //  memset(&si, 0, sizeof(si));
  //  memset(&pi, 0, sizeof(pi));
  //  si.cb = sizeof(si);
  //  si.wShowWindow = SW_SHOWNORMAL;

  //  // --- launch swmm5.exe
  //  exitCode = CreateProcess(NULL, cmdLine, NULL, NULL, 0,
  //			 0, NULL, NULL, &si, &pi);

  //  // --- wait for program to end
  //  exitCode = WaitForSingleObject(pi.hProcess, INFINITE);

  //  // --- retrieve the error code produced by the program
  //  GetExitCodeProcess(pi.hProcess, &exitCode);

  //  // --- release handles
  //  CloseHandle(pi.hProcess);
  //  CloseHandle(pi.hThread);
  //  return exitCode;
}


//-----------------------------------------------------------------------------
int RunSwmmDll(char* inpFile, char* rptFile, char* outFile)
//-----------------------------------------------------------------------------
{
  int err;
  double elapsedTime;
  Project *currentProject =  NULL;

  // --- open a SWMM project

  currentProject = swmm_open(inpFile, rptFile, outFile);

  if (!currentProject->ErrorCode)
  {
    // --- initialize all processing systems
    err = swmm_start(currentProject, TRUE);

    if (err == 0)
    {
      // --- step through the simulation
      do
      {
        // --- allow Windows to process any pending events
        //        ProcessMessages();

        // --- extend the simulation by one routing time step
        err = swmm_step(currentProject, &elapsedTime);

        /////////////////////////////////////////////
        // --- call progress reporting function here,
        //     using elapsedTime as an argument
        /////////////////////////////////////////////

      } while (elapsedTime > 0.0 && err == 0);

      // --- close all processing systems
      swmm_end(currentProject);
    }
  }

  if (currentProject->Fout.mode == SCRATCH_FILE)
    swmm_report(currentProject);

  // --- close the project
  return swmm_close(currentProject);
}


//-----------------------------------------------------------------------------
void ProcessMessages(void)
//-----------------------------------------------------------------------------
{

  /****  Only use this function with a Win32 application *****
  MSG msg;
  while (TRUE)
  {
    if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
    {
      if (msg.message == WM_QUIT) break;
      else
      {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
      }
    }
    else break;
  }
***********************************************************/

}


//-----------------------------------------------------------------------------
IFaceData *OpenSwmmOutFile(char* outFile, int *error)
//-----------------------------------------------------------------------------
{

  IFaceData *faceData = (IFaceData*) calloc(1, sizeof(IFaceData));

  int magic1, magic2, errCode, offset, offset0, version;
  int err;

  // --- open the output file
  faceData->Fout = fopen(outFile, "rb");
  if (faceData->Fout == NULL) return 2;

  // --- check that file contains at least 14 records
  fseek(faceData->Fout, 0L, SEEK_END);
  if (ftell(faceData->Fout) < 14*RECORDSIZE)
  {
    fclose(faceData->Fout);
    free(faceData);
    faceData = NULL;
    *error = 1;
    return faceData;
  }

  // --- read parameters from end of file
  fseek(faceData->Fout, -5*RECORDSIZE, SEEK_END);
  fread(&offset0, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->StartPos, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_Nperiods, RECORDSIZE, 1, faceData->Fout);
  fread(&errCode, RECORDSIZE, 1, faceData->Fout);
  fread(&magic2, RECORDSIZE, 1, faceData->Fout);

  // --- read magic number from beginning of file
  fseek(faceData->Fout, 0L, SEEK_SET);
  fread(&magic1, RECORDSIZE, 1, faceData->Fout);

  // --- perform error checks
  if (magic1 != magic2) err = 1;
  else if (errCode != 0) err = 1;
  else if (faceData->SWMM_Nperiods == 0) err = 1;
  else err = 0;

  // --- quit if errors found
  if (err > 0 )
  {
    fclose(faceData->Fout);
    faceData->Fout = NULL;
    free(faceData);
    faceData = NULL;
    *error = err;
    return faceData;
  }

  // --- otherwise read additional parameters from start of file
  fread(&version, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_FlowUnits, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_Nsubcatch, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_Nnodes, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_Nlinks, RECORDSIZE, 1, faceData->Fout);
  fread(&faceData->SWMM_Npolluts, RECORDSIZE, 1, faceData->Fout);

  // Skip over saved subcatch/node/link input values
  offset = (faceData->SWMM_Nsubcatch+2) * RECORDSIZE  // Subcatchment area
           + (3*faceData->SWMM_Nnodes+4) * RECORDSIZE  // Node type, invert & max depth
           + (5*faceData->SWMM_Nlinks+6) * RECORDSIZE; // Link type, z1, z2, max depth & length
  offset = offset0 + offset;
  fseek(faceData->Fout, offset, SEEK_SET);

  // Read number & codes of computed variables
  fread(&faceData->SubcatchVars, RECORDSIZE, 1, faceData->Fout); // # Subcatch variables
  fseek(faceData->Fout, faceData->SubcatchVars*RECORDSIZE, SEEK_CUR);
  fread(&faceData->NodeVars, RECORDSIZE, 1, faceData->Fout);     // # Node variables
  fseek(faceData->Fout, faceData->NodeVars*RECORDSIZE, SEEK_CUR);
  fread(&faceData->LinkVars, RECORDSIZE, 1, faceData->Fout);     // # Link variables
  fseek(faceData->Fout, faceData->LinkVars*RECORDSIZE, SEEK_CUR);
  fread(&faceData->SysVars, RECORDSIZE, 1, faceData->Fout);     // # System variables

  // --- read data just before start of output results
  offset = faceData->StartPos - 3*RECORDSIZE;
  fseek(faceData->Fout, offset, SEEK_SET);
  fread(&faceData->SWMM_StartDate, sizeof(double), 1, faceData->Fout);
  fread(&faceData->SWMM_ReportStep, RECORDSIZE, 1, faceData->Fout);

  // --- compute number of bytes of results values used per time period
  faceData->BytesPerPeriod = 2*RECORDSIZE +      // date value (a double)
                             (faceData->SWMM_Nsubcatch*faceData->SubcatchVars +
                              faceData->SWMM_Nnodes*faceData->NodeVars+
                              faceData->SWMM_Nlinks*faceData->LinkVars +
                              faceData->SysVars)*RECORDSIZE;

  // --- return with file left open
  return faceData;
}


//-----------------------------------------------------------------------------
int GetSwmmResult(IFaceData *faceData, int iType, int iIndex, int vIndex, int period, float* value)
//-----------------------------------------------------------------------------
{
  int offset;

  // --- compute offset into output file
  *value = 0.0;
  offset = faceData->StartPos + (period-1)*faceData->BytesPerPeriod + 2*RECORDSIZE;
  if ( iType == IFACESUBCATCH )
  {
    offset += RECORDSIZE*(iIndex*faceData->SubcatchVars + vIndex);
  }
  else if (iType == IFACENODE)
  {
    offset += RECORDSIZE*(faceData->SWMM_Nsubcatch*faceData->SubcatchVars +
                          iIndex*faceData->NodeVars + vIndex);
  }
  else if (iType == IFACELINK)
  {
    offset += RECORDSIZE*(faceData->SWMM_Nsubcatch*faceData->SubcatchVars +
                          faceData->SWMM_Nnodes*faceData->NodeVars +
                          iIndex*faceData->LinkVars + vIndex);
  }
  else if (iType == IFACESYS)
  {
    offset += RECORDSIZE*(faceData->SWMM_Nsubcatch*faceData->SubcatchVars +
                          faceData->SWMM_Nnodes*faceData->NodeVars +
                          faceData->SWMM_Nlinks*faceData->LinkVars + vIndex);
  }
  else return 0;

  // --- re-position the file and read the result
  fseek(faceData->Fout, offset, SEEK_SET);
  fread(value, RECORDSIZE, 1, faceData->Fout);
  return 1;
}


//-----------------------------------------------------------------------------
void CloseSwmmOutFile(IFaceData *faceData)
//-----------------------------------------------------------------------------
{


  if (faceData->Fout != NULL)
  {
    fclose(faceData->Fout);
    faceData->Fout = NULL;
  }

  if(faceData)
  {
    free(faceData);
    faceData = NULL;
  }
}

