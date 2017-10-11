//-----------------------------------------------------------------------------
//   swmm5.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.001)
//             03/19/15  (Build 5.1.008)
//   Author:   L. Rossman
//
//   This is the main module of the computational engine for Version 5 of
//   the U.S. Environmental Protection Agency's Storm Water Management Model
//   (SWMM). It contains functions that control the flow of computations.
//
//   Depending on how it is compiled, this engine can be executed either as
//   a command line executable or through a series of calls made to functions
//   in a dynamic link library.
//
//
//   Build 5.1.008:
//   - Support added for the MinGW compiler.
//   - Reporting of project options moved to swmm_start. 
//   - Hot start file now read before routing system opened.
//   - Final routing step adjusted so that total duration not exceeded.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

//**********************************************************
//  Leave only one of the following 3 lines un-commented,
//  depending on the choice of compilation target
//**********************************************************
//#define CLE     /* Compile as a command line executable */
//#define SOL     /* Compile as a shared object library */
#define DLL     /* Compile as a Windows DLL */

// --- define WINDOWS
#undef WINDOWS
#ifdef _WIN32
#define WINDOWS
#endif
#ifdef __WIN32__
#define WINDOWS
#endif

////  ---- following section modified for release 5.1.008.  ////               //(5.1.008)
////
// --- define EXH (MS Windows exception handling)
#undef MINGW       // indicates if MinGW compiler used
#undef EXH         // indicates if exception handling included
#ifdef WINDOWS
#ifndef MINGW
#define EXH
#endif
#endif

// --- include Windows & exception handling headers
#ifdef WINDOWS
#include <windows.h>
#endif
#ifdef EXH
#include <excpt.h>
#endif
////

#ifdef _WINDOWS
#include <direct.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

//-----------------------------------------------------------------------------
//  SWMM's header files
//
//  Note: the directives listed below are also contained in headers.h which
//        is included at the start of most of SWMM's other code modules.
//-----------------------------------------------------------------------------
#include "consts.h"                    // defined constants
#include "macros.h"                    // macros used throughout SWMM
#include "enums.h"                     // enumerated variables
#include "error.h"                     // error message codes
#include "datetime.h"                  // date/time functions
#include "objects.h"                   // definitions of SWMM's data objects
#include "funcs.h"                     // declaration of all global functions
#include "text.h"                      // listing of all text strings
#define  EXTERN                        // defined as 'extern' in headers.h
#include "globals.h"                   // declaration of all global variables

#include "swmm5.h"                     // declaration of exportable functions
//   callable from other programs
#define  MAX_EXCEPTIONS 100            // max. number of exceptions handled

//-----------------------------------------------------------------------------
//  Unit conversion factors
//-----------------------------------------------------------------------------
const double Ucf[10][2] =
{//  US      SI
 { 43200.0, 1097280.0 },         // RAINFALL (in/hr, mm/hr --> ft/sec)
 { 12.0, 304.8 },         // RAINDEPTH (in, mm --> ft)
 { 1036800.0, 26334720.0 },         // EVAPRATE (in/day, mm/day --> ft/sec)
 { 1.0, 0.3048 },         // LENGTH (ft, m --> ft)
 { 2.2956e-5, 0.92903e-5 },         // LANDAREA (ac, ha --> ft2)
 { 1.0, 0.02832 },         // VOLUME (ft3, m3 --> ft3)
 { 1.0, 1.608 },         // WINDSPEED (mph, km/hr --> mph)
 { 1.0, 1.8 },         // TEMPERATURE (deg F, deg C --> deg F)
 { 2.203e-6, 1.0e-6 },         // MASS (lb, kg --> mg)
 { 43560.0, 3048.0 }          // GWFLOW (cfs/ac, cms/ha --> ft/sec)
};
const double Qcf[6] =                  // Flow Conversion Factors:
{ 1.0, 448.831, 0.64632,     // cfs, gpm, mgd --> cfs
      0.02832, 28.317, 2.4466};    // cms, lps, mld --> cfs

////-----------------------------------------------------------------------------
////  Shared variables
////-----------------------------------------------------------------------------
//static int  project->IsOpenFlag;           // TRUE if a project has been opened
//static int  project->IsStartedFlag;        // TRUE if a simulation has been started
//static int  project->SaveResultsFlag;      // TRUE if output to be saved to binary file
//static int  project->ExceptionCount;       // number of exceptions handled
//static int  project->DoRunoff;             // TRUE if runoff is computed
//static int  project->DoRouting;            // TRUE if flow routing is computed

//-----------------------------------------------------------------------------
//  External functions (prototyped in swmm5.h)
//-----------------------------------------------------------------------------
//  swmm_run
//  swmm_open
//  swmm_start
//  swmm_step
//  swmm_end
//  swmm_report
//  swmm_close
//  swmm_getMassBalErr
//  swmm_getVersion

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void execRouting(Project* project, DateTime elapsedTime);

// Exception filtering function
#ifdef WINDOWS
static int  xfilter(Project* project, int xc, DateTime elapsedTime, long step);
#endif

//-----------------------------------------------------------------------------
//  Entry point used to compile a stand-alone executable.
//-----------------------------------------------------------------------------
#ifdef CLE 
int  main(int argc, char *argv[])
//
//  Input:   argc = number of command line arguments
//           argv = array of command line arguments
//  Output:  returns error status
//  Purpose: processes command line arguments.
//
//  Command line for stand-alone operation is: swmm5 f1  f2  f3
//  where f1 = name of input file, f2 = name of report file, and
//  f3 = name of binary output file if saved (or blank if not saved).
//
{
   char *inputFile;
   char *reportFile;
   char *binaryFile;
   char blank[] = "";
   time_t start;
   double runTime;

   // --- initialize flags
   Project* project = NULL;

   //project->IsOpenFlag = FALSE;
   //project->IsStartedFlag = FALSE;
   //project->SaveResultsFlag = TRUE;

   // --- check for proper number of command line arguments
   start = time(0);
   if (argc < 3) writecon(FMT01);
   else
   {
      // --- extract file names from command line arguments
      inputFile = argv[1];
      reportFile = argv[2];
      if (argc > 3) binaryFile = argv[3];
      else          binaryFile = blank;
      writecon(FMT02);

      // --- run SWMM
      //swmm_run(inputFile, reportFile, binaryFile);

      long newHour, oldHour = 0;
      long theDay, theHour;
      DateTime elapsedTime = 0.0;

      // --- open the files & read input data
      // int project->ErrorCode = 0;
      project = swmm_open(inputFile, reportFile, binaryFile);

      project->SaveResultsFlag = TRUE;

      // --- run the simulation if input data OK
      if (!project->ErrorCode)
      {
         // --- initialize values
         swmm_start(project, TRUE);

         // --- execute each time step until elapsed time is re-set to 0
         if (!project->ErrorCode)
         {
            writecon("\n Simulating day: 0     hour:  0");
            do
            {
               swmm_step(project, &elapsedTime);
               newHour = (long)(elapsedTime * 24.0);
               if (newHour > oldHour)
               {
                  theDay = (long)elapsedTime;
                  theHour = (long)((elapsedTime - floor(elapsedTime)) * 24.0);
                  writecon("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                  sprintf(project->Msg, "%-5d hour: %-2d", theDay, theHour);
                  writecon(project->Msg);
                  oldHour = newHour;
               }
            } while (elapsedTime > 0.0 && !project->ErrorCode);
            writecon("\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                     "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            writecon("Simulation complete           ");
         }

         // --- clean up
         swmm_end(project);
      }

      // --- report results
      if (project->Fout.mode == SCRATCH_FILE) swmm_report(project);

      // --- close the system
      //swmm_close(project);
      //return project->ErrorCode;

      // Display closing status on console
      runTime = difftime(time(0), start);
      sprintf(project->Msg, "\n\n... EPA-SWMM completed in %.2f seconds.", runTime);
      writecon(project->Msg);
      if (project->ErrorCode) writecon(FMT03);
      else if (project->WarningCode) writecon(FMT04);
      else                    writecon(FMT05);
   }

   // --- Use the code below if you need to keep the console window visible
   /*
                writecon("    Press Enter to continue...");
                getchar();
                */
   if (project)
      swmm_close(project);


   return 0;
}                                      /* End of main */
#endif

//=============================================================================

int  swmm_run(char* f1, char* f2, char* f3)
//
//  Input:   f1 = name of input file
//           f2 = name of report file
//           f3 = name of binary output file
//  Output:  returns error code
//  Purpose: runs a SWMM simulation.
//
{
   long newHour, oldHour = 0;
   long theDay, theHour;
   DateTime elapsedTime = 0.0;

   // --- open the files & read input data
   // int project->ErrorCode = 0;
   Project* project = swmm_open(f1, f2, f3);

   // --- run the simulation if input data OK
   if (!project->ErrorCode)
   {
      // --- initialize values
      swmm_start(project, TRUE);

      // --- execute each time step until elapsed time is re-set to 0
      if (!project->ErrorCode)
      {
         writecon("\n o  Simulating day: 0     hour:  0");
         do
         {
            swmm_step(project, &elapsedTime);
            newHour = (long)(elapsedTime * 24.0);
            if (newHour > oldHour)
            {
               theDay = (long)elapsedTime;
               theHour = (long)((elapsedTime - floor(elapsedTime)) * 24.0);
               writecon("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
               sprintf(project->Msg, "%-5d hour: %-2d", theDay, theHour);
               writecon(project->Msg);
               oldHour = newHour;
            }
         } while (elapsedTime > 0.0 && !project->ErrorCode);
         writecon("\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                  "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
         writecon("Simulation complete           ");
      }

      // --- clean up
      swmm_end(project);
   }

   // --- report results
   if (project->Fout.mode == SCRATCH_FILE) swmm_report(project);

   // --- close the system
   //swmm_close(project);
   //return project->ErrorCode;
   return swmm_close(project);
}

//=============================================================================

Project* swmm_open(char* f1, char* f2, char* f3)
//
//  Input:   f1 = name of input file
//           f2 = name of report file
//           f3 = name of binary output file
//  Output:  returns error code
//  Purpose: opens a SWMM project.
//
{
#ifdef DLL
#ifdef _WINDOWS
   _fpreset();
#endif
#endif

   Project *project = calloc(1,sizeof(Project));

#ifdef WINDOWS
   // --- begin exception handling here
   __try
      #endif
   {
      // --- initialize error & warning codes
      datetime_setDateFormat(M_D_Y);

      project->ErrorCode = 0;
      project->WarningCode = 0;
      project->IsOpenFlag = FALSE;
      project->IsStartedFlag = FALSE;
      project->ExceptionCount = 0;

      // --- open a SWMM project
      project_open(project, f1, f2, f3);

      if (project->ErrorCode)
         return project;

      project->IsOpenFlag = TRUE;
      report_writeLogo(project);
      writecon(FMT06);

      // --- retrieve project data from input file
      project_readInput(project);
      if (project->ErrorCode)
         return project;

      // --- write project title to report file & validate data
      report_writeTitle(project);
      project_validate(project);

      // --- write input summary to report file if requested
      if (project->RptFlags.input) inputrpt_writeInput(project);
   }
#ifdef WINDOWS
   // --- end of try loop; handle exception here
   __except (xfilter(project, GetExceptionCode(), 0.0, 0))
   {
      project->ErrorCode = ERR_SYSTEM;
   }
#endif

   return project;
}

//=============================================================================

int swmm_start(Project* project, int saveResults)
//
//  Input:   saveResults = TRUE if simulation results saved to binary file 
//  Output:  returns an error code
//  Purpose: starts a SWMM simulation.
//
{
   // --- check that a project is open & no run started
   if (project->ErrorCode) return project->ErrorCode;
   if (!project->IsOpenFlag || project->IsStartedFlag)
   {
      report_writeErrorMsg(project, ERR_NOT_OPEN, "");
      return project->ErrorCode;
   }
   project->ExceptionCount = 0;

#ifdef WINDOWS
   // --- begin exception handling loop here
   __try
      #endif
   {
      // --- initialize runoff, routing & reporting time (in milliseconds)
      project->NewRunoffTime = 0.0;
      project->NewRoutingTime = 0.0;
      project->ReportTime = (double)(1000 * project->ReportStep);
      project->StepCount = 0;
      project->NonConvergeCount = 0;
      project->IsStartedFlag = TRUE;

      // --- initialize global continuity errors
      project->RunoffError = 0.0;
      project->GwaterError = 0.0;
      project->FlowError = 0.0;
      project->QualError = 0.0;

      // --- open rainfall processor (creates/opens a rainfall
      //     interface file and generates any RDII flows)
      if (!project->IgnoreRainfall) rain_open(project);
      if (project->ErrorCode) return project->ErrorCode;

      // --- initialize state of each major system component
      project_init(project);

      // --- see if runoff & routing needs to be computed
      if (project->Nobjects[SUBCATCH] > 0) project->DoRunoff = TRUE;
      else project->DoRunoff = FALSE;
      if (project->Nobjects[NODE] > 0 && !project->IgnoreRouting) project->DoRouting = TRUE;
      else project->DoRouting = FALSE;

      ////  Following section modified for release 5.1.008.  ////                    //(5.1.008)
      ////
      // --- open binary output file
      output_open(project);

      // --- open runoff processor
      if (project->DoRunoff)
         runoff_open(project);

      // --- open & read hot start file if present
      if (!hotstart_open(project)) return project->ErrorCode;

      // --- open routing processor
      if (project->DoRouting) routing_open(project);

      // --- open mass balance and statistics processors
      massbal_open(project);
      stats_open(project);

      // --- write project options to report file
      report_writeOptions(project);
      if (project->RptFlags.controls)
         report_writeControlActionsHeading(project);
      ////
   }

#ifdef WINDOWS
   // --- end of try loop; handle exception here
   __except (xfilter(project, GetExceptionCode(), 0.0, 0))
   {
      project->ErrorCode = ERR_SYSTEM;
   }
#endif

   // --- save saveResults flag to global variable
   project->SaveResultsFlag = saveResults;
   return project->ErrorCode;
}
//=============================================================================

int  swmm_step(Project* project, DateTime* elapsedTime)
//
//  Input:   elapsedTime = current elapsed time in decimal days
//  Output:  updated value of elapsedTime,
//           returns error code
//  Purpose: advances the simulation by one routing time step.
//
{
   // --- check that simulation can proceed
   if (project->ErrorCode)
     return project->ErrorCode;

   if (!project->IsOpenFlag || !project->IsStartedFlag)
   {
      report_writeErrorMsg(project, ERR_NOT_OPEN, "");
      return project->ErrorCode;
   }

#ifdef WINDOWS
   // --- begin exception handling loop here
   __try
      #endif
   {
      // --- if routing time has not exceeded total duration
      if (project->NewRoutingTime < project->TotalDuration)
      {
         // --- route flow & WQ through drainage system
         //     (runoff will be calculated as needed)
         //     (project->NewRoutingTime is updated)
         execRouting(project, *elapsedTime);
      }

      // --- save results at next reporting time
      if (project->NewRoutingTime != project->PreviousReportTime  &&
          project->NewRoutingTime >= project->ReportTime)
      {
         if (project->SaveResultsFlag)
           output_saveResults(project, project->ReportTime);

         project->PreviousReportTime = project->NewRoutingTime;
         project->ReportTime = project->ReportTime + (double)(1000 * project->ReportStep);
      }

      // --- update elapsed time (days)
      if (project->NewRoutingTime < project->TotalDuration)
      {
         *elapsedTime = project->NewRoutingTime / MSECperDAY;
      }

      // --- otherwise end the simulation
      else *elapsedTime = 0.0;
   }

#ifdef WINDOWS
   // --- end of try loop; handle exception here
   __except (xfilter(project, GetExceptionCode(), *elapsedTime, project->StepCount))
   {
      project->ErrorCode = ERR_SYSTEM;
   }
#endif
   return project->ErrorCode;
}

//=============================================================================

void execRouting(Project* project, DateTime elapsedTime)
//
//  Input:   elapsedTime = current elapsed time in decimal days
//  Output:  none
//  Purpose: routes flow & WQ through drainage system over a single time step.
//
{
   double   nextRoutingTime;          // updated elapsed routing time (msec)
   double   routingStep;              // routing time step (sec)

#ifdef WINDOWS
   // --- begin exception handling loop here
   __try
      #endif
   {
      // --- determine when next routing time occurs
      project->StepCount++;

      if (!project->DoRouting)
         routingStep = MIN(project->WetStep, project->ReportStep);
      else
         routingStep = routing_getRoutingStep(project, project->RouteModel, project->RouteStep);

      if (routingStep <= 0.0)
      {
         project->ErrorCode = ERR_TIMESTEP;
         return;
      }

      nextRoutingTime = project->NewRoutingTime + 1000.0 * routingStep;

      ////  Following section added to release 5.1.008.  ////                        //(5.1.008)
      ////
      // --- adjust routing step so that total duration not exceeded
      if (nextRoutingTime > project->TotalDuration)
      {
         routingStep = (project->TotalDuration - project->NewRoutingTime) / 1000.0;
         routingStep = MAX(routingStep, 1. / 1000.0);
         nextRoutingTime = project->TotalDuration;
      }
      ////

      // --- compute runoff until next routing time reached or exceeded
      if (project->DoRunoff)
      {
         while (project->NewRunoffTime < nextRoutingTime)
         {
            runoff_execute(project);

            if (project->ErrorCode)
              return;
         }
      }

      // --- if no runoff analysis, update climate state (for evaporation)
      else
      {
         climate_setState(project, getDateTime(project, project->NewRoutingTime));
      }

      // --- route flows & pollutants through drainage system                //(5.1.008)
      //     (while updating project->NewRoutingTime)                        //(5.1.008)
      if (project->DoRouting)
        routing_execute(project, project->RouteModel, routingStep);
      else
        project->NewRoutingTime = nextRoutingTime;
   }

#ifdef WINDOWS
   // --- end of try loop; handle exception here
   __except (xfilter(project, GetExceptionCode(), elapsedTime, project->StepCount))
   {
      project->ErrorCode = ERR_SYSTEM;
      return;
   }
#endif
}

//=============================================================================

int  swmm_end(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: ends a SWMM simulation.
//
{
   // --- check that project opened and run started
   if (!project->IsOpenFlag)
   {
      report_writeErrorMsg(project, ERR_NOT_OPEN, "");
      return project->ErrorCode;
   }

   if (project->IsStartedFlag)
   {
      // --- write ending records to binary output file
      if (project->Fout.file) output_end(project);

      // --- report mass balance results and system statistics
      if (!project->ErrorCode)
      {
         massbal_report(project);
         stats_report(project);
      }

      // --- close all computing systems
      stats_close(project);
      massbal_close(project);
      if (!project->IgnoreRainfall) rain_close(project);
      if (project->DoRunoff) runoff_close(project);
      if (project->DoRouting) routing_close(project, project->RouteModel);
      hotstart_close(project);
      project->IsStartedFlag = FALSE;
   }
   return project->ErrorCode;
}

//=============================================================================

int  swmm_report(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: writes simulation results to report file.
//
{
   if (project->Fout.mode == SCRATCH_FILE) output_checkFileSize(project);
   if (project->ErrorCode) report_writeErrorCode(project);
   else
   {
      writecon(FMT07);
      report_writeReport(project);
   }
   return project->ErrorCode;
}

//=============================================================================

int  swmm_close(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: closes a SWMM project.
//
{
   if (project->Fout.file) output_close(project);
   if (project->IsOpenFlag) project_close(project);
   report_writeSysTime(project);
   if (project->Finp.file != NULL) fclose(project->Finp.file);
   if (project->Frpt.file != NULL) fclose(project->Frpt.file);
   if (project->Fout.file != NULL)
   {
      fclose(project->Fout.file);
      if (project->Fout.mode == SCRATCH_FILE) remove(project->Fout.name);
   }
   project->IsOpenFlag = FALSE;
   project->IsStartedFlag = FALSE;

   FREE(project);

   return 0;
}

//=============================================================================

int   swmm_getMassBalErr(Project* project, float* runoffErr, float* flowErr,
                                  float* qualErr)
//
//  Input:   none
//  Output:  runoffErr = runoff mass balance error (percent)
//           flowErr   = flow routing mass balance error (percent)
//           qualErr   = quality routing mass balance error (percent)
//           returns an error code
//  Purpose: reports a simulation's mass balance errors.
//
{
   *runoffErr = 0.0;
   *flowErr = 0.0;
   *qualErr = 0.0;

   if (project->IsOpenFlag && !project->IsStartedFlag)
   {
      *runoffErr = (float)project->RunoffError;
      *flowErr = (float)project->FlowError;
      *qualErr = (float)project->QualError;
   }
   return 0;
}

//=============================================================================

int   swmm_getVersion(void)
//
//  Input:   none
//  Output:  returns SWMM engine version number
//  Purpose: retrieves version number of current SWMM engine which
//           uses a format of xyzzz where x = major version number,
//           y = minor version number, and zzz = build number.
//
{
   return VERSION;
}

//=============================================================================
// Additional code added to enable couple
//=============================================================================

int   swmm_getErrorCode(Project* project)
//
//  Input:   none
//  Output:  returns SWMM engine version number
//  Purpose: retrieves version number of current SWMM engine which
//           uses a format of xyzzz where x = major version number,
//           y = minor version number, and zzz = build number.
//
{
   return project->ErrorCode;
}

double  swmm_getDateTime(Project* project, char* beginorend)
{
   if (strcomp(beginorend, "begin"))
   {
      return project->StartDateTime;
   }
   else
   {
      return project->EndDateTime;
   }
}

void  datetime_decodeDateTime(DateTime date, int* y, int* m, int* d, int* h, int* mm, int* s)
{
   datetime_decodeDate(date, y, m, d);
   datetime_decodeTime(date, h, mm, s);
}

char* getErrorMsg(int i)
{
   return error_getMsg(i);
}

int  getObjectTypeCount(Project* project, int type)
{
   return project->Nobjects[type];
}

TNode* getNode(Project* project, int index)
{
   return &project->Node[index];
}

TNode* getNodeById(Project* project, char* id)
{
   int index = project_findObject(project, NODE, id);
   return &project->Node[index];
}

void  setNode(Project* project, char* nodeId, char* propertyName, double value)
{
   int index = project_findObject(project, NODE, nodeId);
   TNode* n1 = &project->Node[index];

   if (strcomp(propertyName, "invertElev"))
   {
      n1->invertElev = value;
   }
   else if (strcomp(propertyName, "crownElev"))
   {
      n1->crownElev = value;
   }
   else if (strcomp(propertyName, "initDepth"))
   {
      n1->initDepth = value;
   }
   else if (strcomp(propertyName, "surDepth"))
   {
      n1->surDepth = value;
   }
   else if (strcomp(propertyName, "newDepth"))
   {
      //to be applied later
      addNodeDepth(project, index, value);
   }
   else if (strcomp(propertyName, "pondedArea"))
   {
      n1->pondedArea = value;
   }
   else if (strcomp(propertyName, "inflow"))
   {
      n1->inflow = value;
   }
   else if (strcomp(propertyName, "outflow"))
   {
      n1->outflow = value;
   }
   else if (strcomp(propertyName, "newLatFlow"))
   {
      //to be applied later
      addNodeLateralInflow(project, index, value);
   }

}

TLink*  getLink(Project* project, int index)
{
   return &project->Link[index];
}

TLink*  getLinkById(Project* project, char* id)
{
   int index = project_findObject(project, LINK, id);
   return &project->Link[index];
}

void  setLink(Project* project, char* linkId, char* propertyName, double value)
{

   TLink* l = &project->Link[project_findObject(project, LINK, linkId)];

   if (strcomp(propertyName, "offset1"))
   {
      l->offset1 = value;
   }
   else if (strcomp(propertyName, "offset2"))
   {
      l->offset2 = value;
   }
   else if (strcomp(propertyName, "q0"))
   {
      l->q0 = value;
   }
   else if (strcomp(propertyName, "cLossInlet"))
   {
      l->cLossInlet = value;
   }
   else if (strcomp(propertyName, "cLossOutlet"))
   {
      l->cLossOutlet = value;
   }
   else if (strcomp(propertyName, "cLossAvg"))
   {
      l->cLossAvg = value;
   }
   else if (strcomp(propertyName, "seepRate"))
   {
      l->seepRate = value;
   }
   else if (strcomp(propertyName, "newFlow"))
   {
      l->newFlow = value;
   }

}

TSubcatch* getSubcatch(Project* project, int index)
{
   return &project->Subcatch[index];
}

TSubcatch* getSubcatchById(Project* project, char* id)
{
   int index = project_findObject(project, SUBCATCH, id);
   return &project->Subcatch[index];
}

void  setSubcatch(Project* project, char* subCatchId, char* propertyName, double value)
{
   int index = project_findObject(project, SUBCATCH, subCatchId);
   TSubcatch* subcatch = &project->Subcatch[index];

   if (strcomp(propertyName, "newRunoff"))
   {
      subcatch->newRunoff = value;
   }
   else if (strcomp(propertyName, "rainfall"))
   {
      addSubcatchRain(project, index, value);
   }
}

//=============================================================================

//=============================================================================
//   General purpose functions
//=============================================================================

double UCF(Project* project, int u)
//
//  Input:   u = integer code of quantity being converted
//  Output:  returns a units conversion factor
//  Purpose: computes a conversion factor from SWMM's internal
//           units to user's units
//
{
   if (u < FLOW)
      return Ucf[u][project->UnitSystem];
   else
      return Qcf[project->FlowUnits];
}

//=============================================================================

char* sstrncpy(char *dest, const char *src, size_t maxlen)
//
//  Input:   dest = string to be copied to
//           src = string to be copied from
//           maxlen = number of characters to copy
//  Output:  returns a pointer to dest
//  Purpose: safe version of standard strncpy function
//
{
   strncpy(dest, src, maxlen);
   dest[maxlen] = '\0';
   return dest;
}

//=============================================================================

int  strcomp(char *s1, char *s2)
//
//  Input:   s1 = a character string
//           s2 = a character string
//  Output:  returns 1 if s1 is same as s2, 0 otherwise
//  Purpose: does a case insensitive comparison of two strings.
//
{
   int i;
   for (i = 0; UCHAR(s1[i]) == UCHAR(s2[i]); i++)
   {
      if (!s1[i + 1] && !s2[i + 1]) return(1);
   }
   return(0);
}

//=============================================================================

char* getTempFileName(Project* project, char* fname)
//
//  Input:   fname = file name string (with max size of MAXFNAME)
//  Output:  returns pointer to file name
//  Purpose: creates a temporary file name with path prepended to it.
//
{
   // For Windows systems:
#ifdef WINDOWS

   char* name = NULL;
   char* dir = NULL;

   // --- set dir to user's choice of a temporary directory
   if (strlen(project->TempDir) > 0)
   {
      _mkdir(project->TempDir);
      dir = project->TempDir;
   }

   // --- use _tempnam to get a pointer to an unused file name
   name = _tempnam(dir, "swmm");
   if (name == NULL) return NULL;

   // --- copy the file name to fname
   if (strlen(name) < MAXFNAME) strncpy(fname, name, MAXFNAME);
   else fname = NULL;

   // --- free the pointer returned by _tempnam
   free(name);

   // --- return the new contents of fname
   return fname;

   // For non-Windows systems:
#else

   // --- use system function mkstemp() to create a temporary file name
   strcpy(fname, "swmmXXXXXX");
   mkstemp(fname);
   return fname;

#endif
}

//=============================================================================

void getElapsedTime(Project* project, DateTime aDate, int* days, int* hrs, int* mins)
//
//  Input:   aDate = simulation calendar date + time
//  Output:  days, hrs, mins = elapsed days, hours & minutes for aDate
//  Purpose: finds elapsed simulation time for a given calendar date
//
{
   DateTime x;
   int secs;
   x = aDate - project->StartDateTime;
   if (x <= 0.0)
   {
      *days = 0;
      *hrs = 0;
      *mins = 0;
   }
   else
   {
      *days = (int)x;
      datetime_decodeTime(x, hrs, mins, &secs);
   }
}

//=============================================================================

DateTime getDateTime(Project* project, double elapsedMsec)
//
//  Input:   elapsedMsec = elapsed milliseconds
//  Output:  returns date/time value
//  Purpose: finds calendar date/time value for elapsed milliseconds of
//           simulation time.
//
{
   return datetime_addSeconds(project->StartDateTime, (elapsedMsec + 1.0) / 1000.0);
}

//=============================================================================

void  writecon(char *s)
//
//  Input:   s = a character string
//  Output:  none
//  Purpose: writes string of characters to the console.
//
{
#ifdef CLE 
   fprintf(stdout, s);
   fflush(stdout);
#endif
}

//=============================================================================

#ifdef WINDOWS
int xfilter(Project* project, int xc, DateTime elapsedTime, long step)
//
//  Input:   xc          = exception code
//           elapsedTime = simulation time when exception occurred (days)
//           step        = step count at time when exception occurred
//  Output:  returns an exception handling code
//  Purpose: exception filtering routine for operating system exceptions
//           under Windows.
//
{
   int  rc;                           // result code
   long hour;                         // current hour of simulation
   char msg[40];                      // exception type text
   char xmsg[120];                    // error message text
   switch (xc)
   {
      case EXCEPTION_ACCESS_VIOLATION:
         sprintf(msg, "\n  Access violation ");
         rc = EXCEPTION_EXECUTE_HANDLER;
         break;
      case EXCEPTION_FLT_DENORMAL_OPERAND:
         sprintf(msg, "\n  Illegal floating point operand ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_FLT_DIVIDE_BY_ZERO:
         sprintf(msg, "\n  Floating point divide by zero ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_FLT_INVALID_OPERATION:
         sprintf(msg, "\n  Illegal floating point operation ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_FLT_OVERFLOW:
         sprintf(msg, "\n  Floating point overflow ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_FLT_STACK_CHECK:
         sprintf(msg, "\n  Floating point stack violation ");
         rc = EXCEPTION_EXECUTE_HANDLER;
         break;
      case EXCEPTION_FLT_UNDERFLOW:
         sprintf(msg, "\n  Floating point underflow ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_INT_DIVIDE_BY_ZERO:
         sprintf(msg, "\n  Integer divide by zero ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      case EXCEPTION_INT_OVERFLOW:
         sprintf(msg, "\n  Integer overflow ");
         rc = EXCEPTION_CONTINUE_EXECUTION;
         break;
      default:
         sprintf(msg, "\n  Exception %d", xc);
         rc = EXCEPTION_EXECUTE_HANDLER;
   }
   hour = (long)(elapsedTime / 1000.0 / 3600.0);
   sprintf(xmsg, "%s at step %d, hour %d", msg, step, hour);
   if (rc == EXCEPTION_EXECUTE_HANDLER ||
       ++project->ExceptionCount >= MAX_EXCEPTIONS)
   {
      strcat(xmsg, " --- execution halted.");
      rc = EXCEPTION_EXECUTE_HANDLER;
   }
   report_writeLine(project, xmsg);
   return rc;
}
#endif

//=============================================================================
