//-----------------------------------------------------------------------------
//   project.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.000)
//             04/14/14  (Build 5.1.004)
//             09/15/14  (Build 5.1.007)
//             03/19/15  (Build 5.1.008)
//             04/30/15  (Build 5.1.009)
//   Author:   L. Rossman
//
//   Project management functions.
//
//   This module provides project-related services such as:
//   o opening a new project and reading its input data
//   o allocating and freeing memory for project objects
//   o setting default values for object properties and options
//   o initializing the internal state of all objects
//   o managing hash tables for identifying objects by ID name
//
//   Build 5.1.004:
//   - Ignore RDII option added.
//
//   Build 5.1.007:
//   - Default monthly adjustments for climate variables included.
//   - User-supplied GW flow equaitions initialized to NULL.
//   - project->Storage node exfiltration object initialized to NULL.
//   - Freeing of memory used for storage node exfiltration included.
//
//   Build 5.1.008:
//   - Constants used for dynamic wave routing moved to dynwave.c.
//   - Input processing of minimum time step & number of
//     parallel threads for dynamic wave routing added.
//   - Default values of hyd. conductivity adjustments added.
//   - Freeing of memory used for outfall pollutant load added.
//
//   Build 5.1.009:
//   - Fixed bug in computing total duration introduced in 5.1.008.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
//#include <malloc.h>
#include <math.h>                                                              //(5.1.008)

#ifdef USE_OPENMP
#include <omp.h>                                                               //(5.1.008)
#endif

#include "headers.h"
#include "lid.h"
#include "hash.h"
#include "mempool.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
////  Constants for DYNWAVE flow routing moved to dynwave.c.  ////             //(5.1.008)

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
//static HTtable* project->Htable[MAX_OBJ_TYPES]; // Hash tables for object ID names
//static char     project->MemPoolAllocated;      // TRUE if memory pool allocated 

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  project_open           (called from swmm_open in swmm5.c)
//  project_close          (called from swmm_close in swmm5.c)
//  project_readInput      (called from swmm_open in swmm5.c)
//  project_readOption     (called from readOption in input.c)
//  project_validate       (called from swmm_open in swmm5.c)
//  project_init           (called from swmm_start in swmm5.c)
//  project_addObject      (called from addObject in input.c)
//  project_createMatrix   (called from openFileForInput in iface.c)
//  project_freeMatrix     (called from iface_closeRoutingFiles)
//  project_findObject
//  project_findID

//-----------------------------------------------------------------------------
//  Function declarations
//-----------------------------------------------------------------------------
static void initPointers(Project* project);
static void setDefaults(Project* project);
static void openFiles(Project* project, char *f1, char *f2, char *f3);
static void createObjects(Project* project);
static void deleteObjects(Project* project);
static void createHashTables(Project* project);
static void deleteHashTables(Project* project);


//=============================================================================

void project_open(Project* project, char *f1, char *f2, char *f3)
//
//  Input:   f1 = pointer to name of input file
//           f2 = pointer to name of report file
//           f3 = pointer to name of binary output file
//  Output:  none
//  Purpose: opens a new SWMM project.
//
{

   initPointers(project);
   setDefaults(project);
   openFiles(project,f1, f2, f3);
}

//=============================================================================

void project_readInput(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: retrieves project data from input file.
//
{
   // --- create hash tables for fast retrieval of objects by ID names
   createHashTables(project);

   // --- count number of objects in input file and create them
   input_countObjects(project);
   createObjects(project);

   // --- read project data from input file
   input_readData(project);

   if ( project->ErrorCode ) return;

   // --- establish starting & ending date/time
   project->StartDateTime = project->StartDate + project->StartTime;
   project->EndDateTime   = project->EndDate + project->EndTime;
   project->ReportStart   = project->ReportStartDate + project->ReportStartTime;
   project->ReportStart   = MAX(project->ReportStart, project->StartDateTime);
   project->PreviousReportTime = -10000000000000000000;

   // --- check for valid starting & ending date/times
   if ( project->EndDateTime <= project->StartDateTime )
   {
      report_writeErrorMsg(project,ERR_START_DATE, "");
   }
   else if ( project->EndDateTime <= project->ReportStart )
   {
      report_writeErrorMsg(project,ERR_REPORT_DATE, "");
   }
   else
   {
      ////  Following code segment was modified for release 5.1.009.  ////           //(5.1.009)
      ////
      // --- compute total duration of simulation in seconds
      project->TotalDuration = floor((project->EndDateTime - project->StartDateTime) * SECperDAY);

      // --- reporting step must be <= total duration
      if ( (double)project->ReportStep > project->TotalDuration )
      {
         project->ReportStep = (int)(project->TotalDuration);
      }

      // --- reporting step can't be < routing step
      if ( (double)project->ReportStep < project->RouteStep )
      {
         report_writeErrorMsg(project,ERR_REPORT_STEP, "");
      }

      // --- convert total duration to milliseconds
      project->TotalDuration *= 1000.0;
   }
   ////
}

//=============================================================================

void project_validate(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: checks validity of project data.
//
{
   int i;
   int j;
   int err;

   // --- validate Curves and TimeSeries
   for ( i=0; i<project->Nobjects[CURVE]; i++ )
   {
      err = table_validate(&project->Curve[i]);
      if ( err ) report_writeErrorMsg(project,ERR_CURVE_SEQUENCE, project->Curve[i].ID);
   }
   for ( i=0; i<project->Nobjects[TSERIES]; i++ )
   {
      err = table_validate(&project->Tseries[i]);
      if ( err ) report_writeTseriesErrorMsg(project,err, &project->Tseries[i]);
   }

   // --- validate hydrology objects
   //     (NOTE: order is important !!!!)
   climate_validate(project);
   lid_validate(project);
   if ( project->Nobjects[SNOWMELT] == 0 ) project->IgnoreSnowmelt = TRUE;
   if ( project->Nobjects[AQUIFER]  == 0 ) project->IgnoreGwater   = TRUE;
   for ( i=0; i<project->Nobjects[GAGE]; i++ )     gage_validate(project,i);
   for (i = 0; i<project->Nobjects[AQUIFER]; i++)  gwater_validateAquifer(project, i);
   for (i = 0; i<project->Nobjects[SUBCATCH]; i++) subcatch_validate(project, i);
   for (i = 0; i<project->Nobjects[SNOWMELT]; i++) snow_validateSnowmelt(project, i);

   // --- compute geometry tables for each shape curve
   j = 0;
   for ( i=0; i<project->Nobjects[CURVE]; i++ )
   {
      if ( project->Curve[i].curveType == SHAPE_CURVE )
      {
         project->Curve[i].refersTo = j;
         project->Shape[j].curve = i;
         if ( !shape_validate(&project->Shape[j], &project->Curve[i]) )
            report_writeErrorMsg(project,ERR_CURVE_SEQUENCE, project->Curve[i].ID);
         j++;
      }
   }

   // --- validate links before nodes, since the latter can
   //     result in adjustment of node depths
   for ( i=0; i<project->Nobjects[NODE]; i++) project->Node[i].oldDepth = project->Node[i].fullDepth;
   for ( i=0; i<project->Nobjects[LINK]; i++) link_validate(project,i);
   for ( i=0; i<project->Nobjects[NODE]; i++) node_validate(project,i);

   // --- adjust time steps if necessary
   if ( project->DryStep < project->WetStep )
   {
      report_writeWarningMsg(project,WARN06, "");
      project->DryStep = project->WetStep;
   }
   if ( project->RouteStep > (double)project->WetStep )
   {
      report_writeWarningMsg(project,WARN07, "");
      project->RouteStep = project->WetStep;
   }

   // --- adjust individual reporting flags to match global reporting flag
   if ( project->RptFlags.subcatchments == ALL )
      for (i=0; i<project->Nobjects[SUBCATCH]; i++) project->Subcatch[i].rptFlag = TRUE;
   if ( project->RptFlags.nodes == ALL )
      for (i=0; i<project->Nobjects[NODE]; i++) project->Node[i].rptFlag = TRUE;
   if ( project->RptFlags.links == ALL )
      for (i=0; i<project->Nobjects[LINK]; i++) project->Link[i].rptFlag = TRUE;

   // --- validate dynamic wave options
   if ( project->RouteModel == DW ) dynwave_validate(project);                                //(5.1.008)

#ifdef USE_OPENMP

#pragma omp parallel                                                           //(5.1.008)
   {
      if ( project->NumThreads == 0 ) project->NumThreads = omp_get_num_threads();                 //(5.1.008)
      else project->NumThreads = MIN(project->NumThreads, omp_get_num_threads());                  //(5.1.008)
   }
#endif

   if ( project->Nobjects[LINK] < 4 * project->NumThreads ) project->NumThreads = 1;                     //(5.1.008)


}

//=============================================================================

void project_close(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: closes a SWMM project.
//
{
   deleteObjects(project);
   deleteHashTables(project);
   clearDataCache(project);
}

//=============================================================================

int  project_init(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: initializes the internal state of all objects.
// 
{
   int j;
   climate_initState(project);
   lid_initState(project);
   for (j=0; j<project->Nobjects[TSERIES]; j++)  table_tseriesInit(&project->Tseries[j]);
   for (j=0; j<project->Nobjects[GAGE]; j++)     gage_initState(project,j);
   for (j = 0; j<project->Nobjects[SUBCATCH]; j++) subcatch_initState(project, j);
   for (j=0; j<project->Nobjects[NODE]; j++)     node_initState(project,j);
   for (j=0; j<project->Nobjects[LINK]; j++)     link_initState(project,j);
   return project->ErrorCode;
}

//=============================================================================

int   project_addObject(Project* project, int type, char *id, int n)
//
//  Input:   type = object type
//           id   = object ID string
//           n    = object index
//  Output:  returns 0 if object already added, 1 if not, -1 if hashing fails
//  Purpose: adds an object ID to a hash table
//
{
   int  result;
   int  len;
   char *newID;

   // --- do nothing if object already placed in hash table
   if ( project_findObject(project,type, id) >= 0 ) return 0;

   // --- use memory from the hash tables' common memory pool to store
   //     a copy of the object's ID string
   len = strlen(id) + 1;
   newID = (char *) Alloc(project,len*sizeof(char));
   strcpy(newID, id);

   // --- insert object's ID into the hash table for that type of object
   result = HTinsert(project->Htable[type], newID, n);
   if ( result == 0 ) result = -1;
   return result;
}

//=============================================================================

int   project_findObject(Project* project, int type, char *id)
//
//  Input:   type = object type
//           id   = object ID
//  Output:  returns index of object with given ID, or -1 if ID not found
//  Purpose: uses hash table to find index of an object with a given ID.
//
{
   return HTfind(project->Htable[type], id);
}

//=============================================================================

char  *project_findID(Project* project, int type, char *id)
//
//  Input:   type = object type
//           id   = ID name being sought
//  Output:  returns pointer to location where object's ID string is stored
//  Purpose: uses hash table to find address of given string entry.
//
{
   return HTfindKey(project->Htable[type], id);
}

//=============================================================================

double ** project_createMatrix(int nrows, int ncols)
//
//  Input:   nrows = number of rows (0-based)
//           ncols = number of columns (0-based)
//  Output:  returns a pointer to a matrix
//  Purpose: allocates memory for a matrix of doubles.
//
{
   int i,j;
   double **a;

   // --- allocate pointers to rows
   a = (double **) malloc(nrows * sizeof(double *));
   if ( !a ) return NULL;

   // --- allocate rows and set pointers to them
   a[0] = (double *) malloc (nrows * ncols * sizeof(double));
   if ( !a[0] ) return NULL;
   for ( i = 1; i < nrows; i++ ) a[i] = a[i-1] + ncols;

   for ( i = 0; i < nrows; i++)
   {
      for ( j = 0; j < ncols; j++) a[i][j] = 0.0;
   }

   // --- return pointer to array of pointers to rows
   return a;
}

//=============================================================================

void project_freeMatrix(double **a)
//
//  Input:   a = matrix of floats
//  Output:  none
//  Purpose: frees memory allocated for a matrix of doubles.
//
{
   if ( a != NULL )
   {
      if ( a[0] != NULL ) free( a[0] );
      free( a );
   }
}

//=============================================================================

int project_readOption(Project* project, char* s1, char* s2)
//
//  Input:   s1 = option keyword
//           s2 = string representation of option's value
//  Output:  returns error code
//  Purpose: reads a project option from a pair of string tokens.
//
//  NOTE:    all project options have default values assigned in setDefaults().
//
{
   int      k, m, h, s;
   double   tStep;
   char     strDate[25];
   DateTime aTime;
   DateTime aDate;

   // --- determine which option is being read
   k = findmatch(s1, OptionWords);
   if ( k < 0 ) return error_setInpError(ERR_KEYWORD, s1);
   switch ( k )
   {
      // --- choice of flow units
      case FLOW_UNITS:
         m = findmatch(s2, FlowUnitWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         project->FlowUnits = m;
         if ( project->FlowUnits <= MGD ) project->UnitSystem = US;
         else                    project->UnitSystem = SI;
         break;

         // --- choice of infiltration modeling method
      case INFIL_MODEL:
         m = findmatch(s2, InfilModelWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         project->InfilModel = m;
         break;

         // --- choice of flow routing method
      case ROUTE_MODEL:
         m = findmatch(s2, RouteModelWords);
         if ( m < 0 ) m = findmatch(s2, OldRouteModelWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         if ( m == NO_ROUTING ) project->IgnoreRouting = TRUE;
         else project->RouteModel = m;
         if ( project->RouteModel == EKW ) project->RouteModel = KW;
         break;

         // --- simulation start date
      case START_DATE:
         if ( !datetime_strToDate(s2, &project->StartDate) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- simulation start time of day
      case START_TIME:
         if ( !datetime_strToTime(s2, &project->StartTime) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- simulation ending date
      case END_DATE:
         if ( !datetime_strToDate(s2, &project->EndDate) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- simulation ending time of day
      case END_TIME:
         if ( !datetime_strToTime(s2, &project->EndTime) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- reporting start date
      case REPORT_START_DATE:
         if ( !datetime_strToDate(s2, &project->ReportStartDate) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- reporting start time of day
      case REPORT_START_TIME:
         if ( !datetime_strToTime(s2, &project->ReportStartTime) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         break;

         // --- day of year when street sweeping begins or when it ends
         //     (year is arbitrarily set to 1947 so that the dayOfYear
         //      function can be applied)
      case SWEEP_START:
      case SWEEP_END:
         strcpy(strDate, s2);
         strcat(strDate, "/1947");
         if ( !datetime_strToDate(strDate, &aDate) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         m = datetime_dayOfYear(aDate);
         if ( k == SWEEP_START ) project->SweepStart = m;
         else project->SweepEnd = m;
         break;

         // --- number of antecedent dry days
      case START_DRY_DAYS:
         project->StartDryDays = atof(s2);
         if ( project->StartDryDays < 0.0 )
         {
            return error_setInpError(ERR_NUMBER, s2);
         }
         break;

         // --- runoff or reporting time steps
         //     (input is in hrs:min:sec format, time step saved as seconds)
      case WET_STEP:
      case DRY_STEP:
      case REPORT_STEP:
         if ( !datetime_strToTime(s2, &aTime) )
         {
            return error_setInpError(ERR_DATETIME, s2);
         }
         datetime_decodeTime(aTime, &h, &m, &s);
         h += 24*(int)aTime;
         s = s + 60*m + 3600*h;
         if ( s <= 0 ) return error_setInpError(ERR_NUMBER, s2);
         switch ( k )
         {
            case WET_STEP:     project->WetStep = s;     break;
            case DRY_STEP:     project->DryStep = s;     break;
            case REPORT_STEP:  project->ReportStep = s;  break;
         }
         break;

         // --- type of damping applied to inertial terms of dynamic wave routing
      case INERT_DAMPING:
         m = findmatch(s2, InertDampingWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         else project->InertDamping = m;
         break;

         // --- Yes/No options (NO = 0, YES = 1)
      case ALLOW_PONDING:
      case SLOPE_WEIGHTING:
      case SKIP_STEADY_STATE:
      case IGNORE_RAINFALL:
      case IGNORE_SNOWMELT:
      case IGNORE_GWATER:
      case IGNORE_ROUTING:
      case IGNORE_QUALITY:
      case IGNORE_RDII:                                                        //(5.1.004)
         m = findmatch(s2, NoYesWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         switch ( k )
         {
            case ALLOW_PONDING:     project->AllowPonding    = m;  break;
            case SLOPE_WEIGHTING:   project->SlopeWeighting  = m;  break;
            case SKIP_STEADY_STATE: project->SkipSteadyState = m;  break;
            case IGNORE_RAINFALL:   project->IgnoreRainfall  = m;  break;
            case IGNORE_SNOWMELT:   project->IgnoreSnowmelt  = m;  break;
            case IGNORE_GWATER:     project->IgnoreGwater    = m;  break;
            case IGNORE_ROUTING:    project->IgnoreRouting   = m;  break;
            case IGNORE_QUALITY:    project->IgnoreQuality   = m;  break;
            case IGNORE_RDII:       project->IgnoreRDII      = m;  break;                 //(5.1.004)
         }
         break;

      case NORMAL_FLOW_LTD:
         m = findmatch(s2, NormalFlowWords);
         if ( m < 0 ) m = findmatch(s2, NoYesWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         project->NormalFlowLtd = m;
         break;

      case FORCE_MAIN_EQN:
         m = findmatch(s2, ForceMainEqnWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         project->ForceMainEqn = m;
         break;

      case LINK_OFFSETS:
         m = findmatch(s2, LinkOffsetWords);
         if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
         project->LinkOffsets = m;
         break;

         // --- compatibility option for selecting solution method for
         //     dynamic wave flow routing (NOT CURRENTLY USED)
      case COMPATIBILITY:
         if      ( strcomp(s2, "3") ) project->Compatibility = SWMM3;
         else if ( strcomp(s2, "4") ) project->Compatibility = SWMM4;
         else if ( strcomp(s2, "5") ) project->Compatibility = SWMM5;
         else return error_setInpError(ERR_KEYWORD, s2);
         break;

         // --- routing or lengthening time step (in decimal seconds)
         //     (lengthening time step is used in Courant stability formula
         //     to artificially lengthen conduits for dynamic wave flow routing
         //     (a value of 0 means that no lengthening is used))
      case ROUTE_STEP:
      case LENGTHENING_STEP:
         if ( !getDouble(s2, &tStep) )
         {
            if ( !datetime_strToTime(s2, &aTime) )
            {
               return error_setInpError(ERR_NUMBER, s2);
            }
            else
            {
               datetime_decodeTime(aTime, &h, &m, &s);
               h += 24*(int)aTime;
               s = s + 60*m + 3600*h;
               tStep = s;
            }
         }
         if ( k == ROUTE_STEP )
         {
            if ( tStep <= 0.0 ) return error_setInpError(ERR_NUMBER, s2);
            project->RouteStep = tStep;
         }
         else project->LengtheningStep = MAX(0.0, tStep);
         break;

         ////  Following code section added to release 5.1.008.  ////                   //(5.1.008)

         // --- minimum variable time step for dynamic wave routing
      case MIN_ROUTE_STEP:
         if ( !getDouble(s2, &project->MinRouteStep) || project->MinRouteStep < 0.0 )
            return error_setInpError(ERR_NUMBER, s2);
         break;

      case NUM_THREADS:
         m = atoi(s2);
         if ( m < 0 ) return error_setInpError(ERR_NUMBER, s2);
         project->NumThreads = m;
         break;
         ////

         // --- safety factor applied to variable time step estimates under
         //     dynamic wave flow routing (value of 0 indicates that variable
         //     time step option not used)
      case VARIABLE_STEP:
         if ( !getDouble(s2, &project->CourantFactor) )
            return error_setInpError(ERR_NUMBER, s2);
         if ( project->CourantFactor < 0.0 || project->CourantFactor > 2.0 )
            return error_setInpError(ERR_NUMBER, s2);
         break;

         // --- minimum surface area (ft2 or sq. meters) associated with nodes
         //     under dynamic wave flow routing
      case MIN_SURFAREA:
         project->MinSurfArea = atof(s2);
         break;

         // --- minimum conduit slope (%)
      case MIN_SLOPE:
         if ( !getDouble(s2, &project->MinSlope) )
            return error_setInpError(ERR_NUMBER, s2);
         if ( project->MinSlope < 0.0 || project->MinSlope >= 100 )
            return error_setInpError(ERR_NUMBER, s2);
         project->MinSlope /= 100.0;
         break;

         // --- maximum trials / time step for dynamic wave routing
      case MAX_TRIALS:
         m = atoi(s2);
         if ( m < 0 ) return error_setInpError(ERR_NUMBER, s2);
         project->MaxTrials = m;
         break;

         // --- head convergence tolerance for dynamic wave routing
      case HEAD_TOL:
         if ( !getDouble(s2, &project->HeadTol) )
         {
            return error_setInpError(ERR_NUMBER, s2);
         }
         break;

         // --- steady state tolerance on system inflow - outflow
      case SYS_FLOW_TOL:
         if ( !getDouble(s2, &project->SysFlowTol) )
         {
            return error_setInpError(ERR_NUMBER, s2);
         }
         project->SysFlowTol /= 100.0;
         break;

         // --- steady state tolerance on nodal lateral inflow
      case LAT_FLOW_TOL:
         if ( !getDouble(s2, &project->LatFlowTol) )
         {
            return error_setInpError(ERR_NUMBER, s2);
         }
         project->LatFlowTol /= 100.0;
         break;

      case TEMPDIR: // Temporary Directory
         sstrncpy(project->TempDir, s2, MAXFNAME);
         break;

   }
   return 0;
}

//=============================================================================

void initPointers(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: assigns NULL to all dynamic arrays for a new project.
//
{
   project->Gage     = NULL;
   project->Subcatch = NULL;
   project->Node     = NULL;
   project->Outfall  = NULL;
   project->Divider  = NULL;
   project->Storage  = NULL;
   project->Link     = NULL;
   project->Conduit  = NULL;
   project->Pump     = NULL;
   project->Orifice  = NULL;
   project->Weir     = NULL;
   project->Outlet   = NULL;
   project->Pollut   = NULL;
   project->Landuse  = NULL;
   project->Pattern  = NULL;
   project->Curve    = NULL;
   project->Tseries  = NULL;
   project->Transect = NULL;
   project->Shape    = NULL;
   project->Aquifer    = NULL;
   project->UnitHyd    = NULL;
   project->Snowmelt   = NULL;

   //add to make sure FREE does not freak out
   project->Rules = NULL;
   project->ActionList = NULL;
   project->Xnode = NULL;
   project->HortInfil = NULL;
   project->GAInfil = NULL;
   project->CNInfil = NULL;
   project->pXsect = NULL;
   project->LidProcs = NULL;
   project->LidGroups = NULL;
   project->LoadingTotals = NULL;
   project->QualTotals = NULL;
   project->StepQualTotals = NULL;
   project->SubcatchStats = NULL;
   project->NodeStats = NULL;
   project->LinkStats = NULL;
   project->StorageStats = NULL;
   project->OutfallStats = NULL;
   project->PumpStats = NULL;
   project->GW = NULL;
   project->LatFlowExpr = NULL;
   project->DeepFlowExpr = NULL;
   project->UHGroup = NULL;
   project->root = NULL;

   project->MemPoolAllocated = FALSE;
}

//=============================================================================

void setDefaults(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: assigns default values to project variables.
//
{
   int i, j;

   // Project title & temp. file path
   for (i = 0; i < MAXTITLE; i++) strcpy(project->Title[i], "");
   strcpy(project->TempDir, "");

   // Interface files
   project->Frain.mode      = SCRATCH_FILE;     // Use scratch rainfall file
   project->Fclimate.mode   = NO_FILE;
   project->Frunoff.mode    = NO_FILE;
   project->Frdii.mode      = NO_FILE;
   project->Fhotstart1.mode = NO_FILE;
   project->Fhotstart2.mode = NO_FILE;
   project->Finflows.mode   = NO_FILE;
   project->Foutflows.mode  = NO_FILE;
   project->Frain.file      = NULL;
   project->Fclimate.file   = NULL;
   project->Frunoff.file    = NULL;
   project->Frdii.file      = NULL;
   project->Fhotstart1.file = NULL;
   project->Fhotstart2.file = NULL;
   project->Finflows.file   = NULL;
   project->Foutflows.file  = NULL;
   project->Fout.file       = NULL;
   project->Fout.mode       = NO_FILE;

   // Analysis options
   project->UnitSystem      = US;               // US unit system
   project->FlowUnits       = CFS;              // CFS flow units
   project->InfilModel      = HORTON;           // Horton infiltration method
   project->RouteModel      = KW;               // Kin. wave flow routing method
   project->AllowPonding    = FALSE;            // No ponding at nodes
   project->InertDamping    = SOME;             // Partial inertial damping
   project->NormalFlowLtd   = BOTH;             // Default normal flow limitation
   project->ForceMainEqn    = H_W;              // Hazen-Williams eqn. for force mains
   project->LinkOffsets     = DEPTH_OFFSET;     // Use depth for link offsets
   project->LengtheningStep = 0;                // No lengthening of conduits
   project->CourantFactor   = 0.0;              // No variable time step
   project->MinSurfArea     = 0.0;              // Force use of default min. surface area
   project->SkipSteadyState = FALSE;            // Do flow routing in steady state periods
   project->IgnoreRainfall  = FALSE;            // Analyze rainfall/runoff
   project->IgnoreRDII      = FALSE;            // Analyze RDII                         //(5.1.004)
   project->IgnoreSnowmelt  = FALSE;            // Analyze snowmelt
   project->IgnoreGwater    = FALSE;            // Analyze groundwater
   project->IgnoreRouting   = FALSE;            // Analyze flow routing
   project->IgnoreQuality   = FALSE;            // Analyze water quality
   project->WetStep         = 300;              // Runoff wet time step (secs)
   project->DryStep         = 3600;             // Runoff dry time step (secs)
   project->RouteStep       = 300.0;            // Routing time step (secs)
   project->MinRouteStep    = 0.5;              // Minimum variable time step (sec)     //(5.1.008)
   project->ReportStep      = 900;              // Reporting time step (secs)
   project->StartDryDays    = 0.0;              // Antecedent dry days
   project->MaxTrials       = 0;                // Force use of default max. trials
   project->HeadTol         = 0.0;              // Force use of default head tolerance
   project->SysFlowTol      = 0.05;             // System flow tolerance for steady state
   project->LatFlowTol      = 0.05;             // Lateral flow tolerance for steady state
   project->NumThreads      = 0;                // Number of parallel threads to use

   // Deprecated options
   project->SlopeWeighting  = TRUE;             // Use slope weighting
   project->Compatibility   = SWMM4;            // Use SWMM 4 up/dn weighting method

   // Starting & ending date/time
   project->StartDate       = datetime_encodeDate(2004, 1, 1);
   project->StartTime       = datetime_encodeTime(0,0,0);
   project->StartDateTime   = project->StartDate + project->StartTime;
   project->EndDate         = project->StartDate;
   project->EndTime         = 0.0;
   project->ReportStartDate = NO_DATE;
   project->ReportStartTime = NO_DATE;
   project->SweepStart      = 1;
   project->SweepEnd        = 365;

   // Reporting options
   project->RptFlags.input         = FALSE;
   project->RptFlags.continuity    = TRUE;
   project->RptFlags.flowStats     = TRUE;
   project->RptFlags.controls      = FALSE;
   project->RptFlags.subcatchments = FALSE;
   project->RptFlags.nodes         = FALSE;
   project->RptFlags.links         = FALSE;
   project->RptFlags.nodeStats     = FALSE;

   // Temperature data
   project->Temp.dataSource  = NO_TEMP;
   project->Temp.tSeries     = -1;
   project->Temp.ta          = 70.0;
   project->Temp.elev        = 0.0;
   project->Temp.anglat      = 40.0;
   project->Temp.dtlong      = 0.0;
   project->Temp.tmax        = MISSING;

   // project->Wind speed data
   project->Wind.type = MONTHLY_WIND;
   for ( i=0; i<12; i++ ) project->Wind.aws[i] = 0.0;

   // project->Snowmelt parameters
   project->Snow.snotmp      = 34.0;
   project->Snow.tipm        = 0.5;
   project->Snow.rnm         = 0.6;

   // project->Snow areal depletion curves for pervious and impervious surfaces
   for ( i=0; i<2; i++ )
   {
      for ( j=0; j<10; j++) project->Snow.adc[i][j] = 1.0;
   }

   // Evaporation rates
   project->Evap.type = CONSTANT_EVAP;
   for (i=0; i<12; i++)
   {
      project->Evap.monthlyEvap[i] = 0.0;
      project->Evap.panCoeff[i]    = 1.0;
   }
   project->Evap.recoveryPattern = -1;
   project->Evap.recoveryFactor  = 1.0;
   project->Evap.tSeries = -1;
   project->Evap.dryOnly = FALSE;

   ////  Following code segment added to release 5.1.007.  ////                   //(5.1.007)
   ////
   // Climate adjustments
   for (i = 0; i < 12; i++)
   {
      project->Adjust.temp[i] = 0.0;   // additive adjustments
      project->Adjust.evap[i] = 0.0;   // additive adjustments
      project->Adjust.rain[i] = 1.0;   // multiplicative adjustments
      project->Adjust.hydcon[i] = 1.0; // hyd. conductivity adjustments                //(5.1.008)
   }
   project->Adjust.rainFactor = 1.0;
   project->Adjust.hydconFactor = 1.0;                                                  //(5.1.008)
   ////
}

//=============================================================================

void openFiles(Project* project, char *f1, char *f2, char *f3)
//
//  Input:   f1 = name of input file
//           f2 = name of report file
//           f3 = name of binary output file
//  Output:  none
//  Purpose: opens a project's input and report files.
//
{
   // --- initialize file pointers to NULL
   project->Finp.file = NULL;
   project->Frpt.file = NULL;
   project->Fout.file = NULL;

   // --- save file names
   sstrncpy(project->Finp.name, f1, MAXFNAME);
   sstrncpy(project->Frpt.name, f2, MAXFNAME);
   sstrncpy(project->Fout.name, f3, MAXFNAME);

   // --- check that file names are not identical
   if (strcomp(f1, f2) || strcomp(f1, f3) || strcomp(f2, f3))
   {
      writecon(FMT11);
      project->ErrorCode = ERR_FILE_NAME;
      return;
   }

   // --- open input and report files
   if ((project->Finp.file = fopen(f1,"rt")) == NULL)
   {
      writecon(FMT12);
      writecon(f1);
      project->ErrorCode = ERR_INP_FILE;
      return;
   }
   if ((project->Frpt.file = fopen(f2,"wt")) == NULL)
   {
      writecon(FMT13);
      project->ErrorCode = ERR_RPT_FILE;
      return;
   }
}

//=============================================================================

void createObjects(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: allocates memory for project's objects.
//
//  NOTE: number of each type of object has already been determined in
//        project_readInput().
//
{
   int j, k;

   // --- allocate memory for each category of object
   if ( project->ErrorCode ) return;
   project->Gage     = (TGage *)     calloc(project->Nobjects[GAGE],     sizeof(TGage));
   project->Subcatch = (TSubcatch *) calloc(project->Nobjects[SUBCATCH], sizeof(TSubcatch));
   project->Node     = (TNode *)     calloc(project->Nobjects[NODE],     sizeof(TNode));
   project->Outfall  = (TOutfall *)  calloc(project->Nnodes[OUTFALL],    sizeof(TOutfall));
   project->Divider  = (TDivider *)  calloc(project->Nnodes[DIVIDER],    sizeof(TDivider));
   project->Storage  = (TStorage *)  calloc(project->Nnodes[STORAGE],    sizeof(TStorage));
   project->Link     = (TLink *)     calloc(project->Nobjects[LINK],     sizeof(TLink));
   project->Conduit  = (TConduit *)  calloc(project->Nlinks[CONDUIT],    sizeof(TConduit));
   project->Pump     = (TPump *)     calloc(project->Nlinks[PUMP],       sizeof(TPump));
   project->Orifice  = (TOrifice *)  calloc(project->Nlinks[ORIFICE],    sizeof(TOrifice));
   project->Weir     = (TWeir *)     calloc(project->Nlinks[WEIR],       sizeof(TWeir));
   project->Outlet   = (TOutlet *)   calloc(project->Nlinks[OUTLET],     sizeof(TOutlet));
   project->Pollut   = (TPollut *)   calloc(project->Nobjects[POLLUT],   sizeof(TPollut));
   project->Landuse  = (TLanduse *)  calloc(project->Nobjects[LANDUSE],  sizeof(TLanduse));
   project->Pattern  = (TPattern *)  calloc(project->Nobjects[TIMEPATTERN],  sizeof(TPattern));
   project->Curve    = (TTable *)    calloc(project->Nobjects[CURVE],    sizeof(TTable));
   project->Tseries  = (TTable *)    calloc(project->Nobjects[TSERIES],  sizeof(TTable));
   project->Aquifer  = (TAquifer *)  calloc(project->Nobjects[AQUIFER],  sizeof(TAquifer));
   project->UnitHyd  = (TUnitHyd *)  calloc(project->Nobjects[UNITHYD],  sizeof(TUnitHyd));
   project->Snowmelt = (TSnowmelt *) calloc(project->Nobjects[SNOWMELT], sizeof(TSnowmelt));
   project->Shape    = (TShape *)    calloc(project->Nobjects[SHAPE],    sizeof(TShape));

   // --- create LID objects
   lid_create(project,project->Nobjects[LID], project->Nobjects[SUBCATCH]);

   // --- create control rules
   project->ErrorCode = controls_create(project,project->Nobjects[CONTROL]);
   if ( project->ErrorCode ) return;

   // --- create cross section transects
   project->ErrorCode = transect_create(project, project->Nobjects[TRANSECT]);
   if ( project->ErrorCode ) return;

   // --- allocate memory for infiltration data
   infil_create(project,project->Nobjects[SUBCATCH], project->InfilModel);

   // --- allocate memory for water quality state variables
   for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
   {
      project->Subcatch[j].initBuildup =
            (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Subcatch[j].oldQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Subcatch[j].newQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Subcatch[j].pondedQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Subcatch[j].totalLoad  = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
   }
   for (j = 0; j < project->Nobjects[NODE]; j++)
   {
      project->Node[j].oldQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Node[j].newQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Node[j].extInflow = NULL;
      project->Node[j].dwfInflow = NULL;
      project->Node[j].rdiiInflow = NULL;
      project->Node[j].treatment = NULL;
   }
   for (j = 0; j < project->Nobjects[LINK]; j++)
   {
      project->Link[j].oldQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Link[j].newQual = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      project->Link[j].totalLoad = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
   }

   // --- allocate memory for land use buildup/washoff functions
   for (j = 0; j < project->Nobjects[LANDUSE]; j++)
   {
      project->Landuse[j].buildupFunc =
            (TBuildup *) calloc(project->Nobjects[POLLUT], sizeof(TBuildup));
      project->Landuse[j].washoffFunc =
            (TWashoff *) calloc(project->Nobjects[POLLUT], sizeof(TWashoff));
   }

   // --- allocate memory for subcatchment landuse factors
   for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
   {
      project->Subcatch[j].landFactor =
            (TLandFactor *) calloc(project->Nobjects[LANDUSE], sizeof(TLandFactor));
      for (k = 0; k < project->Nobjects[LANDUSE]; k++)
      {
         project->Subcatch[j].landFactor[k].buildup =
               (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
      }
   }

   // --- initialize buildup & washoff functions
   for (j = 0; j < project->Nobjects[LANDUSE]; j++)
   {
      for (k = 0; k < project->Nobjects[POLLUT]; k++)
      {
         project->Landuse[j].buildupFunc[k].funcType = NO_BUILDUP;
         project->Landuse[j].buildupFunc[k].normalizer = PER_AREA;
         project->Landuse[j].washoffFunc[k].funcType = NO_WASHOFF;
      }
   }

   // --- initialize rain gage properties
   for (j = 0; j < project->Nobjects[GAGE]; j++)
   {
      project->Gage[j].tSeries = -1;
      strcpy(project->Gage[j].fname, "");
   }

   // --- initialize subcatchment properties
   for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
   {
      project->Subcatch[j].outSubcatch = -1;
      project->Subcatch[j].outNode     = -1;
      project->Subcatch[j].infil       = -1;
      project->Subcatch[j].groundwater = NULL;
      project->Subcatch[j].gwLatFlowExpr = NULL;                                      //(5.1.007)
      project->Subcatch[j].gwDeepFlowExpr = NULL;                                     //(5.1.007)
      project->Subcatch[j].snowpack    = NULL;
      project->Subcatch[j].lidArea     = 0.0;
      for (k = 0; k < project->Nobjects[POLLUT]; k++)
      {
         project->Subcatch[j].initBuildup[k] = 0.0;
      }
   }

   // --- initialize RDII unit hydrograph properties
   for ( j = 0; j < project->Nobjects[UNITHYD]; j++ ) rdii_initUnitHyd(project,j);

   // --- initialize snowmelt properties
   for (j = 0; j < project->Nobjects[SNOWMELT]; j++) snow_initSnowmelt(project, j);

   // --- initialize storage node exfiltration                                //(5.1.007)
   for (j = 0; j < project->Nnodes[STORAGE]; j++) project->Storage[j].exfil = NULL;             //(5.1.007)

   // --- initialize link properties
   for (j = 0; j < project->Nobjects[LINK]; j++)
   {
      project->Link[j].xsect.type   = -1;
      project->Link[j].cLossInlet   = 0.0;
      project->Link[j].cLossOutlet  = 0.0;
      project->Link[j].cLossAvg     = 0.0;
      project->Link[j].hasFlapGate  = FALSE;
   }
   for (j = 0; j < project->Nlinks[PUMP]; j++) project->Pump[j].pumpCurve  = -1;

   // --- initialize reporting flags
   for (j = 0; j < project->Nobjects[SUBCATCH]; j++) project->Subcatch[j].rptFlag = FALSE;
   for (j = 0; j < project->Nobjects[NODE]; j++) project->Node[j].rptFlag = FALSE;
   for (j = 0; j < project->Nobjects[LINK]; j++) project->Link[j].rptFlag = FALSE;

   //  --- initialize curves, time series, and time patterns
   for (j = 0; j < project->Nobjects[CURVE]; j++)   table_init(&project->Curve[j]);
   for (j = 0; j < project->Nobjects[TSERIES]; j++) table_init(&project->Tseries[j]);
   for (j = 0; j < project->Nobjects[TIMEPATTERN]; j++) inflow_initDwfPattern(project,j);
}

//=============================================================================

void deleteObjects(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory allocated for a project's objects.
//
//  NOTE: care is taken to first free objects that are properties of another
//        object before the latter is freed (e.g., we must free a
//        subcatchment's land use factors before freeing the subcatchment).
//
{
   int j, k;

   // --- free memory for landuse factors & groundwater
   if ( project->Subcatch ) for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
   {
      for (k = 0; k < project->Nobjects[LANDUSE]; k++)
      {
         FREE(project->Subcatch[j].landFactor[k].buildup);
      }
      FREE(project->Subcatch[j].landFactor);
      FREE(project->Subcatch[j].groundwater);
      gwater_deleteFlowExpression(project,j);
      FREE(project->Subcatch[j].snowpack);
   }

   // --- free memory for buildup/washoff functions
   if ( project->Landuse ) for (j = 0; j < project->Nobjects[LANDUSE]; j++)
   {
      FREE(project->Landuse[j].buildupFunc);
      FREE(project->Landuse[j].washoffFunc)
   }

   // --- free memory for water quality state variables
   if ( project->Subcatch ) for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
   {
      FREE(project->Subcatch[j].initBuildup);
      FREE(project->Subcatch[j].oldQual);
      FREE(project->Subcatch[j].newQual);
      FREE(project->Subcatch[j].pondedQual);
      FREE(project->Subcatch[j].totalLoad);
   }
   if ( project->Node ) for (j = 0; j < project->Nobjects[NODE]; j++)
   {
      FREE(project->Node[j].oldQual);
      FREE(project->Node[j].newQual);
   }
   if ( project->Link ) for (j = 0; j < project->Nobjects[LINK]; j++)
   {
      FREE(project->Link[j].oldQual);
      FREE(project->Link[j].newQual);
      FREE(project->Link[j].totalLoad);
   }

   // --- free memory used for rainfall infiltration
   infil_delete(project);

   ////  Added for release 5.1.007.  ////                                         //(5.1.007)
   ////
   // --- free memory used for storage exfiltration
   if ( project->Node ) for (j = 0; j < project->Nnodes[STORAGE]; j++)
   {
      if ( project->Storage[j].exfil )
      {
         FREE(project->Storage[j].exfil->btmExfil);
         FREE(project->Storage[j].exfil->bankExfil);
         FREE(project->Storage[j].exfil);
      }
   }
   ////

   // --- free memory used for outfall pollutants loads                       //(5.1.008)
   if ( project->Node ) for (j = 0; j < project->Nnodes[OUTFALL]; j++)                          //(5.1.008)
      FREE(project->Outfall[j].wRouted);                                              //(5.1.008)

   // --- free memory used for nodal inflows & treatment functions
   if ( project->Node ) for (j = 0; j < project->Nobjects[NODE]; j++)
   {
      inflow_deleteExtInflows(project,j);
      inflow_deleteDwfInflows(project,j);
      rdii_deleteRdiiInflow(project,j);
      treatmnt_delete(project,j);
   }

   // --- delete table entries for curves and time series
   if ( project->Tseries ) for (j = 0; j < project->Nobjects[TSERIES]; j++)
      table_deleteEntries(&project->Tseries[j]);
   if ( project->Curve ) for (j = 0; j < project->Nobjects[CURVE]; j++)
      table_deleteEntries(&project->Curve[j]);

   // --- delete cross section transects
   transect_delete(project);

   // --- delete control rules
   controls_delete(project);

   // --- delete LIDs
   lid_delete(project);

   // --- now free each major category of object
   FREE(project->Gage);
   FREE(project->Subcatch);
   FREE(project->Node);
   FREE(project->Outfall);
   FREE(project->Divider);
   FREE(project->Storage);
   FREE(project->Link);
   FREE(project->Conduit);
   FREE(project->Pump);
   FREE(project->Orifice);
   FREE(project->Weir);
   FREE(project->Outlet);
   FREE(project->Pollut);
   FREE(project->Landuse);
   FREE(project->Pattern);
   FREE(project->Curve);
   FREE(project->Tseries);
   FREE(project->Aquifer);
   FREE(project->UnitHyd);
   FREE(project->Snowmelt);
   FREE(project->Shape);
}

//=============================================================================

void createHashTables(Project* project)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: allocates memory for object ID hash tables
//
{   int j;
    project->MemPoolAllocated = FALSE;
     for (j = 0; j < MAX_OBJ_TYPES ; j++)
     {
        project->Htable[j] = HTcreate();
        if ( project->Htable[j] == NULL ) report_writeErrorMsg(project,ERR_MEMORY, "");
     }

      // --- initialize memory pool used to store object ID's
      if (AllocInit(project) == NULL) report_writeErrorMsg(project, ERR_MEMORY, "");
      else project->MemPoolAllocated = TRUE;
}

//=============================================================================

void deleteHashTables(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory allocated for object ID hash tables
//
{
   int j;
   for (j = 0; j < MAX_OBJ_TYPES; j++)
   {
      if ( project->Htable[j] != NULL ) HTfree(project->Htable[j]);
   }

   // --- free object ID memory pool
   if (project->MemPoolAllocated) AllocFreePool(project);
}

//=============================================================================
