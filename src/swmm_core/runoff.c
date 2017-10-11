//-----------------------------------------------------------------------------
//   runoff.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//   Author:   L. Rossman
//
//   Runoff analysis functions.
//
//   Build 5.1.007:
//   - Climate file now opened in climate.c module.
//
//   Build 5.1.008:
//   - Memory for runoff pollutant load now allocated and freed in this module.
//   - Runoff time step chosen so that simulation does not exceed total duration.
//   - State of LIDs considered when choosing wet or dry time step.
//   - More checks added to skip over subcatchments with zero area.
//   - Support added for sending outfall node discharge onto a subcatchment.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include "headers.h"
#include "odesolve.h"

////-----------------------------------------------------------------------------
//// Shared variables
////-----------------------------------------------------------------------------
//static char  project->IsRaining;                // TRUE if precip. falls on study area
//static char  project->HasRunoff;                // TRUE if study area generates runoff
//static char  project->HasSnow;                  // TRUE if any snow cover on study area
//static int   project->Nsteps;                   // number of runoff time steps taken
//static int   project->MaxSteps;                 // final number of runoff time steps
//static long  project->MaxStepsPos;              // position in Runoff interface file
//                                       //    where project->MaxSteps is saved
//
////-----------------------------------------------------------------------------
////  Exportable variables 
////-----------------------------------------------------------------------------
//char    project-> HasWetLids;  // TRUE if any LIDs are wet (used in lidproc.c)           //(5.1.008)
//double* project->OutflowLoad; // exported pollutant mass load (used in surfqual.c)

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
//  extern float* project->SubcatchResults;         // Results vector defined in OUTPUT.C

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
// runoff_open     (called from swmm_start in swmm5.c)
// runoff_execute  (called from swmm_step in swmm5.c)
// runoff_close    (called from swmm_end in swmm5.c)

//-----------------------------------------------------------------------------
// Local functions
//-----------------------------------------------------------------------------
static double runoff_getTimeStep(Project* project, DateTime currentDate);
static void   runoff_initFile(Project* project);
static void   runoff_readFromFile(Project* project);
static void   runoff_saveToFile(Project* project, float tStep);
static void   runoff_getOutfallRunon(Project* project, double tStep);                            //(5.1.008)

//=============================================================================

int runoff_open(Project* project)
//
//  Input:   none
//  Output:  returns the global error code
//  Purpose: opens the runoff analyzer.
//
{
    project->IsRaining = FALSE;
    project->HasRunoff = FALSE;
    project->HasSnow = FALSE;
    project->Nsteps = 0;

    // --- open the Ordinary Differential Equation solver
    if ( !odesolve_open(project,MAXODES) ) report_writeErrorMsg(project,ERR_ODE_SOLVER, "");

    // --- allocate memory for pollutant runoff loads                          //(5.1.008)
    project->OutflowLoad = NULL;
    if ( project->Nobjects[POLLUT] > 0 )
    {
        project->OutflowLoad = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
        if ( !project->OutflowLoad ) report_writeErrorMsg(project,ERR_MEMORY, "");
    }

    // --- see if a runoff interface file should be opened
    switch ( project->Frunoff.mode )
    {
      case USE_FILE:
        if ( (project->Frunoff.file = fopen(project->Frunoff.name, "r+b")) == NULL)
            report_writeErrorMsg(project,ERR_RUNOFF_FILE_OPEN, project->Frunoff.name);
        else runoff_initFile(project);
        break;
      case SAVE_FILE:
        if ( (project->Frunoff.file = fopen(project->Frunoff.name, "w+b")) == NULL)
            report_writeErrorMsg(project,ERR_RUNOFF_FILE_OPEN, project->Frunoff.name);
        else runoff_initFile(project);
        break;
    }

////  Call to climate_openFile() moved to climate_validate().  ////            //(5.1.007)
    return project->ErrorCode;
}

//=============================================================================

void runoff_close(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: closes the runoff analyzer.
//
{
    // --- close the ODE solver
    odesolve_close(project);

    // --- free memory for pollutant runoff loads                              //(5.1.008)
    FREE(project->OutflowLoad);

    // --- close runoff interface file if in use
    if ( project->Frunoff.file )
    {
        // --- write to file number of time steps simulated
        if ( project->Frunoff.mode == SAVE_FILE )
        {
            fseek(project->Frunoff.file, project->MaxStepsPos, SEEK_SET);
            fwrite(&project->Nsteps, sizeof(int), 1, project->Frunoff.file);
        }
        fclose(project->Frunoff.file);
    }

    // --- close climate file if in use
    if ( project->Fclimate.file ) fclose(project->Fclimate.file);
}

//=============================================================================

void runoff_execute(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: computes runoff from each subcatchment at current runoff time.
//
{
    int      j;                        // object index
    int      day;                      // day of calendar year
    double   runoffStep;               // runoff time step (sec)
    double   runoff;                   // subcatchment runoff (ft/sec)
    DateTime currentDate;              // current date/time 
    char     canSweep;                 // TRUE if street sweeping can occur

    if ( project->ErrorCode ) return;

    // --- convert elapsed runoff time in milliseconds to a calendar date
    currentDate = getDateTime(project, project->NewRunoffTime);

    // --- update climatological conditions
    climate_setState(project, currentDate);

    // --- if no subcatchments then simply update runoff elapsed time
    if ( project->Nobjects[SUBCATCH] == 0 )
    {
        project->OldRunoffTime = project->NewRunoffTime;
        project->NewRunoffTime += (double)(1000.0 * project->DryStep);
        project->NewRunoffTime = MIN(project->NewRunoffTime, project->TotalDuration);                     //(5.1.008)
        return;
    }

    // --- update current rainfall at each raingage
    //     NOTE: must examine gages in sequential order due to possible
    //     presence of co-gages (gages that share same rain time series).
    project->IsRaining = FALSE;
    for (j = 0; j < project->Nobjects[GAGE]; j++)
    {
        gage_setState(project, j, currentDate);
        if ( project->Gage[j].rainfall > 0.0 ) project->IsRaining = TRUE;
    }

    // --- read runoff results from interface file if applicable
    if ( project->Frunoff.mode == USE_FILE )
    {
        runoff_readFromFile(project);
        return;
    }

    // --- see if street sweeping can occur on current date
    day = datetime_dayOfYear(currentDate);
    if ( day >= project->SweepStart && day <= project->SweepEnd ) canSweep = TRUE;
    else canSweep = FALSE;

    // --- get runoff time step (in seconds)
    runoffStep = runoff_getTimeStep(project,currentDate);
    if ( runoffStep <= 0.0 )
    {
        project->ErrorCode = ERR_TIMESTEP;
        return;
    }

    // --- update runoff time clock (in milliseconds)
    project->OldRunoffTime = project->NewRunoffTime;
    project->NewRunoffTime += (double)(1000.0 * runoffStep);

////  Following code segment added to release 5.1.008.  ////                   //(5.1.008)
////
    // --- adjust runoff step so that total duration not exceeded
    if ( project->NewRunoffTime > project->TotalDuration )
    {
        runoffStep = (project->TotalDuration - project->OldRunoffTime) / 1000.0;
        project->NewRunoffTime = project->TotalDuration;
    }
////

    // --- update old state of each subcatchment, 
	for (j = 0; j < project->Nobjects[SUBCATCH]; j++) subcatch_setOldState(project, j);

    // --- determine any runon from drainage system outfall nodes              //(5.1.008)
    runoff_getOutfallRunon(project,runoffStep);                                        //(5.1.008)

    // --- determine runon from upstream subcatchments, and implement snow removal
    for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
    {
        if (project-> Subcatch[j].area == 0.0 ) continue;                               //(5.1.008)
		subcatch_getRunon(project, j);
		if (!project->IgnoreSnowmelt) snow_plowSnow(project, j, runoffStep);
    }
    
    // --- determine runoff and pollutant buildup/washoff in each subcatchment
    project->HasSnow = FALSE;
    project->HasRunoff = FALSE;
    project-> HasWetLids = FALSE;                                                        //(5.1.008)
    for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
    {
        // --- find total runoff rate (in ft/sec) over the subcatchment
        //     (the amount that actually leaves the subcatchment (in cfs)
        //     is also computed and is stored inproject-> Subcatch[j].newRunoff)
        if (project-> Subcatch[j].area == 0.0 ) continue;                               //(5.1.008)
		runoff = subcatch_getRunoff(project, j, runoffStep);

        // --- update state of study area surfaces
        if ( runoff > 0.0 ) project->HasRunoff = TRUE;
        if (project-> Subcatch[j].newSnowDepth > 0.0 ) project->HasSnow = TRUE;

        // --- skip pollutant buildup/washoff if quality ignored
        if ( project->IgnoreQuality ) continue;

        // --- add to pollutant buildup if runoff is negligible
		if (runoff < MIN_RUNOFF) surfqual_getBuildup(project, j, runoffStep);

        // --- reduce buildup by street sweeping
        if ( canSweep &&project-> Subcatch[j].rainfall <= MIN_RUNOFF)
			surfqual_sweepBuildup(project, j, currentDate);

        // --- compute pollutant washoff 
		surfqual_getWashoff(project, j, runoff, runoffStep);
    }

    // --- update tracking of system-wide max. runoff rate
    stats_updateMaxRunoff(project);

    // --- save runoff results to interface file if one is used
    project->Nsteps++;
    if ( project->Frunoff.mode == SAVE_FILE )
    {
        runoff_saveToFile(project,(float)runoffStep);
    }

    // --- reset subcatchment runon to 0
    for (j = 0; j < project->Nobjects[SUBCATCH]; j++)project-> Subcatch[j].runon = 0.0;
}

//=============================================================================

double runoff_getTimeStep(Project* project, DateTime currentDate)
//
//  Input:   currentDate = current simulation date/time
//  Output:  time step (sec)
//  Purpose: computes a time step to use for runoff calculations.
//
{
    int  j;
    long timeStep;
    long maxStep = project->DryStep;

    // --- find shortest time until next evaporation or rainfall value
    //     (this represents the maximum possible time step)
    timeStep = datetime_timeDiff(climate_getNextEvapDate(project), currentDate);      //(5.1.008)
    if ( timeStep > 0.0 && timeStep < maxStep ) maxStep = timeStep;            //(5.1.008)
    for (j = 0; j < project->Nobjects[GAGE]; j++)
    {
        timeStep = datetime_timeDiff(gage_getNextRainDate(project, j, currentDate),
                   currentDate);
        if ( timeStep > 0 && timeStep < maxStep ) maxStep = timeStep;
    }

    // --- determine whether wet or dry time step applies
    if ( project->IsRaining || project->HasSnow || project->HasRunoff || project-> HasWetLids ) timeStep = project-> WetStep; //(5.1.008)
    else timeStep = project->DryStep;

    // --- limit time step if necessary
    if ( timeStep > maxStep ) timeStep = maxStep;
    return (double)timeStep;
}

//=============================================================================

void runoff_initFile(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes a Runoff Interface file for saving results.
//
{
    int   nSubcatch;
    int   nPollut;
    int   flowUnits;
    char  fileStamp[] = "SWMM5-RUNOFF";
    char  fStamp[] = "SWMM5-RUNOFF";

    project->MaxSteps = 0;
    if ( project->Frunoff.mode == SAVE_FILE )
    {
        // --- write file stamp, # subcatchments & # pollutants to file
        nSubcatch = project->Nobjects[SUBCATCH];
        nPollut = project->Nobjects[POLLUT];
        flowUnits = project->FlowUnits;
        fwrite(fileStamp, sizeof(char), strlen(fileStamp), project->Frunoff.file);
        fwrite(&nSubcatch, sizeof(int), 1, project->Frunoff.file);
        fwrite(&nPollut, sizeof(int), 1, project->Frunoff.file);
        fwrite(&flowUnits, sizeof(int), 1, project->Frunoff.file);
        project->MaxStepsPos = ftell(project->Frunoff.file); 
        fwrite(&project->MaxSteps, sizeof(int), 1, project->Frunoff.file);
    }

    if ( project->Frunoff.mode == USE_FILE )
    {
        // --- check that interface file contains proper header records
        fread(fStamp, sizeof(char), strlen(fileStamp), project->Frunoff.file);
        if ( strcmp(fStamp, fileStamp) != 0 )
        {
            report_writeErrorMsg(project,ERR_RUNOFF_FILE_FORMAT, "");
            return;
        }
        nSubcatch = -1;
        nPollut = -1;
        flowUnits = -1;
        fread(&nSubcatch, sizeof(int), 1, project->Frunoff.file);
        fread(&nPollut, sizeof(int), 1, project->Frunoff.file);
        fread(&flowUnits, sizeof(int), 1, project->Frunoff.file);
        fread(&project->MaxSteps, sizeof(int), 1, project->Frunoff.file);
        if ( nSubcatch != project->Nobjects[SUBCATCH]
        ||   nPollut   != project->Nobjects[POLLUT]
        ||   flowUnits != project->FlowUnits
        ||   project->MaxSteps  <= 0 )
        {
             report_writeErrorMsg(project,ERR_RUNOFF_FILE_FORMAT, "");
        }
    }
}

//=============================================================================

void  runoff_saveToFile(Project* project, float tStep)
//
//  Input:   tStep = runoff time step (sec)
//  Output:  none
//  Purpose: saves current runoff results to Runoff Interface file.
//
{
    int j;
    int n = MAX_SUBCATCH_RESULTS + project->Nobjects[POLLUT] - 1;
    

    fwrite(&tStep, sizeof(float), 1, project->Frunoff.file);
    for (j=0; j<project->Nobjects[SUBCATCH]; j++)
    {
		subcatch_getResults(project, j, 1.0, project->SubcatchResults);
        fwrite(project->SubcatchResults, sizeof(float), n, project->Frunoff.file);
    }
}

//=============================================================================

void  runoff_readFromFile(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reads runoff results from Runoff Interface file for current time.
//
{
    int    i, j;
    int    nResults;                   // number of results per subcatch.
    int    kount;                      // count of items read from file
    float  tStep;                      // runoff time step (sec)
    TGroundwater* gw;                  // ptr. to Groundwater object

    // --- make sure not past end of file
    if ( project->Nsteps > project->MaxSteps )
    {
         report_writeErrorMsg(project,ERR_RUNOFF_FILE_END, "");
         return;
    }

    // --- replace old state with current one for all subcatchments
	for (j = 0; j < project->Nobjects[SUBCATCH]; j++) subcatch_setOldState(project, j);

    // --- read runoff time step
    kount = 0;
    kount += fread(&tStep, sizeof(float), 1, project->Frunoff.file);

    // --- compute number of results saved for each subcatchment
    nResults = MAX_SUBCATCH_RESULTS + project->Nobjects[POLLUT] - 1;

    // --- for each subcatchment
    for (j = 0; j < project->Nobjects[SUBCATCH]; j++)
    {
        // --- read vector of saved results
        kount += fread(project->SubcatchResults, sizeof(float), nResults, project->Frunoff.file);

        // --- extract hydrologic results, converting units where necessary
        //     (results were saved to file in user's units)
       project-> Subcatch[j].newSnowDepth = project->SubcatchResults[SUBCATCH_SNOWDEPTH] /
                                   UCF(project, RAINDEPTH);
       project-> Subcatch[j].evapLoss     = project->SubcatchResults[SUBCATCH_EVAP] /
		   UCF(project, RAINFALL);
       project-> Subcatch[j].infilLoss    = project->SubcatchResults[SUBCATCH_INFIL] /
		   UCF(project, RAINFALL);
       project-> Subcatch[j].newRunoff    = project->SubcatchResults[SUBCATCH_RUNOFF] /
		   UCF(project, FLOW);
        gw =project-> Subcatch[j].groundwater;
        if ( gw )
        {
			gw->newFlow = project->SubcatchResults[SUBCATCH_GW_FLOW] / UCF(project, FLOW);
            gw->lowerDepth = project->Aquifer[gw->aquifer].bottomElev -
				(project->SubcatchResults[SUBCATCH_GW_ELEV] / UCF(project, LENGTH));
            gw->theta      = project->SubcatchResults[SUBCATCH_SOIL_MOIST];
        }

        // --- extract water quality results
        for (i = 0; i < project->Nobjects[POLLUT]; i++)
        {
           project-> Subcatch[j].newQual[i] = project->SubcatchResults[SUBCATCH_WASHOFF + i];
        }
    }

    // --- report error if not enough values were read
    if ( kount < 1 + project->Nobjects[SUBCATCH] * nResults )
    {
         report_writeErrorMsg(project,ERR_RUNOFF_FILE_READ, "");
         return;
    }

    // --- update runoff time clock
    project->OldRunoffTime = project->NewRunoffTime;
    project->NewRunoffTime = project->OldRunoffTime + (double)(tStep)*1000.0;
    project->NewRunoffTime = MIN(project->NewRunoffTime, project->TotalDuration);                         //(5.1.008)
    project->Nsteps++;
}

//=============================================================================

////  New function added for release 5.1.008.  ////                            //(5.1.008)

void runoff_getOutfallRunon(Project* project, double tStep)
//
//  Input:   tStep = runoff time step (sec)
//  Output:  none
//  Purpose: adds any flow and pollutant loads returned from drainage system
//           outfalls to runoff system subcatchments.
//
{
    int i, k, p;
    double w;

    // --- examine each outfall node
    for (i = 0; i < project->Nnodes[OUTFALL]; i++)
    {
        // --- ignore node if outflow not re-routed onto a subcatchment
        k = project->Outfall[i].routeTo;
        if ( k < 0 ) continue;
        if (project-> Subcatch[k].area == 0.0 ) continue;

        // --- add outfall's flow to subcatchment as runon and re-set routed
        //     flow volume to 0
		subcatch_addRunonFlow(project, k, project->Outfall[i].vRouted / tStep);
        massbal_updateRunoffTotals(project, RUNOFF_RUNON, project->Outfall[i].vRouted);
        project->Outfall[i].vRouted = 0.0;

        // --- add outfall's pollutant load on to subcatchment's wet
        //     deposition load and re-set routed load to 0
        //     (Subcatch.newQual is being used as a temporary load accumulator)
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            w = project->Outfall[i].wRouted[p] * LperFT3;
            massbal_updateLoadingTotals(project, DEPOSITION_LOAD, p, w * project->Pollut[p].mcf);
           project-> Subcatch[k].newQual[p] += w / tStep;
            project->Outfall[i].wRouted[p] = 0.0;
        }
    }
}
