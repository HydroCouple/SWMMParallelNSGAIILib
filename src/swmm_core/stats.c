//-----------------------------------------------------------------------------
//   stats.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//   Author:   L. Rossman (EPA)
//             R. Dickinson (CDM)
//
//   Simulation statistics functions.
//
//   Build 5.1.007:
//   - Exfiltration losses added to storage node statistics.
//
//   Build 5.1.008:
//   - Support for updating groundwater statistics added.
//   - Support for updating maximum reported nodal depths added.
//   - OpenMP parallelization applied to updating node and link flow statistics.
//   - Updating of time that conduit is upstrm/dnstrm full was modified.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <math.h>

#ifdef USE_OPENMP
#include <omp.h>                                                               //(5.1.008)
#endif

#include "headers.h"

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
#define MAX_STATS 5
//static TSysStats       project->SysStats;
//static TMaxStats       project->MaxMassBalErrs[MAX_STATS];
//static TMaxStats       project->MaxCourantCrit[MAX_STATS];
//static TMaxStats       project->MaxFlowTurns[MAX_STATS];
//static double          project->SysOutfallFlow;
//
////-----------------------------------------------------------------------------
////  Exportable variables (shared with statsrpt.c)
////-----------------------------------------------------------------------------
//TSubcatchStats* project->SubcatchStats;
//TNodeStats*     project->NodeStats;
//TLinkStats*     project->LinkStats;
//TStorageStats*  project->StorageStats;
//TOutfallStats*  project->OutfallStats;
//TPumpStats*     project->PumpStats;
//double          project->MaxOutfallFlow;
//double          project->MaxRunoffFlow;
//
////-----------------------------------------------------------------------------
////  Imported variables
////-----------------------------------------------------------------------------
//extern double*         NodeInflow;     // defined in massbal.c
//extern double*         NodeOutflow;    // defined in massbal.c

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  stats_open                    (called from swmm_start in swmm5.c)
//  stats_close                   (called from swmm_end in swmm5.c)
//  stats_report                  (called from swmm_end in swmm5.c)
//  stats_updateSubcatchStats     (called from subcatch_getRunoff)
//  stats_updateGwaterStats       (called from gwater_getGroundwater)          //(5.1.008)
//  stats_updateFlowStats         (called from routing_execute)
//  stats_updateCriticalTimeCount (called from getVariableStep in dynwave.c)
//  stats_updateMaxNodeDepth      (called from output_saveNodeResults)         //(5.1.008)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void stats_updateNodeStats(Project* project, int node, double tStep, DateTime aDate);
static void stats_updateLinkStats(Project* project, int link, double tStep, DateTime aDate);
static void stats_findMaxStats(Project* project);
static void stats_updateMaxStats(TMaxStats maxStats[], int i, int j, double x);

//=============================================================================

int  stats_open(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: opens the simulation statistics system.
//
{
  int j, k;

  // --- set all pointers to NULL
  project->NodeStats = NULL;
  project->LinkStats = NULL;
  project->StorageStats = NULL;
  project->OutfallStats = NULL;
  project->PumpStats = NULL;

  // --- allocate memory for & initialize subcatchment statistics
  project->SubcatchStats = NULL;
  if ( project->Nobjects[SUBCATCH] > 0 )
  {
    project->SubcatchStats = (TSubcatchStats *) calloc(project->Nobjects[SUBCATCH],
                                                       sizeof(TSubcatchStats));
    if ( !project->SubcatchStats )
    {
      report_writeErrorMsg(project, ERR_MEMORY, "");
      return project->ErrorCode;
    }
    for (j=0; j<project->Nobjects[SUBCATCH]; j++)
    {
      project->SubcatchStats[j].precip  = 0.0;
      project->SubcatchStats[j].runon   = 0.0;
      project->SubcatchStats[j].evap    = 0.0;
      project->SubcatchStats[j].infil   = 0.0;
      project->SubcatchStats[j].runoff  = 0.0;
      project->SubcatchStats[j].maxFlow = 0.0;
    }

    ////  Added to release 5.1.008.  ////                                          //(5.1.008)
    ////
    for (j=0; j<project->Nobjects[SUBCATCH]; j++)
    {
      if ( project->Subcatch[j].groundwater == NULL ) continue;
      project->Subcatch[j].groundwater->stats.avgUpperMoist = 0.0;
      project->Subcatch[j].groundwater->stats.avgWaterTable = 0.0;
      project->Subcatch[j].groundwater->stats.infil = 0.0;
      project->Subcatch[j].groundwater->stats.latFlow = 0.0;
      project->Subcatch[j].groundwater->stats.deepFlow = 0.0;
      project->Subcatch[j].groundwater->stats.evap = 0.0;
      project->Subcatch[j].groundwater->stats.maxFlow = 0.0;
    }
    ////
  }

  // --- allocate memory for node & link stats
  if ( project->Nobjects[LINK] > 0 )
  {
    project->NodeStats = (TNodeStats *) calloc(project->Nobjects[NODE], sizeof(TNodeStats));
    project->LinkStats = (TLinkStats *) calloc(project->Nobjects[LINK], sizeof(TLinkStats));
    if ( !project->NodeStats || !project->LinkStats )
    {
      report_writeErrorMsg(project, ERR_MEMORY, "");
      return project->ErrorCode;
    }
  }

  // --- initialize node stats
  if ( project->NodeStats ) for ( j = 0; j < project->Nobjects[NODE]; j++ )
  {
    project->NodeStats[j].avgDepth = 0.0;
    project->NodeStats[j].maxDepth = 0.0;
    project->NodeStats[j].maxDepthDate = project->StartDateTime;
    project->NodeStats[j].maxRptDepth = 0.0;                                        //(5.1.008)
    project->NodeStats[j].volFlooded = 0.0;
    project->NodeStats[j].timeFlooded = 0.0;
    project->NodeStats[j].timeSurcharged = 0.0;
    project->NodeStats[j].timeCourantCritical = 0.0;
    project->NodeStats[j].totLatFlow = 0.0;
    project->NodeStats[j].maxLatFlow = 0.0;
    project->NodeStats[j].maxInflow = 0.0;
    project->NodeStats[j].maxOverflow = 0.0;
    project->NodeStats[j].maxPondedVol = 0.0;
    project->NodeStats[j].maxInflowDate = project->StartDateTime;
    project->NodeStats[j].maxOverflowDate = project->StartDateTime;
  }

  // --- initialize link stats
  if ( project->LinkStats ) for ( j = 0; j < project->Nobjects[LINK]; j++ )
  {
    project->LinkStats[j].maxFlow = 0.0;
    project->LinkStats[j].maxVeloc = 0.0;
    project->LinkStats[j].maxDepth = 0.0;
    project->LinkStats[j].timeSurcharged = 0.0;
    project->LinkStats[j].timeFullUpstream = 0.0;
    project->LinkStats[j].timeFullDnstream = 0.0;
    project->LinkStats[j].timeFullFlow = 0.0;
    project->LinkStats[j].timeCapacityLimited = 0.0;
    project->LinkStats[j].timeCourantCritical = 0.0;
    for (k=0; k<MAX_FLOW_CLASSES; k++)
      project->LinkStats[j].timeInFlowClass[k] = 0.0;
    project->LinkStats[j].flowTurns = 0;
    project->LinkStats[j].flowTurnSign = 0;
  }

  // --- allocate memory for & initialize storage unit statistics
  if ( project->Nnodes[STORAGE] > 0 )
  {
    project->StorageStats = (TStorageStats *) calloc(project->Nnodes[STORAGE],
                                                     sizeof(TStorageStats));
    if ( !project->StorageStats )
    {
      report_writeErrorMsg(project, ERR_MEMORY, "");
      return project->ErrorCode;
    }
    else for ( k = 0; k < project->Nobjects[NODE]; k++ )
    {
      if ( project->Node[k].type != STORAGE ) continue;
      j = project->Node[k].subIndex;
      project->StorageStats[j].initVol = project->Node[k].newVolume;
      project->StorageStats[j].avgVol = 0.0;
      project->StorageStats[j].maxVol = 0.0;
      project->StorageStats[j].maxFlow = 0.0;
      project->StorageStats[j].evapLosses = 0.0;
      project->StorageStats[j].exfilLosses = 0.0;                                 //(5.1.007)
      project->StorageStats[j].maxVolDate = project->StartDateTime;
    }
  }

  // --- allocate memory for & initialize outfall statistics
  if ( project->Nnodes[OUTFALL] > 0 )
  {
    project->OutfallStats = (TOutfallStats *) calloc(project->Nnodes[OUTFALL],
                                                     sizeof(TOutfallStats));
    if ( !project->OutfallStats )
    {
      report_writeErrorMsg(project, ERR_MEMORY, "");
      return project->ErrorCode;
    }
    else for ( j = 0; j < project->Nnodes[OUTFALL]; j++ )
    {
      project->OutfallStats[j].avgFlow = 0.0;
      project->OutfallStats[j].maxFlow = 0.0;
      project->OutfallStats[j].totalPeriods = 0;
      if ( project->Nobjects[POLLUT] > 0 )
      {
        project->OutfallStats[j].totalLoad =
            (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
        if ( !project->OutfallStats[j].totalLoad )
        {
          report_writeErrorMsg(project, ERR_MEMORY, "");
          return project->ErrorCode;
        }
        for (k=0; k<project->Nobjects[POLLUT]; k++)
          project->OutfallStats[j].totalLoad[k] = 0.0;
      }
      else project->OutfallStats[j].totalLoad = NULL;
    }
  }

  // --- allocate memory & initialize pumping statistics
  if ( project->Nlinks[PUMP] > 0 )
  {
    project->PumpStats = (TPumpStats *) calloc(project->Nlinks[PUMP], sizeof(TPumpStats));
    if ( !project->PumpStats )
    {
      report_writeErrorMsg(project, ERR_MEMORY, "");
      return project->ErrorCode;
    }
    else for ( j = 0; j < project->Nlinks[PUMP]; j++ )
    {
      project->PumpStats[j].utilized = 0.0;
      project->PumpStats[j].minFlow  = 0.0;
      project->PumpStats[j].avgFlow  = 0.0;
      project->PumpStats[j].maxFlow  = 0.0;
      project->PumpStats[j].volume   = 0.0;
      project->PumpStats[j].energy   = 0.0;
      project->PumpStats[j].startUps = 0;
      project->PumpStats[j].offCurveLow = 0.0;
      project->PumpStats[j].offCurveHigh = 0.0;
    }
  }

  // --- initialize system stats
  project->MaxRunoffFlow = 0.0;
  project->MaxOutfallFlow = 0.0;
  project->SysStats.maxTimeStep = 0.0;
  project->SysStats.minTimeStep = project->RouteStep;
  project->SysStats.avgTimeStep = 0.0;
  project->SysStats.avgStepCount = 0.0;
  project->SysStats.steadyStateCount = 0.0;
  return 0;
}

//=============================================================================

void  stats_close(Project* project)
//
//  Input:   none
//  Output:  
//  Purpose: closes the simulation statistics system.
//
{
  int j;

  FREE(project->SubcatchStats);
  FREE(project->NodeStats);
  FREE(project->LinkStats);
  FREE(project->StorageStats);
  if ( project->OutfallStats )
  {
    for ( j=0; j<project->Nnodes[OUTFALL]; j++ )
      FREE(project->OutfallStats[j].totalLoad);
    FREE(project->OutfallStats);
  }
  FREE(project->PumpStats);
}

//=============================================================================

void  stats_report(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reports simulation statistics.
//
{
  // --- report flow routing accuracy statistics
  if ( project->Nobjects[LINK] > 0 && project->RouteModel != NO_ROUTING )
  {
    stats_findMaxStats(project);
    report_writeMaxStats(project, project->MaxMassBalErrs, project->MaxCourantCrit, MAX_STATS);
    report_writeMaxFlowTurns(project, project->MaxFlowTurns, MAX_STATS);
    report_writeSysStats(project, &project->SysStats);
  }

  // --- report summary statistics
  statsrpt_writeReport(project);
}

//=============================================================================

void   stats_updateSubcatchStats(Project* project, int j, double rainVol, double runonVol,
                                 double evapVol, double infilVol,
                                 double runoffVol, double runoff)
//
//  Input:   j = subcatchment index
//           rainVol   = rainfall + snowfall volume (ft3)
//           runonVol  = runon volume from other subcatchments (ft3)
//           evapVol   = evaporation volume (ft3)
//           infilVol  = infiltration volume (ft3)
//           runoffVol = runoff volume (ft3)
//           runoff    = runoff rate (cfs)
//  Output:  none
//  Purpose: updates totals of runoff components for a specific subcatchment.
//
{
  project->SubcatchStats[j].precip += rainVol;
  project->SubcatchStats[j].runon  += runonVol;
  project->SubcatchStats[j].evap   += evapVol;
  project->SubcatchStats[j].infil  += infilVol;
  project->SubcatchStats[j].runoff += runoffVol;
  project->SubcatchStats[j].maxFlow = MAX(project->SubcatchStats[j].maxFlow, runoff);
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void  stats_updateGwaterStats(Project* project, int j, double infil, double evap, double latFlow,
                              double deepFlow, double theta, double waterTable,
                              double tStep)
{
  project->Subcatch[j].groundwater->stats.infil += infil * tStep;
  project->Subcatch[j].groundwater->stats.evap += evap * tStep;
  project->Subcatch[j].groundwater->stats.latFlow += latFlow * tStep;
  project->Subcatch[j].groundwater->stats.deepFlow += deepFlow * tStep;
  project->Subcatch[j].groundwater->stats.avgUpperMoist += theta * tStep;
  project->Subcatch[j].groundwater->stats.avgWaterTable += waterTable * tStep;
  project->Subcatch[j].groundwater->stats.finalUpperMoist = theta;
  project->Subcatch[j].groundwater->stats.finalWaterTable = waterTable;
  if ( fabs(latFlow) > fabs(project->Subcatch[j].groundwater->stats.maxFlow) )
  {
    project->Subcatch[j].groundwater->stats.maxFlow = latFlow;
  }
}

//=============================================================================

void  stats_updateMaxRunoff(Project* project)
//
//   Input:   none
//   Output:  updates global variable project->MaxRunoffFlow
//   Purpose: updates value of maximum system runoff rate.
//
{
  int j;
  double sysRunoff = 0.0;

  for (j=0; j<project->Nobjects[SUBCATCH]; j++) sysRunoff += project->Subcatch[j].newRunoff;
  project->MaxRunoffFlow = MAX(project->MaxRunoffFlow, sysRunoff);
}    

//=============================================================================

////  New function added for release 5.1.008.  ////                            //(5.1.008)

void   stats_updateMaxNodeDepth(Project* project, int j, double depth)
//
//   Input:   j = node index
//            depth = water depth at node at current reporting time (ft)
//   Output:  none
//   Purpose: updates a node's maximum depth recorded at reporting times.
//
{
  if ( project->NodeStats != NULL )
    project->NodeStats[j].maxRptDepth = MAX(project->NodeStats[j].maxRptDepth, depth);
}

//=============================================================================

void   stats_updateFlowStats(Project* project, double tStep, DateTime aDate, int stepCount,
                             int steadyState)
//
//  Input:   tStep = routing time step (sec)
//           aDate = current date/time
//           stepCount = # steps required to solve routing at current time period
//           steadyState = TRUE if steady flow conditions exist
//  Output:  none
//  Purpose: updates various flow routing statistics at current time period.
//
{
  int   j;

  // --- update stats only after reporting period begins
  if ( aDate < project->ReportStart ) return;
  project->SysOutfallFlow = 0.0;

  // --- update node & link stats

#ifdef USE_OPENMP
#pragma omp parallel num_threads(project->NumThreads)                                   //(5.1.008)
#endif
  {
#ifdef USE_OPENMP
#pragma omp for                                                            //(5.1.008)
#endif
    for ( j=0; j<project->Nobjects[NODE]; j++ )
      stats_updateNodeStats(project, j, tStep, aDate);

#ifdef USE_OPENMP
#pragma omp for                                                            //(5.1.008)
#endif
    for ( j=0; j<project->Nobjects[LINK]; j++ )
      stats_updateLinkStats(project, j, tStep, aDate);
  }

  // --- update time step stats
  //     (skip initial time step for min. value)
  if ( project->OldRoutingTime > 0 )                                                  //(5.1.008)
  {
    project->SysStats.minTimeStep = MIN(project->SysStats.minTimeStep, tStep);
  }
  project->SysStats.avgTimeStep += tStep;
  project->SysStats.maxTimeStep = MAX(project->SysStats.maxTimeStep, tStep);

  // --- update iteration step count stats
  project->SysStats.avgStepCount += stepCount;

  // --- update count of times in steady state
  project->SysStats.steadyStateCount += steadyState;

  // --- update max. system outfall flow
  project->MaxOutfallFlow = MAX(project->MaxOutfallFlow, project->SysOutfallFlow);
}

//=============================================================================

void stats_updateCriticalTimeCount(Project* project, int node, int link)
//
//  Input:   node = node index
//           link = link index
//  Output:  none
//  Purpose: updates count of times a node or link was time step-critical.
//
{
  if      ( node >= 0 ) project->NodeStats[node].timeCourantCritical += 1.0;
  else if ( link >= 0 ) project->LinkStats[link].timeCourantCritical += 1.0;
}

//=============================================================================

////  Function modified for release 5.1.008.  ////                             //(5.1.008)

void stats_updateNodeStats(Project* project, int j, double tStep, DateTime aDate)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//           aDate = current date/time
//  Output:  none
//  Purpose: updates flow statistics for a node.
//
{
  int    k, p;
  double newVolume = project->Node[j].newVolume;
  double newDepth = project->Node[j].newDepth;
  int    canPond = (project->AllowPonding && project->Node[j].pondedArea > 0.0);

  // --- update depth statistics
  project->NodeStats[j].avgDepth += newDepth;
  if ( newDepth > project->NodeStats[j].maxDepth )
  {
    project->NodeStats[j].maxDepth = newDepth;
    project->NodeStats[j].maxDepthDate = aDate;
  }

  // --- update flooding, ponding, and surcharge statistics
  if ( project->Node[j].type != OUTFALL )
  {
    if ( newVolume > project->Node[j].fullVolume || project->Node[j].overflow > 0.0 )
    {
      project->NodeStats[j].timeFlooded += tStep;
      project->NodeStats[j].volFlooded += project->Node[j].overflow * tStep;
      if ( canPond ) project->NodeStats[j].maxPondedVol =
          MAX(project->NodeStats[j].maxPondedVol,
              (newVolume - project->Node[j].fullVolume));
    }
    if ( newDepth + project->Node[j].invertElev + FUDGE >= project->Node[j].crownElev )
    {
      project->NodeStats[j].timeSurcharged += tStep;
    }
  }

  // --- update storage statistics
  if ( project->Node[j].type == STORAGE )
  {
    k = project->Node[j].subIndex;
    project->StorageStats[k].avgVol += newVolume;
    project->StorageStats[k].evapLosses +=
        project->Storage[project->Node[j].subIndex].evapLoss;
    project->StorageStats[k].exfilLosses +=
        project->Storage[project->Node[j].subIndex].exfilLoss;

    newVolume = MIN(newVolume, project->Node[j].fullVolume);
    if ( newVolume > project->StorageStats[k].maxVol )
    {
      project->StorageStats[k].maxVol = newVolume;
      project->StorageStats[k].maxVolDate = aDate;
    }
    project->StorageStats[k].maxFlow = MAX(project->StorageStats[k].maxFlow, project->Node[j].outflow);
  }

  // --- update outfall statistics
  if ( project->Node[j].type == OUTFALL )
  {
    k = project->Node[j].subIndex;
    if ( project->Node[j].inflow >= MIN_RUNOFF_FLOW )
    {
      project->OutfallStats[k].avgFlow += project->Node[j].inflow;
      project->OutfallStats[k].maxFlow = MAX(project->OutfallStats[k].maxFlow, project->Node[j].inflow);
      project->OutfallStats[k].totalPeriods++;
    }
    for (p=0; p<project->Nobjects[POLLUT]; p++)
    {
      project->OutfallStats[k].totalLoad[p] += project->Node[j].inflow *
                                               project->Node[j].newQual[p] * tStep;
    }
    project->SysOutfallFlow += project->Node[j].inflow;
  }

  // --- update inflow statistics
  project->NodeStats[j].totLatFlow += ( (project->Node[j].oldLatFlow + project->Node[j].newLatFlow) *
                                        0.5 * tStep );
  if ( fabs(project->Node[j].newLatFlow) > fabs(project->NodeStats[j].maxLatFlow) )
    project->NodeStats[j].maxLatFlow = project->Node[j].newLatFlow;
  if ( project->Node[j].inflow > project->NodeStats[j].maxInflow )
  {
    project->NodeStats[j].maxInflow = project->Node[j].inflow;
    project->NodeStats[j].maxInflowDate = aDate;
  }

  // --- update overflow statistics
  if ( project->Node[j].overflow > project->NodeStats[j].maxOverflow )
  {
    project->NodeStats[j].maxOverflow = project->Node[j].overflow;
    project->NodeStats[j].maxOverflowDate = aDate;
  }
}

//=============================================================================

void  stats_updateLinkStats(Project* project, int j, double tStep, DateTime aDate)
//
//  Input:   j = link index
//           tStep = routing time step (sec)
//           aDate = current date/time
//  Output:  none
//  Purpose: updates flow statistics for a link.
//
{
  int    k;
  double q, v;
  double dq;

  // --- update max. flow
  dq = project->Link[j].newFlow - project->Link[j].oldFlow;
  q = fabs(project->Link[j].newFlow);
  if ( q > project->LinkStats[j].maxFlow )
  {
    project->LinkStats[j].maxFlow = q;
    project->LinkStats[j].maxFlowDate = aDate;
  }

  // --- update max. velocity
  v = link_getVelocity(project, j, q, project->Link[j].newDepth);
  if ( v > project->LinkStats[j].maxVeloc )
  {
    project->LinkStats[j].maxVeloc = v;
    //project->LinkStats[j].maxVelocDate = aDate;                                   //(5.1.008)
  }

  // --- update max. depth
  if ( project->Link[j].newDepth > project->LinkStats[j].maxDepth )
  {
    project->LinkStats[j].maxDepth = project->Link[j].newDepth;
  }

  if ( project->Link[j].type == PUMP )
  {
    if ( q >= project->Link[j].qFull )
      project->LinkStats[j].timeFullFlow += tStep;
    if ( q > MIN_RUNOFF_FLOW )
    {
      k = project->Link[j].subIndex;
      project->PumpStats[k].minFlow = MIN(project->PumpStats[k].minFlow, q);
      project->PumpStats[k].maxFlow = project->LinkStats[j].maxFlow;
      project->PumpStats[k].avgFlow += q;
      project->PumpStats[k].volume += q*tStep;
      project->PumpStats[k].utilized += tStep;
      project->PumpStats[k].energy += link_getPower(project, j)*tStep / 3600.0;
      if ( project->Link[j].flowClass == DN_DRY )
        project->PumpStats[k].offCurveLow += tStep;
      if ( project->Link[j].flowClass == UP_DRY )
        project->PumpStats[k].offCurveHigh += tStep;
      if ( project->Link[j].oldFlow < MIN_RUNOFF_FLOW )
        project->PumpStats[k].startUps++;
      project->PumpStats[k].totalPeriods++;
      project->LinkStats[j].timeSurcharged += tStep;
      project->LinkStats[j].timeFullUpstream += tStep;
      project->LinkStats[j].timeFullDnstream += tStep;
    }
  }
  else if ( project->Link[j].type == CONDUIT )
  {

    // --- update time under normal flow & inlet control
    if ( project->Link[j].normalFlow ) project->LinkStats[j].timeNormalFlow += tStep;
    if ( project->Link[j].inletControl ) project->LinkStats[j].timeInletControl += tStep;
    
    // --- update flow classification distribution
    k = project->Link[j].flowClass;
    if ( k >= 0 && k < MAX_FLOW_CLASSES )
    {
      ++project->LinkStats[j].timeInFlowClass[k];
    }

    // --- update time conduit is full
    k = project->Link[j].subIndex;
    if ( q >= project->Link[j].qFull ) project->LinkStats[j].timeFullFlow += tStep;
    if ( project->Conduit[k].capacityLimited )
      project->LinkStats[j].timeCapacityLimited += tStep;

    ////  Following section modified for release 5.1.008.  ////                    //(5.1.008)
    ////
    switch (project->Conduit[k].fullState)
    {
      case ALL_FULL:
        project->LinkStats[j].timeSurcharged += tStep;
        project->LinkStats[j].timeFullUpstream += tStep;
        project->LinkStats[j].timeFullDnstream += tStep;
        break;
      case UP_FULL:
        project->LinkStats[j].timeFullUpstream += tStep;
        break;
      case DN_FULL:
        project->LinkStats[j].timeFullDnstream += tStep;
    }
    ////
  }

  // --- update flow turn count
  k = project->LinkStats[j].flowTurnSign;
  project->LinkStats[j].flowTurnSign = SGN(dq);
  if ( fabs(dq) > 0.001 &&  k * project->LinkStats[j].flowTurnSign < 0 )
    project->LinkStats[j].flowTurns++;
}

//=============================================================================

void  stats_findMaxStats(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: finds nodes & links with highest mass balance errors
//           & highest times Courant time-step critical.
//
{
  int    j;
  double x;

  // --- initialize max. stats arrays
  for (j=0; j<MAX_STATS; j++)
  {
    project->MaxMassBalErrs[j].objType = NODE;
    project->MaxMassBalErrs[j].index   = -1;
    project->MaxMassBalErrs[j].value   = -1.0;
    project->MaxCourantCrit[j].index   = -1;
    project->MaxCourantCrit[j].value   = -1.0;
    project->MaxFlowTurns[j].index     = -1;
    project->MaxFlowTurns[j].value     = -1.0;
  }

  // --- find links with most flow turns
  if ( project->StepCount > 2 )
  {
    for (j=0; j<project->Nobjects[LINK]; j++)
    {
      x = 100.0 * project->LinkStats[j].flowTurns / (2./3.*(project->StepCount-2));
      stats_updateMaxStats(project->MaxFlowTurns, LINK, j, x);
    }
  }

  // --- find nodes with largest mass balance errors
  for (j=0; j<project->Nobjects[NODE]; j++)
  {
    // --- skip terminal nodes and nodes with negligible inflow
    if ( project->Node[j].degree <= 0  ) continue;
    if ( project->NodeInflow[j] <= 0.1 ) continue;

    // --- evaluate mass balance error
    //     (Note: NodeInflow & NodeOutflow include any initial and final
    //            stored volumes, respectively).
    if ( project->NodeInflow[j]  > 0.0 )
      x = 1.0 - project->NodeOutflow[j] / project->NodeInflow[j];
    else if ( project->NodeOutflow[j] > 0.0 ) x = -1.0;
    else                             x = 0.0;
    stats_updateMaxStats(project->MaxMassBalErrs, NODE, j, 100.0*x);
  }

  // --- stop if not using a variable time step
  if ( project->RouteModel != DW || project->CourantFactor == 0.0 ) return;

  // --- find nodes most frequently Courant critical
  if ( project->StepCount == 0 ) return;                                              //(5.1.008)
  for (j=0; j<project->Nobjects[NODE]; j++)
  {
    x = project->NodeStats[j].timeCourantCritical / project->StepCount;
    stats_updateMaxStats(project->MaxCourantCrit, NODE, j, 100.0*x);
  }

  // --- find links most frequently Courant critical
  for (j=0; j<project->Nobjects[LINK]; j++)
  {
    x = project->LinkStats[j].timeCourantCritical / project->StepCount;
    stats_updateMaxStats(project->MaxCourantCrit, LINK, j, 100.0*x);
  }
}

//=============================================================================

void  stats_updateMaxStats(TMaxStats maxStats[], int i, int j, double x)
//
//  Input:   maxStats[] = array of critical statistics values
//           i = object category (NODE or LINK)
//           j = object index
//           x = value of statistic for the object
//  Output:  none
//  Purpose: updates the collection of most critical statistics
//
{
  int   k;
  TMaxStats maxStats1, maxStats2;
  maxStats1.objType = i;
  maxStats1.index   = j;
  maxStats1.value   = x;
  for (k=0; k<MAX_STATS; k++)
  {
    if ( fabs(maxStats1.value) > fabs(maxStats[k].value) )
    {
      maxStats2 = maxStats[k];
      maxStats[k] = maxStats1;
      maxStats1 = maxStats2;
    }
  }
}

//=============================================================================
