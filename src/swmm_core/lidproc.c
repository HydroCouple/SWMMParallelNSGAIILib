//-----------------------------------------------------------------------------
//   lidproc.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/12   (Build 5.1.001)
//             05/19/14   (Build 5.1.006)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//             04/30/15   (Build 5.1.009)
//             08/05/15   (Build 5.1.010)
//   Author:   L. Rossman (US EPA)
//
//   This module computes the hydrologic performance of an LID (Low Impact
//   Development) unit at a given point in time.
//
//   Build 5.1.007:
//   - Euler integration now applied to all LID types except Vegetative
//     Swale which continues to use successive approximation.
//   - LID layer flux routines were re-written to more accurately model
//     flooded conditions.
//
//   Build 5.1.008:
//   - MAX_STATE_VARS replaced with MAX_LAYERS.
//   - Optional soil layer added to Porous Pavement LID.
//   - Rooftop Disconnection added to types of LIDs.
//   - Separate accounting of drain flows added.
//   - Indicator for currently wet LIDs added.
//   - Detailed reporting procedure fixed.
//   - Possibile negative head on Bioretention Cell drain avoided.
//   - Bug in computing flow through Green Roof drainage mat fixed.
//
//   Build 5.1.009:
//   - Fixed typo in net flux rate for vegetative swale LID.
//
//   Build 5.1.010:
//   - New modified version of Green-Ampt used for surface layer infiltration.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lid.h"
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
#define STOPTOL  0.00328     // integration error tolerance in ft (= 1 mm)
#define MINFLOW  2.3e-8      // flow cutoff for dry conditions (= 0.001 in/hr)

//-----------------------------------------------------------------------------
//  Enumerations
//-----------------------------------------------------------------------------
enum LidLayerTypes {
  SURF,                    // surface layer
  SOIL,                    // soil layer
  STOR,                    // storage layer
  PAVE,                    // pavement layer
  DRAIN};                  // underdrain system

enum LidRptVars {
  SURF_INFLOW,             // inflow to surface layer
  TOTAL_EVAP,              // evaporation rate from all layers
  SURF_INFIL,              // infiltration into surface layer
  PAVE_PERC,               // percolation through pavement layer             //(5.1.008)
  SOIL_PERC,               // percolation through soil layer
  STOR_INFIL,              // infiltration from storage layer
  SURF_OUTFLOW,            // outflow from surface layer
  STOR_DRAIN,              // outflow from storage layer
  SURF_DEPTH,              // ponded depth on surface layer
  PAVE_MOIST,              // moisture content of pavement layer             //(5.1.008)
  SOIL_MOIST,              // moisture content of soil layer
  STOR_DEPTH,              // water level in storage layer
  MAX_RPT_VARS};

////  Added to release 5.1.008.  ////                                          //(5.1.008)
//-----------------------------------------------------------------------------
//  Imported variables 
//-----------------------------------------------------------------------------
//extern char project->HasWetLids;      // TRUE if any LIDs are wet (declared in runoff.c)

//-----------------------------------------------------------------------------
//  Local Variables
//-----------------------------------------------------------------------------
//static TLidUnit*   project->theLidUnit;     // ptr. to a subcatchment's LID unit
//static TLidProc*   project->theLidProc;     // ptr. to a LID process

//static double      project->Tsteplp;          // current time step (sec)
////static double     Rainfall;       // current rainfall rate (ft/s)            //(5.1.008)
//static double      project->EvapRatelp;       // evaporation rate (ft/s)
//static double      project->MaxNativeInfillp; // native soil infil. rate limit (ft/s)

//static double      project->SurfaceInflow;  // precip. + runon to LID unit (ft/s)
//static double      project->SurfaceInfil;   // infil. rate from surface layer (ft/s)
//static double      project->SurfaceEvap;    // evap. rate from surface layer (ft/s)
//static double      project->SurfaceOutflow; // outflow from surface layer (ft/s)
//static double      project->SurfaceVolume;  // volume in surface storage (ft)

//static double      project->PaveEvap;       // evap. from pavement layer (ft/s)          //(5.1.008)
//static double      project->PavePerc;       // percolation from pavement layer (ft/s)    //(5.1.008)
//static double      project->PaveVolume;     // volume stored in pavement layer  (ft)     //(5.1.008)

//static double      project->SoilEvap;       // evap. from soil layer (ft/s)
//static double      project->SoilPerc;       // percolation from soil layer (ft/s)
//static double      project->SoilVolume;     // volume in soil/pavement storage (ft)

//static double      project->StorageInflow;  // inflow rate to storage layer (ft/s)
//static double      project->StorageInfil;   // infil. rate from storage layer (ft/s)
//static double      project->StorageEvap;    // evap.rate from storage layer (ft/s)
//static double      project->StorageDrain;   // underdrain flow rate layer (ft/s)
//static double      project->StorageVolume;  // volume in storage layer (ft)

//static double      project->Xold[MAX_LAYERS];  // previous moisture level in LID layers  //(5.1.008)

//-----------------------------------------------------------------------------
//  External Functions (declared in lid.h)
//-----------------------------------------------------------------------------
// lidproc_initWaterBalance  (called by lid_initState)
// lidproc_getOutflow        (called by evalLidUnit in lid.c)
// lidproc_saveResults       (called by evalLidUnit in lid.c)

//-----------------------------------------------------------------------------
// Local Functions
//-----------------------------------------------------------------------------
static void   barrelFluxRates(Project* project, double x[], double f[]);
static void   biocellFluxRates(Project* project, double x[], double f[]);
static void   greenRoofFluxRates(Project* project, double x[], double f[]);
static void   pavementFluxRates(Project* project, double x[], double f[]);
static void   trenchFluxRates(Project* project, double x[], double f[]);
static void   swaleFluxRates(Project* project, double x[], double f[]);
static void   roofFluxRates(Project* project, double x[], double f[]);                           //(5.1.008)

static double getSurfaceOutflowRate(Project* project,double depth);
static double getSurfaceOverflowRate(Project* project,double* surfaceDepth);
static double getPavementPermRate(Project* project);
static double getSoilPercRate(Project* project, double theta);                                   //(5.1.007)
static double getStorageInfilRate(Project* project);
static double getStorageDrainRate(Project* project, double head);
static double getDrainMatOutflow(Project *project, double depth);
static void   getEvapRates(Project* project,double surfaceVol, double paveVol,                  //(5.1.008)
                           double soilVol, double storageVol);


static void   updateWaterBalance(Project *project, TLidUnit *lidUnit, double inflow,
                                 double evap, double infil, double surfFlow,
                                 double drainFlow, double storage);

static int    modpuls_solve(Project* project, int n, double* x, double* xOld, double* xPrev,
                            double* xMin, double* xMax, double* xTol,
                            double* qOld, double* q, double dt, double omega,  //(5.1.007)
                            void (*derivs)(Project*, double*, double*));


//=============================================================================

void lidproc_initWaterBalance(TLidUnit *lidUnit, double initVol)
//
//  Purpose: initializes the water balance components of a LID unit.
//  Input:   lidUnit = a particular LID unit
//           initVol = initial water volume stored in the unit (ft)
//  Output:  none
//
{
  lidUnit->waterBalance.inflow = 0.0;
  lidUnit->waterBalance.evap = 0.0;
  lidUnit->waterBalance.infil = 0.0;
  lidUnit->waterBalance.surfFlow = 0.0;
  lidUnit->waterBalance.drainFlow = 0.0;
  lidUnit->waterBalance.initVol = initVol;
  lidUnit->waterBalance.finalVol = initVol;                                  //(5.1.008)
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double lidproc_getOutflow(Project* project, TLidUnit* lidUnit, TLidProc* lidProc, double inflow,
                          double evap, double infil, double maxInfil,
                          double tStep, double* lidEvap,
                          double* lidInfil, double* lidDrain)
//
//  Purpose: computes runoff outflow from a single LID unit.
//  Input:   lidUnit  = ptr. to specific LID unit being analyzed
//           lidProc  = ptr. to generic LID process of the LID unit
//           inflow   = runoff rate captured by LID unit (ft/s)
//           evap     = potential evaporation rate (ft/s)
//           infil    = infiltration rate of native soil (ft/s)
//           maxInfil = max. infiltration rate to native soil (ft/s)
//           tStep    = time step (sec)
//  Output:  lidEvap  = evaporation rate for LID unit (ft/s)
//           lidInfil = infiltration rate for LID unit (ft/s)
//           lidDrain = drain flow for LID unit (ft/s)
//           returns surface runoff rate from the LID unit (ft/s)
//
{
  int    i;
  double x[MAX_LAYERS];        // layer moisture levels
  double xOld[MAX_LAYERS];     // work vector
  double xPrev[MAX_LAYERS];    // work vector
  double xMin[MAX_LAYERS];     // lower limit on moisture levels
  double xMax[MAX_LAYERS];     // upper limit on moisture levels
  double fOld[MAX_LAYERS];     // previously computed flux rates
  double f[MAX_LAYERS];        // newly computed flux rates

  // convergence tolerance on moisture levels (ft, moisture fraction , ft)
  double xTol[MAX_LAYERS] = {STOPTOL, STOPTOL, STOPTOL, STOPTOL};

  double omega = 0.0;          // integration time weighting

  //... define a pointer to function that computes flux rates through the LID
  void (*fluxRates) (Project*, double *, double *) = NULL;

  //... save references to the LID process and LID unit
  project->theLidProc = lidProc;
  project->theLidUnit = lidUnit;

  //... save evap, max. infil. & time step to shared variables
  project->EvapRatelp = evap;
  project->MaxNativeInfillp = maxInfil;
  project->Tsteplp = tStep;

  //... store current moisture levels in vector x
  x[SURF] =  project->theLidUnit->surfaceDepth;
  x[SOIL] =  project->theLidUnit->soilMoisture;
  x[STOR] =  project->theLidUnit->storageDepth;
  x[PAVE] =  project->theLidUnit->paveMoisture;

  //... initialize layer flux rates and moisture limits
  project->SurfaceInflow  = inflow;
  project->SurfaceInfil   = 0.0;
  project->SurfaceEvap    = 0.0;
  project->SurfaceOutflow = 0.0;
  project->PaveEvap       = 0.0;
  project->PavePerc       = 0.0;
  project->SoilEvap       = 0.0;
  project->SoilPerc       = 0.0;
  project->StorageInflow  = 0.0;
  project->StorageInfil   = 0.0;
  project->StorageEvap    = 0.0;
  project->StorageDrain   = 0.0;

  for (i = 0; i < MAX_LAYERS; i++)
  {
    f[i] = 0.0;
    fOld[i] =  project->theLidUnit->oldFluxRates[i];
    xMin[i] = 0.0;
    xMax[i] = BIG;
    project->Xold[i] = x[i];
  }

  //... find Green-Ampt infiltration from surface layer
  if (  project->theLidProc->lidType == POROUS_PAVEMENT )  project->SurfaceInfil = 0.0;
  else if (  project->theLidUnit->soilInfil.Ks > 0.0 )
  {
    project->SurfaceInfil =
        grnampt_getInfil(project,& project->theLidUnit->soilInfil,  project->Tsteplp,
                         project->SurfaceInflow,  project->theLidUnit->surfaceDepth,
                         MOD_GREEN_AMPT);                                  //(5.1.010)
  }
  else  project->SurfaceInfil = infil;

  //... set moisture limits for soil & storage layers
  if (  project->theLidProc->soil.thickness > 0.0 )
  {
    xMin[SOIL] =  project->theLidProc->soil.wiltPoint;
    xMax[SOIL] =  project->theLidProc->soil.porosity;
  }
  if (  project->theLidProc->pavement.thickness > 0.0 )
  {
    xMax[PAVE] =  project->theLidProc->pavement.voidFrac;
  }
  if (  project->theLidProc->storage.thickness > 0.0 )
  {
    xMax[STOR] =  project->theLidProc->storage.thickness;
  }
  if (  project->theLidProc->lidType == GREEN_ROOF )
  {
    xMax[STOR] =  project->theLidProc->drainMat.thickness;
  }

  //... determine which flux rate function to use
  switch ( project->theLidProc->lidType)
  {
    case BIO_CELL:
    case RAIN_GARDEN:     fluxRates = &biocellFluxRates;  break;
    case GREEN_ROOF:      fluxRates = &greenRoofFluxRates; break;
    case INFIL_TRENCH:    fluxRates = &trenchFluxRates;   break;
    case POROUS_PAVEMENT: fluxRates = &pavementFluxRates; break;
    case RAIN_BARREL:     fluxRates = &barrelFluxRates;   break;
    case ROOF_DISCON:     fluxRates = &roofFluxRates;     break;
    case VEG_SWALE:       fluxRates = &swaleFluxRates;
      omega = 0.5;
      break;
    default:              return 0.0;
  }

  //... update moisture levels and flux rates over the time step
  i = modpuls_solve(project,MAX_LAYERS, x, xOld, xPrev, xMin, xMax, xTol,
                    fOld, f, tStep, omega, fluxRates);

  /** For debugging only ********************************************
    if  (i == 0)
    {
        fprintf(Frpt.file,
        "\n  WARNING 09: integration failed to converge at %s %s",
            theDate, theTime);
        fprintf(Frpt.file,
        "\n              for LID %s placed in subcatchment %s.",
             project->theLidProc->ID, theSubcatch->ID);
    }
*******************************************************************/

  //... add any surface overflow to surface outflow
  if (  project->theLidProc->surface.canOverflow ||  project->theLidUnit->fullWidth == 0.0 )
  {
    project->SurfaceOutflow += getSurfaceOverflowRate(project, &x[SURF]);
  }

  //... save updated results
  project->theLidUnit->surfaceDepth = x[SURF];
  project->theLidUnit->paveMoisture = x[PAVE];
  project->theLidUnit->soilMoisture = x[SOIL];
  project->theLidUnit->storageDepth = x[STOR];
  for (i = 0; i < MAX_LAYERS; i++)  project->theLidUnit->oldFluxRates[i] = f[i];

  //... assign values to LID unit evaporation, infiltration & drain flow
  *lidEvap =  project->SurfaceEvap +  project->PaveEvap +  project->SoilEvap +  project->StorageEvap;
  *lidInfil =  project->StorageInfil;
  *lidDrain =  project->StorageDrain;

  //... return surface outflow (per unit area) from unit
  return  project->SurfaceOutflow;
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void lidproc_saveResults(Project* project, TLidUnit* lidUnit, TLidProc* lidProc, 
                         double ucfRainfall, double ucfRainDepth)
//
//  Purpose: updates the mass balance for an LID unit and saves
//           current flux rates to the LID report file.
//  Input:   lidUnit = ptr. to LID unit
//           lidProc = ptr. to LID process
//           ucfRainfall = units conversion factor for rainfall rate
//           ucfDepth = units conversion factor for rainfall depth
//  Output:  none
//
{
  double ucf;                        // units conversion factor
  double totalEvap;                  // total evaporation rate (ft/s)
  double totalVolume;                // total volume stored in LID (ft)
  double rptVars[MAX_RPT_VARS];      // array of reporting variables
  int    isDry = FALSE;              // true if current state of LID is dry
  double perc,                       // percolation rate (ft/s)
      moist;                      // moisture content

  //... find total evap. rate and stored volume
  totalEvap =  project->SurfaceEvap +  project->PaveEvap +  project->SoilEvap +  project->StorageEvap;
  totalVolume =  project->SurfaceVolume +  project->PaveVolume +  project->SoilVolume +  project->StorageVolume;

  //... update mass balance totals
  updateWaterBalance(project, project->theLidUnit,  project->SurfaceInflow, totalEvap,  project->StorageInfil,
                      project->SurfaceOutflow,  project->StorageDrain, totalVolume);

  //... check if dry-weather conditions hold
  if (  project->SurfaceInflow  < MINFLOW &&
        project->SurfaceOutflow < MINFLOW &&
        project->StorageDrain   < MINFLOW &&
        project->StorageInfil   < MINFLOW &&
        totalEvap      < MINFLOW ) isDry = TRUE;

  //... update status of project->HasWetLids
  if ( !isDry ) project->HasWetLids = TRUE;

  //... write results to LID report file
  if ( lidUnit->rptFile )
  {
    //... convert rate results to original units (in/hr or mm/hr)
    ucf = ucfRainfall;
    rptVars[SURF_INFLOW]  =  project->SurfaceInflow*ucf;
    rptVars[TOTAL_EVAP]   = totalEvap*ucf;
    rptVars[SURF_INFIL]   =  project->SurfaceInfil*ucf;
    rptVars[PAVE_PERC]    =  project->PavePerc*ucf;
    rptVars[SOIL_PERC]    =  project->SoilPerc*ucf;
    rptVars[STOR_INFIL]   =  project->StorageInfil*ucf;
    rptVars[SURF_OUTFLOW] =  project->SurfaceOutflow*ucf;
    rptVars[STOR_DRAIN]   =  project->StorageDrain*ucf;

    //... convert storage results to original units (in or mm)
    ucf = ucfRainDepth;
    rptVars[SURF_DEPTH] =  project->theLidUnit->surfaceDepth*ucf;
    rptVars[PAVE_MOIST] =  project->theLidUnit->paveMoisture;
    rptVars[SOIL_MOIST] =  project->theLidUnit->soilMoisture;
    rptVars[STOR_DEPTH] =  project->theLidUnit->storageDepth*ucf;

    //... if the current LID state is wet but the previous state was dry
    //    then write the saved previous results to the report file thus
    //    marking the end of a dry period

    if ( !isDry &&  project->theLidUnit->rptFile->wasDry )
      fprintf( project->theLidUnit->rptFile->file, "%s",
               project->theLidUnit->rptFile->results);

    //... write the current results to a string which is saved between
    //    reporting periods
    perc = rptVars[SOIL_PERC];
    moist = rptVars[SOIL_MOIST];
    if ( lidProc->lidType == POROUS_PAVEMENT &&
         lidProc->soil.thickness == 0.0 )
    {
      perc = rptVars[PAVE_PERC];
      moist = rptVars[PAVE_MOIST];
    }
    sprintf( project->theLidUnit->rptFile->results,
             "\n%7.3f\t %8.2f\t %8.4f\t %8.2f\t %8.2f\t %8.2f\t %8.2f\t"
             "%8.2f\t %8.2f\t %8.2f\t %8.2f\t",
             project->NewRunoffTime/1000.0/3600.0, rptVars[SURF_INFLOW],
             rptVars[TOTAL_EVAP], rptVars[SURF_INFIL], perc,
             rptVars[STOR_INFIL], rptVars[SURF_OUTFLOW], rptVars[STOR_DRAIN],
             rptVars[SURF_DEPTH], moist, rptVars[STOR_DEPTH]);

    //... if the current LID state is dry
    if ( isDry )
    {
      //... if the previous state was wet then write the current
      //    results to file marking the start of a dry period
      if ( ! project->theLidUnit->rptFile->wasDry )
      {
        fprintf( project->theLidUnit->rptFile->file, "%s",  project->theLidUnit->rptFile->results);
        project->theLidUnit->rptFile->wasDry = TRUE;
      }
    }

    //... if the current LID state is wet
    else
    {
      //... if the previous state was dry then make it wet
      if (  project->theLidUnit->rptFile->wasDry )
      {
        project->theLidUnit->rptFile->wasDry = FALSE;
      }

      //... write the current results to the report file
      fprintf( project->theLidUnit->rptFile->file, "%s",  project->theLidUnit->rptFile->results);
    }
  }
}

//=============================================================================

////  New function for release 5.1.008.  ////                                  //(5.1.008)

void roofFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates for roof disconnection.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double surfaceDepth = x[SURF];

  getEvapRates(project, surfaceDepth, 0.0, 0.0, 0.0);
  project->SurfaceVolume = surfaceDepth;
  project->SurfaceInfil = 0.0;
  if (  project->theLidProc->surface.alpha > 0.0 )
    project->SurfaceOutflow = getSurfaceOutflowRate(project, surfaceDepth);
  else getSurfaceOverflowRate(project,&surfaceDepth);
  project->StorageDrain = MIN( project->theLidProc->drain.coeff/UCF(project,RAINFALL),  project->SurfaceOutflow);
  project->SurfaceOutflow -=  project->StorageDrain;
  f[SURF] = ( project->SurfaceInflow -  project->SurfaceEvap -  project->StorageDrain -  project->SurfaceOutflow);
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void greenRoofFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a green roof.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double surfaceDepth;
  double soilTheta;
  double storageDepth;
  double availVolume;
  double maxRate;

  //... retrieve state variables from work vector
  surfaceDepth = x[SURF];
  soilTheta = x[SOIL];
  storageDepth = x[STOR];

  //... convert state variables to volumes
  project->SurfaceVolume = surfaceDepth *  project->theLidProc->surface.voidFrac;
  project->SoilVolume = soilTheta *  project->theLidProc->soil.thickness;
  project->StorageVolume = storageDepth *  project->theLidProc->storage.voidFrac;

  //... get ET rates
  availVolume =  project->SoilVolume -  project->theLidProc->soil.wiltPoint *
                 project->theLidProc->soil.thickness;
  getEvapRates(project, project->SurfaceVolume, 0.0, availVolume,  project->StorageVolume);              //(5.1.008)

  //... no storage evap if soil layer saturated
  if ( soilTheta >=  project->theLidProc->soil.porosity )  project->StorageEvap = 0.0;

  //... find soil layer perc rate
  project->SoilPerc = getSoilPercRate(project, soilTheta);

  //... find storage (drain mat) outflow rate
  project->StorageInfil = 0.0;
  project->StorageDrain = getDrainMatOutflow(project, storageDepth);

  //... both storage & soil layers are saturated
  if ( storageDepth >=  project->theLidProc->storage.thickness &&
       soilTheta >=  project->theLidProc->soil.porosity )
  {
    //... soil perc can't exceed storage outflow
    if (  project->SoilPerc >  project->StorageDrain )  project->SoilPerc =  project->StorageDrain;

    //... storage outflow can't exceed soil perc
    else  project->StorageDrain = MIN( project->StorageDrain,  project->SoilPerc);
  }

  //... storage and/or soil layers not saturated
  else
  {
    //... limit underdrain flow by volume above drain offset
    if (  project->StorageDrain > 0.0 )
    {
      maxRate = (storageDepth -  project->theLidProc->drain.offset) *
                project->theLidProc->storage.voidFrac /  project->Tsteplp;
      project->StorageDrain = MIN( project->StorageDrain, maxRate);
    }

    //... limit soil perc by available storage volume
    availVolume = ( project->theLidProc->storage.thickness - storageDepth) *
                  project->theLidProc->storage.voidFrac;
    maxRate = availVolume/ project->Tsteplp +  project->StorageEvap +  project->StorageDrain;
    project->SoilPerc = MIN( project->SoilPerc, maxRate);
  }

  //... limit surface infil. by available soil pore volume
  maxRate = ( project->theLidProc->soil.porosity - soilTheta) *
            project->theLidProc->soil.thickness /  project->Tsteplp +  project->SoilPerc;
  project->SurfaceInfil = MIN( project->SurfaceInfil, maxRate);

  // ... find surface outflow rate
  project->SurfaceOutflow = getSurfaceOutflowRate(project, surfaceDepth);

  // ... find net fluxes for each layer
  f[SURF] = ( project->SurfaceInflow -  project->SurfaceEvap -  project->SurfaceInfil -  project->SurfaceOutflow) /
            project->theLidProc->surface.voidFrac;
  f[SOIL] = ( project->SurfaceInfil -  project->SoilEvap -  project->SoilPerc) /
            project->theLidProc->soil.thickness;
  f[STOR] = ( project->SoilPerc -  project->StorageEvap -  project->StorageDrain) /
            project->theLidProc->storage.voidFrac;
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void biocellFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a bio-retention cell LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double surfaceDepth;
  double soilTheta;
  double storageDepth;
  double head;
  double availVolume;
  double maxRate;

  //... retrieve state variables from work vector
  surfaceDepth = x[SURF];
  soilTheta = x[SOIL];
  storageDepth = x[STOR];

  //... convert state variables to volumes
  project->SurfaceVolume = surfaceDepth *  project->theLidProc->surface.voidFrac;
  project->SoilVolume = soilTheta *  project->theLidProc->soil.thickness;
  project->StorageVolume = storageDepth *  project->theLidProc->storage.voidFrac;

  //... get ET rates
  availVolume =  project->SoilVolume -  project->theLidProc->soil.wiltPoint *
                 project->theLidProc->soil.thickness;
  getEvapRates(project, project->SurfaceVolume, 0.0, availVolume,  project->StorageVolume);

  //... no storage evap if soil layer saturated
  if ( soilTheta >=  project->theLidProc->soil.porosity )  project->StorageEvap = 0.0;

  //... find soil layer perc rate
  project->SoilPerc = getSoilPercRate(project, soilTheta);

  //... find infiltration rate out of storage layer
  project->StorageInfil = getStorageInfilRate(project);

  //... find underdrain flow rate
  project->StorageDrain = 0.0;
  head = storageDepth -  project->theLidProc->drain.offset;
  if (  project->theLidProc->drain.coeff > 0.0 && head >= 0.0 )
  {
    if ( storageDepth >=  project->theLidProc->storage.thickness )
    {
      if ( soilTheta >  project->theLidProc->soil.fieldCap )                       //(5.1.008)
      {
        head += (soilTheta -  project->theLidProc->soil.fieldCap) /
                ( project->theLidProc->soil.porosity -  project->theLidProc->soil.fieldCap) *
                project->theLidProc->soil.thickness;
      }
      if ( soilTheta >=  project->theLidProc->soil.porosity ) head += surfaceDepth;
    }
    project->StorageDrain =  getStorageDrainRate(project, head);
  }

  //... special case of no storage layer present
  if (  project->theLidProc->storage.thickness == 0.0 )
  {
    project->StorageEvap = 0.0;
    maxRate = MIN( project->StorageInfil,  project->SoilPerc);
    project->SoilPerc = maxRate;
    project->StorageInfil = maxRate;
  }

  //... both storage & soil layers are saturated
  else if ( storageDepth >=  project->theLidProc->storage.thickness &&
            soilTheta >=  project->theLidProc->soil.porosity )
  {
    //... soil perc can't exceed storage outflow
    maxRate =  project->StorageDrain +  project->StorageInfil;
    if (  project->SoilPerc > maxRate )  project->SoilPerc = maxRate;

    //... storage outflow can't exceed soil perc
    else
    {
      //... use up available drain capacity first
      project->StorageDrain = MIN( project->StorageDrain,  project->SoilPerc);
      project->StorageInfil =  project->SoilPerc -  project->StorageDrain;
    }
  }

  //... layers not saturated
  else
  {
    //... limit underdrain flow by volume above drain offset
    if (  project->StorageDrain > 0.0 )
    {
      maxRate = (storageDepth -  project->theLidProc->drain.offset) *
                project->theLidProc->storage.voidFrac /  project->Tsteplp;
      project->StorageDrain = MIN( project->StorageDrain, maxRate);
    }

    //... limit storage infil. by remaining volume
    maxRate =  project->StorageVolume /  project->Tsteplp -  project->StorageDrain -  project->StorageEvap;
    maxRate = MAX(0.0, maxRate);
    project->StorageInfil = MIN( project->StorageInfil, maxRate);

    //... limit soil perc by available storage volume
    availVolume = ( project->theLidProc->storage.thickness - storageDepth) *
                  project->theLidProc->storage.voidFrac;
    maxRate = availVolume/ project->Tsteplp +  project->StorageEvap +  project->StorageDrain +  project->StorageInfil;
    maxRate = MAX(maxRate, 0.0);
    project->SoilPerc = MIN( project->SoilPerc, maxRate);
  }

  //... limit surface infil. by available soil pore volume
  maxRate = ( project->theLidProc->soil.porosity - soilTheta) *
            project->theLidProc->soil.thickness /  project->Tsteplp +  project->SoilPerc;
  project->SurfaceInfil = MIN( project->SurfaceInfil, maxRate);

  //... find surface layer outflow rate
  project->SurfaceOutflow = getSurfaceOutflowRate(project, surfaceDepth);

  //... compute overall layer flux rates
  f[SURF] = ( project->SurfaceInflow -  project->SurfaceEvap -  project->SurfaceInfil -  project->SurfaceOutflow) /
            project->theLidProc->surface.voidFrac;
  f[SOIL] = ( project->SurfaceInfil -  project->SoilEvap -  project->SoilPerc) /
            project->theLidProc->soil.thickness;
  f[STOR] = ( project->SoilPerc -  project->StorageEvap -  project->StorageInfil -  project->StorageDrain) /
            project->theLidProc->storage.voidFrac;
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void trenchFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of an infiltration trench LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double surfaceDepth;
  double storageDepth;
  double head;
  double availVolume;
  double maxRate;

  //... extract zone depth levels from work vector
  surfaceDepth = x[SURF];
  storageDepth = x[STOR];

  //... convert depths to volumes
  project->SurfaceVolume = surfaceDepth *  project->theLidProc->surface.voidFrac;
  project->SoilVolume = 0.0;
  project->StorageVolume = storageDepth *  project->theLidProc->storage.voidFrac;
  availVolume = ( project->theLidProc->storage.thickness - storageDepth) *
                project->theLidProc->storage.voidFrac;

  //... nominal storage inflow
  project->StorageInflow =  project->SurfaceInflow +  project->SurfaceVolume /  project->Tsteplp;

  //... get ET rate loss for each zone
  getEvapRates(project, project->SurfaceVolume, 0.0, 0.0,  project->StorageVolume);

  //... no storage evap if surface ponded
  if ( surfaceDepth > 0.0 )  project->StorageEvap = 0.0;

  //... find infiltration rate out of storage layer
  project->StorageInfil = getStorageInfilRate(project);

  //... find underdrain flow rate
  project->StorageDrain = 0.0;
  head = storageDepth -  project->theLidProc->drain.offset;
  if (  project->theLidProc->drain.coeff > 0.0 && head >= 0.0 )
  {
    if ( storageDepth >=  project->theLidProc->storage.thickness )
    {
      head += surfaceDepth;
    }
    project->StorageDrain =  getStorageDrainRate(project,head);
  }

  //... limit underdrain flow by volume above drain offset
  if (  project->StorageDrain > 0.0 )
  {
    maxRate = (storageDepth -  project->theLidProc->drain.offset) *
              project->theLidProc->storage.voidFrac /  project->Tsteplp;
    //... add on storage inflow if storage is full
    if ( storageDepth >=  project->theLidProc->storage.thickness )
      maxRate +=  project->StorageInflow;
    project->StorageDrain = MIN( project->StorageDrain, maxRate);
  }

  //... limit storage infil. by remaining volume
  maxRate =  project->StorageVolume /  project->Tsteplp -  project->StorageDrain -  project->StorageEvap;
  maxRate = MAX(0.0, maxRate);
  project->StorageInfil = MIN( project->StorageInfil, maxRate);

  //... limit storage inflow by available storage volume
  availVolume = ( project->theLidProc->storage.thickness - storageDepth) *
                project->theLidProc->storage.voidFrac;
  maxRate = availVolume/ project->Tsteplp +  project->StorageEvap +  project->StorageDrain +  project->StorageInfil;
  maxRate = MAX(maxRate, 0.0);
  project->StorageInflow = MIN( project->StorageInflow, maxRate);

  //... equate surface infil to storage inflow
  project->SurfaceInfil =  project->StorageInflow;

  //... find surface outflow rate
  project->SurfaceOutflow = getSurfaceOutflowRate(project, surfaceDepth);

  // ... find net fluxes for each layer
  f[SURF] =  project->SurfaceInflow -  project->SurfaceEvap -  project->StorageInflow -  project->SurfaceOutflow /
             project->theLidProc->surface.voidFrac;;
  f[STOR] = ( project->StorageInflow -  project->StorageEvap -  project->StorageInfil -  project->StorageDrain) /
            project->theLidProc->storage.voidFrac;
  f[SOIL] = 0.0;
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void pavementFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates for the layers of a porous pavement LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double surfaceDepth;     // depth of water stored on surface (ft)
  double paveTheta;        // moisture content of pavement voids
  double soilTheta;        // moisture content of soil voids
  double storageDepth;     // depth of water in storage layer (ft)
  double pervVolume;       // volume/unit area of pervious pavement (ft)
  double pavePorosity;     // pavement porosity
  double storageInflow;    // inflow rate to storage layer (ft/s)
  double availVolume;
  double maxRate;
  double head;

  //... retrieve state variables from work vector
  surfaceDepth = x[SURF];
  paveTheta = x[PAVE];
  soilTheta = x[SOIL];
  storageDepth = x[STOR];
  pavePorosity =  project->theLidProc->pavement.voidFrac;

  //... convert state variables to volumes
  project->SurfaceVolume = surfaceDepth *  project->theLidProc->surface.voidFrac;
  pervVolume =  project->theLidProc->pavement.thickness *
                (1.0 -  project->theLidProc->pavement.impervFrac);
  project->PaveVolume = paveTheta * pervVolume;
  project->SoilVolume = soilTheta *  project->theLidProc->soil.thickness;
  project->StorageVolume = storageDepth *  project->theLidProc->storage.voidFrac;

  //... get ET rates (arguments are stored volumes in ft)
  availVolume =  project->SoilVolume -  project->theLidProc->soil.wiltPoint *
                 project->theLidProc->soil.thickness;
  getEvapRates(project, project->SurfaceVolume,  project->PaveVolume, availVolume,  project->StorageVolume);

  //... no storage evap if pavement layer saturated
  if ( paveTheta >= pavePorosity ||
       soilTheta >=  project->theLidProc->soil.porosity )  project->StorageEvap = 0.0;

  //... find nominal rate of surface infiltration into pavement
  project->SurfaceInfil =  project->SurfaceInflow + ( project->SurfaceVolume /  project->Tsteplp);

  //... find pavement layer permeability
  project->PavePerc = getPavementPermRate(project);

  //... limit pavement permeability to stored water + surface infil.
  maxRate =  project->PaveVolume/ project->Tsteplp +  project->SurfaceInfil;
  project->PavePerc = MIN( project->PavePerc, maxRate);

  //... find soil layer perc rate
  if (  project->theLidProc->soil.thickness > 0.0 )
    project->SoilPerc = getSoilPercRate(project,soilTheta);
  else
    project->SoilPerc =  project->PavePerc;

  //... find infiltration rate out of storage layer
  project->StorageInfil = getStorageInfilRate(project);

  //... find underdrain flow rate
  project->StorageDrain = 0.0;
  head = storageDepth -  project->theLidProc->drain.offset;
  if (  project->theLidProc->drain.coeff > 0.0 && head >= 0.0 )
  {
    if ( storageDepth >=  project->theLidProc->storage.thickness )
    {
      if (  project->theLidProc->soil.thickness > 0.0 )
      {
        if ( soilTheta >  project->theLidProc->soil.fieldCap )
        {
          head += (soilTheta -  project->theLidProc->soil.fieldCap) /
                  ( project->theLidProc->soil.porosity -
                    project->theLidProc->soil.fieldCap) *
                  project->theLidProc->soil.thickness;
          if ( soilTheta >=  project->theLidProc->soil.porosity )
          {
            head += paveTheta / pavePorosity *
                    project->theLidProc->pavement.thickness;
          }
        }
      }
      else head += paveTheta / pavePorosity *
                   project->theLidProc->pavement.thickness;
      if ( paveTheta >= pavePorosity ) head += surfaceDepth;
    }
    project->StorageDrain =  getStorageDrainRate(project, head);
  }

  //... storage layer is saturated
  if ( storageDepth >=  project->theLidProc->storage.thickness )
  {
    //... if soil layer present and is saturated
    if (  project->theLidProc->soil.thickness > 0.0 &&
          soilTheta >=  project->theLidProc->soil.porosity )
    {
      //... soil perc can't exceed storage outflow
      maxRate =  project->StorageDrain +  project->StorageInfil;
      if (  project->SoilPerc > maxRate )  project->SoilPerc = maxRate;

      //... storage outflow can't exceed soil perc
      else
      {
        //... use up available drain capacity first
        project->StorageDrain = MIN( project->StorageDrain,  project->SoilPerc);
        project->StorageInfil =  project->SoilPerc -  project->StorageDrain;
      }
    }

    //... pavement layer is saturated
    if ( paveTheta >= pavePorosity &&  project->SurfaceInfil > MIN_RUNOFF )
    {
      //... pavement outflow can't exceed surface infil or soil perc.
      project->PavePerc = MIN( project->SurfaceInfil,  project->PavePerc);
      project->PavePerc = MIN( project->PavePerc,  project->SoilPerc);

      //... pavement outflow can't exceed storage outflow
      maxRate =  project->StorageEvap +  project->StorageDrain +  project->StorageInfil;
      if (  project->PavePerc > maxRate )
      {
        project->PavePerc = maxRate;
        project->SurfaceInfil =  project->PavePerc;
      }

      //... storage outflow can't exceed pavement perm.
      else
      {
        project->StorageDrain = MIN( project->StorageDrain,  project->PavePerc);
        project->StorageInfil =  project->PavePerc -  project->StorageDrain;
      }

      //... soil perc must equal pavement perc
      project->SoilPerc =  project->PavePerc;
    }
  }

  //... storage layer not full
  else
  {
    //... limit underdrain flow by volume above drain offset
    if (  project->StorageDrain > 0.0 )
    {
      maxRate = (storageDepth -  project->theLidProc->drain.offset) *
                project->theLidProc->storage.voidFrac /  project->Tsteplp;
      project->StorageDrain = MIN( project->StorageDrain, maxRate);
    }

    //... limit storage infil. by remaining volume
    maxRate =  project->StorageVolume /  project->Tsteplp -  project->StorageDrain -  project->StorageEvap;
    maxRate = MAX(0.0, maxRate);
    project->StorageInfil = MIN( project->StorageInfil, maxRate);

    //... limit soil/pavement outflow by available storage volume
    availVolume = ( project->theLidProc->storage.thickness - storageDepth) *
                  project->theLidProc->storage.voidFrac;
    maxRate = availVolume/ project->Tsteplp +  project->StorageEvap +  project->StorageDrain +  project->StorageInfil;
    maxRate = MAX(maxRate, 0.0);
    if (  project->theLidProc->soil.thickness > 0.0 )
    {
      project->SoilPerc = MIN( project->SoilPerc, maxRate);
      maxRate = ( project->theLidProc->soil.porosity - soilTheta) *
                project->theLidProc->soil.thickness /  project->Tsteplp +  project->SoilPerc;
    }
    project->PavePerc = MIN( project->PavePerc, maxRate);

    //... limit pavement inflow by available pavement volume
    availVolume = (pavePorosity - paveTheta) * pervVolume;
    maxRate = availVolume /  project->Tsteplp +  project->PavePerc;
    project->SurfaceInfil = MIN( project->SurfaceInfil, maxRate);
  }

  //... surface outflow
  project->SurfaceOutflow = getSurfaceOutflowRate(project,surfaceDepth);

  //... compute overall layer flux rates
  f[SURF] =  project->SurfaceInflow -  project->SurfaceEvap -  project->SurfaceInfil -  project->SurfaceOutflow;
  f[PAVE] = ( project->SurfaceInfil -  project->PaveEvap -  project->PavePerc) / pervVolume;
  if (  project->theLidProc->soil.thickness > 0.0)
  {
    f[SOIL] = ( project->PavePerc -  project->SoilEvap -  project->SoilPerc) /  project->theLidProc->soil.thickness;
    storageInflow =  project->SoilPerc;
  }
  else
  {
    f[SOIL] = 0.0;
    storageInflow =  project->PavePerc;
  }
  f[STOR] = (storageInflow -  project->StorageEvap -  project->StorageInfil -  project->StorageDrain) /
            project->theLidProc->storage.voidFrac;
}

//=============================================================================

void swaleFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates from a vegetative swale LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double depth;            // depth of surface water in swale (ft)
  double topWidth;         // top width of full swale (ft)
  double botWidth;         // bottom width of swale (ft)
  double length;           // length of swale (ft)
  double surfInflow;       // inflow rate to swale (cfs)
  double surfWidth;        // top width at current water depth (ft)
  double surfArea;         // surface area of current water depth (ft2)
  double flowArea;         // x-section flow area (ft2)
  double lidArea;          // surface area of full swale (ft2)
  double hydRadius;        // hydraulic radius for current depth (ft)
  double slope;            // slope of swale side wall (run/rise)
  double volume;           // swale volume at current water depth (ft3)
  double dVdT;             // change in volume w.r.t. time (cfs)
  double dStore;           // depression storage depth (ft)
  double xDepth;           // depth above depression storage (ft)

  //... retrieve state variable from work vector
  depth = x[SURF];
  depth = MIN(depth,  project->theLidProc->surface.thickness);

  //... depression storage depth
  dStore = 0.0;

  //... get swale's bottom width
  //    (0.5 ft minimum to avoid numerical problems)
  slope =  project->theLidProc->surface.sideSlope;
  topWidth =  project->theLidUnit->fullWidth;
  topWidth = MAX(topWidth, 0.5);
  botWidth = topWidth - 2.0 * slope *  project->theLidProc->surface.thickness;
  if ( botWidth < 0.5 )
  {
    botWidth = 0.5;
    slope = 0.5 * (topWidth - 0.5) /  project->theLidProc->surface.thickness;
  }

  //... swale's length
  lidArea =  project->theLidUnit->area;
  length = lidArea / topWidth;

  //... top width, surface area and flow area of current ponded depth
  surfWidth = botWidth + 2.0 * slope * depth;
  surfArea = length * surfWidth;
  flowArea = (depth * (botWidth + slope * depth)) *
             project->theLidProc->surface.voidFrac;

  //... wet volume and effective depth
  volume = length * flowArea;

  //... surface inflow into swale (cfs)
  surfInflow =  project->SurfaceInflow * lidArea;

  //... ET rate in cfs
  project->SurfaceEvap =  project->EvapRatelp * surfArea;
  project->SurfaceEvap = MIN( project->SurfaceEvap, volume/ project->Tsteplp);

  //... infiltration rate to native soil in cfs
  project->StorageInfil =  project->SurfaceInfil * surfArea;

  //... no surface outflow if depth below depression storage
  xDepth = depth - dStore;
  if ( xDepth <= ZERO )  project->SurfaceOutflow = 0.0;

  //... otherwise compute a surface outflow
  else
  {
    //... modify flow area to remove depression storage,
    flowArea -= (dStore * (botWidth + slope * dStore)) *
                project->theLidProc->surface.voidFrac;
    if ( flowArea < ZERO )  project->SurfaceOutflow = 0.0;
    else
    {
      //... compute hydraulic radius
      botWidth = botWidth + 2.0 * dStore * slope;
      hydRadius = botWidth + 2.0 * xDepth * sqrt(1.0 + slope*slope);
      hydRadius = flowArea / hydRadius;

      //... use Manning Eqn. to find outflow rate in cfs
      project->SurfaceOutflow =  project->theLidProc->surface.alpha * flowArea *
                                 pow(hydRadius, 2./3.);
    }
  }

  //... net flux rate (dV/dt) in cfs
  dVdT = surfInflow -  project->SurfaceEvap -  project->StorageInfil -  project->SurfaceOutflow;           //(5.1.009)

  //... when full, any net positive inflow becomes spillage
  if ( depth ==  project->theLidProc->surface.thickness && dVdT > 0.0 )
  {
    project->SurfaceOutflow += dVdT;
    dVdT = 0.0;
  }

  //... convert flux rates to ft/s
  project->SurfaceEvap /= lidArea;
  project->StorageInfil /= lidArea;
  project->SurfaceOutflow /= lidArea;
  f[SURF] = dVdT / surfArea;
  f[SOIL] = 0.0;
  f[STOR] = 0.0;

  //... assign values to layer volumes
  project->SurfaceVolume = volume / lidArea;
  project->SoilVolume = 0.0;
  project->StorageVolume = 0.0;
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void barrelFluxRates(Project* project, double x[], double f[])
//
//  Purpose: computes flux rates for a rain barrel LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
  double storageDepth = x[STOR];
  double head;
  double maxValue;

  //... assign values to layer volumes
  project->SurfaceVolume = 0.0;
  project->SoilVolume = 0.0;
  project->StorageVolume = storageDepth;

  //... initialize flows
  project->SurfaceInfil = 0.0;
  project->SurfaceOutflow = 0.0;
  project->StorageDrain = 0.0;

  //... compute outflow if time since last rain exceeds drain delay
  //    (dryTime is updated in lid.evalLidUnit at each time step)
  if (  project->theLidProc->drain.delay == 0.0 ||
        project->theLidUnit->dryTime >=  project->theLidProc->drain.delay )
  {
    head = storageDepth -  project->theLidProc->drain.offset;
    if ( head > 0.0 )
    {
      project->StorageDrain = getStorageDrainRate(project,head);
      maxValue = (head/ project->Tsteplp);
      project->StorageDrain = MIN( project->StorageDrain, maxValue);
    }
  }

  //... limit inflow to available storage
  project->StorageInflow =  project->SurfaceInflow;
  maxValue = ( project->theLidProc->storage.thickness - storageDepth) /  project->Tsteplp +
             project->StorageDrain;
  project->StorageInflow = MIN( project->StorageInflow, maxValue);
  project->SurfaceInfil =  project->StorageInflow;

  //... assign values to layer flux rates
  f[SURF] =  project->SurfaceInflow -  project->StorageInflow;
  f[STOR] =  project->StorageInflow -  project->StorageDrain;
  f[SOIL] = 0.0;
}

//=============================================================================

double getSurfaceOutflowRate(Project* project, double depth)
//
//  Purpose: computes outflow rate from a LID's surface layer.
//  Input:   depth = depth of ponded water on surface layer (ft)
//  Output:  returns outflow from surface layer (ft/s)
//
//  Note: this function should not be applied to swales or rain barrels.
//
{
  double delta;
  double outflow;

  //... no outflow if ponded depth below storage depth
  delta = depth -  project->theLidProc->surface.thickness;
  if ( delta < 0.0 ) return 0.0;

  //... compute outflow from overland flow Manning equation
  outflow =  project->theLidProc->surface.alpha * pow(delta, 5.0/3.0) *
             project->theLidUnit->fullWidth /  project->theLidUnit->area;
  outflow = MIN(outflow, delta /  project->Tsteplp);
  return outflow;
}

//=============================================================================

double getPavementPermRate(Project* project)
//
//  Purpose: computes reduced permeability of a pavement layer due to
//           clogging.
//  Input:   none
//  Output:  returns the reduced permeability of the pavement layer (ft/s).
//
{
  double permRate;
  double permReduction;

  permReduction =  project->theLidProc->pavement.clogFactor;
  if ( permReduction > 0.0 )
  {
    permReduction =  project->theLidUnit->waterBalance.inflow / permReduction;
    permReduction = MIN(permReduction, 1.0);
  }
  permRate =  project->theLidProc->pavement.kSat * (1.0 - permReduction);
  return permRate;
}

//=============================================================================

////  This function was modified for release 5.1.007.  ////                    //(5.1.007)

double getSoilPercRate(Project* project, double theta)
//
//  Purpose: computes percolation rate of water through a LID's soil layer.
//  Input:   theta = moisture content (fraction)
//  Output:  returns percolation rate within soil layer (ft/s)
//
{
  double percRate;         // percolation rate (ft/s)
  double delta;            // moisture deficit
  double maxValue;         // max. allowable perc. rate (ft/s)

  // ... max. drainable soil moisture
  maxValue = (theta -  project->theLidProc->soil.fieldCap) *
             project->theLidProc->soil.thickness /  project->Tsteplp;
  if ( maxValue <= 0.0 ) return 0.0;

  // ... perc rate = unsaturated hydraulic conductivity
  delta =  project->theLidProc->soil.porosity - theta;
  percRate =  project->theLidProc->soil.kSat * exp(-delta *  project->theLidProc->soil.kSlope);

  //... rate limited by drainable moisture content
  percRate = MIN(percRate, maxValue);
  return percRate;
}

//=============================================================================

double getStorageInfilRate(Project* project)
//
//  Purpose: computes infiltration rate between storage zone and
//           native soil beneath a LID.
//  Input:   depth = depth of water storage zone (ft)
//  Output:  returns infiltration rate (ft/s)
//
{
  double infil = 0.0;
  double clogFactor = 0.0;

  if (  project->theLidProc->storage.kSat == 0.0 ) return 0.0;
  if (  project->MaxNativeInfillp == 0.0 ) return 0.0;

  //... reduction due to clogging
  clogFactor =  project->theLidProc->storage.clogFactor;
  if ( clogFactor > 0.0 )
  {
    clogFactor =  project->theLidUnit->waterBalance.inflow / clogFactor;
    clogFactor = MIN(clogFactor, 1.0);
  }

  //... infiltration rate = storage Ksat reduced by any clogging
  infil =  project->theLidProc->storage.kSat * (1.0 - clogFactor);

  //... limit infiltration rate by any groundwater-imposed limit
  return MIN(infil,  project->MaxNativeInfillp);
}

//=============================================================================

////  This function was modified for release 5.1.007.  ////                    //(5.1.007)

double  getStorageDrainRate(Project* project, double head)
//
//  Purpose: computes underdrain flow rate in a LID's storage layer.
//  Input:   head = head of water above underdrain (ft)
//           inflow = rate of inflow to storage zone (ft/s)
//  Output:  returns flow in underdrain (ft/s)
//
//  Note:    drain eqn. is evaluated in user's units.
{
  double outflow = 0.0;

  if (  project->theLidProc->drain.coeff > 0.0 && head > ZERO )
  {
    // ... evaluate underdrain flow rate equation
    head *= UCF(project,RAINDEPTH);
    outflow =  project->theLidProc->drain.coeff *
               pow(head,  project->theLidProc->drain.expon);
    outflow /= UCF(project,RAINFALL);
  }
  return outflow;
}

//=============================================================================

////  This function was modified for release 5.1.007.  ////                    //(5.1.007)

double getDrainMatOutflow(Project *project, double depth)
{
  double result =  project->SoilPerc;
  if (  project->theLidProc->drainMat.alpha > 0.0 )
  {
    result =  project->theLidProc->drainMat.alpha * pow(depth, 5.0/3.0) *
              project->theLidUnit->fullWidth /  project->theLidUnit->area *
              project->theLidProc->drainMat.voidFrac;                                //(5.1.008)
  }
  return result;
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void getEvapRates(Project *project, double surfaceVol, double paveVol, double soilVol,
                  double storageVol)
//
//  Purpose: computes surface, pavement, soil, and storage evaporation rates.
//  Input:   surfaceVol = volume/area of ponded water on surface layer (ft)
//           paveVol    = volume/area of water in pavement pores (ft)
//           soilVol    = volume/area of water in soil (or pavement) pores (ft)
//           storageVol = volume/area of water in storage layer (ft)
//  Output:  none
//
{
  double availEvap;

  //... surface evaporation flux
  availEvap =  project->EvapRatelp;
  project->SurfaceEvap = MIN(availEvap, surfaceVol/ project->Tsteplp);
  project->SurfaceEvap = MAX(0.0,  project->SurfaceEvap);
  availEvap = MAX(0.0, (availEvap -  project->SurfaceEvap));

  //... no subsurface evap if water is infiltrating
  if (  project->SurfaceInfil > 0.0 )
  {
    project->PaveEvap = 0.0;
    project->SoilEvap = 0.0;
    project->StorageEvap = 0.0;
  }
  else
  {
    //... pavement evaporation flux
    project->PaveEvap = MIN(availEvap, paveVol /  project->Tsteplp);
    availEvap = MAX(0.0, (availEvap -  project->PaveEvap));

    //... soil evaporation flux
    project->SoilEvap = MIN(availEvap, soilVol /  project->Tsteplp);
    availEvap = MAX(0.0, (availEvap -  project->SoilEvap));

    //... storage evaporation flux
    project->StorageEvap = MIN(availEvap, storageVol /  project->Tsteplp);
  }
}

//=============================================================================

double getSurfaceOverflowRate(Project* project, double* surfaceDepth)
//
//  Purpose: finds surface overflow rate from a LID unit.
//  Input:   surfaceDepth = depth of water stored in surface layer (ft)
//  Output:  returns the overflow rate (ft/s)
//
{
  double delta = *surfaceDepth -  project->theLidProc->surface.thickness;
  if (  delta <= 0.0 ) return 0.0;
  *surfaceDepth =  project->theLidProc->surface.thickness;
  return delta *  project->theLidProc->surface.voidFrac /  project->Tsteplp;
}

//=============================================================================

void updateWaterBalance(Project* project, TLidUnit *lidUnit, double inflow, double evap,
                        double infil, double surfFlow, double drainFlow, double storage)
//
//  Purpose: updates components of the water mass balance for a LID unit
//           over the current time step.
//  Input:   lidUnit   = a particular LID unit
//           inflow    = runon + rainfall to the LID unit (ft/s)
//           evap      = evaporation rate from the unit (ft/s)
//           infil     = infiltration out the bottom of the unit (ft/s)
//           surfFlow  = surface runoff from the unit (ft/s)
//           drainFlow = underdrain flow from the unit
//           storage   = volume of water stored in the unit (ft)
//  Output:  none
//
{
  lidUnit->waterBalance.inflow += inflow *  project->Tsteplp;
  lidUnit->waterBalance.evap += evap *  project->Tsteplp;
  lidUnit->waterBalance.infil += infil *  project->Tsteplp;
  lidUnit->waterBalance.surfFlow += surfFlow *  project->Tsteplp;
  lidUnit->waterBalance.drainFlow += drainFlow *  project->Tsteplp;
  lidUnit->waterBalance.finalVol = storage;
}

//=============================================================================

int modpuls_solve(Project* project, int n, double* x, double* xOld, double* xPrev,
                  double* xMin, double* xMax, double* xTol,
                  double* qOld, double* q, double dt, double omega,            //(5.1.007)
                  void (*derivs)(Project*, double*, double*))
//
//  Purpose: solves system of equations dx/dt = q(x) for x at end of time step
//           dt using a modified Puls method.
//  Input:   n = number of state variables
//           x = vector of state variables
//           xOld = state variable values at start of time step
//           xPrev = state variable values from previous iteration
//           xMin = lower limits on state variables
//           xMax = upper limits on state variables
//           xTol = convergence tolerances on state variables
//           qOld = flux rates at start of time step
//           q = flux rates at end of time step
//           dt = time step (sec)
//           omega = time weighting parameter (use 0 for Euler method          //(5.1.007)
//                   or 0.5 for modified Puls method)                          //(5.1.007)
//           derivs = pointer to function that computes flux rates q as a
//                    function of state variables x
//  Output:  returns number of steps required for convergence (or 0 if
//           process doesn't converge)
//
{
  int i;
  int canStop;
  int steps = 1;
  int maxSteps = 20;

  //... initialize state variable values
  for (i=0; i<n; i++)
  {
    xOld[i] = x[i];
    xPrev[i] = x[i];
  }

  //... repeat until convergence achieved
  while (steps < maxSteps)
  {
    //... compute flux rates for current state levels
    canStop = 1;
    derivs(project,x, q);

    //... update state levels based on current flux rates
    for (i=0; i<n; i++)
    {
      x[i] = xOld[i] + (omega*qOld[i] + (1.0 - omega)*q[i]) * dt;
      x[i] = MIN(x[i], xMax[i]);
      x[i] = MAX(x[i], xMin[i]);

      if ( omega > 0.0 &&                                                //(5.1.007)
           fabs(x[i] - xPrev[i]) > xTol[i] ) canStop = 0;
      xPrev[i] = x[i];
    }

    //... return if process converges
    if (canStop) return steps;
    steps++;
  }

  //... no convergence so return 0
  return 0;
}
