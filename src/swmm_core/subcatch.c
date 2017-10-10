//-----------------------------------------------------------------------------
//   subcatch.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.000)
//             04/19/14  (Build 5.1.006)
//             03/19/15  (Build 5.1.008)
//             04/30/15  (Build 5.1.009)
//             08/05/15  (Build 5.1.010)
//   Author:   L. Rossman
//
//   Subcatchment runoff functions.
//
//   Build 5.1.008:
//   - Support added for keeping separate track of drain outflows from LIDs.
//   - Processing of inflow/outflow volumes over a time step was refactored. 
//   - Reported subcatchment runoff includes both surface runoff and LID
//     drain flows, even though latter can be routed elsewhere.
//   - Runon now distributed only over non-LID area of a subcatchment, unless
//     LID covers full area.
//   - Pollutant buildup and washoff functions were moved to surfqual.c.
//
//   Build 5.1.009:
//   - Runon for full LID subcatchment added to statistical summary.
//
//   Build 5.1.010:
//   - Fixed a bug introduced in 5.1.008 that forgot to include LID
//     exfiltration as inflow sent to GW routine.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <string.h>
#include "headers.h"
#include "lid.h"
#include "odesolve.h"

//-----------------------------------------------------------------------------
// Constants 
//-----------------------------------------------------------------------------
const double MCOEFF    = 1.49;              // constant in Manning Eq.
const double MEXP      = 1.6666667;         // exponent in Manning Eq.
const double ODETOL    = 0.0001;            // acceptable error for ODE solver

////-----------------------------------------------------------------------------
//// Globally shared variables   
////-----------------------------------------------------------------------------
//// Volumes (ft3) for a subcatchment over a time step                           //(5.1.008)
//double     project->Vevap;         // evaporation
//double     project->Vpevap;        // pervious area evaporation
//double     project->Vinfil;        // non-LID infiltration
//double     project->Vinflow;       // non-LID precip + snowmelt + runon + ponded water
//double     project->Voutflow;      // non-LID runoff to subcatchment's outlet
//double     project->VlidIn;        // impervious area flow to LID units
//double     project->VlidInfil;     // infiltration from LID units
//double     project->VlidOut;       // surface outflow from LID units
//double     project->VlidDrain;     // drain outflow from LID units
//double     project->VlidReturn;    // LID outflow returned to pervious area

//-----------------------------------------------------------------------------
// Locally shared variables   
//-----------------------------------------------------------------------------
static  TSubarea* theSubarea;     // subarea to which getDdDt() is applied
static  char *RunoffRoutingWords[] = { w_OUTLET,  w_IMPERV, w_PERV, NULL};

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)   
//-----------------------------------------------------------------------------
//  subcatch_readParams        (called from parseLine in input.c)
//  subcatch_readSubareaParams (called from parseLine in input.c)
//  subcatch_readLanduseParams (called from parseLine in input.c)
//  subcatch_readInitBuildup   (called from parseLine in input.c)

//  subcatch_validate          (called from project_validate)
//  subcatch_initState         (called from project_init)

//  subcatch_setOldState       (called from runoff_execute)
//  subcatch_getRunon          (called from runoff_execute)
//  subcatch_addRunon          (called from subcatch_getRunon,
//                              lid_addDrainRunon, & runoff_getOutfallRunon)   //(5.1.008)
//  subcatch_getRunoff         (called from runoff_execute)
//  subcatch_hadRunoff         (called from runoff_execute)

//  subcatch_getFracPerv       (called from gwater_initState)
//  subcatch_getStorage        (called from massbal_getRunoffError)
//  subcatch_getDepth          (called from findPondedLoads in surfqual.c)

//  subcatch_getWtdOutflow     (called from addWetWeatherInflows in routing.c)
//  subcatch_getResults        (called from output_saveSubcatchResults)

//-----------------------------------------------------------------------------
// Function declarations
//-----------------------------------------------------------------------------
static void   getNetPrecip(Project* project, int j, double* netPrecip, double tStep);
static double getSubareaRunoff(Project* project, int subcatch, int subarea, double area,         //(5.1.008)
              double rainfall, double evap, double tStep);
static double getSubareaInfil(Project* project, int j, TSubarea* subarea, double precip,
              double tStep);
static double findSubareaRunoff(TSubarea* subarea, double tRunoff);            //(5.1.008)
static void   updatePondedDepth(Project* project, TSubarea* subarea, double* tx);
static void   getDdDt(Project* project, double t, double* d, double* dddt);

//=============================================================================

int  subcatch_readParams(Project* project, int j, char* tok[], int ntoks)
//
//  Input:   j = subcatchment index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads subcatchment parameters from a tokenized  line of input data.
//
//  Data has format:
//    Name  RainGage  Outlet  Area  %Imperv  Width  Slope CurbLength  Snowpack  
//
{
    int    i, k, m;
    char*  id;
    double x[9];

    // --- check for enough tokens
    if ( ntoks < 8 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that named subcatch exists
    id = project_findID(project,SUBCATCH, tok[0]);
    if ( id == NULL ) return error_setInpError(ERR_NAME, tok[0]);

    // --- check that rain gage exists
    k = project_findObject(project,GAGE, tok[1]);
    if ( k < 0 ) return error_setInpError(ERR_NAME, tok[1]);
    x[0] = k;

    // --- check that outlet node or subcatch exists
    m = project_findObject(project,NODE, tok[2]);
    x[1] = m;
    m = project_findObject(project,SUBCATCH, tok[2]);
    x[2] = m;
    if ( x[1] < 0.0 && x[2] < 0.0 )
        return error_setInpError(ERR_NAME, tok[2]);

    // --- read area, %imperv, width, slope, & curb length
    for ( i = 3; i < 8; i++)
    {
        if ( ! getDouble(tok[i], &x[i]) || x[i] < 0.0 )
            return error_setInpError(ERR_NUMBER, tok[i]);
    }

    // --- if snowmelt object named, check that it exists
    x[8] = -1;
    if ( ntoks > 8 )
    {
        k = project_findObject(project,SNOWMELT, tok[8]);
        if ( k < 0 ) return error_setInpError(ERR_NAME, tok[8]);
        x[8] = k;
    }

    // --- assign input values to subcatch's properties
    project->Subcatch[j].ID = id;
    project->Subcatch[j].gage        = (int)x[0];
    project->Subcatch[j].outNode     = (int)x[1];
    project->Subcatch[j].outSubcatch = (int)x[2];
	project->Subcatch[j].area = x[3] / UCF(project, LANDAREA);
    project->Subcatch[j].fracImperv  = x[4] / 100.0;
	project->Subcatch[j].width = x[5] / UCF(project, LENGTH);
    project->Subcatch[j].slope       = x[6] / 100.0;
    project->Subcatch[j].curbLength  = x[7];

    // --- create the snow pack object if it hasn't already been created
    if ( x[8] >= 0 )
    {
        if ( !snow_createSnowpack(project,j, (int)x[8]) )
            return error_setInpError(ERR_MEMORY, "");
    }
    return 0;
}

//=============================================================================

int subcatch_readSubareaParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads subcatchment's subarea parameters from a tokenized 
//           line of input data.
//
//  Data has format:
//    Subcatch  Imperv_N  Perv_N  Imperv_S  Perv_S  PctZero  RouteTo (PctRouted)
//
{
    int    i, j, k, m;
    double x[7];

    // --- check for enough tokens
    if ( ntoks < 7 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that named subcatch exists
    j = project_findObject(project,SUBCATCH, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);

    // --- read in Mannings n, depression storage, & PctZero values
    for (i = 0; i < 5; i++)
    {
        if ( ! getDouble(tok[i+1], &x[i])  || x[i] < 0.0 )
            return error_setInpError(ERR_NAME, tok[i+1]);
    }

    // --- check for valid runoff routing keyword
    m = findmatch(tok[6], RunoffRoutingWords);
    if ( m < 0 ) return error_setInpError(ERR_KEYWORD, tok[6]);

    // --- get percent routed parameter if present (default is 100)
    x[5] = m;
    x[6] = 1.0;
    if ( ntoks >= 8 )
    {
        if ( ! getDouble(tok[7], &x[6]) || x[6] < 0.0 || x[6] > 100.0 )
            return error_setInpError(ERR_NUMBER, tok[7]);
        x[6] /= 100.0;
    }

    // --- assign input values to each type of subarea
    project->Subcatch[j].subArea[IMPERV0].N = x[0];
    project->Subcatch[j].subArea[IMPERV1].N = x[0];
    project->Subcatch[j].subArea[PERV].N    = x[1];

    project->Subcatch[j].subArea[IMPERV0].dStore = 0.0;
	project->Subcatch[j].subArea[IMPERV1].dStore = x[2] / UCF(project, RAINDEPTH);
	project->Subcatch[j].subArea[PERV].dStore = x[3] / UCF(project, RAINDEPTH);

    project->Subcatch[j].subArea[IMPERV0].fArea  = project->Subcatch[j].fracImperv * x[4] / 100.0;
    project->Subcatch[j].subArea[IMPERV1].fArea  = project->Subcatch[j].fracImperv * (1.0 - x[4] / 100.0);
    project->Subcatch[j].subArea[PERV].fArea     = (1.0 - project->Subcatch[j].fracImperv);

    // --- assume that all runoff from each subarea goes to subcatch outlet
    for (i = IMPERV0; i <= PERV; i++)
    {
        project->Subcatch[j].subArea[i].routeTo = TO_OUTLET;
        project->Subcatch[j].subArea[i].fOutlet = 1.0;
    }

    // --- modify routing if pervious runoff routed to impervious area
    //     (fOutlet is the fraction of runoff not routed)
    
    k = (int)x[5];
    if ( project->Subcatch[j].fracImperv == 0.0
    ||   project->Subcatch[j].fracImperv == 1.0 ) k = TO_OUTLET;
    if ( k == TO_IMPERV && project->Subcatch[j].fracImperv )
    {
        project->Subcatch[j].subArea[PERV].routeTo = k;
        project->Subcatch[j].subArea[PERV].fOutlet = 1.0 - x[6];
    }

    // --- modify routing if impervious runoff routed to pervious area
    if ( k == TO_PERV )
    {
        project->Subcatch[j].subArea[IMPERV0].routeTo = k;
        project->Subcatch[j].subArea[IMPERV1].routeTo = k;
        project->Subcatch[j].subArea[IMPERV0].fOutlet = 1.0 - x[6];
        project->Subcatch[j].subArea[IMPERV1].fOutlet = 1.0 - x[6];
    }
    return 0;
}

//=============================================================================

int subcatch_readLanduseParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads assignment of landuses to subcatchment from a tokenized 
//           line of input data.
//
//  Data has format:
//    Subcatch  landuse  percent .... landuse  percent
//
{
    int     j, k, m;
    double  f;

    // --- check for enough tokens
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that named subcatch exists
    j = project_findObject(project,SUBCATCH, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);

    // --- process each pair of landuse - percent items
    for ( k = 2; k <= ntoks; k = k+2)
    {
        // --- check that named land use exists and is followed by a percent
        m = project_findObject(project,LANDUSE, tok[k-1]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[k-1]);
        if ( k+1 > ntoks ) return error_setInpError(ERR_ITEMS, "");
        if ( ! getDouble(tok[k], &f) )
            return error_setInpError(ERR_NUMBER, tok[k]);

        // --- store land use fraction in subcatch's landFactor property
        project->Subcatch[j].landFactor[m].fraction = f/100.0;
    }
    return 0;
}

//=============================================================================

int subcatch_readInitBuildup(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads initial pollutant buildup on subcatchment from 
//           tokenized line of input data.
//
//  Data has format:
//    Subcatch  pollut  initLoad .... pollut  initLoad
//
{
    int     j, k, m;
    double  x;

    // --- check for enough tokens
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that named subcatch exists
    j = project_findObject(project,SUBCATCH, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);

    // --- process each pair of pollutant - init. load items
    for ( k = 2; k <= ntoks; k = k+2)
    {
        // --- check for valid pollutant name and loading value
        m = project_findObject(project,POLLUT, tok[k-1]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[k-1]);
        if ( k+1 > ntoks ) return error_setInpError(ERR_ITEMS, "");
        if ( ! getDouble(tok[k], &x) )
            return error_setInpError(ERR_NUMBER, tok[k]);

        // --- store loading in subcatch's initBuildup property
        project->Subcatch[j].initBuildup[m] = x;
    }
    return 0;
}

//=============================================================================

void  subcatch_validate(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: checks for valid subcatchment input parameters.
//
{
    int     i;
    double  area;
    double  nonLidArea = project->Subcatch[j].area;

    // --- check for ambiguous outlet name
    if ( project->Subcatch[j].outNode >= 0 && project->Subcatch[j].outSubcatch >= 0 )
        report_writeErrorMsg(project,ERR_SUBCATCH_OUTLET, project->Subcatch[j].ID);

    // --- validate subcatchment's groundwater component 
    gwater_validate(project,j);

    // --- validate placement of LIDs in the subcatchment
    nonLidArea -= project->Subcatch[j].lidArea;

    // --- compute alpha (i.e. WCON in old SWMM) for overland flow
    //     NOTE: the area which contributes to alpha for both imperv
    //     subareas w/ and w/o depression storage is the total imperv area.
    for (i = IMPERV0; i <= PERV; i++)
    {
        if ( i == PERV )
        {
            area = (1.0 - project->Subcatch[j].fracImperv) * nonLidArea;
        }
        else
        {
             area = project->Subcatch[j].fracImperv * nonLidArea;
        }
        project->Subcatch[j].subArea[i].alpha = 0.0;
        if ( area > 0.0 && project->Subcatch[j].subArea[i].N > 0.0 )
        {
            project->Subcatch[j].subArea[i].alpha = MCOEFF * project->Subcatch[j].width / area *
                sqrt(project->Subcatch[j].slope) / project->Subcatch[j].subArea[i].N;
        }
    }
}

//=============================================================================

void  subcatch_initState(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: Initializes the state of a subcatchment.
//
{
    int    i;

    // --- initialize rainfall, runoff, & snow depth
    project->Subcatch[j].rainfall = 0.0;
    project->Subcatch[j].oldRunoff = 0.0;
    project->Subcatch[j].newRunoff = 0.0;
    project->Subcatch[j].oldSnowDepth = 0.0;
    project->Subcatch[j].newSnowDepth = 0.0;
    project->Subcatch[j].runon = 0.0;
    project->Subcatch[j].evapLoss = 0.0;                                                //(5.1.008)
    project->Subcatch[j].infilLoss = 0.0;                                               //(5.1.008)

    // --- set isUsed property of subcatchment's rain gage
    i = project->Subcatch[j].gage;
    if ( i >= 0 )
    {
        project->Gage[i].isUsed = TRUE;
        if ( project->Gage[i].coGage >= 0 ) project->Gage[project->Gage[i].coGage].isUsed = TRUE;
    }

    // --- initialize state of infiltration, groundwater, & snow pack objects
    if ( project->Subcatch[j].infil == j )  infil_initState(project,j, project->InfilModel);
    if ( project->Subcatch[j].groundwater ) gwater_initState(project,j);
    if ( project->Subcatch[j].snowpack )    snow_initSnowpack(project,j);

    // --- initialize state of sub-areas
    for (i = IMPERV0; i <= PERV; i++)
    {
        project->Subcatch[j].subArea[i].depth  = 0.0;
        project->Subcatch[j].subArea[i].inflow = 0.0;
        project->Subcatch[j].subArea[i].runoff = 0.0;
    }

    // --- initialize runoff quality
	surfqual_initState(project, j);
}

//=============================================================================

void subcatch_setOldState(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: replaces old state of subcatchment with new state.
//
{
    int i;
    project->Subcatch[j].oldRunoff = project->Subcatch[j].newRunoff;
    project->Subcatch[j].oldSnowDepth = project->Subcatch[j].newSnowDepth;
    for (i = IMPERV0; i <= PERV; i++)
    {
        project->Subcatch[j].subArea[i].inflow = 0.0;
    }
    for (i = 0; i < project->Nobjects[POLLUT]; i++)
    {
        project->Subcatch[j].oldQual[i] = project->Subcatch[j].newQual[i];
        project->Subcatch[j].newQual[i] = 0.0;
    }
    lid_setOldGroupState(project,j);                                                   //(5.1.008)
}

//=============================================================================

double subcatch_getFracPerv(Project* project, int j)
//
//  Purpose: determines what fraction of subcatchment area, including any LID
//           area, is pervious.
//  Input:   j = subcatchment index
//  Output:  returns fraction of area with pervious cover
//
{
    double fracPerv = 1.0 - project->Subcatch[j].fracImperv;

    if ( project->Subcatch[j].lidArea > 0.0 )
    {
        fracPerv = (fracPerv * (project->Subcatch[j].area - project->Subcatch[j].lidArea) + 
                    lid_getPervArea(project,j)) / project->Subcatch[j].area;
        fracPerv = MIN(fracPerv, 1.0);
    }
    return fracPerv;
}

//=============================================================================

double subcatch_getStorage(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  returns total volume of stored water (ft3)
//  Purpose: finds total volume of water stored on a subcatchment's surface
//           and its LIDs at the current time.
//
{
    int    i;
    double v = 0.0;

    for ( i = IMPERV0; i <= PERV; i++)
    {
        v += project->Subcatch[j].subArea[i].depth * project->Subcatch[j].subArea[i].fArea;
    }
    return v * (project->Subcatch[j].area - project->Subcatch[j].lidArea) +
           lid_getStoredVolume(project,j);
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

void subcatch_getRunon(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: Routes runoff from a subcatchment to its outlet subcatchment
//           or between its subareas.
//
{
    int    k;                          // outlet subcatchment index
    int    p;                          // pollutant index
    double q;                          // runon to outlet subcatchment (ft/sec)
    double q1, q2;                     // runoff from imperv. areas (ft/sec)
    double pervArea;                   // subcatchment pervious area (ft2)

    // --- add previous period's runoff from this subcatchment to the
    //     runon of the outflow subcatchment, if it exists
    k = project->Subcatch[j].outSubcatch;
    q = project->Subcatch[j].oldRunoff;
    if ( k >= 0 && k != j )
    {
        subcatch_addRunonFlow(project,k, q);
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            project->Subcatch[k].newQual[p] += q * project->Subcatch[j].oldQual[p] * LperFT3;
        }
    }

    // --- add any LID underdrain flow sent from this subcatchment to
    //     other subcatchments
    if ( project->Subcatch[j].lidArea > 0.0 ) lid_addDrainRunon(project,j);

    // --- add to sub-area inflow any outflow from other subarea in previous period
    //     (NOTE: no transfer of runoff pollutant load, since runoff loads are
    //     based on runoff flow from entire subcatchment.)

    // --- Case 1: imperv --> perv
    if ( project->Subcatch[j].fracImperv < 1.0 &&
         project->Subcatch[j].subArea[IMPERV0].routeTo == TO_PERV )
    {
        // --- add area-wtd. outflow from imperv1 subarea to perv area inflow
        q1 = project->Subcatch[j].subArea[IMPERV0].runoff *
             project->Subcatch[j].subArea[IMPERV0].fArea;
        q2 = project->Subcatch[j].subArea[IMPERV1].runoff *
             project->Subcatch[j].subArea[IMPERV1].fArea;
        q = q1 + q2;
        project->Subcatch[j].subArea[PERV].inflow += q *
             (1.0 - project->Subcatch[j].subArea[IMPERV0].fOutlet) /
             project->Subcatch[j].subArea[PERV].fArea;
    }

    // --- Case 2: perv --> imperv
    if ( project->Subcatch[j].fracImperv > 0.0 &&
         project->Subcatch[j].subArea[PERV].routeTo == TO_IMPERV &&
         project->Subcatch[j].subArea[IMPERV1].fArea > 0.0 )
    {
        q = project->Subcatch[j].subArea[PERV].runoff;
        project->Subcatch[j].subArea[IMPERV1].inflow +=
            q * (1.0 - project->Subcatch[j].subArea[PERV].fOutlet) *
            project->Subcatch[j].subArea[PERV].fArea /
            project->Subcatch[j].subArea[IMPERV1].fArea;
    }

    // --- Add any return flow from LID units to pervious subarea
    if ( project->Subcatch[j].lidArea > 0.0 && project->Subcatch[j].fracImperv < 1.0 )
    {
        pervArea = project->Subcatch[j].subArea[PERV].fArea *
                   (project->Subcatch[j].area - project->Subcatch[j].lidArea);
        q = lid_getFlowToPerv(project,j);
        if ( pervArea > 0.0 )
        {
            project->Subcatch[j].subArea[PERV].inflow += q / pervArea;
        }
    }
}

//=============================================================================

////  New function added to release 5.1.008.  ////                            //(5.1.008)

void  subcatch_addRunonFlow(Project* project, int k, double q)
//
//  Input:   k = subcatchment index
//           q = runon flow rate (cfs) to subcatchment k
//  Output:  none
//  Purpose: Updates the total runon flow (ft/s) seen by a subcatchment that
//           receives runon flow from an upstream subcatchment.
//
{
    int i;
    double nonLidArea;

    // --- distribute runoff from upstream subcatchment (in cfs)
    //     uniformly over the non-LID area of current subcatchment (ft/sec)
    if ( project->Subcatch[k].area <= 0.0 ) return;
    nonLidArea = project->Subcatch[k].area - project->Subcatch[k].lidArea; 
    if ( nonLidArea > 0.0 ) q = q / nonLidArea;
    else                    q = q / project->Subcatch[k].area;
    project->Subcatch[k].runon += q;

    // --- assign this flow to the 3 types of subareas
    for (i = IMPERV0; i <= PERV; i++)
    {
        project->Subcatch[k].subArea[i].inflow += q;
    }
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double subcatch_getRunoff(Project* project, int j, double tStep)
//
//  Input:   j = subcatchment index
//           tStep = time step (sec)
//  Output:  returns total runoff produced by subcatchment (ft/sec)
//  Purpose: Computes runoff & new storage depth for subcatchment.
//
//  The 'runoff' value returned by this function is the total runoff
//  generated (in ft/sec) by the subcatchment before any internal
//  re-routing is applied. It is used to compute pollutant washoff.
//
//  The 'outflow' value computed here (in cfs) is the surface runoff
//  that actually leaves the subcatchment after any LID controls are
//  applied and is saved to project->Subcatch[j].newRunoff. 
//
{
    int    i;                          // subarea index
    double nonLidArea;                 // non-LID portion of subcatch area (ft2)
    double area;                       // sub-area or subcatchment area (ft2)
    double netPrecip[3];               // subarea net precipitation (ft/sec)
    double vRain;                      // rainfall (+ snowfall) volume (ft3)
    double vRunon    = 0.0;            // runon volume from other areas (ft3)
    double vOutflow  = 0.0;            // runoff volume leaving subcatch (ft3)
    double runoff    = 0.0;            // total runoff flow on subcatch (cfs)
    double evapRate  = 0.0;            // potential evaporation rate (ft/sec)

    // --- initialize shared water balance variables
    project->Vevap     = 0.0;
    project->Vpevap    = 0.0;
    project->Vinfil    = 0.0;
    project->Voutflow  = 0.0;
    project->VlidIn    = 0.0;
    project->VlidInfil = 0.0;
    project->VlidOut   = 0.0;
    project->VlidDrain = 0.0;
    project->VlidReturn = 0.0;

    // --- find volume of inflow to non-LID portion of subcatchment as existing
    //     ponded water + any runon volume from upstream areas;
    //     rainfall and snowmelt will be added as each sub-area is analyzed
    nonLidArea = project->Subcatch[j].area - project->Subcatch[j].lidArea;
    vRunon = project->Subcatch[j].runon * tStep * nonLidArea;
    project->Vinflow = vRunon + subcatch_getDepth(project,j) * nonLidArea;

////  Added to release 5.1.009.  ////                                          //(5.1.009)
    // --- find LID runon only if LID occupies full subcatchment
    if ( nonLidArea == 0.0 )
        vRunon = project->Subcatch[j].runon * tStep * project->Subcatch[j].area;
////

    // --- get net precip. (rainfall + snowfall + snowmelt) on the 3 types
    //     of subcatchment sub-areas and update project->Vinflow with it
    getNetPrecip(project,j, netPrecip, tStep);

    // --- find potential evaporation rate
    if ( project->Evap.dryOnly && project->Subcatch[j].rainfall > 0.0 ) evapRate = 0.0;
    else evapRate = project->Evap.rate;

    // --- examine each type of sub-area (impervious w/o depression storage,
    //     impervious w/ depression storage, and pervious)
    if ( nonLidArea > 0.0 ) for (i = IMPERV0; i <= PERV; i++)
    {
        // --- get runoff from sub-area updating project->Vevap, project->Vpevap,
        //     project->Vinfil & project->Voutflow)
        area = nonLidArea * project->Subcatch[j].subArea[i].fArea;
        project->Subcatch[j].subArea[i].runoff =
            getSubareaRunoff(project,j, i, area, netPrecip[i], evapRate, tStep);
        runoff += project->Subcatch[j].subArea[i].runoff * area;
    }

    // --- evaluate any LID treatment provided (updating project->Vevap,
    //     project->Vpevap, project->VlidInfil, project->VlidIn, project->VlidOut, & project->VlidDrain)
    if ( project->Subcatch[j].lidArea > 0.0 )
    {
        lid_getRunoff(project,j, tStep);
    }

    // --- update groundwater levels & flows if applicable
    if ( !project->IgnoreGwater && project->Subcatch[j].groundwater )
    {
        gwater_getGroundwater(project,j, project->Vpevap, project->Vinfil+project->VlidInfil, tStep);             //(5.1.010)
    }

    // --- save subcatchment's total loss rates (ft/s)
    area = project->Subcatch[j].area;
    project->Subcatch[j].evapLoss = project->Vevap / tStep / area;
    project->Subcatch[j].infilLoss = (project->Vinfil + project->VlidInfil) / tStep / area;

    // --- find net surface runoff volume
    //     (project->VlidDrain accounts for LID drain flows)
    vOutflow = project->Voutflow      // runoff from all non-LID areas
               - project->VlidIn      // runoff treated by LID units
               + project->VlidOut;    // runoff from LID units
    project->Subcatch[j].newRunoff = vOutflow / tStep;

    // --- obtain external precip. volume (without any snowmelt)
    vRain = project->Subcatch[j].rainfall * tStep * area;

    // --- update the cumulative stats for this subcatchment
    stats_updateSubcatchStats(project,j, vRain, vRunon, project->Vevap, project->Vinfil + project->VlidInfil,
                              vOutflow + project->VlidDrain,
                              project->Subcatch[j].newRunoff + project->VlidDrain/tStep);

    // --- include this subcatchment's contribution to overall flow balance
    //     only if its outlet is a drainage system node
    if ( project->Subcatch[j].outNode == -1 && project->Subcatch[j].outSubcatch != j )
    {
        vOutflow = 0.0;
    }

    // --- update mass balances
    massbal_updateRunoffTotals(project,RUNOFF_RAINFALL, vRain);
    massbal_updateRunoffTotals(project,RUNOFF_EVAP, project->Vevap);
    massbal_updateRunoffTotals(project,RUNOFF_INFIL, project->Vinfil+project->VlidInfil);
    massbal_updateRunoffTotals(project,RUNOFF_RUNOFF, vOutflow);

    // --- return area-averaged runoff (ft/s)
    return runoff / area;
}

//=============================================================================

void getNetPrecip(Project* project, int j, double* netPrecip, double tStep)
{
//
//  Purpose: Finds combined rainfall + snowmelt on a subcatchment.
//  Input:   j = subcatchment index
//           tStep = time step (sec)
//  Output:  netPrecip = rainfall + snowmelt over each type of subarea (ft/s)
//
    int    i, k;
    double rainfall = 0.0;             // rainfall (ft/sec)
    double snowfall = 0.0;             // snowfall (ft/sec)

    // --- get current rainfall or snowfall from rain gage (in ft/sec)
    k = project->Subcatch[j].gage;
    if ( k >= 0 )
    {
        gage_getPrecip(project, k, &rainfall, &snowfall);
    }

	double value = 0;
	if (containsSubcatchRain(project, j, &value))
	{
		rainfall = value;
	}

    // --- assign total precip. rate to subcatch's rainfall property
    project->Subcatch[j].rainfall = rainfall + snowfall;

    // --- determine net precipitation input (netPrecip) to each sub-area

    // --- if subcatch has a snowpack, then base netPrecip on possible snow melt
    if ( project->Subcatch[j].snowpack && !project->IgnoreSnowmelt )
    {
        project->Subcatch[j].newSnowDepth = 
            snow_getSnowMelt(project,j, rainfall, snowfall, tStep, netPrecip);
    }

    // --- otherwise netPrecip is just sum of rainfall & snowfall
    else
    {
        for (i=IMPERV0; i<=PERV; i++) netPrecip[i] = rainfall + snowfall;
    }
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double subcatch_getDepth(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  returns average depth of ponded water (ft)
//  Purpose: finds average depth of water over the non-LID portion of a
//           subcatchment
//
{
    int    i;
    double fArea;
    double depth = 0.0;

    for (i = IMPERV0; i <= PERV; i++)
    {
        fArea = project->Subcatch[j].subArea[i].fArea;
        if ( fArea > 0.0 ) depth += project->Subcatch[j].subArea[i].depth * fArea;
    }
    return depth;
}

//=============================================================================

double subcatch_getWtdOutflow(Project* project, int j, double f)
//
//  Input:   j = subcatchment index
//           f = weighting factor.
//  Output:  returns weighted runoff value
//  Purpose: computes wtd. combination of old and new subcatchment runoff.
//
{
    if ( project->Subcatch[j].area == 0.0 ) return 0.0;
    return (1.0 - f) * project->Subcatch[j].oldRunoff + f * project->Subcatch[j].newRunoff;
}

//=============================================================================

void  subcatch_getResults(Project* project, int j, double f, float x[])
//
//  Input:   j = subcatchment index
//           f = weighting factor
//  Output:  x = array of results
//  Purpose: computes wtd. combination of old and new subcatchment results.
//
{
    int    p;                          // pollutant index
    int    k;                          // rain gage index
    double f1 = 1.0 - f;
    double z;
    double runoff;
    TGroundwater* gw;                  // ptr. to groundwater object

    // --- retrieve rainfall for current report period
    k = project->Subcatch[j].gage;
    if ( k >= 0 ) x[SUBCATCH_RAINFALL] = (float)project->Gage[k].reportRainfall;
    else          x[SUBCATCH_RAINFALL] = 0.0f;

    // --- retrieve snow depth
    z = ( f1 * project->Subcatch[j].oldSnowDepth +
          f * project->Subcatch[j].newSnowDepth ) * UCF(project,RAINDEPTH);
    x[SUBCATCH_SNOWDEPTH] = (float)z;

    // --- retrieve runoff and losses
	x[SUBCATCH_EVAP] = (float)(project->Subcatch[j].evapLoss * UCF(project, EVAPRATE));
	x[SUBCATCH_INFIL] = (float)(project->Subcatch[j].infilLoss * UCF(project, RAINFALL));
    runoff = f1 * project->Subcatch[j].oldRunoff + f * project->Subcatch[j].newRunoff;

////  Following code segement added to release 5.1.008.  ////                  //(5.1.008)
////
    // --- add any LID drain flow to reported runoff
    if ( project->Subcatch[j].lidArea > 0.0 )
    {
		runoff += f1 * lid_getDrainFlow(project, j, PREVIOUS) +
			f * lid_getDrainFlow(project, j, CURRENT);
    }
////

    // --- if runoff is really small, report it as zero
    if ( runoff < MIN_RUNOFF * project->Subcatch[j].area ) runoff = 0.0;                //(5.1.008)
	x[SUBCATCH_RUNOFF] = (float)(runoff * UCF(project, FLOW));

    // --- retrieve groundwater results
    gw = project->Subcatch[j].groundwater;
    if ( gw )
    {
		z = (f1 * gw->oldFlow + f * gw->newFlow) * project->Subcatch[j].area * UCF(project, FLOW);
        x[SUBCATCH_GW_FLOW] = (float)z;
		z = (project->Aquifer[gw->aquifer].bottomElev + gw->lowerDepth) * UCF(project, LENGTH);
        x[SUBCATCH_GW_ELEV] = (float)z;
        z = gw->theta;
        x[SUBCATCH_SOIL_MOIST] = (float)z;
    }
    else
    {
        x[SUBCATCH_GW_FLOW] = 0.0f;
        x[SUBCATCH_GW_ELEV] = 0.0f;
        x[SUBCATCH_SOIL_MOIST]  = 0.0f;
    }

    // --- retrieve pollutant washoff
    if ( !project->IgnoreQuality ) for (p = 0; p < project->Nobjects[POLLUT]; p++ )
    {
        if ( runoff == 0.0 ) z = 0.0;
        else z = f1 * project->Subcatch[j].oldQual[p] + f * project->Subcatch[j].newQual[p];
        x[SUBCATCH_WASHOFF+p] = (float)z;
    }
}


//=============================================================================
//                              SUB-AREA METHODS
//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double getSubareaRunoff(Project* project, int j, int i, double area, double precip, double evap,
    double tStep)
//
//  Purpose: computes runoff & losses from a subarea over the current time step.
//  Input:   j = subcatchment index
//           i = subarea index
//           area = sub-area area (ft2)
//           precip = rainfall + snowmelt over subarea (ft/sec)
//           evap = evaporation (ft/sec)
//           tStep = time step (sec)
//  Output:  returns runoff rate from the sub-area (cfs);
//           updates shared variables project->Vinflow, project->Vevap, project->Vpevap, project->Vinfil & project->Voutflow.
//
{
    double    tRunoff;                 // time over which runoff occurs (sec)
    double    surfMoisture;            // surface water available (ft/sec)
    double    surfEvap;                // evap. used for surface water (ft/sec)
    double    infil = 0.0;             // infiltration rate (ft/sec)
    double    runoff = 0.0;            // runoff rate (ft/sec)
    TSubarea* subarea;                 // pointer to subarea being analyzed

    // --- no runoff if no area
    if ( area == 0.0 ) return 0.0;

    // --- assign pointer to current subarea
    subarea = &project->Subcatch[j].subArea[i];

    // --- assume runoff occurs over entire time step
    tRunoff = tStep;

    // --- determine evaporation loss rate
    surfMoisture = subarea->depth / tStep;
    surfEvap = MIN(surfMoisture, evap);

    // --- compute infiltration loss rate
	if (i == PERV) infil = getSubareaInfil(project,j, subarea, precip, tStep);

    // --- add precip to other subarea inflows
    subarea->inflow += precip;
    surfMoisture += subarea->inflow;

    // --- update total inflow, evaporation & infiltration volumes
    project->Vinflow += precip * area * tStep;
    project->Vevap += surfEvap * area * tStep;
    if ( i == PERV ) project->Vpevap += project->Vevap;
    project->Vinfil += infil * area * tStep;

    // --- if losses exceed available moisture then no ponded water remains
    if ( surfEvap + infil >= surfMoisture )
    {
        subarea->depth = 0.0;
    }

    // --- otherwise reduce inflow by losses and update depth
    //     of ponded water and time over which runoff occurs
    else
    {
        subarea->inflow -= surfEvap + infil;
        updatePondedDepth(project,subarea, &tRunoff);
    }

    // --- compute runoff based on updated ponded depth
    runoff = findSubareaRunoff(subarea, tRunoff);

    // --- compute runoff volume leaving subcatchment for mass balance purposes
    //     (fOutlet is the fraction of this subarea's runoff that goes to the
    //     subcatchment outlet as opposed to another subarea of the subcatchment)
    project->Voutflow += subarea->fOutlet * runoff * area * tStep;
    return runoff;
}

//=============================================================================

double getSubareaInfil(Project* project, int j, TSubarea* subarea, double precip, double tStep)
//
//  Purpose: computes infiltration rate at current time step.
//  Input:   j = subcatchment index
//           subarea = ptr. to a subarea
//           precip = rainfall + snowmelt over subarea (ft/sec)
//           tStep = time step (sec)
//  Output:  returns infiltration rate (ft/s)
//
{
    double infil = 0.0;                     // actual infiltration rate (ft/sec)

    // --- compute infiltration rate 
    infil = infil_getInfil(project, j, project->InfilModel, tStep, precip,
                           subarea->inflow, subarea->depth);

    // --- limit infiltration rate by available void space in unsaturated
    //     zone of any groundwater aquifer
    if ( !project->IgnoreGwater && project->Subcatch[j].groundwater )
    {
        infil = MIN(infil,
                    project->Subcatch[j].groundwater->maxInfilVol/tStep);
    }
    return infil;
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double findSubareaRunoff(TSubarea* subarea, double tRunoff)
//
//  Purpose: computes runoff (ft/s) from subarea after current time step.
//  Input:   subarea = ptr. to a subarea
//           tRunoff = time step over which runoff occurs (sec)
//  Output:  returns runoff rate (ft/s)
//
{
    double xDepth = subarea->depth - subarea->dStore;
    double runoff = 0.0;

    if ( xDepth > ZERO )
    {
        // --- case where nonlinear routing is used
        if ( subarea->N > 0.0 )
        {
            runoff = subarea->alpha * pow(xDepth, MEXP);
        }

        // --- case where no routing is used (Mannings N = 0)
        else
        {
            runoff = xDepth / tRunoff;
            subarea->depth = subarea->dStore;
        }
    }
    else
    {    
        runoff = 0.0;
    }
    return runoff;
}

//=============================================================================

void updatePondedDepth(Project* project, TSubarea* subarea, double* dt)
//
//  Input:   subarea = ptr. to a subarea,
//           dt = time step (sec)
//  Output:  dt = time ponded depth is above depression storage (sec)
//  Purpose: computes new ponded depth over subarea after current time step.
//
{
    double ix = subarea->inflow;       // excess inflow to subarea (ft/sec)    //(5.1.008)
    double dx;                         // depth above depression storage (ft)
    double tx = *dt;                   // time over which dx > 0 (sec)

    // --- see if not enough inflow to fill depression storage (dStore)
    if ( subarea->depth + ix*tx <= subarea->dStore )
    {
        subarea->depth += ix * tx;
    }

    // --- otherwise use the ODE solver to integrate flow depth
    else
    {
        // --- if depth < dStore then fill up dStore & reduce time step
        dx = subarea->dStore - subarea->depth;
        if ( dx > 0.0 && ix > 0.0 )
        {
            tx -= dx / ix;
            subarea->depth = subarea->dStore;
        }

        // --- now integrate depth over remaining time step tx
        if ( subarea->alpha > 0.0 && tx > 0.0 )
        {
            theSubarea = subarea;
            odesolve_integrate(project,&(subarea->depth), 1, 0, tx, ODETOL, tx,
                               getDdDt);
        }
        else
        {
            if ( tx < 0.0 ) tx = 0.0;
            subarea->depth += ix * tx;
        }
    }

    // --- do not allow ponded depth to go negative
    if ( subarea->depth < 0.0 ) subarea->depth = 0.0;

    // --- replace original time step with time ponded depth
    //     is above depression storage
    *dt = tx;
}

//=============================================================================

void  getDdDt(Project* project, double t, double* d, double* dddt)
//
//  Input:   t = current time (not used)
//           d = stored depth (ft)
//  Output   dddt = derivative of d with respect to time
//  Purpose: evaluates derivative of stored depth w.r.t. time
//           for the subarea whose runoff is being computed.
//
{
    double ix = theSubarea->inflow;                                            //(5.1.008)
    double rx = *d - theSubarea->dStore;
    if ( rx < 0.0 )
    {
        rx = 0.0;
    }
    else
    {
        rx = theSubarea->alpha * pow(rx, MEXP);
    }
    *dddt = ix - rx;
}

//=============================================================================
