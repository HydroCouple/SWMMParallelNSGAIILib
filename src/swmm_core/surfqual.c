//-----------------------------------------------------------------------------
//   surfqual.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/15  (Build 5.1.008)
//   Author:   L. Rossman
//
//   Subcatchment water quality functions.
//
//   Build 5.1.008:
//   - Pollutant surface buildup and washoff functions were moved here from
//     subcatch.c.
//   - Support for separate accounting of LID drain flows included. 
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <string.h>
#include "headers.h"
#include "lid.h"

//-----------------------------------------------------------------------------
//  Imported variables 
//-----------------------------------------------------------------------------
// Declared in RUNOFF.C
//extern  double*    OutflowLoad;   // exported pollutant mass load
//
//// Volumes (ft3) for a subcatchment over a time step declared in SUBCATCH.C
//extern double      project->Vinfil;        // non-LID infiltration
//extern double      project->Vinflow;       // non-LID precip + snowmelt + runon + ponded water
//extern double      project->Voutflow;      // non-LID runoff to subcatchment's outlet
//extern double      VlidIn;        // inflow to LID units
//extern double      VlidInfil;     // infiltration from LID units
//extern double      VlidOut;       // surface outflow from LID units
//extern double      project->VlidDrain;     // drain outflow from LID units
//extern double      project->VlidReturn;    // LID outflow returned to pervious area

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)   
//-----------------------------------------------------------------------------
//  surfqual_initState         (called from subcatch_initState)
//  surfqual_getWashoff        (called from runoff_execute)
//  surfqual_getBuildup        (called from runoff_execute)
//  surfqual_sweepBuildup      (called from runoff_execute)
//  surfqual_getWtdWashoff     (called from addWetWeatherInflows in routing.c)

//-----------------------------------------------------------------------------
// Function declarations
//-----------------------------------------------------------------------------
static void  findWashoffLoads(Project* project, int j, double runoff);
static void  findPondedLoads(Project* project, int j, double tStep);
static void  findLidLoads(Project* project, int j, double tStep);

//=============================================================================

void surfqual_initState(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: initializes pollutant buildup, ponded mass, and washoff.
//
{
    int p;

    // --- initialize washoff quality
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        project->Subcatch[j].oldQual[p] = 0.0;
        project->Subcatch[j].newQual[p] = 0.0;
        project->Subcatch[j].pondedQual[p] = 0.0;
    }

    // --- initialize pollutant buildup
	landuse_getInitBuildup(project,project->Subcatch[j].landFactor,  project->Subcatch[j].initBuildup,
		project->Subcatch[j].area, project->Subcatch[j].curbLength);
}

//=============================================================================

void surfqual_getBuildup(Project* project, int j, double tStep)
//
//  Input:   j = subcatchment index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: adds to pollutant buildup on subcatchment surface.
//
{
    int     i;                         // land use index
    int     p;                         // pollutant index
    double  f;                         // land use fraction
    double  area;                      // land use area (acres or hectares)
    double  curb;                      // land use curb length (user units)
    double  oldBuildup;                // buildup at start of time step
    double  newBuildup;                // buildup at end of time step

    // --- consider each landuse
    for (i = 0; i < project->Nobjects[LANDUSE]; i++)
    {
        // --- skip landuse if not in subcatch
        f = project->Subcatch[j].landFactor[i].fraction;
        if ( f == 0.0 ) continue;

        // --- get land area (in acres or hectares) & curb length
		area = f * project->Subcatch[j].area * UCF(project, LANDAREA);
        curb = f * project->Subcatch[j].curbLength;

        // --- examine each pollutant
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            // --- see if snow-only buildup is in effect
            if (project->Pollut[p].snowOnly 
            && project->Subcatch[j].newSnowDepth < 0.001/12.0) continue;

            // --- use land use's buildup function to update buildup amount
            oldBuildup = project->Subcatch[j].landFactor[i].buildup[p];        
            newBuildup = landuse_getBuildup(project,i, p, area, curb, oldBuildup,
                         tStep);
            newBuildup = MAX(newBuildup, oldBuildup);
            project->Subcatch[j].landFactor[i].buildup[p] = newBuildup;
            massbal_updateLoadingTotals(project,BUILDUP_LOAD, p, 
                                       (newBuildup - oldBuildup));
       }
    }
}

//=============================================================================

void surfqual_sweepBuildup(Project* project, int j, DateTime aDate)
//
//  Input:   j = subcatchment index
//           aDate = current date/time
//  Output:  none
//  Purpose: reduces pollutant buildup over a subcatchment if sweeping occurs.
//
{
    int     i;                         // land use index
    int     p;                         // pollutant index
    double  oldBuildup;                // buildup before sweeping (lbs or kg)
    double  newBuildup;                // buildup after sweeping (lbs or kg)

    // --- no sweeping if there is snow on plowable impervious area
    if ( project->Subcatch[j].snowpack != NULL &&
         project->Subcatch[j].snowpack->wsnow[IMPERV0] > MIN_TOTAL_DEPTH ) return;

    // --- consider each land use
    for (i = 0; i < project->Nobjects[LANDUSE]; i++)
    {
        // --- skip land use if not in subcatchment 
        if ( project->Subcatch[j].landFactor[i].fraction == 0.0 ) continue;

        // --- see if land use is subject to sweeping
        if ( project->Landuse[i].sweepInterval == 0.0 ) continue;

        // --- see if sweep interval has been reached
        if ( aDate - project->Subcatch[j].landFactor[i].lastSwept >=
            project->Landuse[i].sweepInterval )
        {
            // --- update time when last swept
            project->Subcatch[j].landFactor[i].lastSwept = aDate;

            // --- examine each pollutant
            for (p = 0; p < project->Nobjects[POLLUT]; p++)
            {
                // --- reduce buildup by the fraction available
                //     times the sweeping effic.
                oldBuildup = project->Subcatch[j].landFactor[i].buildup[p];
                newBuildup = oldBuildup * (1.0 - project->Landuse[i].sweepRemoval *
                             project->Landuse[i].washoffFunc[p].sweepEffic);
                newBuildup = MIN(oldBuildup, newBuildup);
                newBuildup = MAX(0.0, newBuildup);
                project->Subcatch[j].landFactor[i].buildup[p] = newBuildup;

                // --- update mass balance totals
                massbal_updateLoadingTotals(project,SWEEPING_LOAD, p,
                                            oldBuildup - newBuildup);
            }
        }
    }
}

//=============================================================================

void  surfqual_getWashoff(Project* project, int j, double runoff, double tStep)
//
//  Input:   j = subcatchment index
//           runoff = total subcatchment runoff before internal re-routing or
//                    LID controls (ft/sec)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: computes new runoff quality for a subcatchment.
//
//  Considers three pollutant generating streams that are combined together:
//  1. washoff of pollutant buildup as described by the project's land use
//     washoff functions,
//  2. complete mix mass balance of pollutants in surface ponding on
//     non-LID area due to runon, wet deposition, infiltration, & evaporation,
//  3. wet deposition and runon over LID areas.
//
{
    int    p;                // pollutant index
    int    hasOutflow;       // TRUE if subcatchment has outflow
    double cOut;             // final washoff concentration (mass/ft3)
    double massLoad;         // pollut. mass load (mass)
    double vLidRain;         // rainfall volume on LID area (ft3)
    double vLidRunon;        // external runon volume to LID area (ft3)
    double vSurfOut;         // surface runoff volume leaving subcatchment (ft3)
    double vOut1;            // runoff volume prior to LID treatment (ft3)
    double vOut2;            // runoff volume after LID treatment (ft3)
    double area;             // subcatchment area (ft2)

    // --- return if there is no area or no pollutants
    area = project->Subcatch[j].area;
    if ( project->Nobjects[POLLUT] == 0 || area == 0.0 ) return;

    // --- find contributions from washoff, runon and wet precip. to OutflowLoad
    for (p = 0; p < project->Nobjects[POLLUT]; p++) project->OutflowLoad[p] = 0.0;
    findWashoffLoads(project,j, runoff);
    findPondedLoads(project,j, tStep);
    findLidLoads(project,j, tStep);

    // --- contribution from direct rainfall on LID areas
    vLidRain = project->Subcatch[j].rainfall * project->Subcatch[j].lidArea * tStep;

    // --- contribution from upstream runon onto LID areas
    //     (only if LIDs occupy full subcatchment)
    vLidRunon = 0.0;
    if ( area == project->Subcatch[j].lidArea )
    {
        vLidRunon = project->Subcatch[j].runon * area * tStep;
    }

    // --- runoff volume before LID treatment (ft3)
    //     (project->Voutflow, computed in subcatch_getRunoff, is subcatchment
    //      runoff volume before LID treatment)
    vOut1 = project->Voutflow + vLidRain + vLidRunon;             

    // --- surface runoff + LID drain flow volume leaving the subcatchment
    //     (Subcatch.newRunoff, computed in subcatch_getRunoff, includes
    //      any surface runoff reduction from LID treatment)
    vSurfOut = project->Subcatch[j].newRunoff * tStep;
    vOut2 = vSurfOut + project->VlidDrain;

    // --- determine if subcatchment outflow is below a small cutoff
    hasOutflow = (vOut2 > MIN_RUNOFF * area * tStep);

    // --- for each pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- convert washoff load to a concentration
        cOut = 0.0;
        if ( vOut1 > 0.0 ) cOut = project->OutflowLoad[p] / vOut1;

        // --- assign any difference between pre- and post-LID
        //     loads (with LID return flow included) to BMP removal
        if ( project->Subcatch[j].lidArea > 0.0 )
        {    
            massLoad = cOut * (vOut1 - vOut2 - project->VlidReturn) * project->Pollut[p].mcf;
            massbal_updateLoadingTotals(project,BMP_REMOVAL_LOAD, p, massLoad);
        }

        // --- update subcatchment's cumulative runoff load in lbs (or kg)
        massLoad = cOut * vOut2 * project->Pollut[p].mcf;
        project->Subcatch[j].totalLoad[p] += massLoad;

        // --- update mass balance for surface runoff load routed to a
        //     conveyance system node
        //     (loads from LID drains are accounted for below since they
        //     can go to different outlets than parent subcatchment)
        if ( (project->Subcatch[j].outNode >= 0 || project->Subcatch[j].outSubcatch == j) )
        {
            massLoad = cOut * vSurfOut * project->Pollut[p].mcf;
            massbal_updateLoadingTotals(project,RUNOFF_LOAD, p, massLoad);
        }
        
        // --- save new washoff concentration
        if ( !hasOutflow ) cOut = 0.0;
        project->Subcatch[j].newQual[p] = cOut / LperFT3;
    }

    // --- add contribution of LID drain flows to mass balance 
    if ( project->Subcatch[j].lidArea > 0.0 )
    {
        lid_addDrainLoads(project,j, project->Subcatch[j].newQual, tStep);
    }
}

//=============================================================================

double surfqual_getWtdWashoff(Project* project, int j, int p, double f)
//
//  Input:   j = subcatchment index
//           p = pollutant index
//           f = weighting factor
//  Output:  returns pollutant washoff value
//  Purpose: finds wtd. combination of old and new washoff for a pollutant.
//
{
    return (1.0 - f) * project->Subcatch[j].oldRunoff * project->Subcatch[j].oldQual[p] +
           f * project->Subcatch[j].newRunoff *project->Subcatch[j].newQual[p];
}

//=============================================================================

void findPondedLoads(Project* project, int j, double tStep)
//
//  Input:   j = subcatchment index
//           tStep = time step (sec)
//  Output:  updates pondedQual and OutflowLoad 
//  Purpose: mixes wet deposition and runon pollutant loading with existing
//           ponded pollutant mass to compute an ouflow loading.
//
{
    int    p;           // pollutant index
    double cPonded,     // ponded concentration (mass/ft3)
           wPonded,     // pollutant mass in ponded water (mass)
           bmpRemoval,  // load reduction by best mgt. practice (mass)
           vRain,       // volume of direct precipitation (ft3)
           wRain,       // wet deposition pollutant load (mass)
           wRunon,      // external runon pollutant load (mass)
           wInfil,      // pollutant load lost to infiltration (mass)
           wOutflow,    // ponded water contribution to runoff load (mass)
           fullArea,    // full subcatchment area (ft2)
           nonLidArea;  // non-LID area (ft2)

    // --- subcatchment and non-LID areas
    if ( project->Subcatch[j].area == project->Subcatch[j].lidArea ) return;
    fullArea = project->Subcatch[j].area;
    nonLidArea = fullArea - project->Subcatch[j].lidArea;

    // --- compute precip. volume over time step (ft3)
    vRain = project->Subcatch[j].rainfall * nonLidArea * tStep;

    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- update mass balance for wet deposition
        wRain = project->Pollut[p].pptConcen * LperFT3 * vRain;
        massbal_updateLoadingTotals(project,DEPOSITION_LOAD, p, wRain * project->Pollut[p].mcf);

        // --- surface is dry and has no runon -- add any remaining mass
        //     to overall mass balance's FINAL_LOAD category
        if ( project->Vinflow == 0.0 )
        {
            massbal_updateLoadingTotals(project,FINAL_LOAD, p,
                project->Subcatch[j].pondedQual[p] * project->Pollut[p].mcf);
            project->Subcatch[j].pondedQual[p] = 0.0;
        }
        else
        {
            // --- find concen. of ponded water
            //     (newQual[] temporarily holds runon mass loading)
            wRunon = project->Subcatch[j].newQual[p] * tStep;
            wPonded = project->Subcatch[j].pondedQual[p] + wRain + wRunon;
            cPonded = wPonded / project->Vinflow;

            // --- mass lost to infiltration
            wInfil = cPonded * project->Vinfil;
            wInfil = MIN(wInfil, wPonded);
            massbal_updateLoadingTotals(project,INFIL_LOAD, p, wInfil * project->Pollut[p].mcf);
            wPonded -= wInfil;

            // --- mass lost to runoff
            wOutflow = cPonded * project->Voutflow;
            wOutflow = MIN(wOutflow, wPonded);
            wPonded -= wOutflow;

            // --- reduce outflow load by average BMP removal
            bmpRemoval = landuse_getAvgBmpEffic(project,j, p) * wOutflow;
            massbal_updateLoadingTotals(project,BMP_REMOVAL_LOAD, p,
                bmpRemoval*project->Pollut[p].mcf);
            wOutflow -= bmpRemoval;

            // --- update ponded mass (using newly computed ponded depth)
			project->Subcatch[j].pondedQual[p] = cPonded * subcatch_getDepth(project, j) * nonLidArea;
            project->OutflowLoad[p] += wOutflow;
        }
    }
}

//=============================================================================

void  findWashoffLoads(Project* project, int j, double runoff)
//
//  Input:   j = subcatchment index
//           runoff = subcatchment runoff before internal re-routing or
//                    LID controls (ft/sec)
//  Output:  updates OutflowLoad array
//  Purpose: computes pollutant washoff loads for each land use and adds these
//           to the subcatchment's total outflow loads.
//
{
    int    i,                          // land use index
           p,                          // pollutant index
           k;                          // co-pollutant index
    double w,                          // co-pollutant load (mass)
           area = project->Subcatch[j].area;    // subcatchment area (ft2)
    
    // --- examine each land use
    if ( runoff < MIN_RUNOFF ) return;
    for (i = 0; i < project->Nobjects[LANDUSE]; i++)
    {
        if ( project->Subcatch[j].landFactor[i].fraction > 0.0 )
        {
            // --- compute load generated by washoff function
            for (p = 0; p < project->Nobjects[POLLUT]; p++)
            {
                project->OutflowLoad[p] += landuse_getWashoffLoad(project,
                    i, p, area, project->Subcatch[j].landFactor, runoff, project->Voutflow);
            }
        }
    }

    // --- compute contribution from any co-pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- check if pollutant p has a co-pollutant k
        k = project->Pollut[p].coPollut;
        if ( k >= 0 )
        {
            // --- compute addition to washoff from co-pollutant
            w = project->Pollut[p].coFraction * project->OutflowLoad[k];

            // --- add this washoff to buildup mass balance totals
            //     so that things will balance
            massbal_updateLoadingTotals(project,BUILDUP_LOAD, p, w * project->Pollut[p].mcf);

            // --- then also add it to the total washoff load
            project->OutflowLoad[p] += w;
        }
    }
}

//=============================================================================

void  findLidLoads(Project* project, int j, double tStep)
//
//  Input:   j = subcatchment index
//           tStep = time step (sec)
//  Output:  updates OutflowLoad array
//  Purpose: finds addition to subcatchment pollutant loads from wet deposition 
//           and upstream runon to LID areas.
//
{
    int    p;                // pollutant index
    int    useRunon;         // = 1 if LIDs receive upstream runon loads
    double lidArea,          // area occupied by LID units (ft2)
           vLidRain,         // direct precip. falling on LID areas (ft3)
           wLidRain,         // wet deposition pollut. load on LID areas (mass)
           wLidRunon;        // runon pollut. load seen by LID areas (mass)

    // --- find rainfall volume seen by LIDs
    lidArea = project->Subcatch[j].lidArea;
    if ( lidArea == 0.0 ) return;
    vLidRain = project->Subcatch[j].rainfall * lidArea * tStep;

    // --- use upstream runon load only if LIDs occupy full subcatchment
    //     (for partial LID coverage, runon loads were directed onto non-LID area)
    useRunon = (lidArea == project->Subcatch[j].area);

    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- wet deposition load on LID area
        wLidRain = project->Pollut[p].pptConcen * vLidRain * LperFT3;
        massbal_updateLoadingTotals(project,DEPOSITION_LOAD, p, wLidRain * project->Pollut[p].mcf);

        // --- runon load to LID area from other subcatchments
        if ( useRunon ) wLidRunon = project->Subcatch[j].newQual[p] * tStep;
        else            wLidRunon = 0.0;

        // --- update total outflow pollutant load (mass)
        project->OutflowLoad[p] += wLidRain + wLidRunon;
    }
}
