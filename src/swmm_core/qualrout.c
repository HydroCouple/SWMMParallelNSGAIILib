//-----------------------------------------------------------------------------
//   qualrout.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             04/02/15   (Build 5.1.008)
//             04/30/15   (Build 5.1.009)
//             08/05/15   (Build 5.1.010)
//   Author:   L. Rossman
//
//   Water quality routing functions.
//
//   Build 5.1.008:
//   - Pollutant mass lost to seepage flow added to mass balance totals.
//   - Pollutant concen. increased when evaporation occurs.
//
//   Build 5.1.009:
//   - Criterion for dry link/storage node changed to avoid concen. blowup.
//
//   Build 5.1.010:
//   - Entire module re-written to be more compact and easier to follow.
//   - Neglible depth limit replaced with a negligible volume limit.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const double ZeroVolume = 0.0353147; // 1 liter in ft3

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  qualrout_init            (called by swmm_start)
//  qualrout_execute         (called by routing_execute)

//-----------------------------------------------------------------------------
//  Function declarations
//-----------------------------------------------------------------------------
static void  findLinkMassFlow(Project* project, int i, double tStep);
static void  findNodeQual(Project* project, int j);
static void  findLinkQual(Project* project, int i, double tStep);
static void  findSFLinkQual(Project* project, int i, double qSeep, double fEvap, double tStep);
static void  findStorageQual(Project* project, int j, double tStep);
static void  updateHRT(Project* project, int j, double v, double q, double tStep);
static double getReactedQual(Project* project, int p, double c, double v1, double tStep);
static double getMixedQual(double c, double v1, double wIn, double qIn,
              double tStep);
//=============================================================================

void    qualrout_init(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes water quality concentrations in all nodes and links.
//
{
    int     i, p, isWet;
    double  c;

    for (i = 0; i < project->Nobjects[NODE]; i++)
    {
        isWet = ( project->Node[i].newDepth > FUDGE );
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            c = 0.0;
            if ( isWet ) c = project->Pollut[p].initConcen;
            project->Node[i].oldQual[p] = c;
            project->Node[i].newQual[p] = c;
        }
    }

    for (i = 0; i < project->Nobjects[LINK]; i++)
    {
        isWet = ( project->Link[i].newDepth > FUDGE );
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            c = 0.0;
            if ( isWet ) c = project->Pollut[p].initConcen;
            project->Link[i].oldQual[p] = c;
            project->Link[i].newQual[p] = c;
        }
    }
}

//=============================================================================

void qualrout_execute(Project* project, double tStep)
//
//  Input:   tStep = routing time step (sec)
//  Output:  none
//  Purpose: routes water quality constituents through the drainage
//           network over the current time step.
//
{
    int    i, j;
    double qIn, vAvg;

    // --- find mass flow each link contributes to its downstream node
    for ( i = 0; i < project->Nobjects[LINK]; i++ ) findLinkMassFlow(project,i, tStep);

    // --- find new water quality concentration at each node  
    for (j = 0; j < project->Nobjects[NODE]; j++)
    {
        // --- get node inflow and average volume
        qIn = project->Node[j].inflow;
        vAvg = (project->Node[j].oldVolume + project->Node[j].newVolume) / 2.0;
        
        // --- save inflow concentrations if treatment applied
        if ( project->Node[j].treatment )
        {
            if ( qIn < ZERO ) qIn = 0.0;
            treatmnt_setInflow(project, qIn, project->Node[j].newQual);
        }
       
        // --- find new quality at the node 
        if ( project->Node[j].type == STORAGE || project->Node[j].oldVolume > FUDGE )
        {
            findStorageQual(project,j, tStep);
        }
        else findNodeQual(project,j);

        // --- apply treatment to new quality values
		if (project->Node[j].treatment) treatmnt_treat(project, j, qIn, vAvg, tStep);
    }

    // --- find new water quality in each link
    for ( i = 0; i < project->Nobjects[LINK]; i++ ) findLinkQual(project,i, tStep);
}

//=============================================================================

double getMixedQual(double c, double v1, double wIn, double qIn, double tStep)
//
//  Input:   c = concentration in reactor at start of time step (mass/ft3)
//           v1 = volume in reactor at start of time step (ft3)
//           wIn = mass inflow rate (mass/sec)
//           qIn = flow inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  returns pollutant concentration at end of time step (mass/ft3)
//  Purpose: finds pollutant concentration within a completely mixed reactor.
//
{
    double vIn, cIn, cMax;

    // --- if no inflow then reactor concentration is unchanged
    if ( qIn <= ZERO ) return c;

    // --- compute concentration of any inflow
    vIn = qIn * tStep;
    cIn = wIn * tStep / vIn;

    // --- mixture concen. can't exceed either original or inflow concen.
    cMax = MAX(c, cIn);

    // --- mix inflow with current reactor contents
    c = (c*v1 + wIn*tStep) / (v1 + vIn);
    c = MIN(c, cMax);
    c = MAX(c, 0.0);
    return c;
}


//=============================================================================

void findLinkMassFlow(Project* project, int i, double tStep)
//
//  Input:   i = link index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: adds constituent mass flow out of link to the total
//           accumulation at the link's downstream node.
//
//  Note:    project->Node[].newQual[], the accumulator variable, already contains
//           contributions from runoff and other external inflows from
//           calculations made in routing_execute().
{
    int    j, p;
    double qLink, w;

    // --- find inflow to downstream node
    qLink = project->Link[i].newFlow;

    // --- identify index of downstream node
    j = project->Link[i].node2;
    if ( qLink < 0.0 ) j = project->Link[i].node1;
    qLink = fabs(qLink);

    // --- examine each pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- temporarily accumulate inflow load in project->Node[j].newQual
        w = qLink * project->Link[i].oldQual[p];
        project->Node[j].newQual[p] += w;

        // --- update total load transported by link
        project->Link[i].totalLoad[p] += w * tStep;
    }
}

//=============================================================================

void findNodeQual(Project* project, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: finds new quality in a node with no storage volume.
//
{
    int    p;
    double qNode;

    // --- if there is flow into node then concen. = mass inflow/node flow
    qNode = project->Node[j].inflow;
    if ( qNode > ZERO )
    {
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            project->Node[j].newQual[p] /= qNode;
        }
    }

    // --- otherwise concen. is 0
    else for (p = 0; p < project->Nobjects[POLLUT]; p++) project->Node[j].newQual[p] = 0.0;
}

//=============================================================================

void findLinkQual(Project* project, int i, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step.
//
{
    int    j,                // upstream node index
           k,                // conduit index
           p;                // pollutant index
    double wIn,              // pollutant mass inflow rate (mass/sec)
           qIn,              // inflow rate (cfs)
           qSeep,            // rate of seepage loss (cfs)
           v1,               // link volume at start of time step (ft3)
           v2,               // link volume at end of time step (ft3)
           c1,               // current concentration within link (mass/ft3)
           c2,               // new concentration within link (mass/ft3)
           vEvap,            // volume lost to evaporation (ft3)
           vLosses,          // evap. + seepage volume loss (ft3)
           fEvap,            // evaporation concentration factor
           barrels;          // number of barrels in conduit

    // --- identify index of upstream node
    j = project->Link[i].node1;
    if ( project->Link[i].newFlow < 0.0 ) j = project->Link[i].node2;

    // --- link quality is that of upstream node when
    //     link is not a conduit or is a dummy link
    if ( project->Link[i].type != CONDUIT || project->Link[i].xsect.type == DUMMY )
    {
        for (p = 0; p < project->Nobjects[POLLUT]; p++)
        {
            project->Link[i].newQual[p] = project->Node[j].newQual[p];
        }
        return;
    }

    // --- get flow rates and evaporation loss
    k = project->Link[i].subIndex;
    barrels = project->Conduit[k].barrels;
    qIn  = fabs(project->Conduit[k].q1) * barrels;
    qSeep = project->Conduit[k].seepLossRate * barrels;
    vEvap = project->Conduit[k].evapLossRate * barrels * tStep;

    // --- get starting and ending volumes
    v1 = project->Link[i].oldVolume;
    v2 = project->Link[i].newVolume;
    vLosses = qSeep*tStep + vEvap;

    // --- compute factor by which concentrations are increased due to
    //     evaporation loss 
    fEvap = 1.0;
    if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;

    // --- Steady Flow routing requires special treatment
    if ( project->RouteModel == SF )
    {
        findSFLinkQual(project,i, qSeep, fEvap, tStep);
        return;
    }

    // --- adjust inflow to compensate for volume change under Dynamic
    //     Wave routing (which produces just a single (out)flow rate
    //     for a conduit)
    if ( project->RouteModel == DW )
    {
        qIn = qIn + (v2 + vLosses - v1) / tStep; 
        qIn = MAX(qIn, 0.0);
    }

    // --- examine each pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step
        c1 = project->Link[i].oldQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(project,p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- reduce concen. by 1st-order reaction
        c2 = getReactedQual(project,p, c1, v1, tStep);

        // --- mix resulting contents with inflow from upstream node
        wIn = project->Node[j].newQual[p]*qIn;
        c2 = getMixedQual(c2, v1, wIn, qIn, tStep);

        // --- set concen. to zero if remaining volume is negligible
        if ( v2 < ZeroVolume )
        {
            massbal_addToFinalStorage(project,p, c2 * v2);
            c2 = 0.0;
        }

        // --- assign new concen. to link
        project->Link[i].newQual[p] = c2;
    }
}

//=============================================================================

void  findSFLinkQual(Project* project, int i, double qSeep, double fEvap, double tStep)
//
//  Input:   i = link index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a link at end of the current time step for
//           Steady Flow routing.
//
{
    int j = project->Link[i].node1;
    int p;
    double c1, c2;
    double lossRate;

    // --- examine each pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- conduit's quality equals upstream node quality
        c1 = project->Node[j].newQual[p];

        // --- update mass balance accounting for seepage loss
        massbal_addSeepageLoss(project,p, qSeep*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- apply first-order decay over travel time
        c2 = c1;
        if ( project->Pollut[p].kDecay > 0.0 )
        {
            c2 = c1 * exp(-project->Pollut[p].kDecay * tStep);
            c2 = MAX(0.0, c2);
            lossRate = (c1 - c2) * project->Link[i].newFlow;
            massbal_addReactedMass(project,p, lossRate);
        }
        project->Link[i].newQual[p] = c2;
    }
}

//=============================================================================

void  findStorageQual(Project* project, int j, double tStep)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: finds new quality in a node with storage volume.
//  
{
    int    p,                // pollutant index
           k;                // storage unit index
    double qIn,              // inflow rate (cfs)
           wIn,              // pollutant mass inflow rate (mass)
           v1,               // volume at start of time step (ft3)
           c1,               // initial pollutant concentration (mass/ft3)
           c2,               // final pollutant concentration (mass/ft3)
           qExfil = 0.0,     // exfiltration rate from storage unit (cfs)
           vEvap = 0.0,      // evaporation loss from storage unit (ft3)
           fEvap = 1.0;      // evaporation concentration factor

    // --- get inflow rate & initial volume
    qIn = project->Node[j].inflow;
    v1 = project->Node[j].oldVolume;

    // -- for storage nodes
    if ( project->Node[j].type == STORAGE )
    {    
        // --- update hydraulic residence time
        //     (HRT can be used in treatment functions)
        updateHRT(project,j, project->Node[j].oldVolume, qIn, tStep);

        // --- get exfiltration rate and evaporation loss
        k = project->Node[j].subIndex;
        qExfil = project->Storage[k].exfilLoss / tStep;
        vEvap = project->Storage[k].evapLoss;

        // --- compute factor by which concentrations are increased due to
        //     evaporation loss (avoiding huge factors as storage unit
        //     dries out completely)
        if ( vEvap > 0.0 && v1 > ZeroVolume ) fEvap += vEvap / v1;
    }

    // --- for each pollutant
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- start with concen. at start of time step 
        c1 = project->Node[j].oldQual[p];

        // --- update mass balance accounting for exfiltration loss
        massbal_addSeepageLoss(project,p, qExfil*c1);

        // --- increase concen. by evaporation factor
        c1 *= fEvap;

        // --- apply first order reaction only if no separate treatment function
        if ( project->Node[j].treatment == NULL ||
             project->Node[j].treatment[p].equation == NULL )
        {
            c1 = getReactedQual(project,p, c1, v1, tStep);
        }

        // --- mix resulting contents with inflow from all sources
        //     (temporarily accumulated in project->Node[j].newQual)
        wIn = project->Node[j].newQual[p];
        c2 = getMixedQual(c1, v1, wIn, qIn, tStep);

        // --- set concen. to zero if remaining volume is negligible
        if ( project->Node[j].newVolume <= ZeroVolume )
        {
            massbal_addToFinalStorage(project,p, c2 * project->Node[j].newVolume);
            c2 = 0.0;
        }

        // --- assign new concen. to node
        project->Node[j].newQual[p] = c2;
    }
}

//=============================================================================

void updateHRT(Project* project, int j, double v, double q, double tStep)
//
//  Input:   j = node index
//           v = storage volume (ft3)
//           q = inflow rate (cfs)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: updates hydraulic residence time (i.e., water age) at a 
//           storage node.
//
{
    int    k = project->Node[j].subIndex;
    double hrt = project->Storage[k].hrt;
    if ( v < ZERO ) hrt = 0.0;
    else hrt = (hrt + tStep) * v / (v + q*tStep);
    project->Storage[k].hrt = MAX(hrt, 0.0);
}

//=============================================================================

double getReactedQual(Project* project, int p, double c, double v1, double tStep)
//
//  Input:   p = pollutant index
//           c = initial concentration (mass/ft3)
//           v1 = initial volume (ft3)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: applies a first order reaction to a pollutant over a given
//           time step.
//
{
    double c2, lossRate;
    double kDecay = project->Pollut[p].kDecay;

    if ( kDecay == 0.0 ) return c;
    c2 = c * (1.0 - kDecay * tStep);
    c2 = MAX(0.0, c2);
    lossRate = (c - c2) * v1 / tStep;
    massbal_addReactedMass(project,p, lossRate);
    return c2;
}
 
