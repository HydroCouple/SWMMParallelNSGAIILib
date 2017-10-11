//-----------------------------------------------------------------------------
//   gwater.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.000)
//             09/15/14  (Build 5.1.007)
//             03/19/15  (Build 5.1.008)
//             08/05/15  (Build 5.1.010)
//   Author:   L. Rossman
//
//   Groundwater functions.
//
//   Build 5.1.007:
//   - User-supplied function for deep project->GW seepage flow added.
//   - New variable names for use in user-supplied project->GW flow equations added.
//
//   Build 5.1.008:
//   - More variable names for user-supplied project->GW flow equations added.
//   - Subcatchment area made into a shared variable.
//   - Evaporation loss initialized to 0.
//   - Support for collecting project->GW statistics added.
//
//   Build 5.1.010:
//   - Unsaturated hydraulic conductivity added to project->GW flow equation variables.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "odesolve.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const double GWTOL = 0.0001;    // ODE solver tolerance
static const double XTOL  = 0.001;     // tolerance for moisture & depth

enum   GWstates {THETA,                // moisture content of upper project->GW zone
                 LOWERDEPTH};          // depth of lower saturated project->GW zone

enum   GWvariables {
	     gwvHGW,                       // water table height (ft)
         gwvHSW,                       // surface water height (ft)
         gwvHCB,                       // channel bottom height (ft)           //(5.1.007)
         gwvHGS,                       // ground surface height (ft)           //(5.1.007)
         gwvKS,                        // sat. hyd. condutivity (ft/s)         //(5.1.007)
         gwvK,                         // unsat. hyd. conductivity (ft/s)      //(5.1.008)
         gwvTHETA,                     // upper zone moisture content          //(5.1.008)
         gwvPHI,                       // soil porosity                        //(5.1.008)
         gwvFI,                        // surface infiltration (ft/s)          //(5.1.008)
         gwvFU,                        // uper zone percolation rate (ft/s)    //(5.1.008)
         gwvA,                         // subcatchment area (ft2)              //(5.1.008)
         gwvMAX};

// Names of project->GW variables that can be used in project->GW outflow expression
static char* GWVarWords[] = {"HGW", "HSW", "HCB", "HGS", "KS", "K",            //(5.1.010)
                             "THETA", "PHI", "FI", "FU", "A", NULL};           //(5.1.008)

////-----------------------------------------------------------------------------
////  Shared variables
////-----------------------------------------------------------------------------
////  NOTE: all flux rates are in ft/sec, all depths are in ft.
//static double    project->Area;            // subcatchment area (ft2)                   //(5.1.008)
//static double    project->Infil;           // infiltration rate from surface
//static double    project->MaxEvap;         // max. evaporation rate
//static double    project->AvailEvap;       // available evaporation rate
//static double    project->UpperEvap;       // evaporation rate from upper project->GW zone
//static double    project->LowerEvap;       // evaporation rate from lower project->GW zone
//static double    project->UpperPerc;       // percolation rate from upper to lower zone
//static double    project->LowerLoss;       // loss rate from lower project->GW zone
//static double    project->GWFlow;          // flow rate from lower zone to conveyance node
//static double    project->MaxUpperPerc;    // upper limit on project->UpperPerc
//static double    project->MaxGWFlowPos;    // upper limit on project->GWFlow when its positve
//static double    project->MaxGWFlowNeg;    // upper limit on project->GWFlow when its negative
//static double    project->FracPerv;        // fraction of surface that is pervious
//static double    project->TotalDepth;      // total depth of project->GW aquifer
//static double    project->Theta;           // moisture content of upper zone
//static double    project->HydCon;          // unsaturated hydraulic conductivity (ft/s) //(5.1.010)
//static double    project->Hgw;             // ht. of saturated zone
//static double    project->Hstar;           // ht. from aquifer bottom to node invert
//static double    project->Hsw;             // ht. from aquifer bottom to water surface
//static double    project->Tstep;           // current time step (sec)
//static TAquifer  project->A;               // aquifer being analyzed
//static TGroundwater* project->GW;          // groundwater object being analyzed
//static MathExpr* project->LatFlowExpr;     // user-supplied lateral project->GW flow expression  //(5.1.007)
//static MathExpr* project->DeepFlowExpr;    // user-supplied deep project->GW flow expression     //(5.1.007)

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  gwater_readAquiferParams     (called by input_readLine)
//  gwater_readGroundwaterParams (called by input_readLine)
//  gwater_readFlowExpression    (called by input_readLine)
//  gwater_deleteFlowExpression  (called by deleteObjects in project.c)
//  gwater_validateAquifer       (called by swmm_open)
//  gwater_validate              (called by subcatch_validate) 
//  gwater_initState             (called by subcatch_initState)
//  gwater_getVolume             (called by massbal_open & massbal_getGwaterError)
//  gwater_getGroundwater        (called by getSubareaRunoff in subcatch.c)
//  gwater_getState              (called by saveRunoff in hotstart.c)
//  gwater_setState              (called by readRunoff in hotstart.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void   getDxDt(Project* project, double t, double* x, double* dxdt);
static void   getFluxes(Project* project, double upperVolume, double lowerDepth);
static void   getEvapRates(Project* project, double theta, double upperDepth);
static double getUpperPerc(Project* project, double theta, double upperDepth);
static double getGWFlow(Project* project, double lowerDepth);
static void   updateMassBal(Project* project, double area,  double tStep);

// Used to process custom project->GW outflow equations
static int    getVariableIndex(Project* project, char* s);
static double getVariableValue(Project* project, int varIndex);

//=============================================================================

int gwater_readAquiferParams(Project* project, int j, char* tok[], int ntoks)
//
//  Input:   j = aquifer index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error message
//  Purpose: reads aquifer parameter values from line of input data
//
//  Data line contains following parameters:
//    ID, porosity, wiltingPoint, fieldCapacity,     conductivity,
//    conductSlope, tensionSlope, upperEvapFraction, lowerEvapDepth,
//    gwRecession,  bottomElev,   waterTableElev,    upperMoisture
//    (evapPattern)
//
{
    int   i, p;
    double x[12];
    char *id;

    // --- check that aquifer exists
    if ( ntoks < 13 ) return error_setInpError(project,ERR_ITEMS, "");
	id = project_findID(project, AQUIFER, tok[0]);
    if ( id == NULL ) return error_setInpError(project,ERR_NAME, tok[0]);

    // --- read remaining tokens as numbers
    for (i = 0; i < 11; i++) x[i] = 0.0;
    for (i = 1; i < 13; i++)
    {
        if ( ! getDouble(tok[i], &x[i-1]) )
            return error_setInpError(project,ERR_NUMBER, tok[i]);
    }

    // --- read upper evap pattern if present
    p = -1;
    if ( ntoks > 13 )
    {
		p = project_findObject(project, TIMEPATTERN, tok[13]);
	if ( p < 0 ) return error_setInpError(project,ERR_NAME, tok[13]);
    }

    // --- assign parameters to aquifer object
    project->Aquifer[j].ID = id;
    project->Aquifer[j].porosity       = x[0];
    project->Aquifer[j].wiltingPoint   = x[1];
    project->Aquifer[j].fieldCapacity  = x[2];
	project->Aquifer[j].conductivity = x[3] / UCF(project, RAINFALL);
    project->Aquifer[j].conductSlope   = x[4];
	project->Aquifer[j].tensionSlope = x[5] / UCF(project, LENGTH);
    project->Aquifer[j].upperEvapFrac  = x[6];
	project->Aquifer[j].lowerEvapDepth = x[7] / UCF(project, LENGTH);
	project->Aquifer[j].lowerLossCoeff = x[8] / UCF(project, RAINFALL);
	project->Aquifer[j].bottomElev = x[9] / UCF(project, LENGTH);
	project->Aquifer[j].waterTableElev = x[10] / UCF(project, LENGTH);
    project->Aquifer[j].upperMoisture  = x[11];
    project->Aquifer[j].upperEvapPat   = p;
    return 0;
}

//=============================================================================

int gwater_readGroundwaterParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads groundwater inflow parameters for a subcatchment from
//           a line of input data.
//
//  Data format is:
//  subcatch  aquifer  node  surfElev  a1  b1  a2  b2  a3  fixedDepth +
//            (nodeElev  bottomElev  waterTableElev  upperMoisture )
//
{
    int    i, j, k, m, n;
    double x[11];
    TGroundwater* gw;

    // --- check that specified subcatchment, aquifer & node exist
    if ( ntoks < 3 ) return error_setInpError(project,ERR_ITEMS, "");
	j = project_findObject(project, SUBCATCH, tok[0]);
    if ( j < 0 ) return error_setInpError(project,ERR_NAME, tok[0]);

    // --- check for enough tokens
    if ( ntoks < 11 ) return error_setInpError(project,ERR_ITEMS, "");

    // --- check that specified aquifer and node exists
	k = project_findObject(project, AQUIFER, tok[1]);
    if ( k < 0 ) return error_setInpError(project,ERR_NAME, tok[1]);
	n = project_findObject(project, NODE, tok[2]);
    if ( n < 0 ) return error_setInpError(project,ERR_NAME, tok[2]);

    // -- read in the flow parameters
    for ( i = 0; i < 7; i++ )
    {
        if ( ! getDouble(tok[i+3], &x[i]) ) 
            return error_setInpError(project,ERR_NUMBER, tok[i+3]);
    }

    // --- read in optional depth parameters
    for ( i = 7; i < 11; i++)
    {
        x[i] = MISSING;
        m = i + 3;
        if ( ntoks > m && *tok[m] != '*' )
        {    
            if (! getDouble(tok[m], &x[i]) ) 
                return error_setInpError(project,ERR_NUMBER, tok[m]);
            if ( i < 10 ) x[i] /= UCF(project,LENGTH);
        }
    }

    // --- create a groundwater flow object
    if ( !project->Subcatch[j].groundwater )
    {
        gw = (TGroundwater *) malloc(sizeof(TGroundwater));
        if ( !gw ) return error_setInpError(project,ERR_MEMORY, "");
        project->Subcatch[j].groundwater = gw;
    }
    else gw = project->Subcatch[j].groundwater;

    // --- populate the groundwater flow object with its parameters
    gw->aquifer    = k;
    gw->node       = n;
    gw->surfElev   = x[0] / UCF(project,LENGTH);
    gw->a1         = x[1];
    gw->b1         = x[2];
    gw->a2         = x[3];
    gw->b2         = x[4];
    gw->a3         = x[5];
    gw->fixedDepth = x[6] / UCF(project,LENGTH);
    gw->nodeElev   = x[7];                       //already converted to ft.
    gw->bottomElev     = x[8];
    gw->waterTableElev = x[9];
    gw->upperMoisture  = x[10];
    return 0;
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

int gwater_readFlowExpression(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads mathematical expression for lateral or deep groundwater
//           flow for a subcatchment from a line of input data.
//
//  Format is: subcatch LATERAL/DEEP <expr>
//     where subcatch is the ID of the subcatchment, LATERAL is for lateral
//     project->GW flow, DEEP is for deep project->GW flow and <expr> is any well-formed math
//     expression. 
//
{
    int   i, j, k;
    char  exprStr[MAXLINE+1];
    MathExpr* expr;

    // --- return if too few tokens
    if ( ntoks < 3 ) return error_setInpError(project,ERR_ITEMS, "");

    // --- check that subcatchment exists
	j = project_findObject(project, SUBCATCH, tok[0]);
    if ( j < 0 ) return error_setInpError(project,ERR_NAME, tok[0]);

    // --- check if expression is for lateral or deep project->GW flow
    k = 1;
    if ( match(tok[1], "LAT") ) k = 1;
    else if ( match(tok[1], "DEEP") ) k = 2;
    else return error_setInpError(project,ERR_KEYWORD, tok[1]);

    // --- concatenate remaining tokens into a single string
    strcpy(exprStr, tok[2]);
    for ( i = 3; i < ntoks; i++)
    {
        strcat(exprStr, " ");
        strcat(exprStr, tok[i]);
    }

    // --- delete any previous flow eqn.
    if ( k == 1 ) mathexpr_delete(project->Subcatch[j].gwLatFlowExpr);
    else          mathexpr_delete(project->Subcatch[j].gwDeepFlowExpr);

    // --- create a parsed expression tree from the string expr
    //     (getVariableIndex is the function that converts a project->GW
    //      variable's name into an index number) 
    expr = mathexpr_create(project,exprStr, getVariableIndex);
    if ( expr == NULL ) return error_setInpError(project,ERR_TREATMENT_EXPR, "");

    // --- save expression tree with the subcatchment
    if ( k == 1 ) project->Subcatch[j].gwLatFlowExpr = expr;
    else          project->Subcatch[j].gwDeepFlowExpr = expr;
    return 0;
}

//=============================================================================

void gwater_deleteFlowExpression(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: deletes a subcatchment's custom groundwater flow expressions.     //(5.1.007)
//
{
    mathexpr_delete(project->Subcatch[j].gwLatFlowExpr);
    mathexpr_delete(project->Subcatch[j].gwDeepFlowExpr);
}

//=============================================================================

void  gwater_validateAquifer(Project* project, int j)
//
//  Input:   j = aquifer index
//  Output:  none
//  Purpose: validates groundwater aquifer properties .
//
{
	int p;

    if ( project->Aquifer[j].porosity          <= 0.0 
    ||   project->Aquifer[j].fieldCapacity     >= project->Aquifer[j].porosity
    ||   project->Aquifer[j].wiltingPoint      >= project->Aquifer[j].fieldCapacity
    ||   project->Aquifer[j].conductivity      <= 0.0
    ||   project->Aquifer[j].conductSlope      <  0.0
    ||   project->Aquifer[j].tensionSlope      <  0.0
    ||   project->Aquifer[j].upperEvapFrac     <  0.0
    ||   project->Aquifer[j].lowerEvapDepth    <  0.0
    ||   project->Aquifer[j].waterTableElev    <  project->Aquifer[j].bottomElev
    ||   project->Aquifer[j].upperMoisture     >  project->Aquifer[j].porosity 
    ||   project->Aquifer[j].upperMoisture     <  project->Aquifer[j].wiltingPoint )
        report_writeErrorMsg(project,ERR_AQUIFER_PARAMS, project->Aquifer[j].ID);

    p = project->Aquifer[j].upperEvapPat;
    if ( p >= 0 && project->Pattern[p].type != MONTHLY_PATTERN )
    {
        report_writeErrorMsg(project,ERR_AQUIFER_PARAMS, project->Aquifer[j].ID);
    }
}

//=============================================================================

void  gwater_validate(Project* project, int j)
{
    TAquifer a;         // Aquifer data structure
    TGroundwater* gw;   // Groundwater data structure
    
    gw = project->Subcatch[j].groundwater;
    if ( gw )
    {
        a = project->Aquifer[gw->aquifer];

        // ... use aquifer values for missing groundwater parameters
        if ( gw->bottomElev == MISSING ) gw->bottomElev = a.bottomElev;
        if ( gw->waterTableElev == MISSING ) gw->waterTableElev = a.waterTableElev;
        if ( gw->upperMoisture == MISSING ) gw->upperMoisture = a.upperMoisture;

        // ... ground elevation can't be below water table elevation
        if ( gw->surfElev < gw->waterTableElev )
            report_writeErrorMsg(project,ERR_GROUND_ELEV, project->Subcatch[j].ID);
    }
}

//=============================================================================

void  gwater_initState(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  none
//  Purpose: initializes state of subcatchment's groundwater.
//
{
    TAquifer a;         // Aquifer data structure
    TGroundwater* gw;   // Groundwater data structure
    
    gw = project->Subcatch[j].groundwater;
    if ( gw )
    {
        a = project->Aquifer[gw->aquifer];

        // ... initial moisture content
        gw->theta = gw->upperMoisture;
        if ( gw->theta >= a.porosity )
        {
            gw->theta = a.porosity - XTOL;
        }

        // ... initial depth of lower (saturated) zone
        gw->lowerDepth = gw->waterTableElev - gw->bottomElev;
        if ( gw->lowerDepth >= gw->surfElev - gw->bottomElev )
        {
            gw->lowerDepth = gw->surfElev - gw->bottomElev - XTOL;
        }

        // ... initial lateral groundwater outflow
        gw->oldFlow = 0.0;
        gw->newFlow = 0.0;
        gw->evapLoss = 0.0;                                                    //(5.1.008)

        // ... initial available infiltration volume into upper zone
        gw->maxInfilVol = (gw->surfElev - gw->waterTableElev) *
                          (a.porosity - gw->theta) /
                          subcatch_getFracPerv(project,j);
    }
}

//=============================================================================

void gwater_getState(Project* project, int j, double x[])
//
//  Input:   j = subcatchment index
//  Output:  x[] = array of groundwater state variables
//  Purpose: retrieves state of subcatchment's groundwater.
//
{
    TGroundwater* gw = project->Subcatch[j].groundwater;
    x[0] = gw->theta;
    x[1] = gw->bottomElev + gw->lowerDepth;
    x[2] = gw->newFlow;
    x[3] = gw->maxInfilVol;
}

//=============================================================================

void gwater_setState(Project* project, int j, double x[])
//
//  Input:   j = subcatchment index
//           x[] = array of groundwater state variables
//  Purpose: assigns values to a subcatchment's groundwater state.
//
{
    TGroundwater* gw = project->Subcatch[j].groundwater;
    if ( gw == NULL ) return;
    gw->theta = x[0];
    gw->lowerDepth = x[1] - gw->bottomElev;
    gw->oldFlow = x[2];
    if ( x[3] != MISSING ) gw->maxInfilVol = x[3];
}

//=============================================================================

double gwater_getVolume(Project* project, int j)
//
//  Input:   j = subcatchment index
//  Output:  returns total volume of groundwater in ft/ft2
//  Purpose: finds volume of groundwater stored in upper & lower zones
//
{
    TAquifer a;
    TGroundwater* gw;
    double upperDepth;
    gw = project->Subcatch[j].groundwater;
    if ( gw == NULL ) return 0.0;
    a = project->Aquifer[gw->aquifer];
    upperDepth = gw->surfElev - gw->bottomElev - gw->lowerDepth;
    return (upperDepth * gw->theta) + (gw->lowerDepth * a.porosity);
}

//=============================================================================

void gwater_getGroundwater(Project* project, int j, double evap, double infil, double tStep)
//
//  Purpose: computes groundwater flow from subcatchment during current time step.
//  Input:   j     = subcatchment index
//           evap  = pervious surface evaporation volume consumed (ft3)
//           infil = surface infiltration volume (ft3)
//           tStep = time step (sec)
//  Output:  none
//

//  Note: local "area" variable was replaced with shared variable "project->Area". //   //(5.1.008)

{
    int    n;                          // node exchanging groundwater
    double x[2];                       // upper moisture content & lower depth 
    double vUpper;                     // upper vol. available for percolation
    double nodeFlow;                   // max. possible project->GW flow from node

    // --- save subcatchment's groundwater and aquifer objects to 
    //     shared variables
    project->GW = project->Subcatch[j].groundwater;
    if ( project->GW == NULL ) return;
    project->LatFlowExpr = project->Subcatch[j].gwLatFlowExpr;                                   //(5.1.007)
    project->DeepFlowExpr = project->Subcatch[j].gwDeepFlowExpr;                                 //(5.1.007)
    project->A = project->Aquifer[project->GW->aquifer];

    // --- get fraction of total area that is pervious
    project->FracPerv = subcatch_getFracPerv(project,j);
    if ( project->FracPerv <= 0.0 ) return;
    project->Area = project->Subcatch[j].area;

    // --- convert infiltration volume (ft3) to equivalent rate
    //     over entire project->GW (subcatchment) area
    infil = infil / project->Area / tStep;
    project->Infil = infil;
    project->Tstep = tStep;

    // --- convert pervious surface evaporation already exerted (ft3)
    //     to equivalent rate over entire project->GW (subcatchment) area
    evap = evap / project->Area / tStep;

    // --- convert max. surface evap rate (ft/sec) to a rate
    //     that applies to project->GW evap (project->GW evap can only occur
    //     through the pervious land surface area)
    project->MaxEvap = project->Evap.rate * project->FracPerv;

    // --- available subsurface evaporation is difference between max.
    //     rate and pervious surface evap already exerted
    project->AvailEvap = MAX((project->MaxEvap - evap), 0.0);

    // --- save total depth & outlet node properties to shared variables
    project->TotalDepth = project->GW->surfElev - project->GW->bottomElev;
    if ( project->TotalDepth <= 0.0 ) return;
    n = project->GW->node;

    // --- establish min. water table height above aquifer bottom at which
    //     project->GW flow can occur (override node's invert if a value was provided
    //     in the project->GW object)
    if ( project->GW->nodeElev != MISSING ) project->Hstar = project->GW->nodeElev - project->GW->bottomElev;
    else project->Hstar = project->Node[n].invertElev - project->GW->bottomElev;
    
    // --- establish surface water height (relative to aquifer bottom)
    //     for drainage system node connected to the project->GW aquifer
    if ( project->GW->fixedDepth > 0.0 )
    {
        project->Hsw = project->GW->fixedDepth + project->Node[n].invertElev - project->GW->bottomElev;
    }
    else project->Hsw = project->Node[n].newDepth + project->Node[n].invertElev - project->GW->bottomElev;

    // --- store state variables (upper zone moisture content, lower zone
    //     depth) in work vector x
    x[THETA] = project->GW->theta;
    x[LOWERDEPTH] = project->GW->lowerDepth;

    // --- set limit on percolation rate from upper to lower project->GW zone
    vUpper = (project->TotalDepth - x[LOWERDEPTH]) * (x[THETA] - project->A.fieldCapacity);
    vUpper = MAX(0.0, vUpper); 
    project->MaxUpperPerc = vUpper / tStep;

    // --- set limit on project->GW flow out of aquifer based on volume of lower zone
    project->MaxGWFlowPos = x[LOWERDEPTH]*project->A.porosity / tStep;

    // --- set limit on project->GW flow into aquifer from drainage system node
    //     based on min. of capacity of upper zone and drainage system
    //     inflow to the node
    project->MaxGWFlowNeg = (project->TotalDepth - x[LOWERDEPTH]) * (project->A.porosity - x[THETA])
                   / tStep;
    nodeFlow = (project->Node[n].inflow + project->Node[n].newVolume/tStep) / project->Area;
    project->MaxGWFlowNeg = -MIN(project->MaxGWFlowNeg, nodeFlow);
    
    // --- integrate eqns. for d(project->Theta)/dt and d(LowerDepth)/dt
    //     NOTE: ODE solver must have been initialized previously
    odesolve_integrate(project,x, 2, 0, tStep, GWTOL, tStep, getDxDt);
    
    // --- keep state variables within allowable bounds
    x[THETA] = MAX(x[THETA], project->A.wiltingPoint);
    if ( x[THETA] >= project->A.porosity )
    {
        x[THETA] = project->A.porosity - XTOL;
        x[LOWERDEPTH] = project->TotalDepth - XTOL;
    }
    x[LOWERDEPTH] = MAX(x[LOWERDEPTH],  0.0);
    if ( x[LOWERDEPTH] >= project->TotalDepth )
    {
        x[LOWERDEPTH] = project->TotalDepth - XTOL;
    }

    // --- save new values of state values
    project->GW->theta = x[THETA];
    project->GW->lowerDepth  = x[LOWERDEPTH];
	getFluxes(project, project->GW->theta, project->GW->lowerDepth);
    project->GW->oldFlow = project->GW->newFlow;
    project->GW->newFlow = project->GWFlow;
    project->GW->evapLoss = project->UpperEvap + project->LowerEvap;

    //--- find max. infiltration volume (as depth over
    //    the pervious portion of the subcatchment)
    //    that upper zone can support in next time step
    project->GW->maxInfilVol = (project->TotalDepth - x[LOWERDEPTH]) *
                      (project->A.porosity - x[THETA]) / project->FracPerv;

    // --- update project->GW mass balance
    updateMassBal(project,project->Area, tStep);

    // --- update project->GW statistics                                                //(5.1.008)
    stats_updateGwaterStats(project, j, infil, project->GW->evapLoss, project->GWFlow, project->LowerLoss,         //(5.1.008)
        project->GW->theta, project->GW->lowerDepth + project->GW->bottomElev, tStep);                    //(5.1.008)
}

//=============================================================================

void updateMassBal(Project* project, double area, double tStep)
//
//  Input:   area  = subcatchment area (ft2)
//           tStep = time step (sec)
//  Output:  none
//  Purpose: updates project->GW mass balance with volumes of water fluxes.
//
{
    double vInfil;                     // infiltration volume
    double vUpperEvap;                 // upper zone evap. volume
    double vLowerEvap;                 // lower zone evap. volume
    double vLowerPerc;                 // lower zone deep perc. volume
    double vGwater;                    // volume of exchanged groundwater
    double ft2sec = area * tStep;

    vInfil     = project->Infil * ft2sec;
    vUpperEvap = project->UpperEvap * ft2sec;
    vLowerEvap = project->LowerEvap * ft2sec;
    vLowerPerc = project->LowerLoss * ft2sec;
    vGwater    = 0.5 * (project->GW->oldFlow + project->GW->newFlow) * ft2sec;
    massbal_updateGwaterTotals(project,vInfil, vUpperEvap, vLowerEvap, vLowerPerc,
                               vGwater);
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void  getFluxes(Project* project, double theta, double lowerDepth)
//
//  Input:   upperVolume = vol. depth of upper zone (ft)
//           upperDepth  = depth of upper zone (ft)
//  Output:  none
//  Purpose: computes water fluxes into/out of upper/lower project->GW zones.
//
{
    double upperDepth;

    // --- find upper zone depth
    lowerDepth = MAX(lowerDepth, 0.0);
    lowerDepth = MIN(lowerDepth, project->TotalDepth);
    upperDepth = project->TotalDepth - lowerDepth;

    // --- save lower depth and theta to global variables
    project->Hgw = lowerDepth;
    project->Theta = theta;

    // --- find evaporation rate from both zones
    getEvapRates(project,theta, upperDepth);

    // --- find percolation rate from upper to lower zone
    project->UpperPerc = getUpperPerc(project,theta, upperDepth);
    project->UpperPerc = MIN(project->UpperPerc, project->MaxUpperPerc);

    // --- find loss rate to deep project->GW
    if ( project->DeepFlowExpr != NULL )
        project->LowerLoss = mathexpr_eval(project, project->DeepFlowExpr, getVariableValue) /
                    UCF(project,RAINFALL);
    else
        project->LowerLoss = project->A.lowerLossCoeff * lowerDepth / project->TotalDepth;
    project->LowerLoss = MIN(project->LowerLoss, lowerDepth/project->Tstep);

    // --- find project->GW flow rate from lower zone to drainage system node
    project->GWFlow = getGWFlow(project,lowerDepth);
    if ( project->LatFlowExpr != NULL )
    {
		project->GWFlow += mathexpr_eval(project, project->LatFlowExpr, getVariableValue) / UCF(project, GWFLOW);
    }
    if ( project->GWFlow >= 0.0 ) project->GWFlow = MIN(project->GWFlow, project->MaxGWFlowPos);
    else project->GWFlow = MAX(project->GWFlow, project->MaxGWFlowNeg);
}

//=============================================================================

void  getDxDt(Project* project, double t, double* x, double* dxdt)
//
//  Input:   t    = current time (not used)
//           x    = array of state variables
//  Output:  dxdt = array of time derivatives of state variables
//  Purpose: computes time derivatives of upper moisture content 
//           and lower depth.
//
{
    double qUpper;    // inflow - outflow for upper zone (ft/sec)
    double qLower;    // inflow - outflow for lower zone (ft/sec)
    double denom;

    getFluxes(project,x[THETA], x[LOWERDEPTH]);
    qUpper = project->Infil - project->UpperEvap - project->UpperPerc;
    qLower = project->UpperPerc - project->LowerLoss - project->LowerEvap - project->GWFlow;

    // --- d(upper zone moisture)/dt = (net upper zone flow) /
    //                                 (upper zone depth)
    denom = project->TotalDepth - x[LOWERDEPTH];
    if (denom > 0.0)
        dxdt[THETA] = qUpper / denom;
    else
        dxdt[THETA] = 0.0;

    // --- d(lower zone depth)/dt = (net lower zone flow) /
    //                              (upper zone moisture deficit)
    denom = project->A.porosity - x[THETA];
    if (denom > 0.0)
        dxdt[LOWERDEPTH] = qLower / denom;
    else
        dxdt[LOWERDEPTH] = 0.0;
}

//=============================================================================

void getEvapRates(Project* project, double theta, double upperDepth)
//
//  Input:   theta      = moisture content of upper zone
//           upperDepth = depth of upper zone (ft)
//  Output:  none
//  Purpose: computes evapotranspiration out of upper & lower zones.
//
{
    int    p, month;
    double f;
    double lowerFrac, upperFrac;

    // --- no project->GW evaporation when infiltration is occurring
    project->UpperEvap = 0.0;
    project->LowerEvap = 0.0;
    if ( project->Infil > 0.0 ) return;

    // --- get monthly-adjusted upper zone evap fraction
    upperFrac = project->A.upperEvapFrac;
    f = 1.0;
    p = project->A.upperEvapPat;
    if ( p >= 0 )
    {
        month = datetime_monthOfYear(getDateTime(project,project->NewRunoffTime));
        f = project->Pattern[p].factor[month-1];
    }
    upperFrac *= f;

    // --- upper zone evaporation requires that soil moisture
    //     be above the wilting point
    if ( theta > project->A.wiltingPoint )
    {
        // --- actual evap is upper zone fraction applied to max. potential
        //     rate, limited by the available rate after any surface evap 
        project->UpperEvap = upperFrac * project->MaxEvap;
        project->UpperEvap = MIN(project->UpperEvap, project->AvailEvap);
    }

    // --- check if lower zone evaporation is possible
    if ( project->A.lowerEvapDepth > 0.0 )
    {
        // --- find the fraction of the lower evaporation depth that
        //     extends into the saturated lower zone
        lowerFrac = (project->A.lowerEvapDepth - upperDepth) / project->A.lowerEvapDepth;
        lowerFrac = MAX(0.0, lowerFrac);
        lowerFrac = MIN(lowerFrac, 1.0);

        // --- make the lower zone evap rate proportional to this fraction
        //     and the evap not used in the upper zone
        project->LowerEvap = lowerFrac * (1.0 - upperFrac) * project->MaxEvap;
        project->LowerEvap = MIN(project->LowerEvap, (project->AvailEvap - project->UpperEvap));
    }
}

//=============================================================================

double getUpperPerc(Project* project, double theta, double upperDepth)
//
//  Input:   theta      = moisture content of upper zone
//           upperDepth = depth of upper zone (ft)
//  Output:  returns percolation rate (ft/sec)
//  Purpose: finds percolation rate from upper to lower zone.
//
{
    double delta;                       // unfilled water content of upper zone
    double dhdz;                        // avg. change in head with depth
    double hydcon;                      // unsaturated hydraulic conductivity

    // --- no perc. from upper zone if no depth or moisture content too low    
    if ( upperDepth <= 0.0 || theta <= project->A.fieldCapacity ) return 0.0;

    // --- compute hyd. conductivity as function of moisture content
    delta = theta - project->A.porosity;
    hydcon = project->A.conductivity * exp(delta * project->A.conductSlope);

    // --- compute integral of dh/dz term
    delta = theta - project->A.fieldCapacity;
    dhdz = 1.0 + project->A.tensionSlope * 2.0 * delta / upperDepth;

    // --- compute upper zone percolation rate
    project->HydCon = hydcon;                                                           //(5.1.010)
    return hydcon * dhdz;
}

//=============================================================================

double getGWFlow(Project* project, double lowerDepth)
//
//  Input:   lowerDepth = depth of lower zone (ft)
//  Output:  returns groundwater flow rate (ft/sec)
//  Purpose: finds groundwater outflow from lower saturated zone.
//
{
    double q, t1, t2, t3;

    // --- water table must be above project->Hstar for flow to occur
    if ( lowerDepth <= project->Hstar ) return 0.0;

    // --- compute groundwater component of flow
    if ( project->GW->b1 == 0.0 ) t1 = project->GW->a1;
	else t1 = project->GW->a1 * pow((lowerDepth - project->Hstar)*UCF(project, LENGTH), project->GW->b1);

    // --- compute surface water component of flow
    if ( project->GW->b2 == 0.0 ) t2 = project->GW->a2;
    else if (project->Hsw > project->Hstar)
    {
		t2 = project->GW->a2 * pow((project->Hsw - project->Hstar)*UCF(project, LENGTH), project->GW->b2);
    }
    else t2 = 0.0;

    // --- compute groundwater/surface water interaction term
	t3 = project->GW->a3 * lowerDepth * project->Hsw * UCF(project, LENGTH) * UCF(project, LENGTH);

    // --- compute total groundwater flow
	q = (t1 - t2 + t3) / UCF(project, GWFLOW);
    if ( q < 0.0 && project->GW->a3 != 0.0 ) q = 0.0;
    return q;
}

//=============================================================================

int  getVariableIndex(Project* project, char* s)
//
//  Input:   s = name of a groundwater variable
//  Output:  returns index of groundwater variable
//  Purpose: finds position of project->GW variable in list of project->GW variable names.
//
{
    int k;

    k = findmatch(s, GWVarWords);
    if ( k >= 0 ) return k;
    return -1;
}

//=============================================================================

double getVariableValue(Project* project, int varIndex)
//
//  Input:   varIndex = index of a project->GW variable
//  Output:  returns current value of project->GW variable
//  Purpose: finds current value of a project->GW variable.
//
{
    switch (varIndex)
    {
    case gwvHGW:  return project->Hgw * UCF(project,LENGTH);
	case gwvHSW:  return project->Hsw * UCF(project, LENGTH);
	case gwvHCB:  return project->Hstar * UCF(project, LENGTH);
	case gwvHGS:  return project->TotalDepth * UCF(project, LENGTH);
	case gwvKS:   return project->A.conductivity * UCF(project, RAINFALL);
	case gwvK:    return project->HydCon * UCF(project, RAINFALL);                               //(5.1.010)
    case gwvTHETA:return project->Theta;                                                //(5.1.008)
    case gwvPHI:  return project->A.porosity;                                           //(5.1.008)
	case gwvFI:   return project->Infil * UCF(project, RAINFALL);                                //(5.1.008)
	case gwvFU:   return project->UpperPerc * UCF(project, RAINFALL);                            //(5.1.008)
	case gwvA:    return project->Area * UCF(project, LANDAREA);                                 //(5.1.008)
    default:      return 0.0;
    }
}
