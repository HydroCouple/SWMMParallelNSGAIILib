//-----------------------------------------------------------------------------
//   treatmnt.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             03/19/15   (Build 5.1.008)
//   Author:   L. Rossman
//
//   Pollutant treatment functions.
//
//   Build 5.1.008:
//   - A bug in evaluating recursive calls to treatment functions was fixed. 
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const int PVMAX = 5;            // number of process variables
enum   ProcessVarType {pvHRT,          // hydraulic residence time
                       pvDT,           // time step duration
                       pvFLOW,         // flow rate
                       pvDEPTH,        // water height above invert
                       pvAREA};        // storage surface area

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
//static int     project->ErrCode;                // treatment error code
//static int     project->J;                      // index of node being analyzed
//static double  project->Dt;                     // curent time step (sec)
//static double  project->Q;                      // node inflow (cfs)
//static double  project->V;                      // node volume (ft3)
//static double* project->R;                      // array of pollut. removals
//static double* project->Cin;                    // node inflow concentrations
//static TTreatment* Treatment; // defined locally in treatmnt_treat()         //(5.1.008)

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  treatmnt_open           (called from routing_open)
//  treatment_close         (called from routing_close)
//  treatmnt_readExpression (called from parseLine in input.c)
//  treatmnt_delete         (called from deleteObjects in project.c)
//  treatmnt_setInflow      (called from qualrout_execute)
//  treatmnt_treat          (called from findNodeQual in qualrout.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    createTreatment(Project* project, int node);
static double getRemoval(Project* project, int pollut);
static int    getVariableIndex(Project* project, char* s);
static double getVariableValue(Project* project, int varCode);


//=============================================================================

int  treatmnt_open(Project* project)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: allocates memory for computing pollutant removals by treatment.
//
{
    project->R = NULL;
    project->Cin = NULL;
    if ( project->Nobjects[POLLUT] > 0 )
    {
        project->R = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
        project->Cin = (double *) calloc(project->Nobjects[POLLUT], sizeof(double));
        if ( project->R == NULL || project->Cin == NULL)
        {
            report_writeErrorMsg(project,ERR_MEMORY, "");
            return FALSE;
        }
    }
    return TRUE;
}

//=============================================================================

void treatmnt_close(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: frees memory used for computing pollutant removals by treatment.
//
{
    FREE(project->R);
    FREE(project->Cin);
}

//=============================================================================

int  treatmnt_readExpression(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads a treatment expression from a tokenized line of input.
//
{
    char  s[MAXLINE+1];
    char* expr;
    int   i, j, k, p;
    MathExpr* equation;                // ptr. to a math. expression

    // --- retrieve node & pollutant
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(project,NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    p = project_findObject(project,POLLUT, tok[1]);
    if ( p < 0 ) return error_setInpError(ERR_NAME, tok[1]);

    // --- concatenate remaining tokens into a single string
    strcpy(s, tok[2]);
    for ( i=3; i<ntoks; i++)
    {
        strcat(s, " ");
        strcat(s, tok[i]);
    }

    // --- check treatment type
    if      ( UCHAR(s[0]) == 'R' ) k = 0;
    else if ( UCHAR(s[0]) == 'C' ) k = 1;
    else return error_setInpError(ERR_KEYWORD, tok[2]);

    // --- start treatment expression after equals sign
    expr = strchr(s, '=');
    if ( expr == NULL ) return error_setInpError(ERR_KEYWORD, "");
    else expr++;

    // --- create treatment objects at node j if they don't already exist
    if ( project->Node[j].treatment == NULL )
    {
        if ( !createTreatment(project,j) ) return error_setInpError(ERR_MEMORY, "");
    }

    // --- create a parsed expression tree from the string expr
    //     (getVariableIndex is the function that converts a treatment
    //      variable's name into an index number) 
    equation = mathexpr_create(project,expr, getVariableIndex);
    if ( equation == NULL )
        return error_setInpError(ERR_TREATMENT_EXPR, "");

    // --- save the treatment parameters in the node's treatment object
    project->Node[j].treatment[p].treatType = k;
    project->Node[j].treatment[p].equation = equation;
    return 0;
}

//=============================================================================

void treatmnt_delete(Project* project, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: deletes the treatment objects for each pollutant at a node.
//
{
    int p;
    if ( project->Node[j].treatment )
    {
        for (p=0; p<project->Nobjects[POLLUT]; p++)
            mathexpr_delete(project->Node[j].treatment[p].equation);
        free(project->Node[j].treatment);
    }
    project->Node[j].treatment = NULL;
}

//=============================================================================

void  treatmnt_setInflow(Project* project, double qIn, double wIn[])
//
//  Input:   j = node index
//           qIn = flow inflow rate (cfs)
//           wIn = pollutant mass inflow rate (mass/sec)
//  Output:  none
//  Purpose: computes and saves array of inflow concentrations to a node.
//
{
    int    p;
    if ( qIn > 0.0 )
        for (p = 0; p < project->Nobjects[POLLUT]; p++) project->Cin[p] = wIn[p]/qIn;
    else
        for (p = 0; p < project->Nobjects[POLLUT]; p++) project->Cin[p] = 0.0;
}

//=============================================================================

void  treatmnt_treat(Project* project, int j, double q, double v, double tStep)
//
//  Input:   j     = node index
//           q     = inflow to node (cfs)
//           v     = volume of node (ft3)
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: updates pollutant concentrations at a node after treatment.
//
{
    int    p;                          // pollutant index
    double cOut;                       // concentration after treatment
    double massLost;                   // mass lost by treatment per time step
    TTreatment* treatment;             // pointer to treatment object          //(5.1.008)

    // --- set locally shared variables for node j
    if ( project->Node[j].treatment == NULL ) return;
    project->ErrCode = 0;
    project->J  = j;                            // current node
    project->Dt = tStep;                        // current time step
    project->Q  = q;                            // current inflow rate
    project->V  = v;                            // current node volume

    // --- initialze each removal to indicate no value 
    for ( p = 0; p < project->Nobjects[POLLUT]; p++) project->R[p] = -1.0;

    // --- determine removal of each pollutant
    for ( p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        // --- removal is zero if there is no treatment equation
        treatment = &project->Node[j].treatment[p];                                     //(5.1.008)
        if ( treatment->equation == NULL ) project->R[p] = 0.0;                         //(5.1.008)

        // --- no removal for removal-type expression when there is no inflow 
	    else if ( treatment->treatType == REMOVAL && q <= ZERO ) project->R[p] = 0.0;   //(5.1.008)

        // --- otherwise evaluate the treatment expression to find project->R[p]
        else getRemoval(project,p);
    }

    // --- check for error condition
    if ( project->ErrCode == ERR_CYCLIC_TREATMENT )
    {
         report_writeErrorMsg(project,ERR_CYCLIC_TREATMENT, project->Node[project->J].ID);
    }

    // --- update nodal concentrations and mass balances
    else for ( p = 0; p < project->Nobjects[POLLUT]; p++ )
    {
        if ( project->R[p] == 0.0 ) continue;
        treatment = &project->Node[j].treatment[p];                                     //(5.1.008)

        // --- removal-type treatment equations get applied to inflow stream

        if ( treatment->treatType == REMOVAL )                                 //(5.1.008)
        {
            // --- if no pollutant in inflow then cOut is current nodal concen.
            if ( project->Cin[p] == 0.0 ) cOut = project->Node[j].newQual[p];

            // ---  otherwise apply removal to influent concen.
            else cOut = (1.0 - project->R[p]) * project->Cin[p];

            // --- cOut can't be greater than mixture concen. at node
            //     (i.e., in case node is a storage unit) 
            cOut = MIN(cOut, project->Node[j].newQual[p]);
        }

        // --- concentration-type equations get applied to nodal concentration
        else
        {
            cOut = (1.0 - project->R[p]) * project->Node[j].newQual[p];
        }

        // --- mass lost must account for any initial mass in storage 
        massLost = (project->Cin[p]*q*tStep + project->Node[j].oldQual[p]*project->Node[j].oldVolume - 
                   cOut*(q*tStep + project->Node[j].oldVolume)) / tStep; 
        massLost = MAX(0.0, massLost); 

        // --- add mass loss to mass balance totals and revise nodal concentration
        massbal_addReactedMass(project,p, massLost);
        project->Node[j].newQual[p] = cOut;
    }
}

//=============================================================================

int  createTreatment(Project* project, int j)
//
//  Input:   j = node index
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: creates a treatment object for each pollutant at a node.
//
{
    int p;
    project->Node[j].treatment = (TTreatment *) calloc(project->Nobjects[POLLUT],
                                              sizeof(TTreatment));
    if ( project->Node[j].treatment == NULL )
    {
        return FALSE;
    }
    for (p = 0; p < project->Nobjects[POLLUT]; p++)
    {
        project->Node[j].treatment[p].equation = NULL;
    }
    return TRUE;
}

//=============================================================================

int  getVariableIndex(Project* project, char* s)
//
//  Input:   s = name of a process variable or pollutant
//  Output:  returns index of process variable or pollutant
//  Purpose: finds position of process variable/pollutant in list of names.
//
{
    // --- check for a process variable first
    int k;
    int m = PVMAX;                     // PVMAX is number of process variables

    k = findmatch(s, ProcessVarWords);
    if ( k >= 0 ) return k;

    // --- then check for a pollutant concentration
    k = project_findObject(project,POLLUT, s);
    if ( k >= 0 ) return (k + m);

    // --- finally check for a pollutant removal
    if ( UCHAR(s[0]) == 'R' && s[1] == '_')
    {
        k = project_findObject(project,POLLUT, s+2);
        if ( k >= 0 ) return (project->Nobjects[POLLUT] + k + m);
    }
    return -1;
}

//=============================================================================

double getVariableValue(Project* project, int varCode)
//
//  Input:   varCode = code number of process variable or pollutant
//  Output:  returns current value of variable
//  Purpose: finds current value of a process variable or pollutant concen.,
//           making reference to the node being evaluated which is stored in
//           shared variable project->J.
//
{
    int    p;
    double a1, a2, y;
    TTreatment* treatment;                                                     //(5.1.008)

    // --- variable is a process variable
    if ( varCode < PVMAX )
    {
        switch ( varCode )
        {
          case pvHRT:                                 // HRT in hours
            if ( project->Node[project->J].type == STORAGE )
            {
                return project->Storage[project->Node[project->J].subIndex].hrt / 3600.0;
            }
            else return 0.0;

          case pvDT:
            return project->Dt;                                // time step in seconds

          case pvFLOW:
            return project->Q * UCF(project, FLOW);                     // flow in user's units

          case pvDEPTH:
            y = (project->Node[project->J].oldDepth + project->Node[project->J].newDepth) / 2.0;
			return y * UCF(project, LENGTH);                   // depth in ft or m

          case pvAREA:
			  a1 = node_getSurfArea(project, project->J, project->Node[project->J].oldDepth);
			  a2 = node_getSurfArea(project, project->J, project->Node[project->J].newDepth);
			return (a1 + a2) / 2.0 * UCF(project, LENGTH) * UCF(project, LENGTH);
            
          default: return 0.0;
        }
    }

    // --- variable is a pollutant concentration
    else if ( varCode < PVMAX + project->Nobjects[POLLUT] )
    {
        p = varCode - PVMAX;
        treatment = &project->Node[project->J].treatment[p];                                     //(5.1.008)
        if ( treatment->treatType == REMOVAL ) return project->Cin[p];                  //(5.1.008)
        return project->Node[project->J].newQual[p];
    }

    // --- variable is a pollutant removal
    else
    {
        p = varCode - PVMAX - project->Nobjects[POLLUT];
        if ( p >= project->Nobjects[POLLUT] ) return 0.0;
        return getRemoval(project,p);
    }
}

//=============================================================================

double  getRemoval(Project* project, int p)
//
//  Input:   p = pollutant index
//  Output:  returns fractional removal of pollutant
//  Purpose: computes removal of a specific pollutant
//
{
    double c0 = project->Node[project->J].newQual[p];    // initial node concentration
    double r;                          // removal value
    TTreatment* treatment;                                                     //(5.1.008)

    // --- case where removal already being computed for another pollutant
    if ( project->R[p] > 1.0 || project->ErrCode )
    {
        project->ErrCode = 1;
        return 0.0;
    }

    // --- case where removal already computed
    if ( project->R[p] >= 0.0 && project->R[p] <= 1.0 ) return project->R[p];

    // --- set project->R[p] to value > 1 to show that value is being sought
    //     (prevents infinite recursive calls in case two removals
    //     depend on each other)
    project->R[p] = 10.0;

    // --- case where current concen. is zero
    if ( c0 == 0.0 )
    {
        project->R[p] = 0.0;
        return 0.0;
    }

    // --- apply treatment eqn.
    treatment = &project->Node[project->J].treatment[p];                                         //(5.1.008)
    r = mathexpr_eval(project, treatment->equation, getVariableValue);                  //(5.1.008)
    r = MAX(0.0, r);

    // --- case where treatment eqn. is for removal
    if ( treatment->treatType == REMOVAL )                                     //(5.1.008)
    {
        r = MIN(1.0, r);
        project->R[p] = r;
    }

    // --- case where treatment eqn. is for effluent concen.
    else
    {
        r = MIN(c0, r);
        project->R[p] = 1.0 - r/c0;
    }
    return project->R[p];
}

//=============================================================================
