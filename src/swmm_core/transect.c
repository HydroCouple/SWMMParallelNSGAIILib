//-----------------------------------------------------------------------------
//   transect.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//   Author:   L. Rossman
//
//   Geometry processing for irregular cross-section transects.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
// #define MAXSTATION 1500                // max. number of stations in a transect

////-----------------------------------------------------------------------------
////  Shared variables
////-----------------------------------------------------------------------------
//static int    project->Ntransects;              // total number of transects
//static int    project->Nstations;               // number of stations in current transect
//static double  project->Station[MAXSTATION+1];  // x-coordinate of each station
//static double  project->Elev[MAXSTATION+1];     // elevation of each station
//static double  project->Nleft;                  // Manning's n for left overbank
//static double  project->Nright;                 // Manning's n for right overbank
//static double  project->Nchannel;               // Manning's n for main channel
//static double  project->Xleftbank;              // station where left overbank ends
//static double  project->Xrightbank;             // station where right overbank begins
//static double  project->Xfactor;                // multiplier for station spacing
//static double  project->Yfactor;                // factor added to station elevations
//static double  project->Lfactor;                // main channel/flood plain length

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)   
//-----------------------------------------------------------------------------
//  transect_create      (called by createObjects in project.c)
//  transect_delete      (called by deleteObjects in project.c)
//  transect_readParams  (called by parseLine in input.c)
//  transect_validate    (called by input_readData)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    setParams(Project* project, int transect, char* id, double x[]);
static int    setManning(Project* project, double n[]);
static int    addStation(Project* project, double x, double y);
static double getFlow(Project* project, int k, double a, double wp, int findFlow);
static void   getGeometry(Project* project, int i, int j, double y);
static void   getSliceGeom(Project* project, int k, double y, double yu, double yd, double *w,
              double *a, double *wp);
static void   setMaxSectionFactor(Project* project, int transect);

//=============================================================================

int transect_create(Project* project, int n)
//
//  Input:   n = number of transect objects to create
//  Output:  returns an error code
//  Purpose: creates an array of cross-section transects.
//
{
    project->Ntransects = n;
    if ( n == 0 ) return 0;
    project->Transect = (TTransect *) calloc(project->Ntransects, sizeof(TTransect));
    if ( project->Transect == NULL ) return ERR_MEMORY;
    project->Nchannel = 0.0;
    project->Nleft = 0.0;
    project->Nright = 0.0;
    project->Nstations = 0;
    return 0;
}

//=============================================================================

void transect_delete(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: deletes memory allocated for all transects.
//
{
    if ( project->Ntransects == 0 ) return;
    FREE(project->Transect);
    project->Ntransects = 0;
}

//=============================================================================

int transect_readParams(Project* project, int* count, char* tok[], int ntoks)
//
//  Input:   count = transect index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  updated value of count,
//           returns an error code
//  Purpose: read parameters of a transect from a tokenized line of input data.
//
//  Format of transect data follows that used for HEC-2 program:
//    NC  nLeft  nRight  nChannel
//    X1  name  nSta  xLeftBank  xRightBank  0  0  0  xFactor  yFactor
//    GR  Elevation  Station  ... 
//
{
    int    i, k;
    int    index = *count;             // transect index
    int    errcode;                    // error code
    double x[10];                      // parameter values
    char*  id;                         // transect ID name

    // --- match first token to a transect keyword
    k = findmatch(tok[0], TransectKeyWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);

    // --- read parameters associated with keyword
    switch ( k )
    {
      // --- NC line: Manning n values
      case 0:

        // --- finish processing the previous transect
		  transect_validate(project, index - 1);

        // --- read Manning's n values
        if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i <= 3; i++)
        {
            if ( ! getDouble(tok[i], &x[i]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
		return setManning(project, x);

      // --- X1 line: identifies start of next transect
      case 1:

        // --- check that transect was already added to project
        //     (by input_countObjects)
        if ( ntoks < 10 ) return error_setInpError(ERR_ITEMS, "");
        id = project_findID(project,TRANSECT, tok[1]);
        if ( id == NULL ) return error_setInpError(ERR_NAME, tok[1]);

        // --- read in rest of numerical values on data line
        for ( i = 2; i < 10; i++ )
        {
            if ( ! getDouble(tok[i], &x[i]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }

        // --- update total transect count
        *count = index + 1;

        // --- transfer parameter values to transect's properties
		return setParams(project, index, id, x);

      // --- GR line: station elevation & location data
      case 2:

        // --- check that line contains pairs of data values
        if ( (ntoks - 1) % 2 > 0 ) return error_setInpError(ERR_ITEMS, "");

        // --- parse each pair of Elevation-Station values
        i = 1;
        while ( i < ntoks )
        {
            if ( ! getDouble(tok[i], &x[1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
            if ( ! getDouble(tok[i+1], &x[2]) )
                return error_setInpError(ERR_NUMBER, tok[i+1]);
			errcode = addStation(project, x[1], x[2]);
            if ( errcode ) return errcode;
            i += 2;
        }
        return 0;
    }
    return 0;
}

//=============================================================================

void  transect_validate(Project* project, int j)
//
//  Input:   j = transect index
//  Output:  none
//  Purpose: validates transect data and creates its geometry tables.
//
{
    int    i, nLast;
    double dy, y, ymin, ymax;
    double oldNchannel = project->Nchannel;

    // --- check for valid transect data
    if ( j < 0 || j >= project->Ntransects ) return;
    if ( project->Nstations < 2 ) 
    {
        report_writeErrorMsg(project,ERR_TRANSECT_TOO_FEW, project->Transect[j].ID);
        return;
    }
    if ( project->Nstations >= MAXSTATION )
    {
        report_writeErrorMsg(project,ERR_TRANSECT_TOO_MANY, project->Transect[j].ID);
        return;
    }
    if ( project->Nchannel <= 0.0 )
    {
        report_writeErrorMsg(project,ERR_TRANSECT_MANNING, project->Transect[j].ID);
        return;
    }
    if ( project->Xleftbank > project->Xrightbank )
    {
        report_writeErrorMsg(project,ERR_TRANSECT_OVERBANK, project->Transect[j].ID);
        return;
    }

    // --- adjust main channel's Mannings n to make its equivalent
    //     length equal to that of entire flood plain
    project->Nchannel = project->Nchannel * sqrt(project->Lfactor);
    project->Transect[j].lengthFactor = project->Lfactor;

    // --- find max. depth across transect
    ymax = project->Elev[1];
    ymin = project->Elev[1];
    for (i = 2; i <= project->Nstations; i++)
    {
        ymax = MAX(project->Elev[i], ymax);
        ymin = MIN(project->Elev[i], ymin);
    }
    if ( ymin >= ymax )
    {
        report_writeErrorMsg(project,ERR_TRANSECT_NO_DEPTH, project->Transect[j].ID);
        return;
    }
    project->Transect[j].yFull = ymax - ymin;

    // --- add vertical sides to transect to reach full ht. on both ends
    project->Station[0] = project->Station[1];
    project->Elev[0] = ymax;
    project->Nstations++;
    project->Station[project->Nstations] = project->Station[project->Nstations-1];
    project->Elev[project->Nstations] = project->Elev[0];

    // --- determine size & depth increment for geometry tables
    project->Transect[j].nTbl = N_TRANSECT_TBL;
    dy = (ymax - ymin) / (double)(project->Transect[j].nTbl - 1);

    // --- set 1st table entries to zero
    project->Transect[j].areaTbl[0] = 0.0;
    project->Transect[j].hradTbl[0] = 0.0;
    project->Transect[j].widthTbl[0] = 0.0;

    // --- compute geometry for each depth increment
    y = ymin;
    project->Transect[j].wMax = 0.0;
    for (i = 1; i < project->Transect[j].nTbl; i++)
    {
        y += dy;
        project->Transect[j].areaTbl[i] = 0.0;
        project->Transect[j].hradTbl[i] = 0.0;
        project->Transect[j].widthTbl[i] = 0.0;
		getGeometry(project, i, j, y);
    }

    // --- determine max. section factor 
	setMaxSectionFactor(project, j);

    // --- normalize geometry table entries
    //     (full cross-section values are last table entries)
    nLast = project->Transect[j].nTbl - 1;
    project->Transect[j].aFull = project->Transect[j].areaTbl[nLast];
    project->Transect[j].rFull = project->Transect[j].hradTbl[nLast];
    project->Transect[j].wMax = project->Transect[j].widthTbl[nLast];

    for (i = 1; i <= nLast; i++)
    {
        project->Transect[j].areaTbl[i] /= project->Transect[j].aFull;
        project->Transect[j].hradTbl[i] /= project->Transect[j].rFull;
        project->Transect[j].widthTbl[i] /= project->Transect[j].wMax;
    }

    // --- set width at 0 height equal to width at 4% of max. height
    project->Transect[j].widthTbl[0] = project->Transect[j].widthTbl[1];

    // --- save unadjusted main channel roughness 
    project->Transect[j].roughness = oldNchannel;
}

//=============================================================================

int  setManning(Project* project, double n[])
//
//  Input:   n[] = array of Manning's n values
//  Output:  returns an error code
//  Purpose: sets Manning's n for overbanks and main channel of a transect.
//
{
    int i;
    for (i=1; i<=3; i++)
    {
        if ( n[i] < 0.0 ) return ERR_NUMBER;
    }
    if ( n[1] > 0.0 ) project->Nleft = n[1];
    if ( n[2] > 0.0 ) project->Nright = n[2];
    if ( n[3] > 0.0 ) project->Nchannel = n[3];
    if ( project->Nleft == 0.0  ) project->Nleft = project->Nchannel;
    if ( project->Nright == 0.0 ) project->Nright = project->Nchannel;
    return 0;
}

//=============================================================================

int  setParams(Project* project, int j, char* id, double x[])
//
//  Input:   j = transect index
//           id = transect ID name
//           x[] = array of parameter values
//  Output:  returns an error code
//  Purpose: assigns parameter values to current transect being processed.
//
{
    if ( j < 0 || j >= project->Ntransects ) return ERR_NUMBER;
    project->Transect[j].ID = id;                         // ID name
    project->Xleftbank = x[3] / UCF(project,LENGTH);              // left overbank location
    project->Xrightbank = x[4] / UCF(project,LENGTH);             // right overbank location
    project->Lfactor = x[7];                              // channel/bank length
    if ( project->Lfactor == 0.0 ) project->Lfactor = 1.0;
    project->Xfactor = x[8];                              // station location multiplier
    if ( project->Xfactor == 0.0 ) project->Xfactor = 1.0;
    project->Xleftbank *= project->Xfactor;                        // adjusted left bank
    project->Xrightbank *= project->Xfactor;                       // adjusted right bank
    project->Yfactor = x[9] / UCF(project,LENGTH);                // elevation offset
    project->Nstations = 0;
    return 0;
}

//=============================================================================

int  addStation(Project* project, double y, double x)
//
//  Input:   y = station elevation value
//           x = station distance value
//  Output:  returns an error code
//  Purpose: adds a new station to the transect currently being processed.
//
{
    // --- check for valid number of stations
    if ( project->Nstations < 0 ) return ERR_TRANSECT_UNKNOWN;
    project->Nstations++;
    if ( project->Nstations >= MAXSTATION ) return 0;

    // --- add station distance, modified by distance multiplier
	project->Station[project->Nstations] = x * project->Xfactor / UCF(project, LENGTH);

    // --- add station elevation, modified by offset elevation
	project->Elev[project->Nstations] = (y + project->Yfactor) / UCF(project, LENGTH);

    // --- check if station distances are non-increasing
    if ( project->Nstations > 1
        && project->Station[project->Nstations] < project->Station[project->Nstations-1] )
        return ERR_TRANSECT_SEQUENCE;
    return 0;    
}

//=============================================================================

void  getGeometry(Project* project, int i, int j, double y)
//
//  Input:   i = index of current entry in geometry tables
//           j = transect index
//           y = depth of current entry in geometry tables
//  Output:  none
//  Purpose: computes entries in a transect's geometry tables at a given depth. 
//
{
    int    k;                // station index
    double ylo,              // lower elev. of transect slice
           yhi,              // higher elev. of transect slice
           w,                // top width of transect slice
           wp,               // wetted perimeter of transect slice
           wpSum,            // total wetted perimeter across transect
           a,                // area of transect slice
           aSum,             // total area across transect
           q,                // flow across transect slices with same roughness
           qSum;             // total flow across transect
    int   findFlow;          // true if flow thru area slice needs updating

    // --- initialize
    wpSum = 0.0;
    aSum = 0.0;
    qSum = 0.0;

    // --- examine each horizontal station from left to right
    for (k = 1; k <= project->Nstations; k++)
    {
        // --- determine low & high elevations for transect sub-section
        if ( project->Elev[k-1] >= project->Elev[k] )
        {
            yhi = project->Elev[k-1];
            ylo = project->Elev[k];
        }
        else
        {
            yhi = project->Elev[k];
            ylo = project->Elev[k-1];
        }

        // --- skip station if its totally dry
        if ( ylo >= y ) continue;

        // --- get top width, area & wetted perimeter values for transect
        //     slice between station k and k-1
		getSliceGeom(project, k, y, ylo, yhi, &w, &a, &wp);

        // --- update total transect values
        wpSum += wp;
        aSum += a;
        project->Transect[j].areaTbl[i] += a;
        project->Transect[j].widthTbl[i] += w;

        // --- must update flow if station elevation is above water level
        if ( project->Elev[k] >= y ) findFlow = TRUE;
        else findFlow = FALSE;

        // --- update flow across transect if called for
		q = getFlow(project, k, aSum, wpSum, findFlow);
        if ( q > 0.0 )
        {
            qSum += q;
            aSum = 0.0;
            wpSum = 0.0;
        }

    }   // next station k 

    // --- find hyd. radius table entry solving Manning eq. with
    //     total flow, total area, and main channel n
    aSum = project->Transect[j].areaTbl[i];
    if ( aSum == 0.0 ) project->Transect[j].hradTbl[i] = project->Transect[j].hradTbl[i-1];
    else project->Transect[j].hradTbl[i] = pow(qSum * project->Nchannel / 1.49 / aSum, 1.5);
}

//=============================================================================

void getSliceGeom(Project* project, int k, double y, double ylo, double yhi, double *w,
                  double *a, double *wp)
//
//  Input:   k = station index
//           y = water elevation
//           ylo = transect elevation on low side of slice
//           yhi = transect elevation on high side of slice
//  Output   w = width of transect slice
//           a = area of transect slice
//           wp = wetted perimeter of transect slice
//  Purpose: finds area, width & wetted perim. for slice of transect that
//           is covered by given water depth.
//
//      yhi  |           
//           |
//        y  |**********
//           |********** --> slice of transect being analyzed
//      ylo  |**********|
//           |**********|
//           |**********|
//         Station    Station
//           k-1        k
//
{
    double width, ratio;

    // --- compute width & wetted perimeter of transect slice
    width = fabs(project->Station[k] - project->Station[k-1]);
    (*w) = width;
    (*wp) = sqrt(width * width + (yhi - ylo) * (yhi - ylo));
    (*a)  = 0.0;

    // --- find area for completely submerged slice
    if ( y > yhi )
    {
        (*a) = width * ( (y - yhi) + (y - ylo) ) / 2.0;
    }

    // --- otherwise find area and adjust width & wetted perim. for
    //     partly submerged slice
    else if ( yhi > ylo )
    {
         ratio = (y - ylo) / (yhi - ylo);
         (*a) = width * (yhi - ylo) / 2.0 * ratio * ratio;
         (*w) *= ratio;
         (*wp) *= ratio;
     }
}

//=============================================================================

double getFlow(Project* project, int k, double a, double wp, int findFlow)
//
//  Input:   k = index of station at end of transect sub-section
//           a = flow area of sub-section
//           wp = wetted perimeter of flow area of sub-section
//           findFlow = TRUE if flow needs updating 
//  Output:  returns normal flow (per unit of slope)
//  Purpose: finds flow through a sub-section of a transect.
//
{
    double n;                          // Manning's n

    if ( findFlow == FALSE)
    {
        // --- flow needs updating if we are at last station
        if ( k == project->Nstations - 1 ) findFlow = TRUE;

        // --- flow needs updating if we are at end of left overbank and
        //     there is a change in Manning's n and section not vertical
        else if ( project->Station[k] == project->Xleftbank )
        {
            if ( project->Nleft != project->Nchannel &&
                project->Station[k] != project->Station[k-1] ) findFlow = TRUE;
        }

        // --- flow needs updating if we are at start of right overbank and
        //     there is a change in Manning's n and section not vertical
        else if ( project->Station[k] == project->Xrightbank )
        {
            if ( project->Nright != project->Nchannel &&
                project->Station[k] != project->Station[k+1] ) findFlow = TRUE;
        }
    }

    // --- if flow needs updating
    if ( findFlow )
    {
        // --- find value of Manning's n to use
        n = project->Nchannel;
        if ( project->Station[k-1] < project->Xleftbank ) n = project->Nleft;
        if ( project->Station[k] > project->Xrightbank )  n = project->Nright;

        // --- compute flow through flow area
        return PHI / n * a * pow(a/wp, 2./3.);
    }
    return 0.0;
}

//=============================================================================

void setMaxSectionFactor(Project* project, int j)
//
//  Input:   j = transect index
//  Output:  none
//  Purpose: determines the maximum section factor for a transect and the
//           area where this maxumum occurs.
//
{
    int    i;
    double sf;

    project->Transect[j].aMax = 0.0;
    project->Transect[j].sMax = 0.0;
    for (i=1; i<project->Transect[j].nTbl; i++)
    {
        sf = project->Transect[j].areaTbl[i] * pow(project->Transect[j].hradTbl[i], 2./3.);
        if ( sf > project->Transect[j].sMax )
        {
            project->Transect[j].sMax = sf;
            project->Transect[j].aMax = project->Transect[j].areaTbl[i];
        }
    }
}

//=============================================================================
