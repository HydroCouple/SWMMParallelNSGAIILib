//-----------------------------------------------------------------------------
//   rdii.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             04/04/14   (Build 5.1.003)
//             04/14/14   (Build 5.1.004)
//             09/15/14   (Build 5.1.007)
//   Author:   L. Rossman (EPA)
//             R. Dickinson (CDM)
//
//   RDII processing functions.
//
//   Note: RDII means rainfall dependent infiltration/inflow,
//         UH means unit hydrograph.
//
//   Build 5.1.007:
//   - Ignore RDII option implemented.
//   - Rainfall climate adjustment implemented.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include "headers.h"

//-----------------------------------------------------------------------------
// Definition of 4-byte integer, 4-byte real and 8-byte real types
//-----------------------------------------------------------------------------
//#define INT4  int
//#define REAL4 float
//#define REAL8 double
#define FILE_STAMP "SWMM5-RDII"

//-----------------------------------------------------------------------------
// Constants
//-----------------------------------------------------------------------------
const double ZERO_RDII = 0.0001;       // Minimum non-zero RDII inflow (cfs)
const char   FileStamp[] = FILE_STAMP;

////-----------------------------------------------------------------------------
//// Data Structures
////----------------------------------------------------------------------------
enum FileTypes {BINARY, TEXT};         // File mode types
//
//typedef struct                         // Data for a single unit hydrograph
//{                                      // -------------------------------------
//   double*   pastRain;                 // array of past rainfall values
//   char*     pastMonth;                // month in which past rainfall occurred
//   int       period;                   // current UH time period
//   int       hasPastRain;              // true if > 0 past periods with rain
//   int       maxPeriods;               // max. past rainfall periods
//   long      drySeconds;               // time since last nonzero rainfall
//   double    iaUsed;                   // initial abstraction used (in or mm)
//}  TUHData;
//
//typedef struct                         // Data for a unit hydrograph group
//{                                      //---------------------------------
//   int       isUsed;                   // true if UH group used by any nodes
//   int       rainInterval;             // time interval for RDII processing (sec)
//   double    area;                     // sewered area covered by UH's gage (ft2)
//   double    rdii;                     // rdii flow (in rainfall units)
//   DateTime  gageDate;                 // calendar date of rain gage period
//   DateTime  lastDate;                 // date of last rdii computed
//   TUHData   uh[3];                    // data for each unit hydrograph
//}  TUHGroup;

////-----------------------------------------------------------------------------
//// Shared Variables
////-----------------------------------------------------------------------------
//static TUHGroup*  project->UHGroup;             // processing data for each UH group
//static int        project->RdiiStep;            // RDII time step (sec)
//static int        project->NumRdiiNodes;        // number of nodes w/ RDII data
//static int*       project->RdiiNodeIndex;       // indexes of nodes w/ RDII data
//static REAL4*     project->RdiiNodeFlow;        // inflows for nodes with RDII          //(5.1.003)
//static int        project->RdiiFlowUnits;       // RDII flow units code
//static DateTime   project->RdiiStartDate;       // start date of RDII inflow period
//static DateTime   project->RdiiEndDate;         // end date of RDII inflow period
//static double     project->TotalRainVol;        // total rainfall volume (ft3)
//static double     project->TotalRdiiVol;        // total RDII volume (ft3)
//static int        project->RdiiFileType;        // type (binary/text) of RDII file

//-----------------------------------------------------------------------------
// Imported Variables
//-----------------------------------------------------------------------------
extern double     Qcf[];               // flow units conversion factors
                                       // (see swmm5.c)
//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  rdii_readRdiiInflow     (called from parseLine in input.c)
//  rdii_deleteRdiiInflow   (called from deleteObjects in project.c)
//  rdii_initUnitHyd        (called from createObjects in project.c)
//  rdii_readUnitHydParams  (called from parseLine in input.c)
//  rdii_openRdii           (called from rain_open)
//  rdii_closeRdii          (called from rain_close)
//  rdii_getNumRdiiFlows    (called from addRdiiInflows in routing.c)
//  rdii_getRdiiFlow        (called from addRdiiInflows in routing.c)

//-----------------------------------------------------------------------------
// Function Declarations
//-----------------------------------------------------------------------------
// --- functions used to create a RDII file
static int    readOldUHFormat(Project* project, int j, int m, char* tok[], int ntoks);
static void   setUnitHydParams(Project* project, int j, int i, int m, double x[]);
static void   createRdiiFile(Project* project);
static int    getNumRdiiNodes(Project* project);
static void   validateRdii(Project* project);

static void   openRdiiProcessor(Project* project);
static int    allocRdiiMemory(Project* project);
static int    getRainInterval(Project* project, int i);
static int    getMaxPeriods(Project* project, int i, int k);
static void   initGageData(Project* project);
static void   initUnitHydData(Project* project);
static int    openNewRdiiFile(Project* project);
static void   getRainfall(Project* project, DateTime currentDate);

static double applyIA(Project* project, int j, int k, DateTime aDate, double dt,
              double rainDepth);
static void   updateDryPeriod(Project* project, int j, int k, double rain, int gageInterval);
static void   getUnitHydRdii(Project* project, DateTime currentDate);
static double getUnitHydConvol(Project* project, int j, int k, int gageInterval);
static double getUnitHydOrd(Project* project, int j, int m, int k, double t);

static int    getNodeRdii(Project* project);
static void   saveRdiiFlows(Project* project, DateTime currentDate);
static void   closeRdiiProcessor(Project* project);
static void   freeRdiiMemory(Project* project);

// --- functions used to read an existing RDII file
static int   readRdiiFileHeader(Project* project);
static void  readRdiiFlows(Project* project);

static void  openRdiiTextFile(Project* project);
static int   readRdiiTextFileHeader(Project* project);
static void  readRdiiTextFlows(Project* project);

//=============================================================================
//                   Management of RDII-Related Data
//=============================================================================

int rdii_readRdiiInflow(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads properties of an RDII inflow from a line of input.
//
{
    int    j, k;
    double a;
    TRdiiInflow* inflow;

    // --- check for proper number of items
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");

    // --- check that node receiving RDII exists
    j = project_findObject(project, NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);

    // --- check that RDII unit hydrograph exists
    k = project_findObject(project,UNITHYD, tok[1]);
    if ( k < 0 ) return error_setInpError(ERR_NAME, tok[1]);

    // --- read in sewer area value
    if ( !getDouble(tok[2], &a) || a < 0.0 )
        return error_setInpError(ERR_NUMBER, tok[2]);

    // --- create the RDII inflow object if it doesn't already exist
    inflow = project->Node[j].rdiiInflow;
    if ( inflow == NULL )
    {
        inflow = (TRdiiInflow *) malloc(sizeof(TRdiiInflow));
        if ( !inflow ) return error_setInpError(ERR_MEMORY, "");
    }

    // --- assign UH & area to inflow object
    inflow->unitHyd = k;
    inflow->area = a / UCF(project,LANDAREA);

    // --- assign inflow object to node
    project->Node[j].rdiiInflow = inflow;
    return 0;
}

//=============================================================================

void rdii_initUnitHyd(Project* project, int j)
//
//  Input:   j = UH group index
//  Output:  none
//  Purpose: initializes properties of a unit hydrograph group.
//
{
    int i;                             // individual UH index
    int m;                             // month index

    for ( m=0; m<12; m++)
    {
        for (i=0; i<3; i++)
        {
            project->UnitHyd[j].iaMax[m][i]   = 0.0;
            project->UnitHyd[j].iaRecov[m][i] = 0.0;
            project->UnitHyd[j].iaInit[m][i]  = 0.0;
            project->UnitHyd[j].r[m][i]       = 0.0;
            project->UnitHyd[j].tPeak[m][i]   = 0;
            project->UnitHyd[j].tBase[m][i]   = 0;
        }
    }
}

//=============================================================================

int rdii_readUnitHydParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads parameters of an RDII unit hydrograph from a line of input.
//
{
    int i, j, k, m, g;
    double x[6];

    // --- check that RDII UH object exists in database
    j = project_findObject(project,UNITHYD, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);

    // --- assign UH ID to name in hash table
    if ( project->UnitHyd[j].ID == NULL )
        project->UnitHyd[j].ID = project_findID(project,UNITHYD, tok[0]);

    // --- line has 2 tokens; assign rain gage to UH object
    if ( ntoks == 2 )
    {
        g = project_findObject(project, GAGE, tok[1]);
        if ( g < 0 ) return error_setInpError(ERR_NAME, tok[1]);
        project->UnitHyd[j].rainGage = g;
        return 0;
    }
    else if ( ntoks < 6 ) return error_setInpError(ERR_ITEMS, "");

    // --- find which month UH params apply to
    m = datetime_findMonth(tok[1]);
    if ( m == 0 )
    {
        if ( !match(tok[1], w_ALL) )
            return error_setInpError(ERR_KEYWORD, tok[1]);
    }

    // --- find type of UH being specified
    k = findmatch(tok[2], UHTypeWords);

    // --- if no type match, try using older UH line format
    if ( k < 0 ) return readOldUHFormat(project,j, m, tok, ntoks);

    // --- read the R-T-K parameters
    for ( i = 0; i < 3; i++ )
    {
        if ( ! getDouble(tok[i+3], &x[i]) )
            return error_setInpError(ERR_NUMBER, tok[i+3]);
    }

    // --- read the IA parameters if present
    for (i = 3; i < 6; i++)
    {
        x[i] = 0.0;
        if ( ntoks > i+3 )
        {
            if ( ! getDouble(tok[i+3], &x[i]) )
                return error_setInpError(ERR_NUMBER, tok[i+2]);
        }
    }

    // --- save UH params
    setUnitHydParams(project,j, k, m, x);
    return 0;
}

//=============================================================================

int readOldUHFormat(Project* project, int j, int m, char* tok[], int ntoks)
//
//  Input:   j = unit hydrograph index
//           m = month of year (0 = all months)
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads parameters of a set of RDII unit hydrographs from a line of
//           input.
//
{
    int    i, k;
    double p[9], x[6];

    // --- check for proper number of tokens
    if ( ntoks < 11 ) return error_setInpError(ERR_ITEMS, "");

    // --- read 3 sets of r-t-k values
    for ( i = 0; i < 9; i++ )
    {
        if ( ! getDouble(tok[i+2], &p[i]) )
            return error_setInpError(ERR_NUMBER, tok[i+2]);
    }

    // --- read initial abstraction parameters
    for (i = 0; i < 3; i++)
    {
        x[i+3] = 0.0;
        if ( ntoks > i+11 )
        {
            if ( ! getDouble(tok[i+11], &x[i+3]) )
                return error_setInpError(ERR_NUMBER, tok[i+11]);
        }
    }

    // --- save UH parameters
    for ( k = 0; k < 3; k++)
    {
        for ( i = 0; i < 3; i++)
        {
            x[i] = p[3*k + i];
            setUnitHydParams(project,j, k, m, x);
        }
    }
    return 0;
}

//=============================================================================

void setUnitHydParams(Project* project, int j, int i, int m, double x[])
//
//  Input:   j = unit hydrograph index
//           i = type of UH response (short, medium or long term)
//           m = month of year (0 = all months)
//           x = array of UH parameters
//  Output:  none
//  Purpose: assigns parameters to a unit hydrograph for a specified month of year.
//
{
    int    m1, m2;                     // start/end month indexes
    double t,                          // UH time to peak (hrs)
           k,                          // UH k-value
           tBase;                      // UH base time (hrs)

    // --- find range of months that share same parameter values
    if ( m == 0 )
    {
        m1 = 0;
        m2 = 11;
    }
    else
    {
        m1 = m-1;
        m2 = m1;
    }

    // --- for each month in the range
    for (m=m1; m<=m2; m++)
    {
        // --- set UH response ratio, time to peak, & base time
        project->UnitHyd[j].r[m][i] = x[0];
        t = x[1];
        k = x[2];
        tBase = t * (1.0 + k);                              // hours
        project->UnitHyd[j].tPeak[m][i] = (long)(t * 3600.);         // seconds
        project->UnitHyd[j].tBase[m][i] = (long)(tBase * 3600.);     // seconds

        // -- set initial abstraction parameters
        project->UnitHyd[j].iaMax[m][i]   = x[3];
        project->UnitHyd[j].iaRecov[m][i] = x[4];
        project->UnitHyd[j].iaInit[m][i]  = x[5];
    }
}

//=============================================================================

void rdii_deleteRdiiInflow(Project* project, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: deletes the RDII inflow object for a node.
//
{
    if ( project->Node[j].rdiiInflow )
    {
        free(project->Node[j].rdiiInflow);
        project->Node[j].rdiiInflow = NULL;
    }
}


//=============================================================================
//                 Reading Inflow Data From a RDII File
//=============================================================================

void rdii_openRdii(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens an exisiting RDII interface file or creates a new one.
//
{
    char  fStamp[] = FILE_STAMP;

    project->RdiiNodeIndex = NULL;
    project->RdiiNodeFlow = NULL;
    project->NumRdiiNodes = 0;
    project->RdiiStartDate = NO_DATE;

    // --- create the RDII file if existing file not being used
    if ( project->IgnoreRDII ) return;                                                  //(5.1.004)
    if ( project->Frdii.mode != USE_FILE ) createRdiiFile(project);
    if ( project->Frdii.mode == NO_FILE || project->ErrorCode ) return;

    // --- try to open the RDII file in binary mode
    project->Frdii.file = fopen(project->Frdii.name, "rb");
    if ( project->Frdii.file == NULL)
    {
        if ( project->Frdii.mode == SCRATCH_FILE )
        {
            report_writeErrorMsg(project,ERR_RDII_FILE_SCRATCH, "");
        }
        else
        {
            report_writeErrorMsg(project,ERR_RDII_FILE_OPEN, project->Frdii.name);
        }
        return;
    }

    // --- check for valid file stamp
    fread(fStamp, sizeof(char), strlen(FileStamp), project->Frdii.file);
    if ( strcmp(fStamp, FileStamp) == 0 )
    {
        project->RdiiFileType = BINARY;
        project->ErrorCode = readRdiiFileHeader(project);
    }

    // --- if stamp invalid try to open the file in text mode
    else
    {
        fclose(project->Frdii.file);
        project->RdiiFileType = TEXT;
        openRdiiTextFile(project);
    }

    // --- catch any error
    if ( project->ErrorCode )
    {
        report_writeErrorMsg(project,project->ErrorCode, project->Frdii.name);
    }

    // --- read the first set of RDII flows form the file
    else readRdiiFlows(project);
}

//=============================================================================

void openRdiiTextFile(Project* project)
{
    // --- try to open the RDII file in text mode
    project->Frdii.file = fopen(project->Frdii.name, "rt");
    if ( project->Frdii.file == NULL)
    {
        if ( project->Frdii.mode == SCRATCH_FILE )
        {
            report_writeErrorMsg(project,ERR_RDII_FILE_SCRATCH, "");
        }
        else
        {
            report_writeErrorMsg(project,ERR_RDII_FILE_OPEN, project->Frdii.name);
        }
        return;
    }

    // --- read header records from file
    project->ErrorCode = readRdiiTextFileHeader(project);
    if ( project->ErrorCode )
    {
        report_writeErrorMsg(project,project->ErrorCode, project->Frdii.name);
    }
}

//=============================================================================

void rdii_closeRdii(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: closes the RDII interface file.
//
{
    if ( project->Frdii.file ) fclose(project->Frdii.file);
    if ( project->Frdii.mode == SCRATCH_FILE ) remove(project->Frdii.name);
    FREE(project->RdiiNodeIndex);
    FREE(project->RdiiNodeFlow);
}

//=============================================================================

int rdii_getNumRdiiFlows(Project* project, DateTime aDate)
//
//  Input:   aDate = current date/time
//  Output:  returns 0 if no RDII flow or number of nodes with RDII inflows
//  Purpose: finds number of RDII inflows at a specified date.
//
{
    // --- default result is 0 indicating no RDII inflow at specified date
    if ( project->NumRdiiNodes == 0 ) return 0;
    if ( !project->Frdii.file ) return 0;

    // --- keep reading RDII file as need be
    while ( !feof(project->Frdii.file) )
    {
        // --- return if date of current RDII inflow not reached yet
        if ( project->RdiiStartDate == NO_DATE ) return 0;
        if ( aDate < project->RdiiStartDate ) return 0;

        // --- return RDII node count if specified date falls
        //     within time interval of current RDII inflow
        if ( aDate < project->RdiiEndDate ) return project->NumRdiiNodes;

        // --- otherwise get next date and RDII flow values from file
        else readRdiiFlows(project);
    }
    return 0;
}

//=============================================================================

void rdii_getRdiiFlow(Project* project, int i, int* j, double* q)
//
//  Input:   i = RDII node index
//           j = pointer to project node index
//           q = pointer to RDII flow rate
//  Output:  sets node index and RDII inflow for node
//  Purpose: finds index and current RDII inflow for an RDII node.
//
{
    if ( i >= 0 && i < project->NumRdiiNodes )
    {
        *j = project->RdiiNodeIndex[i];
        *q = project->RdiiNodeFlow[i];
    }
}

//=============================================================================

int readRdiiFileHeader(Project* project)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads header information from a binary RDII file.
//
{
    int i, j;

    // --- extract time step and number of RDII nodes
    fread(&project->RdiiStep, sizeof(INT4), 1, project->Frdii.file);
    if ( project->RdiiStep <= 0 ) return ERR_RDII_FILE_FORMAT;
    fread(&project->NumRdiiNodes, sizeof(INT4), 1, project->Frdii.file);
    if ( project->NumRdiiNodes <= 0 ) return ERR_RDII_FILE_FORMAT;

    // --- allocate memory for project->RdiiNodeIndex & project->RdiiNodeFlow arrays
    project->RdiiNodeIndex = (int *) calloc(project->NumRdiiNodes, sizeof(int));
    if ( !project->RdiiNodeIndex ) return ERR_MEMORY;
    project->RdiiNodeFlow = (REAL4 *) calloc(project->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !project->RdiiNodeFlow ) return ERR_MEMORY;

    // --- read indexes of RDII nodes
    if ( feof(project->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    fread(project->RdiiNodeIndex, sizeof(INT4), project->NumRdiiNodes, project->Frdii.file);
    for ( i=0; i<project->NumRdiiNodes; i++ )
    {
        j = project->RdiiNodeIndex[i];
        if ( project->Node[j].rdiiInflow == NULL ) return ERR_RDII_FILE_FORMAT;
    }
    if ( feof(project->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    return 0;
}

//=============================================================================

int readRdiiTextFileHeader(Project* project)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads header information from a text RDII file.
//
{
    int   i;
    char  line[MAXLINE+1];             // line from RDII data file
    char  s1[MAXLINE+1];               // general string variable
    char  s2[MAXLINE+1];

    // --- check for correct file type
    fgets(line, MAXLINE, project->Frdii.file);
    sscanf(line, "%s", s1);
    if ( strcmp(s1, "SWMM5") != 0 ) return ERR_RDII_FILE_FORMAT;

    // --- skip title line
    fgets(line, MAXLINE, project->Frdii.file);

    // --- read RDII UH time step interval (sec)
    project->RdiiStep = 0;
    fgets(line, MAXLINE, project->Frdii.file);
    sscanf(line, "%d", &project->RdiiStep);
    if ( project->RdiiStep <= 0 ) return ERR_RDII_FILE_FORMAT;

    // --- skip over line with number of constituents (= 1 for RDII)
    fgets(line, MAXLINE, project->Frdii.file);

    // --- read flow units
    fgets(line, MAXLINE, project->Frdii.file);
    sscanf(line, "%s %s", s1, s2);
    project->RdiiFlowUnits = findmatch(s2, FlowUnitWords);
    if ( project->RdiiFlowUnits < 0 ) return ERR_RDII_FILE_FORMAT;

    // --- read number of RDII nodes
    fgets(line, MAXLINE, project->Frdii.file);
    if ( sscanf(line, "%d", &project->NumRdiiNodes) < 1 ) return ERR_RDII_FILE_FORMAT;

    // --- allocate memory for project->RdiiNodeIndex & project->RdiiNodeFlow arrays
    project->RdiiNodeIndex = (int *) calloc(project->NumRdiiNodes, sizeof(int));
    if ( !project->RdiiNodeIndex ) return ERR_MEMORY;
    project->RdiiNodeFlow = (REAL4 *) calloc(project->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !project->RdiiNodeFlow ) return ERR_MEMORY;

    // --- read names of RDII nodes from file & save their indexes
    for ( i=0; i<project->NumRdiiNodes; i++ )
    {
        if ( feof(project->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
        fgets(line, MAXLINE, project->Frdii.file);
        sscanf(line, "%s", s1);
        project->RdiiNodeIndex[i] = project_findObject(project,NODE, s1);
    }

    // --- skip column heading line
    if ( feof(project->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    fgets(line, MAXLINE, project->Frdii.file);
    return 0;
}

//=============================================================================

void readRdiiFlows(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reads date and flow values of next RDII inflows from RDII file.
//
{
    if ( project->RdiiFileType == TEXT ) readRdiiTextFlows(project);
    else
    {
        project->RdiiStartDate = NO_DATE;
        project->RdiiEndDate = NO_DATE;
        if ( feof(project->Frdii.file) ) return;
        fread(&project->RdiiStartDate, sizeof(DateTime), 1, project->Frdii.file);
        if ( project->RdiiStartDate == NO_DATE ) return;
        if ( fread(project->RdiiNodeFlow, sizeof(REAL4), project->NumRdiiNodes, project->Frdii.file)      //(5.1.003)
            < (size_t)project->NumRdiiNodes ) project->RdiiStartDate = NO_DATE;
        else project->RdiiEndDate = datetime_addSeconds(project->RdiiStartDate, project->RdiiStep);
    }
}

//=============================================================================

void readRdiiTextFlows(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reads date and flow values of next RDII inflows from RDII file.
//
{
    int    i, n;
    int    yr = 0, mon = 0, day = 0,
		   hr = 0, min = 0, sec = 0;   // year, month, day, hour, minute, second
    double x;                          // RDII flow in original units          //(5.1.003)
    char   line[MAXLINE+1];            // line from RDII data file
    char   s[MAXLINE+1];               // node ID label (not used)

    project->RdiiStartDate = NO_DATE;
    for (i=0; i<project->NumRdiiNodes; i++)
    {
        if ( feof(project->Frdii.file) ) return;
        fgets(line, MAXLINE, project->Frdii.file);
        n = sscanf(line, "%s %d %d %d %d %d %d %f",
            s, &yr, &mon, &day, &hr, &min, &sec, &x);
        if ( n < 8 ) return;
        project->RdiiNodeFlow[i] = (REAL4)(x / Qcf[project->RdiiFlowUnits]);                     //(5.1.003)
    }
    project->RdiiStartDate = datetime_encodeDate(yr, mon, day) +
                    datetime_encodeTime(hr, min, sec);
    project->RdiiEndDate = datetime_addSeconds(project->RdiiStartDate, project->RdiiStep);
}


//=============================================================================
//                   Creation of a RDII Interface File
//=============================================================================

void createRdiiFile(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: computes time history of RDII inflows and saves them to file.
//
{
    int      hasRdii;                  // true when total RDII > 0
    double   elapsedTime;              // current elapsed time (sec)
    double   duration;                 // duration being analyzed (sec)
    DateTime currentDate;              // current calendar date/time

    // --- set RDII reporting time step to Runoff wet step
    project->RdiiStep = project->WetStep;

    // --- count nodes with RDII data
    project->NumRdiiNodes = getNumRdiiNodes(project);

    // --- if no RDII nodes then re-set RDII file usage to NO_FILE
    if ( project->NumRdiiNodes == 0 )
    {
        project->Frdii.mode = NO_FILE;
        return;
    }

    // --- otherwise set file usage to SCRATCH if originally set to NO_FILE
    else if ( project->Frdii.mode == NO_FILE ) project->Frdii.mode = SCRATCH_FILE;

    // --- validate RDII data
    validateRdii(project);
    initGageData(project);
    if ( project->ErrorCode ) return;

    // --- open RDII processing system
	openRdiiProcessor(project);
    if ( !project->ErrorCode )
    {
        // --- initialize rain gage & UH processing data
		initUnitHydData(project);

        // --- convert total simulation duration from millisec to sec
        duration = project->TotalDuration / 1000.0;

        // --- examine rainfall record over each project->RdiiStep time step
        elapsedTime = 0.0;
        while ( elapsedTime <= duration && !project->ErrorCode )
        {
            // --- compute current calendar date/time
            currentDate = project->StartDateTime + elapsedTime / SECperDAY;

            // --- update rainfall at all rain gages
            getRainfall(project,currentDate);

            // --- compute convolutions of past rainfall with UH's
            getUnitHydRdii(project,currentDate);

            // --- find RDII at all nodes
			hasRdii = getNodeRdii(project);

            // --- save RDII at all nodes to file for current date
            if ( hasRdii ) saveRdiiFlows(project,currentDate);

            // --- advance one time step
            elapsedTime += project->RdiiStep;
        }
    }

    // --- close RDII processing system
	closeRdiiProcessor(project);
}

//=============================================================================

int  getNumRdiiNodes(Project* project)
//
//  Input:   none
//  Output:  returns node count
//  Purpose: counts number of nodes that receive RDII inflow.
//
{
    int j,                             // node index
        n;                             // node count

    n = 0;
    for (j=0; j<project->Nobjects[NODE]; j++)
    {
        if ( project->Node[j].rdiiInflow ) n++;
    }
    return n;
}

//=============================================================================

void validateRdii(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: validates UH and RDII inflow object data.
//
{
    int    i,                          // node index
           j,                          // UH group index
           k,                          // individual UH index
           m;                          // month index
    double rsum;                       // sum of UH r-values
//  long   gageInterval;               // rain gage time interval

    // --- check each unit hydrograph for consistency
    for (j=0; j<project->Nobjects[UNITHYD]; j++)
    {
        for (m=0; m<12; m++)
        {
            rsum = 0.0;
            for (k=0; k<3; k++)
            {
                // --- if no base time then UH doesn't exist
                if ( project->UnitHyd[j].tBase[m][k] == 0 ) continue;

                // --- restriction on time to peak being less than the
                //     rain gage's recording interval no longer applies

                // --- can't have negative UH parameters
                if ( project->UnitHyd[j].tPeak[m][k] < 0.0 )
                {
                    report_writeErrorMsg(project,ERR_UNITHYD_TIMES, project->UnitHyd[j].ID);
                }

                // --- can't have negative UH response ratio
                if ( project->UnitHyd[j].r[m][k] < 0.0 )
                {
                    report_writeErrorMsg(project,ERR_UNITHYD_RATIOS, project->UnitHyd[j].ID);
                }
                else rsum += project->UnitHyd[j].r[m][k];
            }
            if ( rsum > 1.01 )
            {
                report_writeErrorMsg(project,ERR_UNITHYD_RATIOS, project->UnitHyd[j].ID);
            }
        }
    }

    // --- check each node's RDII inflow object
    for (i=0; i<project->Nobjects[NODE]; i++)
    {
        if ( project->Node[i].rdiiInflow )
        {
            // --- check that sewer area is non-negative
            if ( project->Node[i].rdiiInflow->area < 0.0 )
            {
                report_writeErrorMsg(project,ERR_RDII_AREA, project->Node[i].ID);
            }
        }
    }
}

//=============================================================================

void openRdiiProcessor(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens RDII processing system.
//
{
    int j;                             // object index
    int n;                             // RDII node count

    // --- set RDII processing arrays to NULL
    project->UHGroup = NULL;
    project->RdiiNodeIndex = NULL;
    project->RdiiNodeFlow = NULL;
    project->TotalRainVol = 0.0;
    project->TotalRdiiVol = 0.0;

    // --- allocate memory used for RDII processing
	if (!allocRdiiMemory(project))
    {
        report_writeErrorMsg(project,ERR_MEMORY, "");
        return;
    }

    // --- open & initialize RDII file
	if (!openNewRdiiFile(project))
    {
        report_writeErrorMsg(project,ERR_RDII_FILE_SCRATCH, "");
        return;
    }

    // --- identify index of each node with RDII inflow
    n = 0;
    for (j=0; j<project->Nobjects[NODE]; j++)
    {
        if ( project->Node[j].rdiiInflow )
        {
            project->RdiiNodeIndex[n] = j;
            n++;
        }
    }
}

//=============================================================================

int  allocRdiiMemory(Project* project)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: allocates memory used for RDII processing .
//
//
{
    int i;                             // UH group index
    int k;                             // UH index
    int n;                             // number of past rain periods

    // --- allocate memory for RDII processing data for UH groups
    project->UHGroup = (TUHGroup *) calloc(project->Nobjects[UNITHYD], sizeof(TUHGroup));
    if ( !project->UHGroup ) return FALSE;

    // --- allocate memory for past rainfall data for each UH in each group
    for (i=0; i<project->Nobjects[UNITHYD]; i++)
    {
        project->UHGroup[i].rainInterval = getRainInterval(project,i);
        for (k=0; k<3; k++)
        {
            project->UHGroup[i].uh[k].pastRain = NULL;
            project->UHGroup[i].uh[k].pastMonth = NULL;
            project->UHGroup[i].uh[k].maxPeriods = getMaxPeriods(project,i, k);
            n = project->UHGroup[i].uh[k].maxPeriods;
            if ( n > 0 )
            {
                project->UHGroup[i].uh[k].pastRain =
                    (double *) calloc(n, sizeof(double));
                if ( !project->UHGroup[i].uh[k].pastRain ) return FALSE;
                project->UHGroup[i].uh[k].pastMonth =
                    (char *) calloc(n, sizeof(char));
                if ( !project->UHGroup[i].uh[k].pastMonth ) return FALSE;
            }
        }
    }

    // --- allocate memory for RDII indexes & inflow at each node w/ RDII data
    project->RdiiNodeIndex = (int *) calloc(project->NumRdiiNodes, sizeof(int));
    if ( !project->RdiiNodeIndex ) return FALSE;
    project->RdiiNodeFlow = (REAL4 *) calloc(project->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !project->RdiiNodeFlow ) return FALSE;
    return TRUE;
}

//=============================================================================

int  getRainInterval(Project* project, int i)
//
//  Input:   i = UH group index
//  Output:  returns a time interval (sec)
//  Purpose: finds rainfall processing time interval for a unit hydrograph group.
//
{
    int ri;        // rainfal processing time interval for the UH group
    int tLimb;     // duration of a UH's rising & falling limbs
    int k, m;

    // --- begin with UH group time step equal to wet runoff step
    ri = project->WetStep;

    // --- examine each UH in the group
    for (m=0; m<12; m++)
    {
        for (k=0; k<3; k++)
        {
            // --- make sure the UH exists
            if ( project->UnitHyd[i].tPeak[m][k] > 0 )
            {
                // --- reduce time step if rising/falling limb is smaller
                tLimb = project->UnitHyd[i].tPeak[m][k];
                ri = MIN(ri, tLimb);
                tLimb = project->UnitHyd[i].tBase[m][k] - tLimb;
                if ( tLimb > 0 ) ri = MIN(ri, tLimb);
            }
        }
    }
    return ri;
}

//=============================================================================

int  getMaxPeriods(Project* project, int i, int k)
//
//  Input:   i = UH group index
//           k = UH index
//  Output:  returns number of past rainfall values
//  Purpose: finds number of past rainfall values to save for a UH.
//
{
    int   m,                           // month index
          n,                           // number of time periods
          nMax,                        // maximum number of time periods
          rainInterval;                // rainfall processing interval (sec)

    // --- examine each monthly set of UHs
    rainInterval = project->UHGroup[i].rainInterval;
    nMax = 0;
    for (m=0; m<12; m++)
    {
        // --- compute number of time periods in UH base
        n = (project->UnitHyd[i].tBase[m][k] / rainInterval) + 1;

        // --- update number of time periods to be saved
        nMax = MAX(n, nMax);
    }
    return nMax;
}

//=============================================================================

void initGageData(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes state of Unit Hydrograph rain gages.
//
{
    int i;                             // unit hyd. index
    int g;                             // rain gage index

    // --- first initialize the state of each rain gage
    for (g=0; g<project->Nobjects[GAGE]; g++)
    {
        if ( project->Gage[g].tSeries >= 0 )
        {
            table_tseriesInit(&project->Tseries[project->Gage[g].tSeries]);
        }
        gage_initState(project,g);
    }

    // --- then flag each gage that is used by a Unit Hydrograph set
    for (i=0; i<project->Nobjects[UNITHYD]; i++)
    {
        g = project->UnitHyd[i].rainGage;
        if ( g >= 0 )
        {
            project->Gage[g].isUsed = TRUE;

            // --- if UH's gage uses same time series as a previous gage,
            //     then assign the latter gage to the UH
            if ( project->Gage[g].coGage >= 0 )
            {
                project->UnitHyd[i].rainGage = project->Gage[g].coGage;
                project->Gage[project->Gage[g].coGage].isUsed = TRUE;
            }
        }
    }
}

//=============================================================================

void initUnitHydData(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes unit hydrograph processing data.
//
{
    int i,                             // UH group index
        j,                             // node index
        k,                             // UH index
        n;                             // RDII node index
//  int g,                             // rain gage index
    int month;                         // month index

    // --- initialize project->UHGroup entries for each Unit Hydrograph
    month = datetime_monthOfYear(project->StartDateTime) - 1;
    for (i=0; i<project->Nobjects[UNITHYD]; i++)
    {
        for (k=0; k<3; k++)
        {
            // --- make the first recorded rainfall begin a new RDII event
            // --- (new RDII event occurs when dry period > base of longest UH)
            project->UHGroup[i].uh[k].drySeconds =
                (project->UHGroup[i].uh[k].maxPeriods * project->UHGroup[i].rainInterval) + 1;
            project->UHGroup[i].uh[k].period = project->UHGroup[i].uh[k].maxPeriods + 1;
            project->UHGroup[i].uh[k].hasPastRain = FALSE;

            // --- assign initial abstraction used
            project->UHGroup[i].uh[k].iaUsed = project->UnitHyd[i].iaInit[month][k];
        }

        // --- initialize gage date to simulation start date
        project->UHGroup[i].gageDate = project->StartDateTime;
        project->UHGroup[i].area = 0.0;
        project->UHGroup[i].rdii = 0.0;
    }

    // --- assume each UH group is not used
    for (i=0; i<project->Nobjects[UNITHYD]; i++) project->UHGroup[i].isUsed = FALSE;

    // --- look at each node with RDII inflow
    for (n=0; n<project->NumRdiiNodes; n++)
    {
        // --- mark as used the UH group associated with the node
        j = project->RdiiNodeIndex[n];
        i = project->Node[j].rdiiInflow->unitHyd;
        project->UHGroup[i].isUsed = TRUE;

        // --- add node's sewer area to UH group's area
        project->UHGroup[i].lastDate = project->StartDateTime;
        project->UHGroup[i].area += project->Node[j].rdiiInflow->area;
    }
}

//=============================================================================

int openNewRdiiFile(Project* project)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: opens a new RDII interface file.
//
{
    int j;                             // node index

    // --- create a temporary file name if scratch file being used
    if ( project->Frdii.mode == SCRATCH_FILE ) getTempFileName(project, project->Frdii.name);

    // --- open the RDII file as a formatted text file
    project->Frdii.file = fopen(project->Frdii.name, "w+b");
    if ( project->Frdii.file == NULL )
    {
        return FALSE;
    }

    // --- write file stamp to RDII file
    fwrite(FileStamp, sizeof(char), strlen(FileStamp), project->Frdii.file);

    // --- initialize the contents of the file with RDII time step (sec),
    //     number of RDII nodes, and index of each node
    fwrite(&project->RdiiStep, sizeof(INT4), 1, project->Frdii.file);
    fwrite(&project->NumRdiiNodes, sizeof(INT4), 1, project->Frdii.file);
    for (j=0; j<project->Nobjects[NODE]; j++)
    {
        if ( project->Node[j].rdiiInflow ) fwrite(&j, sizeof(INT4), 1, project->Frdii.file);
    }
    return TRUE;
}

//=============================================================================

void getRainfall(Project* project, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: determines rainfall at current RDII processing date.
//
//
{
    int      j;                        // UH group index
    int      k;                        // UH index
    int      g;                        // rain gage index
    int      i;                        // past rainfall index
    int      month;                    // month of current date
    int      rainInterval;             // rainfall interval (sec)
    double   rainDepth;                // rainfall depth (inches or mm)
    double   excessDepth;              // excess rainfall depth (inches or mm))
    DateTime gageDate;                 // calendar date for rain gage

    // --- examine each UH group
    month = datetime_monthOfYear(currentDate) - 1;
    for (g = 0; g < project->Nobjects[GAGE]; g++) project->Gage[g].isCurrent = FALSE;
    for (j = 0; j < project->Nobjects[UNITHYD]; j++)
    {
        // --- repeat until gage's date reaches or exceeds current date
        g = project->UnitHyd[j].rainGage;
        rainInterval = project->UHGroup[j].rainInterval;
        while ( project->UHGroup[j].gageDate < currentDate )
        {
            // --- get rainfall volume over gage's recording interval
            //     at gage'a current date (in original depth units)
            gageDate = project->UHGroup[j].gageDate;
            project->Adjust.rainFactor = project->Adjust.rain[datetime_monthOfYear(gageDate)-1]; //(5.1.007)
            if (!project->Gage[g].isCurrent)
            {
                gage_setState(project, g, gageDate);
                project->Gage[g].isCurrent = TRUE;
            }
            rainDepth = project->Gage[g].rainfall * (double)rainInterval / 3600.0;

            // --- update amount of total rainfall volume (ft3)
            project->TotalRainVol += rainDepth / UCF(project,RAINDEPTH) * project->UHGroup[j].area;

            // --- compute rainfall excess for each UH in the group
            for (k=0; k<3; k++)
            {
                // --- adjust rainfall volume for any initial abstraction
                excessDepth = applyIA(project,j, k, gageDate, rainInterval, rainDepth);

                // --- adjust extent of dry period for the UH
                updateDryPeriod(project,j, k, excessDepth, rainInterval);

                // --- add rainfall to list of past values,
                //     wrapping array index if necessary
                i = project->UHGroup[j].uh[k].period;
                if ( i >= project->UHGroup[j].uh[k].maxPeriods ) i = 0;
                project->UHGroup[j].uh[k].pastRain[i] = excessDepth;
                project->UHGroup[j].uh[k].pastMonth[i] = (char)month;
                project->UHGroup[j].uh[k].period = i + 1;
            }

            // --- advance rain date by gage recording interval
            project->UHGroup[j].gageDate = datetime_addSeconds(gageDate, rainInterval);
        }
    }
}

//=============================================================================

double  applyIA(Project* project, int j, int k, DateTime aDate, double dt, double rainDepth)
//
//  Input:   j = UH group index
//           k = unit hydrograph index
//           aDate = current date/time
//           dt = time interval (sec)
//           rainDepth = unadjusted rain depth (in or mm)
//  Output:  returns rainfall adjusted for initial abstraction (IA)
//  Purpose: adjusts rainfall for any initial abstraction and updates the
//           amount of available initial abstraction actually used.
//
{
    int m;
    double ia, netRainDepth;

    // --- determine amount of unused IA
    m = datetime_monthOfYear(aDate) - 1;
    ia = project->UnitHyd[j].iaMax[m][k] - project->UHGroup[j].uh[k].iaUsed;
    ia = MAX(ia, 0.0);

    // --- case where there's some rainfall
    if ( rainDepth > 0.0 )
    {
        // --- reduce rain depth by unused IA
        netRainDepth = rainDepth - ia;
        netRainDepth = MAX(netRainDepth, 0.0);

        // --- update amount of IA used up
        ia = rainDepth - netRainDepth;
        project->UHGroup[j].uh[k].iaUsed += ia;
    }

    // --- case where there's no rainfall
    else
    {
        // --- recover a portion of the IA already used
        project->UHGroup[j].uh[k].iaUsed -= dt / 86400. * project->UnitHyd[j].iaRecov[m][k];
        project->UHGroup[j].uh[k].iaUsed = MAX(project->UHGroup[j].uh[k].iaUsed, 0.0);
        netRainDepth = 0.0;
    }
    return netRainDepth;
}

//=============================================================================

void updateDryPeriod(Project* project, int j, int k, double rainDepth, int rainInterval)
//
//  Input:   j = UH group index
//           k = unit hydrograph index
//           rainDepth = excess rain depth (in or mm)
//           rainInterval = rainfall time interval (sec)
//  Output:  none
//  Purpose: adjusts the length of the dry period between rainfall events.
//
{
    int i;

    // --- if rainfall occurs
    if ( rainDepth > 0.0 )
    {
        // --- if previous dry period long enough then begin
        //     new RDII event with time period index set to 0
        if ( project->UHGroup[j].uh[k].drySeconds >= rainInterval *
            project->UHGroup[j].uh[k].maxPeriods )
        {
            for (i=0; i<project->UHGroup[j].uh[k].maxPeriods; i++)
            {
                project->UHGroup[j].uh[k].pastRain[i] = 0.0;
            }
            project->UHGroup[j].uh[k].period = 0;
        }
        project->UHGroup[j].uh[k].drySeconds = 0;
        project->UHGroup[j].uh[k].hasPastRain = TRUE;
    }

    // --- if no rainfall, update duration of dry period
    else
    {
        project->UHGroup[j].uh[k].drySeconds += rainInterval;
        if ( project->UHGroup[j].uh[k].drySeconds >=
            rainInterval * project->UHGroup[j].uh[k].maxPeriods )
        {
            project->UHGroup[j].uh[k].hasPastRain = FALSE;
        }
        else project->UHGroup[j].uh[k].hasPastRain = TRUE;
    }
}

//=============================================================================

void getUnitHydRdii(Project* project, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: computes RDII generated by past rainfall for each UH group.
//
{
    int   j;                           // UH group index
    int   k;                           // UH index
    int   rainInterval;                // rainfall time interval (sec)

    // --- examine each UH group
    for (j=0; j<project->Nobjects[UNITHYD]; j++)
    {
        // --- skip calculation if group not used by any RDII node or if
        //     current date hasn't reached last date RDII was computed
        if ( !project->UHGroup[j].isUsed ) continue;
        if ( currentDate < project->UHGroup[j].lastDate ) continue;

        // --- update date RDII last computed
        project->UHGroup[j].lastDate = project->UHGroup[j].gageDate;

        // --- perform convolution for each UH in the group
        rainInterval = project->UHGroup[j].rainInterval;
        project->UHGroup[j].rdii = 0.0;
        for (k=0; k<3; k++)
        {
            if ( project->UHGroup[j].uh[k].hasPastRain )
            {
                project->UHGroup[j].rdii += getUnitHydConvol(project,j, k, rainInterval);
            }
        }
    }
}

//=============================================================================

double getUnitHydConvol(Project* project, int j, int k, int rainInterval)
//
//  Input:   j = UH group index
//           k = UH index
//           rainInterval = rainfall time interval (sec)
//  Output:  returns a RDII flow value
//  Purpose: computes convolution of Unit Hydrographs with past rainfall.
//
{
    int    i;                          // previous rainfall period index
    int    m;                          // month of year index
    int    p;                          // UH time period index
    int    pMax;                       // max. number of periods
    double t;                          // UH time value (sec)
    double u;                          // UH ordinate
    double v;                          // rainfall volume
    double rdii;                       // RDII flow
    TUHData* uh;                       // UH data

    // --- initialize RDII, rain period index and UH period index
    rdii = 0.0;
    uh = &project->UHGroup[j].uh[k];
    i = uh->period - 1;
    if ( i < 0 ) i = uh->maxPeriods - 1;
    pMax = uh->maxPeriods;
    p = 1;

    // --- evaluate each time period of UH's
    while ( p < pMax )
    {
        // --- if rain period has rainfall
        v = uh->pastRain[i];
        m = uh->pastMonth[i];
        if ( v > 0.0 )
        {
            // --- find mid-point time of UH period in seconds
            t = ((double)(p) - 0.5) * (double)rainInterval;

            // --- convolute rain volume with UH ordinate
            u = getUnitHydOrd(project,j, m, k, t) * project->UnitHyd[j].r[m][k];
            rdii += u * v;
        }

        // --- move to next UH period & previous rainfall period
        p = p + 1;
        i = i - 1;
        if ( i < 0 ) i = uh->maxPeriods - 1;
    }
    return rdii;
}

//=============================================================================

double getUnitHydOrd(Project* project, int h, int m, int k, double t)
//
//  Input:   h = index of UH group
//           m = month index
//           k = individual UH index
//           t = UH time (sec)
//  Output:  returns ordinate of a unit hydrograph
//  Purpose: gets ordinate of a particular unit hydrograph at specified time.
//
{
    double qPeak;                      // peak flow of unit hydrograph
    double f;                          // fraction of time to/from peak on UH
    double t1;                         // time to peak on UH (sec)
    double t2;                         // time after peak on UH (sec)
    double tBase;                      // base time of UH (sec)

    // --- return 0 if past end of UH time base
    tBase = project->UnitHyd[h].tBase[m][k];
    if ( t >= tBase ) return 0.0;

    // --- compute peak value of UH in original rainfall units (in/hr or mm/hr)
    qPeak = 2. / tBase * 3600.0;

    // --- break UH base into times before & after peak flow
    t1 = project->UnitHyd[h].tPeak[m][k];
    t2 = tBase - t1;

    // --- find UH flow at time t
    if ( t <= t1 ) f = t / t1;
    else           f = 1.0 - (t - t1) / t2;
    return MAX(f, 0.0) * qPeak;
}

//=============================================================================

int getNodeRdii(Project* project)
//
//  Input:   none
//  Output:  returns TRUE if any node has RDII inflow, FALSE if not
//  Purpose: computes current RDII inflow at each node.
//
{
    int   hasRdii = FALSE;             // true if any node has some RDII
    int   i;                           // UH group index
    int   j;                           // node index
    int   n;                           // number of nodes w/ RDII
    double rdii;                       // RDII flow (cfs)

    // --- examine each node w/ RDII data
    for (n = 0; n < project->NumRdiiNodes; n++)
    {
        // --- identify node's index in project's data base
        j = project->RdiiNodeIndex[n];

        // --- apply node's sewer area to UH RDII to get node RDII in CFS
        i = project->Node[j].rdiiInflow->unitHyd;
        rdii = project->UHGroup[i].rdii * project->Node[j].rdiiInflow->area / UCF(project,RAINFALL);
        if ( rdii < ZERO_RDII ) rdii = 0.0;
        else hasRdii = TRUE;

        // --- update total RDII volume
        project->RdiiNodeFlow[n] = (REAL4)rdii;
        if ( rdii > 0.0 )
        {
            project->TotalRdiiVol += rdii * (double)project->RdiiStep;
        }
    }
    return hasRdii;
}

//=============================================================================

void saveRdiiFlows(Project* project, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: saves current set of RDII inflows in current flow units to file.
//
{
    fwrite(&currentDate, sizeof(DateTime), 1, project->Frdii.file);
    fwrite(project->RdiiNodeFlow, sizeof(REAL4), project->NumRdiiNodes, project->Frdii.file);             //(5.1.003)
}

//=============================================================================

void  closeRdiiProcessor(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: closes RDII processing system.
//
{
    // --- write rainfall & RDII totals to report file
    if ( !project->ErrorCode )
    {
        report_writeRdiiStats(project,project->TotalRainVol, project->TotalRdiiVol);
    }

    // --- free allocated memory and close RDII file
    freeRdiiMemory(project);
    if ( project->Frdii.file ) fclose(project->Frdii.file);
}

//=============================================================================

void freeRdiiMemory(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory used for RDII processing.
//
{
    int i;
    int k;
    if ( project->UHGroup )
    {
        for (i = 0; i < project->Nobjects[UNITHYD]; i++)
        {
            for (k=0; k<3; k++)
            {
                FREE(project->UHGroup[i].uh[k].pastRain);
                FREE(project->UHGroup[i].uh[k].pastMonth);
            }
        }
        FREE(project->UHGroup);
    }
    FREE(project->RdiiNodeIndex);
    FREE(project->RdiiNodeFlow);
}
