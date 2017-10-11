//-----------------------------------------------------------------------------
//   globals.h
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/19/14  (Build 5.1.000)
//            04/14/14  (Build 5.1.004)
//            09/15/14  (Build 5.1.007)
//            03/19/15  (Build 5.1.008)
//   Author:  L. Rossman
//
//   Global Variables
//
//   Build 5.1.004:
//   - Ignore RDII option added.
//
//   Build 5.1.007:
//   - Monthly climate variable adjustments added.
//
//   Build 5.1.008:
//   - Number of parallel threads for dynamic wave routing added.
//   - Minimum dynamic wave routing variable time step added.
//
//-----------------------------------------------------------------------------

#ifndef GLOBALS_H
#define GLOBALS_H

#include "hash.h"
#include "enums.h"
#include "consts.h"
#include "objects.h"
#include <time.h>

#define USE_OPENMP

#ifdef __cplusplus
extern "C" {
#endif

//typedef struct alloc_root_s alloc_root_t;
//-----------------------------------------------------------------------------
//  Constants - moved from transect.c
//-----------------------------------------------------------------------------
#define MAXSTATION 1500                // max. number of stations in a transect
#define MAX_STATS 5
#define REAL4 float
#define INT4  int
#define REAL4 float
#define REAL8 double


struct Project
{
    //EXTERN TFile
    TFile	Finp;                     // Input file
    TFile	Fout;                     // Output file
    TFile	Frpt;                     // Report file
    TFile	Fclimate;                 // Climate file
    TFile	Frain;                    // Rainfall file
    TFile	Frunoff;                  // Runoff file
    TFile	Frdii;                    // RDII inflow file
    TFile	Fhotstart1;               // Hot start input file
    TFile	Fhotstart2;               // Hot start output file
    TFile	Finflows;                 // Inflows routing file
    TFile	Foutflows;                // Outflows routing file

    //EXTERN long
    long	Nperiods;                 // Number of reporting periods
    long	StepCount;                // Number of routing steps used
    long	NonConvergeCount;         // Number of non-converging steps

    //EXTERN char
    char	Msg[MAXMSG + 1];            // Text of output message
    char	Title[MAXTITLE][MAXMSG + 1];// Project title
    char	TempDir[MAXFNAME + 1];      // Temporary file directory

    //EXTERN TRptFlags
    TRptFlags	RptFlags;                 // Reporting options

    //EXTERN int
    int	Nobjects[MAX_OBJ_TYPES];  // Number of each object type
    int	Nnodes[MAX_NODE_TYPES];   // Number of each node sub-type
    int	Nlinks[MAX_LINK_TYPES];   // Number of each link sub-type
    int	UnitSystem;               // Unit system
    int	FlowUnits;                // Flow units
    int	InfilModel;               // Infiltration method
    int	RouteModel;               // Flow routing method
    int	ForceMainEqn;             // Flow equation for force mains
    int	LinkOffsets;              // Link offset convention
    int	AllowPonding;             // Allow water to pond at nodes
    int     InertDamping;             // Degree of inertial damping
    int	NormalFlowLtd;            // Normal flow limited
    int	SlopeWeighting;           // Use slope weighting
    int	Compatibility;            // SWMM 5/3/4 compatibility
    int	SkipSteadyState;          // Skip over steady state periods
    int	IgnoreRainfall;           // Ignore rainfall/runoff
    int	IgnoreRDII;               // Ignore RDII                     //(5.1.004)
    int	IgnoreSnowmelt;           // Ignore snowmelt
    int	IgnoreGwater;             // Ignore groundwater
    int	IgnoreRouting;            // Ignore flow routing
    int	IgnoreQuality;            // Ignore water quality
    int	ErrorCode;                // Error code number
    int	WarningCode;              // Warning code number
    int	WetStep;                  // Runoff wet time step (sec)
    int	DryStep;                  // Runoff dry time step (sec)
    int	ReportStep;               // Reporting time step (sec)
    int	SweepStart;               // Day of year when sweeping starts
    int	SweepEnd;                 // Day of year when sweeping ends
    int	MaxTrials;                // Max. trials for DW routing
    int	NumThreads;               // Number of parallel threads used //(5.1.008)

    //EXTERN double
    double	RouteStep;                // Routing time step (sec)
    double	MinRouteStep;             // Minimum variable time step (sec) //(5.1.008)
    double	LengtheningStep;          // Time step for lengthening (sec)
    double	StartDryDays;             // Antecedent dry days
    double	CourantFactor;            // Courant time step factor
    double	MinSurfArea;              // Minimum nodal surface area
    double	MinSlope;                 // Minimum conduit slope
    double	RunoffError;              // Runoff continuity error
    double	GwaterError;              // Groundwater continuity error
    double	FlowError;                // Flow routing error
    double	QualError;                // Quality routing error
    double	HeadTol;                  // DW routing head tolerance (ft)
    double	SysFlowTol;               // Tolerance for steady system flow
    double	LatFlowTol;               // Tolerance for steady nodal inflow

    //EXTERN DateTime
    DateTime	StartDate;                // Starting date
    DateTime	StartTime;                // Starting time
    DateTime	StartDateTime;            // Starting Date+Time
    DateTime	EndDate;                  // Ending date
    DateTime	EndTime;                  // Ending time
    DateTime	EndDateTime;              // Ending Date+Time
    DateTime	ReportStartDate;          // Report start date
    DateTime	ReportStartTime;          // Report start time
    DateTime	ReportStart;              // Report start Date+Time
    //added to support iteration without writing
    DateTime	PreviousReportTime;              // Report start Date+Time

    //EXTERN double
    double	ReportTime;               // Current reporting time (msec)
    double	OldRunoffTime;            // Previous runoff time (msec)
    double	NewRunoffTime;            // Current runoff time (msec)
    double	OldRoutingTime;           // Previous routing time (msec)
    double	NewRoutingTime;           // Current routing time (msec)
    double	TotalDuration;            // Simulation duration (msec)
    double  RoutingStep;
    //		double  PreviousReportTime;

    TTemp      Temp;                     // Temperature data
    TEvap      Evap;                     // Evaporation data
    TWind      Wind;                     // Wind speed data
    TSnow      Snow;                     // Snow melt data
    TAdjust    Adjust;                   // Climate adjustments             //(5.1.007)

    TSnowmelt* Snowmelt;                 // Array of snow melt objects
    TGage*     Gage;                     // Array of rain gages
    TSubcatch* Subcatch;                 // Array of subcatchments
    TAquifer*  Aquifer;                  // Array of groundwater aquifers
    TUnitHyd*  UnitHyd;                  // Array of unit hydrographs
    TNode*     Node;                     // Array of nodes
    TOutfall*  Outfall;                  // Array of outfall nodes
    TDivider*  Divider;                  // Array of divider nodes
    TStorage*  Storage;                  // Array of storage nodes
    TLink*     Link;                     // Array of links
    TConduit*  Conduit;                  // Array of conduit links
    TPump*     Pump;                     // Array of pump links
    TOrifice*  Orifice;                  // Array of orifice links
    TWeir*     Weir;                     // Array of weir links
    TOutlet*   Outlet;                   // Array of outlet device links
    TPollut*   Pollut;                   // Array of pollutants
    TLanduse*  Landuse;                  // Array of landuses
    TPattern*  Pattern;                  // Array of time patterns
    TTable*    Curve;                    // Array of curve tables
    TTable*    Tseries;                  // Array of time series tables
    TTransect* Transect;                 // Array of transect data
    TShape*    Shape;                    // Array of custom conduit shapes


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from swmm.c
    //-----------------------------------------------------------------------------
    int  IsOpenFlag;           // TRUE if a project has been opened
    int  IsStartedFlag;        // TRUE if a simulation has been started
    int  SaveResultsFlag;      // TRUE if output to be saved to binary file
    int  ExceptionCount;       // number of exceptions handled
    int  DoRunoff;             // TRUE if runoff is computed
    int  DoRouting;            // TRUE if flow routing is computed

    //-----------------------------------------------------------------------------
    //  Shared variables
    //-----------------------------------------------------------------------------
    HTtable* Htable[MAX_OBJ_TYPES]; // Hash tables for object ID names
    char MemPoolAllocated;      // TRUE if memory pool allocated
    alloc_root_t * root;


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from climate.c
    //-----------------------------------------------------------------------------
    // Temperature variables
    double    Tmin;                 // min. daily temperature (deg F)
    double    Tmax;                 // max. daily temperature (deg F)
    double    Trng;                 // 1/2 range of daily temperatures
    double    Trng1;                // prev. max - current min. temp.
    double    Tave;                 // average daily temperature (deg F)
    double    Hrsr;                 // time of min. temp. (hrs)
    double    Hrss;                 // time of max. temp (hrs)
    double    Hrday;                // avg. of min/max temp times
    double    Dhrdy;                // hrs. between min. & max. temp. times
    double    Dydif;                // hrs. between max. & min. temp. times
    DateTime  LastDay;              // date of last day with temp. data
    TMovAve   Tma;                  // moving average of daily temperatures //(5.1.010)

    // Evaporation variables
    DateTime  NextEvapDate;         // next date when evap. rate changes
    double    NextEvapRate;         // next evaporation rate (user units)

    // Climate file variables
    int      FileFormat;            // file format (see ClimateFileFormats)
    int      FileYear;              // current year of file data
    int      FileMonth;             // current month of year of file data
    int      FileDay;               // current day of month of file data
    int      FileLastDay;           // last day of current month of file data
    int      FileElapsedDays;       // number of days read from file
    double   FileValue[4];          // current day's values of climate data
    double   FileData[4][32];       // month's worth of daily climate data
    char     FileLine[MAXLINE + 1];   // line from climate data file

    int      FileFieldPos[4];       // start of data fields for file record //(5.1.007)
    int      FileDateFieldPos;      // start of date field for file record  //(5.1.007)
    int      FileWindType;          // wind speed type;                     //(5.1.007)

    //-----------------------------------------------------------------------------
    //  Shared variables
    //-----------------------------------------------------------------------------
    TRule*       Rules;           // array of control rules
    TActionList* ActionList;      // linked list of control actions
    int      InputState;                   // state of rule interpreter
    int      RuleCount;                    // total number of rules
    double   ControlValue;                 // value of controller variable
    double   SetPoint;                     // value of controller setpoint
    DateTime CurrentDate;                  // current date in whole days
    DateTime CurrentTime;                  // current time of day (decimal)
    DateTime ElapsedTime;                  // elasped simulation time (decimal days)


    //-----------------------------------------------------------------------------
    //  Shared Variables - moved from dynwave.c
    //-----------------------------------------------------------------------------
    double  VariableStep;           // size of variable time step (sec)
    TXnode* Xnode;                  // extended nodal information

    double  Omega;                  // actual under-relaxation parameter
    int     Steps;                  // number of Picard iterations


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from iface.c
    //-----------------------------------------------------------------------------
    int      IfaceFlowUnits;        // flow units for routing interface file
    int      IfaceStep;             // interface file time step (sec)
    int      NumIfacePolluts;       // number of pollutants in interface file
    int*     IfacePolluts;          // indexes of interface file pollutants
    int      NumIfaceNodes;         // number of nodes on interface file
    int*     IfaceNodes;            // indexes of nodes on interface file
    double** OldIfaceValues;        // interface flows & WQ at previous time
    double** NewIfaceValues;        // interface flows & WQ at next time
    double   IfaceFrac;             // fraction of interface file time step
    DateTime OldIfaceDate;          // previous date of interface values
    DateTime NewIfaceDate;          // next date of interface values


    //-----------------------------------------------------------------------------
    //  Local Variables - moved from infil.c
    //-----------------------------------------------------------------------------
    THorton*   HortInfil;
    TGrnAmpt*  GAInfil;
    TCurveNum* CNInfil;
    double Fumax;   // saturated water volume in upper soil zone (ft)

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from input.c
    //-----------------------------------------------------------------------------
    char *Tok[MAXTOKS];             // String tokens from line of input
    int  Ntokens;                   // Number of tokens in line of input
    int  Mobjects[MAX_OBJ_TYPES];   // Working number of objects of each type
    int  Mnodes[MAX_NODE_TYPES];    // Working number of node objects
    int  Mlinks[MAX_LINK_TYPES];    // Working number of link objects

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from input.c
    //-----------------------------------------------------------------------------
    double   Beta1;
    double   C1;
    double   C2;
    double   Afull;
    double   Qfull;
    TXsect*  pXsect;


    //-----------------------------------------------------------------------------
    // Globally shared variables   - moved from subcatch.c
    //-----------------------------------------------------------------------------
    // Volumes (ft3) for a subcatchment over a time step                           //(5.1.008)
    double     Vevap;         // evaporation
    double     Vpevap;        // pervious area evaporation
    double     Vinfil;        // non-LID infiltration
    double     Vinflow;       // non-LID precip + snowmelt + runon + ponded water
    double     Voutflow;      // non-LID runoff to subcatchment's outlet
    double     VlidIn;        // impervious area flow to LID units
    double     VlidInfil;     // infiltration from LID units
    double     VlidOut;       // surface outflow from LID units
    double     VlidDrain;     // drain outflow from LID units
    double     VlidReturn;    // LID outflow returned to pervious area
    char HasWetLids;
    TSubarea* theSubarea;     // subarea to which getDdDt() is applied

    //-----------------------------------------------------------------------------
    //  Shared Variables - moved from lid.c
    //-----------------------------------------------------------------------------
    TLidProc*  LidProcs;            // array of LID processes
    int        LidCount;            // number of LID processes
    TLidGroup* LidGroups;           // array of LID process groups
    int        GroupCount;          // number of LID groups (subcatchments)

    double     EvapRate;            // evaporation rate (ft/s)
    double     NativeInfil;         // native soil infil. rate (ft/s)
    double     MaxNativeInfil;      // native soil infil. rate limit (ft/s)

    //-----------------------------------------------------------------------------
    //  Shared variables   - moved from massbal.c
    //-----------------------------------------------------------------------------
    TRunoffTotals    RunoffTotals;    // overall surface runoff continuity totals
    TLoadingTotals*  LoadingTotals;   // overall WQ washoff continuity totals
    TGwaterTotals    GwaterTotals;    // overall groundwater continuity totals
    TRoutingTotals   FlowTotals;      // overall routed flow continuity totals
    TRoutingTotals*  QualTotals;      // overall routed WQ continuity totals
    TRoutingTotals   StepFlowTotals;  // routed flow totals over time step
    TRoutingTotals   OldStepFlowTotals;
    TRoutingTotals*  StepQualTotals;  // routed WQ totals over time step

    //-----------------------------------------------------------------------------
    //  Shared variables
    //-----------------------------------------------------------------------------

    TSysStats       SysStats;
    TMaxStats       MaxMassBalErrs[MAX_STATS];
    TMaxStats       MaxCourantCrit[MAX_STATS];
    TMaxStats       MaxFlowTurns[MAX_STATS];
    double          SysOutfallFlow;

    //-----------------------------------------------------------------------------
    //  Exportable variables (shared with statsrpt.c)
    //-----------------------------------------------------------------------------
    TSubcatchStats* SubcatchStats;
    TNodeStats*     NodeStats;
    TLinkStats*     LinkStats;
    TStorageStats*  StorageStats;
    TOutfallStats*  OutfallStats;
    TPumpStats*     PumpStats;
    double          MaxOutfallFlow;
    double          MaxRunoffFlow;

    //-----------------------------------------------------------------------------
    //  Imported variables
    //-----------------------------------------------------------------------------
    double*         NodeInflow;     // defined in massbal.c
    double*         NodeOutflow;    // defined in massbal.c


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from gwater.c
    //-----------------------------------------------------------------------------
    //  NOTE: all flux rates are in ft/sec, all depths are in ft.
    double    Area;            // subcatchment area (ft2)                   //(5.1.008)
    double    Infil;           // infiltration rate from surface
    double    MaxEvap;         // max. evaporation rate
    double    AvailEvap;       // available evaporation rate
    double    UpperEvap;       // evaporation rate from upper GW zone
    double    LowerEvap;       // evaporation rate from lower GW zone
    double    UpperPerc;       // percolation rate from upper to lower zone
    double    LowerLoss;       // loss rate from lower GW zone
    double    GWFlow;          // flow rate from lower zone to conveyance node
    double    MaxUpperPerc;    // upper limit on UpperPerc
    double    MaxGWFlowPos;    // upper limit on GWFlow when its positve
    double    MaxGWFlowNeg;    // upper limit on GWFlow when its negative
    double    FracPerv;        // fraction of surface that is pervious
    double    TotalDepth;      // total depth of GW aquifer
    double    Theta;           // moisture content of upper zone
    double    HydCon;          // unsaturated hydraulic conductivity (ft/s) //(5.1.010)
    double    Hgw;             // ht. of saturated zone
    double    Hstar;           // ht. from aquifer bottom to node invert
    double    Hsw;             // ht. from aquifer bottom to water surface
    double    Tstep;           // current time step (sec)
    TAquifer  A;               // aquifer being analyzed
    TGroundwater* GW;          // groundwater object being analyzed
    MathExpr* LatFlowExpr;     // user-supplied lateral GW flow expression  //(5.1.007)
    MathExpr* DeepFlowExpr;    // user-supplied deep GW flow expression     //(5.1.007)

    //-----------------------------------------------------------------------------
    //  Imported variables - moved from statsrpt.c
    //-----------------------------------------------------------------------------
    //double          MaxOutfallFlow;
    //double          MaxRunoffFlow;
    //double*         NodeInflow;             // defined in MASSBAL.C
    //double*         NodeOutflow;            // defined in massbal.c
    double TotalArea;

    //-----------------------------------------------------------------------------
    //  Imported variables - moved from report.c
    //-----------------------------------------------------------------------------
    REAL4* SubcatchResults;         // Results vectors defined in OUTPUT.C
    REAL4* NodeResults;             //  "
    REAL4* LinkResults;             //  "
    //		char   ErrString[81];           // defined in ERROR.C

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from report.c
    //-----------------------------------------------------------------------------
    time_t SysTime;

    //-----------------------------------------------------------------------------
    //  Shared variables
    //-----------------------------------------------------------------------------
    INT4      IDStartPos;           // starting file position of ID names
    INT4      InputStartPos;        // starting file position of input data
    INT4      OutputStartPos;       // starting file position of output data
    INT4      BytesPerPeriod;       // bytes saved per simulation time period
    INT4      NsubcatchResults;     // number of subcatchment output variables
    INT4      NnodeResults;         // number of node output variables
    INT4      NlinkResults;         // number of link output variables
    INT4      NumSubcatch;          // number of subcatchments reported on
    INT4      NumNodes;             // number of nodes reported on
    INT4      NumLinks;             // number of links reported on
    INT4      NumPolluts;           // number of pollutants reported on
    REAL4     SysResults[MAX_SYS_RESULTS];    // values of system output vars.


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from rain.c
    //-----------------------------------------------------------------------------
    TRainStats RainStats;                  // see objects.h for definition
    int        Condition;                  // rainfall condition code
    int        TimeOffset;                 // time offset of rainfall reading (sec)
    int        DataOffset;                 // start of data on line of input
    int        ValueOffset;                // start of rain value on input line
    int        RainType;                   // rain measurement type code
    int        Interval;                   // rain measurement interval (sec)
    double     UnitsFactor;                // units conversion factor
    float      RainAccum;                  // rainfall depth accumulation
    char       *StationID;                 // station ID appearing in rain file
    DateTime   AccumStartDate;             // date when accumulation begins
    DateTime   PreviousDate;               // date of previous rainfall record
    int        GageIndex;                  // index of rain gage analyzed
    int        hasStationName;             // true if data contains station name


    //-----------------------------------------------------------------------------
    // Imported Variables - moved from rdii.c
    //-----------------------------------------------------------------------------
    //double     Qcf[];               // flow units conversion factors
    // (see swmm5.c)

    //-----------------------------------------------------------------------------
    // Shared Variables - moved from rdii.c
    //-----------------------------------------------------------------------------
    TUHGroup*  UHGroup;             // processing data for each UH group
    int        RdiiStep;            // RDII time step (sec)
    int        NumRdiiNodes;        // number of nodes w/ RDII data
    int*       RdiiNodeIndex;       // indexes of nodes w/ RDII data
    REAL4*     RdiiNodeFlow;        // inflows for nodes with RDII          //(5.1.003)
    int        RdiiFlowUnits;       // RDII flow units code
    DateTime   RdiiStartDate;       // start date of RDII inflow period
    DateTime   RdiiEndDate;         // end date of RDII inflow period
    double     TotalRainVol;        // total rainfall volume (ft3)
    double     TotalRdiiVol;        // total RDII volume (ft3)
    int        RdiiFileType;        // type (binary/text) of RDII file

    //-----------------------------------------------------------------------------
    // Shared variables
    //-----------------------------------------------------------------------------
    char  IsRaining;                // TRUE if precip. falls on study area
    char  HasRunoff;                // TRUE if study area generates runoff
    char  HasSnow;                  // TRUE if any snow cover on study area
    int   Nsteps;                   // number of runoff time steps taken
    int   MaxSteps;                 // final number of runoff time steps
    long  MaxStepsPos;              // position in Runoff interface file
    //    where MaxSteps is saved

    //-----------------------------------------------------------------------------
    //  Exportable variables
    //-----------------------------------------------------------------------------
    //char    HasWetLids;  // TRUE if any LIDs are wet (used in lidproc.c)           //(5.1.008)
    double* OutflowLoad; // exported pollutant mass load (used in surfqual.c)
    TLidUnit*  theLidUnit;     // ptr. to a subcatchment's LID unit
    TLidProc*  theLidProc;     // ptr. to a LID process

    double     Tsteplp;          // current time step (sec)
    //static double     Rainfall;       // current rainfall rate (ft/s)            //(5.1.008)
    double     EvapRatelp;       // evaporation rate (ft/s)
    double     MaxNativeInfillp; // native soil infil. rate limit (ft/s)

    double     SurfaceInflow;  // precip. + runon to LID unit (ft/s)
    double     SurfaceInfil;   // infil. rate from surface layer (ft/s)
    double     SurfaceEvap;    // evap. rate from surface layer (ft/s)
    double     SurfaceOutflow; // outflow from surface layer (ft/s)
    double     SurfaceVolume;  // volume in surface storage (ft)

    double     PaveEvap;       // evap. from pavement layer (ft/s)          //(5.1.008)
    double     PavePerc;       // percolation from pavement layer (ft/s)    //(5.1.008)
    double     PaveVolume;     // volume stored in pavement layer  (ft)     //(5.1.008)

    double     SoilEvap;       // evap. from soil layer (ft/s)
    double     SoilPerc;       // percolation from soil layer (ft/s)
    double     SoilVolume;     // volume in soil/pavement storage (ft)

    double     StorageInflow;  // inflow rate to storage layer (ft/s)
    double     StorageInfil;   // infil. rate from storage layer (ft/s)
    double     StorageEvap;    // evap.rate from storage layer (ft/s)
    double     StorageDrain;   // underdrain flow rate layer (ft/s)
    double     StorageVolume;  // volume in storage layer (ft)

    double     Xold[MAX_LAYERS];  // previous moisture level in LID layers  //(5.1.008)

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from toposort.c
    //-----------------------------------------------------------------------------
    int* InDegree;                  // number of incoming links to each node
    int* StartPos;                  // start of a node's outlinks in AdjList
    int* AdjList;                   // list of outlink indexes for each node
    int* Stack;                     // array of nodes "reached" during sorting
    int  First;                     // position of first node in stack
    int  Last;                      // position of last node added to stack

    char* Examined;                 // TRUE if node included in spanning tree
    char* InTree;                   // state of each link in spanning tree:
    // 0 = unexamined,
    // 1 = in spanning tree,
    // 2 = chord of spanning tree
    int*  LoopLinks;                // list of links which forms a loop
    int   LoopLinksLast;            // number of links in a loop

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from tansect.c
    //-----------------------------------------------------------------------------
    int    Ntransects;              // total number of transects
    int    Nstations;               // number of stations in current transect
    double  Station[MAXSTATION+1];  // x-coordinate of each station
    double  Elev[MAXSTATION+1];     // elevation of each station
    double  Nleft;                  // Manning's n for left overbank
    double  Nright;                 // Manning's n for right overbank
    double  Nchannel;               // Manning's n for main channel
    double  Xleftbank;              // station where left overbank ends
    double  Xrightbank;             // station where right overbank begins
    double  Xfactor;                // multiplier for station spacing
    double  Yfactor;                // factor added to station elevations
    double  Lfactor;                // main channel/flood plain length


    //-----------------------------------------------------------------------------
    //  Shared variables - moved from treatment.c
    //-----------------------------------------------------------------------------
    int     ErrCode;                // treatment error code
    int     J;                      // index of node being analyzed
    double  Dt;                     // curent time step (sec)
    double  Q;                      // node inflow (cfs)
    double  V;                      // node volume (ft3)
    double* R;                      // array of pollut. removals
    double* Cin;                    // node inflow concentrations

    char  ErrString[256];

    //-----------------------------------------------------------------------------
    //  Shared variables - moved from mathexpr.c
    //-----------------------------------------------------------------------------
     int    Err;
     int    Bc;
     int    PrevLex, CurLex;
     int    Len, Pos;
     char   *S;
     char   Token[255];
     int    Ivar;
     double Fvalue;


     //-----------------------------------------------------------------------------
     //  Shared variables - moved from odesolve.c
     //-----------------------------------------------------------------------------
     int      nmax;      // max. number of equations
     double*  y;         // dependent variable
     double*  yscal;     // scaling factors
     double*  yerr;      // integration errors
     double*  ytemp;     // temporary values of y
     double*  dydx;      // derivatives of y
     double*  ak;        // derivatives at intermediate points

     //-----------------------------------------------------------------------------
     //  Shared variables - moved from routing.c
     //-----------------------------------------------------------------------------
     int* SortedLinks;


     //-----------------------------------------------------------------------------
     //  Shared variables - moved from shape.c
     //-----------------------------------------------------------------------------
     double Atotal;
     double Ptotal;
};

//typedef struct Project Project;

#ifdef __cplusplus
}
#endif

#endif
