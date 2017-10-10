//-----------------------------------------------------------------------------
//   climate.c
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/20/10 (Build 5.1.001)
//            09/15/14 (Build 5.1.007)
//            03/19/15 (Build 5.1.008)
//            08/05/15 (Build 5.1.010)
//   Author:  L. Rossman
//
//   Climate related functions.
//
//   Build 5.1.007:
//   - NCDC GHCN climate file format added.
//   - Monthly adjustments for temperature, evaporation & rainfall added.
//
//   Build 5.1.008:
//   - Monthly adjustments for hyd. conductivity added.
//   - Time series evaporation rates can now vary within a day.
//   - Evaporation rates are now properly updated when only flow routing
//     is being simulated.
//
//   Build 5.1.010:
//   - Hargreaves evaporation now computed using 7-day average temperatures.
//             
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
enum ClimateFileFormats {UNKNOWN_FORMAT,
                         USER_PREPARED,     // SWMM 5's own user format
                         GHCND,             // NCDC GHCN Daily format          //(5.1.007)
                         TD3200,            // NCDC TD3200 format
                         DLY0204};          // Canadian DLY02 or DLY04 format
static const int    MAXCLIMATEVARS  = 4;
static const int    MAXDAYSPERMONTH = 32;

// These variables are used when processing climate files.
enum   ClimateVarType {TMIN, TMAX, EVAP, WIND};
enum   WindSpeedType  {WDMV, AWND};                                            //(5.1.007)
static char* ClimateVarWords[] = {"TMIN", "TMAX", "EVAP", "WDMV", "AWND",      //(5.1.007)
                                  NULL};

////  Added for release 5.1.010.  ////                                         //(5.1.010)
////-----------------------------------------------------------------------------
////  Data Structures
////-----------------------------------------------------------------------------
//typedef struct
//{
//    double    tAve;          // moving avg. for daily temperature (deg F)
//    double    tRng;          // moving avg. for daily temp. range (deg F)
//    double    ta[7];         // data window for tAve
//    double    tr[7];         // data window for tRng
//    int       count;         // length of moving average window
//    int       maxCount;      // maximum length of moving average window
//    int       front;         // index of front of moving average window
//} TMovAve;
////

////-----------------------------------------------------------------------------
////  Shared variables
////-----------------------------------------------------------------------------
//// Temperature variables
//static double    project->Tmin;                 // min. daily temperature (deg F)
//static double    project->Tmax;                 // max. daily temperature (deg F)
//static double    project->Trng;                 // 1/2 range of daily temperatures
//static double    project->Trng1;                // prev. max - current min. temp.
//static double    project->Tave;                 // average daily temperature (deg F)
//static double    project->Hrsr;                 // time of min. temp. (hrs)
//static double    project->Hrss;                 // time of max. temp (hrs)
//static double    project->Hrday;                // avg. of min/max temp times
//static double    project->Dhrdy;                // hrs. between min. & max. temp. times
//static double    project->Dydif;                // hrs. between max. & min. temp. times
//static DateTime  project->LastDay;              // date of last day with temp. data
//static TMovAve   Tma;                  // moving average of daily temperatures //(5.1.010)
//
//// Evaporation variables
//static DateTime  project->NextEvapDate;         // next date when evap. rate changes
//static double    project->NextEvapRate;         // next evaporation rate (user units)
//
//// Climate file variables
//static int      project->FileFormat;            // file format (see ClimateFileFormats)
//static int      project->FileYear;              // current year of file data
//static int      project->FileMonth;             // current month of year of file data
//static int      project->FileDay;               // current day of month of file data
//static int      project->FileLastDay;           // last day of current month of file data
//static int      project->FileElapsedDays;       // number of days read from file
//static double   project->FileValue[4];          // current day's values of climate data
//static double   project->FileData[4][32];       // month's worth of daily climate data
//static char     project->FileLine[MAXLINE+1];   // line from climate data file
//
//static int      project->FileFieldPos[4];       // start of data fields for file record //(5.1.007)
//static int      project->FileDateFieldPos;      // start of date field for file record  //(5.1.007)
//static int      project->FileWindType;          // wind speed type;                     //(5.1.007)

//-----------------------------------------------------------------------------
//  External functions (defined in funcs.h)
//-----------------------------------------------------------------------------
//  climate_readParams                 // called by input_parseLine
//  climate_readEvapParams             // called by input_parseLine
//  climate_validate                   // called by project_validate
//  climate_openFile                   // called by runoff_open
//  climate_initState                  // called by project_init
//  climate_setState                   // called by runoff_execute
//  climate_getNextEvapDate            // called by runoff_getTimeStep         //(5.1.008)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int  getFileFormat(Project* project);
static void readFileLine(Project* project, int *year, int *month);
static void readUserFileLine(Project* project, int *year, int *month);
static void readTD3200FileLine(Project* project, int *year, int *month);
static void readDLY0204FileLine(Project* project, int *year, int *month);
static void readFileValues(Project* project);

static void setNextEvapDate(Project* project, DateTime thedate);                                 //(5.1.008)
static void setEvap(Project* project, DateTime theDate);
static void setTemp(Project* project, DateTime theDate);
static void setWind(Project* project, DateTime theDate);
static void updateTempTimes(Project* project, int day);
static void updateTempMoveAve(Project* project, double tmin, double tmax);                       //(5.1.010)
static double getTempEvap(Project* project, int day, double ta, double tr);                      //(5.1.010)

static void updateFileValues(Project* project, DateTime theDate);
static void parseUserFileLine(Project* project);
static void parseTD3200FileLine(Project* project);
static void parseDLY0204FileLine(Project* project);
static void setTD3200FileValues(Project* project,int param);

static int  isGhcndFormat(Project* project, char* line);                                         //(5.1.007)
static void readGhcndFileLine(Project* project, int *year, int *month);                          //(5.1.007)
static void parseGhcndFileLine(Project* project);                                          //(5.1.007)

//=============================================================================

int  climate_readParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads climate/temperature parameters from input line of data
//
//  Format of data can be
//    TIMESERIES  name
//    FILE        name
//    WINDSPEED   MONTHLY  v1  v2  ...  v12
//    WINDSPEED   FILE
//    SNOWMELT    v1  v2  ...  v6
//    ADC         IMPERV/PERV  v1  v2  ...  v10
//
{
    int      i, j, k;
    double   x[6], y;
    DateTime aDate;

    // --- identify keyword
    k = findmatch(tok[0], TempKeyWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);
    switch (k)
    {
      case 0: // Time series name
        // --- check that time series name exists
        if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
        i = project_findObject(project,TSERIES, tok[1]);
        if ( i < 0 ) return error_setInpError(ERR_NAME, tok[1]);

        // --- record the time series as being the data source for temperature
        project->Temp.dataSource = TSERIES_TEMP;
        project->Temp.tSeries = i;
        project->Tseries[i].refersTo = TSERIES_TEMP;
        break;

      case 1: // Climate file
        // --- record file as being source of temperature data
        if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
        project->Temp.dataSource = FILE_TEMP;

        // --- save name and usage mode of external climate file
        project->Fclimate.mode = USE_FILE;
        sstrncpy(project->Fclimate.name, tok[1], MAXFNAME);

        // --- save starting date to read from file if one is provided
        project->Temp.fileStartDate = NO_DATE;
        if ( ntoks > 2 )
        {
            if ( *tok[2] != '*')
            {
                if ( !datetime_strToDate(tok[2], &aDate) )
                    return error_setInpError(ERR_DATETIME, tok[2]);
                project->Temp.fileStartDate = aDate;
            }
        }
        break;

      case 2: // project->Wind speeds
        // --- check if wind speeds will be supplied from climate file
        if ( strcomp(tok[1], w_FILE) )
        {
            project->Wind.type = FILE_WIND;
        }

        // --- otherwise read 12 monthly avg. wind speed values
        else
        {
            if ( ntoks < 14 ) return error_setInpError(ERR_ITEMS, "");
            project->Wind.type = MONTHLY_WIND;
            for (i=0; i<12; i++)
            {
                if ( !getDouble(tok[i+2], &y) )
                    return error_setInpError(ERR_NUMBER, tok[i+2]);
                project->Wind.aws[i] = y;
            }
        }
        break;

      case 3: // Snowmelt params
        if ( ntoks < 7 ) return error_setInpError(ERR_ITEMS, "");
        for (i=1; i<7; i++)
        {
            if ( !getDouble(tok[i], &x[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        // --- convert deg. C to deg. F for snowfall temperature
        if ( project->UnitSystem == SI ) x[0] = 9./5.*x[0] + 32.0;
        project->Snow.snotmp = x[0];
        project->Snow.tipm   = x[1];
        project->Snow.rnm    = x[2];
		project->Temp.elev = x[3] / UCF(project, LENGTH);
        project->Temp.anglat = x[4];
        project->Temp.dtlong = x[5] / 60.0;
        break;

      case 4:  // Areal Depletion Curve data
        // --- check if data is for impervious or pervious areas
        if ( ntoks < 12 ) return error_setInpError(ERR_ITEMS, "");
        if      ( match(tok[1], w_IMPERV) ) i = 0;
        else if ( match(tok[1], w_PERV)   ) i = 1;
        else return error_setInpError(ERR_KEYWORD, tok[1]);

        // --- read 10 fractional values
        for (j=0; j<10; j++)
        {
            if ( !getDouble(tok[j+2], &y) || y < 0.0 || y > 1.0 )
                return error_setInpError(ERR_NUMBER, tok[j+2]);
            project->Snow.adc[i][j] = y;
        }
        break;
    }
    return 0;
}

//=============================================================================

int climate_readEvapParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads evaporation parameters from input line of data.
//
//  Data formats are:
//    CONSTANT  value
//    MONTHLY   v1 ... v12
//    TIMESERIES name
//    TEMPERATURE
//    FILE      (v1 ... v12)
//    RECOVERY   name
//    DRY_ONLY   YES/NO
//
{
    int i, k;
    double x;

    // --- find keyword indicating what form the evaporation data is in
    k = findmatch(tok[0], EvapTypeWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);

    // --- check for RECOVERY pattern data
    if ( k == RECOVERY )
    {
        if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
        i = project_findObject(project,TIMEPATTERN, tok[1]);
        if ( i < 0 ) return error_setInpError(ERR_NAME, tok[1]);
        project->Evap.recoveryPattern = i;
        return 0;
    }

    // --- check for no evaporation in wet periods
    if ( k == DRYONLY )
    {
        if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
        if      ( strcomp(tok[1], w_NO ) )  project->Evap.dryOnly = FALSE;
        else if ( strcomp(tok[1], w_YES ) ) project->Evap.dryOnly = TRUE;
        else return error_setInpError(ERR_KEYWORD, tok[1]);
        return 0;
    }

    // --- process data depending on its form
    project->Evap.type = k;
    if ( k != TEMPERATURE_EVAP && ntoks < 2 )
        return error_setInpError(ERR_ITEMS, "");
    switch ( k )
    {
      case CONSTANT_EVAP:
        // --- for constant evap., fill monthly avg. values with same number
        if ( !getDouble(tok[1], &x) )
            return error_setInpError(ERR_NUMBER, tok[1]);
        for (i=0; i<12; i++) project->Evap.monthlyEvap[i] = x;
        break;

      case MONTHLY_EVAP:
        // --- for monthly evap., read a value for each month of year
        if ( ntoks < 13 ) return error_setInpError(ERR_ITEMS, "");
        for ( i=0; i<12; i++)
            if ( !getDouble(tok[i+1], &project->Evap.monthlyEvap[i]) )
                return error_setInpError(ERR_NUMBER, tok[i+1]);
        break;

      case TIMESERIES_EVAP:
        // --- for time series evap., read name of time series
		  i = project_findObject(project, TSERIES, tok[1]);
        if ( i < 0 ) return error_setInpError(ERR_NAME, tok[1]);
        project->Evap.tSeries = i;
        project->Tseries[i].refersTo = TIMESERIES_EVAP;
        break;

      case FILE_EVAP:
        // --- for evap. from climate file, read monthly pan coeffs.
        //     if they are provided (default values are 1.0)
        if ( ntoks > 1 )
        {
            if ( ntoks < 13 ) return error_setInpError(ERR_ITEMS, "");
            for (i=0; i<12; i++)
            {
                if ( !getDouble(tok[i+1], &project->Evap.panCoeff[i]) )
                    return error_setInpError(ERR_NUMBER, tok[i+1]);
            }
        }
        break;
    }
    return 0;
}

//=============================================================================

////  New function added to release 5.1.007.  ////                             //(5.1.007)

int climate_readAdjustments(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads adjustments to monthly evaporation or rainfall
//           from input line of data.
//
//  Data formats are:
//    TEMPERATURE   v1 ... v12
//    EVAPORATION   v1 ... v12
//    RAINFALL      v1 ... v12
//    CONDUCTIVITY  v1 ... v12                                                 //(5.1.008)
{
    int i;
    if (ntoks == 1) return 0;

    if ( match(tok[0], "TEMP") )
    {
        if ( ntoks < 13 )  return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &project->Adjust.temp[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        return 0;
    }

    if ( match(tok[0], "EVAP") )
    {
        if ( ntoks < 13 )  return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &project->Adjust.evap[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        return 0;
    }

    if ( match(tok[0], "RAIN") )
    {
        if ( ntoks < 13 )  return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &project->Adjust.rain[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        return 0;
    }

////  Following code segment added for release 5.1.008.  ////                  //(5.1.008)
////
    if ( match(tok[0], "CONDUCT") )
    {
        if ( ntoks < 13 )  return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &project->Adjust.hydcon[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        return 0;
    }
////
    return error_setInpError(ERR_KEYWORD, tok[0]);
}

//=============================================================================

void climate_validate(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: validates climatological variables
//
{
    int       i;                                                               //(5.1.007)
    double    a, z, pa;

    // --- check if climate data comes from external data file                 //(5.1.007)
    if ( project->Wind.type == FILE_WIND || project->Evap.type == FILE_EVAP ||
         project->Evap.type == TEMPERATURE_EVAP )
    {
        if ( project->Fclimate.mode == NO_FILE )
        {
            report_writeErrorMsg(project,ERR_NO_CLIMATE_FILE, "");
        }
    }

    // --- open the climate data file                                          //(5.1.007)
    if ( project->Fclimate.mode == USE_FILE ) climate_openFile(project);                       //(5.1.007)

    // --- snow melt parameters tipm & rnm must be fractions
    if ( project->Snow.tipm < 0.0 ||
         project->Snow.tipm > 1.0 ||
         project->Snow.rnm  < 0.0 ||
         project->Snow.rnm  > 1.0 ) report_writeErrorMsg(project,ERR_SNOWMELT_PARAMS, "");

    // --- latitude should be between -90 & 90 degrees
    a = project->Temp.anglat;
    if ( a <= -89.99 ||
         a >= 89.99  ) report_writeErrorMsg(project,ERR_SNOWMELT_PARAMS, "");
    else project->Temp.tanAnglat = tan(a * PI / 180.0);

    // --- compute psychrometric constant
    z = project->Temp.elev / 1000.0;
    if ( z <= 0.0 ) pa = 29.9;
    else  pa = 29.9 - 1.02*z + 0.0032*pow(z, 2.4); // atmos. pressure
    project->Temp.gamma = 0.000359 * pa;

    // --- convert units of monthly temperature & evap adjustments             //(5.1.007)
    for (i = 0; i < 12; i++)
    {
        if (project->UnitSystem == SI) project->Adjust.temp[i] *= 9.0/5.0;
        project->Adjust.evap[i] /= UCF(project,EVAPRATE);
    }
}

//=============================================================================

void climate_openFile(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a climate file and reads in first set of values.
//
{
    int i, m, y;

    // --- open the file
    if ( (project->Fclimate.file = fopen(project->Fclimate.name, "rt")) == NULL )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_OPEN, project->Fclimate.name);
        return;
    }

    // --- initialize values of file's climate variables
    //     (project->Temp.ta was previously initialized in project.c)
    project->FileValue[TMIN] = project->Temp.ta;
    project->FileValue[TMAX] = project->Temp.ta;
    project->FileValue[EVAP] = 0.0;
    project->FileValue[WIND] = 0.0;

    // --- find climate file's format
    project->FileFormat = getFileFormat(project);
    if ( project->FileFormat == UNKNOWN_FORMAT )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_READ, project->Fclimate.name);
        return;
    }

    // --- position file to begin reading climate file at either user-specified
    //     month/year or at start of simulation period.
    rewind(project->Fclimate.file);
    strcpy(project->FileLine, "");
    if ( project->Temp.fileStartDate == NO_DATE )
        datetime_decodeDate(project->StartDate, &project->FileYear, &project->FileMonth, &project->FileDay);
    else
        datetime_decodeDate(project->Temp.fileStartDate, &project->FileYear, &project->FileMonth, &project->FileDay);
    while ( !feof(project->Fclimate.file) )
    {
        strcpy(project->FileLine, "");
        readFileLine(project,&y, &m);
        if ( y == project->FileYear && m == project->FileMonth ) break;
    }
    if ( feof(project->Fclimate.file) )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_END_OF_FILE, project->Fclimate.name);
        return;
    }

    // --- initialize file dates and current climate variable values
    if ( !project->ErrorCode )
    {
        project->FileElapsedDays = 0;
        project->FileLastDay = datetime_daysPerMonth(project->FileYear, project->FileMonth);
        readFileValues(project);
        for (i=TMIN; i<=WIND; i++)
        {
            if ( project->FileData[i][project->FileDay] == MISSING ) continue;
            project->FileValue[i] = project->FileData[i][project->FileDay];
        }
    }
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void climate_initState(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes climate state variables.
//
{
    project->LastDay = NO_DATE;
    project->Temp.tmax = MISSING;
    project->Snow.removed = 0.0;
    project->NextEvapDate = project->StartDate;
    project->NextEvapRate = 0.0;

    // --- initialize variables for time series evaporation
    if ( project->Evap.type == TIMESERIES_EVAP && project->Evap.tSeries >= 0  )
    {
        // --- initialize project->NextEvapDate & project->NextEvapRate to first entry of
        //     time series whose date <= the simulation start date
        table_getFirstEntry(&project->Tseries[project->Evap.tSeries],
                            &project->NextEvapDate, &project->NextEvapRate);
        if ( project->NextEvapDate < project->StartDate )
        {  
            setNextEvapDate(project,project->StartDate);
        }
        project->Evap.rate = project->NextEvapRate / UCF(project,EVAPRATE);

        // --- find the next time evaporation rates change after this
		setNextEvapDate(project, project->NextEvapDate);
    }

////  Following section added to release 5.1.010.  ////                        //(5.1.010)
    // --- initialize variables for temperature evaporation
    if ( project->Evap.type == TEMPERATURE_EVAP )
    {
        project->Tma.maxCount = sizeof(project->Tma.ta) / sizeof(double);
        project->Tma.count = 0;
        project->Tma.front = 0;
        project->Tma.tAve = 0.0;
        project->Tma.tRng = 0.0;
    }
////
}

//=============================================================================

void climate_setState(Project* project, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets climate variables for current date.
//
{
	if (project->Fclimate.mode == USE_FILE) updateFileValues(project, theDate);
	if (project->Temp.dataSource != NO_TEMP) setTemp(project, theDate);
	setEvap(project, theDate);
	setWind(project, theDate);
    project->Adjust.rainFactor = project->Adjust.rain[datetime_monthOfYear(theDate)-1];          //(5.1.007)
    project->Adjust.hydconFactor = project->Adjust.hydcon[datetime_monthOfYear(theDate)-1];      //(5.1.008)
	setNextEvapDate(project, theDate);                                                  //(5.1.008)
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

DateTime climate_getNextEvapDate(Project* project)
//
//  Input:   none
//  Output:  returns the current value of project->NextEvapDate
//  Purpose: gets the next date when evaporation rate changes.
//
{
    return project->NextEvapDate;
}

//=============================================================================

////  Modified from what was previously named climate_getNextEvap.  ////       //(5.1.008)

void setNextEvapDate(Project* project, DateTime theDate)
//
//  Input:   theDate = current simulation date
//  Output:  sets a new value for project->NextEvapDate
//  Purpose: finds date for next change in evaporation after the current date.
//
{
    int    yr, mon, day, k;
    double d, e;

    // --- do nothing if current date hasn't reached the current next date
    if ( project->NextEvapDate > theDate ) return;

    switch ( project->Evap.type )
    {
      // --- for constant evaporation, use a next date far in the future
      case CONSTANT_EVAP:
         project->NextEvapDate = theDate + 365.;
         break;

      // --- for monthly evaporation, use the start of the next month
      case MONTHLY_EVAP:
        datetime_decodeDate(theDate, &yr, &mon, &day);
        if ( mon == 12 )
        {
            mon = 1;
            yr++;
        }
        else mon++;
        project->NextEvapDate = datetime_encodeDate(yr, mon, 1);
        break;

      // --- for time series evaporation, find the next entry in the
      //     series on or after the current date
      case TIMESERIES_EVAP:
        k = project->Evap.tSeries;
        if ( k >= 0 )
        {
            project->NextEvapDate = theDate + 365.;
            while ( table_getNextEntry(&project->Tseries[k], &d, &e) &&
                    d <= project->EndDateTime )
            {
                if ( d >= theDate )
                {
                    project->NextEvapDate = d;
                    project->NextEvapRate = e;
                    break;
                }
            }
        }
        break;

      // --- for climate file daily evaporation, use the next day
      case FILE_EVAP:
        project->NextEvapDate = floor(theDate) + 1.0;
        break;

      default: project->NextEvapDate = theDate + 365.;
    }
}

//=============================================================================

void updateFileValues(Project* project, DateTime theDate)
//
//  Input:   theDate = current simulation date
//  Output:  none
//  Purpose: updates daily climate variables for new day or reads in
//           another month worth of values if a new month begins.
//
//  NOTE:    counters project->FileElapsedDays, project->FileDay, project->FileMonth, project->FileYear and
//           project->FileLastDay were initialized in climate_openFile().
//
{
    int i;
    int deltaDays;

    // --- see if a new day has begun
    deltaDays = (int)(floor(theDate) - floor(project->StartDateTime));
    if ( deltaDays > project->FileElapsedDays )
    {
        // --- advance day counters
        project->FileElapsedDays++;
        project->FileDay++;

        // --- see if new month of data needs to be read from file
        if ( project->FileDay > project->FileLastDay )
        {
            project->FileMonth++;
            if ( project->FileMonth > 12 )
            {
                project->FileMonth = 1;
                project->FileYear++;
            }
			readFileValues(project);
            project->FileDay = 1;
            project->FileLastDay = datetime_daysPerMonth(project->FileYear, project->FileMonth);
        }

        // --- set climate variables for new day
        for (i=TMIN; i<=WIND; i++)
        {
            // --- no change in current value if its missing
            if ( project->FileData[i][project->FileDay] == MISSING ) continue;
            project->FileValue[i] = project->FileData[i][project->FileDay];
        }
    }
}

//=============================================================================

void setTemp(Project* project, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: updates temperatures for new simulation date.
//
{
    int      j;                        // snow data object index
    int      k;                        // time series index
    int      mon;                      // month of year                        //(5.1.007)
    int      day;                      // day of year
    DateTime theDay;                   // calendar day
    double   hour;                     // hour of day
    double   tmp;                      // temporary temperature

    // --- see if a new day has started
    mon = datetime_monthOfYear(theDate);                                       //(5.1.007)
    theDay = floor(theDate);
    if ( theDay > project->LastDay )
    {
        // --- update min. & max. temps & their time of day
        day = datetime_dayOfYear(theDate);
        if ( project->Temp.dataSource == FILE_TEMP )
        {
            project->Tmin = project->FileValue[TMIN] + project->Adjust.temp[mon-1];                       //(5.1.007)
            project->Tmax = project->FileValue[TMAX] + project->Adjust.temp[mon-1];                       //(5.1.007)
            if ( project->Tmin > project->Tmax )
            {
                tmp = project->Tmin;
                project->Tmin = project->Tmax;
                project->Tmax = tmp;
            }
            updateTempTimes(project, day);
            if ( project->Evap.type == TEMPERATURE_EVAP )
            {
				updateTempMoveAve(project, project->Tmin, project->Tmax);                                 //(5.1.010)
				project->FileValue[EVAP] = getTempEvap(project, day, project->Tma.tAve, project->Tma.tRng);        //(5.1.010)
            }
        }

        // --- compute snow melt coefficients based on day of year
        project->Snow.season = sin(0.0172615*(day-81.0));
        for (j=0; j<project->Nobjects[SNOWMELT]; j++)
        {
            snow_setMeltCoeffs(project, j, project->Snow.season);
        }

        // --- update date of last day analyzed
        project->LastDay = theDate;
    }

    // --- for min/max daily temps. from climate file,
    //     compute hourly temp. by sinusoidal interp.
    if ( project->Temp.dataSource == FILE_TEMP )
    {
        hour = (theDate - theDay) * 24.0;
        if ( hour < project->Hrsr )
            project->Temp.ta = project->Tmin + project->Trng1/2.0 * sin(PI/project->Dydif * (project->Hrsr - hour));
        else if ( hour >= project->Hrsr && hour <= project->Hrss )
            project->Temp.ta = project->Tave + project->Trng * sin(PI/project->Dhrdy * (project->Hrday - hour));
        else
            project->Temp.ta = project->Tmax - project->Trng * sin(PI/project->Dydif * (hour - project->Hrss));
    }

    // --- for user-supplied temperature time series,
    //     get temperature value from time series
    if ( project->Temp.dataSource == TSERIES_TEMP )
    {
        k = project->Temp.tSeries;
        if ( k >= 0)
        {
            project->Temp.ta = table_tseriesLookup(&project->Tseries[k], theDate, TRUE);

            // --- convert from deg. C to deg. F if need be
            if ( project->UnitSystem == SI )
            {
                project->Temp.ta = (9./5.) * project->Temp.ta + 32.0;
            }

            // --- apply climate change adjustment factor                      //(5.1.007)
            project->Temp.ta += project->Adjust.temp[mon-1];                                     //(5.1.007)
        }
    }

    // --- compute saturation vapor pressure
    project->Temp.ea = 8.1175e6 * exp(-7701.544 / (project->Temp.ta + 405.0265) );
}

//=============================================================================

void setEvap(Project* project, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets evaporation rate (ft/sec) for a specified date.
//
{
    int k;
    int mon = datetime_monthOfYear(theDate);                                   //(5.1.007)

    switch ( project->Evap.type )
    {
      case CONSTANT_EVAP:
		  project->Evap.rate = project->Evap.monthlyEvap[0] / UCF(project, EVAPRATE);
        break;

      case MONTHLY_EVAP:
		  project->Evap.rate = project->Evap.monthlyEvap[mon - 1] / UCF(project, EVAPRATE);
        break;

      case TIMESERIES_EVAP:
        if ( theDate >= project->NextEvapDate )
			project->Evap.rate = project->NextEvapRate / UCF(project, EVAPRATE);
        break;

      case FILE_EVAP:
		  project->Evap.rate = project->FileValue[EVAP] / UCF(project, EVAPRATE);
        project->Evap.rate *= project->Evap.panCoeff[mon-1];
        break;

      case TEMPERATURE_EVAP:
		  project->Evap.rate = project->FileValue[EVAP] / UCF(project, EVAPRATE);
        break;

      default: project->Evap.rate = 0.0;
    }

    // --- apply climate change adjustment                                     //(5.1.007)
    project->Evap.rate += project->Adjust.evap[mon-1];                                           //(5.1.007)

    // --- set soil recovery factor
    project->Evap.recoveryFactor = 1.0;
    k = project->Evap.recoveryPattern;
    if ( k >= 0 && project->Pattern[k].type == MONTHLY_PATTERN )
    {
        project->Evap.recoveryFactor = project->Pattern[k].factor[mon-1];                        //(5.1.007)
    }
}

//=============================================================================

void setWind(Project* project, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets wind speed (mph) for a specified date.
//
{
    int yr, mon, day;

    switch ( project->Wind.type )
    {
      case MONTHLY_WIND:
        datetime_decodeDate(theDate, &yr, &mon, &day);
        project->Wind.ws = project->Wind.aws[mon-1] / UCF(project, WINDSPEED);
        break;

      case FILE_WIND:
        project->Wind.ws = project->FileValue[WIND];
        break;

      default: project->Wind.ws = 0.0;
    }
}

//=============================================================================

void updateTempTimes(Project* project, int day)
//
//  Input:   day = day of year
//  Output:  none
//  Purpose: computes time of day when min/max temperatures occur.
//           (min. temp occurs at sunrise, max. temp. at 3 hrs. < sunset)
//
{
    double decl;                       // earth's declination
    double hrang;                      // hour angle of sunrise/sunset
    double arg;

    decl  = 0.40928*cos(0.017202*(172.0-day));
    arg = -tan(decl)*project->Temp.tanAnglat;
    if      ( arg <= -1.0 ) arg = PI;
    else if ( arg >= 1.0 )  arg = 0.0;
    else                    arg = acos(arg);
    hrang = 3.8197 * arg;
    project->Hrsr  = 12.0 - hrang + project->Temp.dtlong;
    project->Hrss  = 12.0 + hrang + project->Temp.dtlong - 3.0;
    project->Dhrdy = project->Hrsr - project->Hrss;
    project->Dydif = 24.0 + project->Hrsr - project->Hrss;
    project->Hrday = (project->Hrsr + project->Hrss) / 2.0;
    project->Tave  = (project->Tmin + project->Tmax) / 2.0;
    project->Trng  = (project->Tmax - project->Tmin) / 2.0;
    if ( project->Temp.tmax == MISSING ) project->Trng1 = project->Tmax - project->Tmin;
    else                        project->Trng1 = project->Temp.tmax - project->Tmin;
    project->Temp.tmax = project->Tmax;
}

//=============================================================================

////  This function was modified for release 5.1.010.  ////                    //(5.1.010)

double getTempEvap(Project* project, int day, double tave, double trng)
//
//  Input:   day = day of year
//           tave = 7-day average temperature (deg F)
//           trng = 7-day average daily temperature range (deg F)
//  Output:  returns evaporation rate in user's units (US:in/day, SI:mm/day)
//  Purpose: uses Hargreaves method to compute daily evaporation rate
//           from daily average temperatures and Julian day.
//
{
    double a = 2.0*PI/365.0;
    double ta = (tave - 32.0)*5.0/9.0;           //average temperature (deg C)
    double tr = trng*5.0/9.0;                    //temperature range (deg C)
    double lamda = 2.50 - 0.002361 * ta;         //latent heat of vaporization
    double dr = 1.0 + 0.033*cos(a*day);          //relative earth-sun distance
    double phi = project->Temp.anglat*2.0*PI/360.0;       //latitude angle (rad)
    double del = 0.4093*sin(a*(284+day));        //solar declination angle (rad)
    double omega = acos(-tan(phi)*tan(del));     //sunset hour angle (rad)
    double ra = 37.6*dr*                         //extraterrestrial radiation
                (omega*sin(phi)*sin(del) +
                 cos(phi)*cos(del)*sin(omega));
    double e = 0.0023*ra/lamda*sqrt(tr)*(ta+17.8);    //evap. rate (mm/day)
    if ( e < 0.0 ) e = 0.0;
    if ( project->UnitSystem == US ) e /= MMperINCH;           //evap rate (in/day)
    return e;
}

//=============================================================================

int  getFileFormat(Project* project)
//
//  Input:   none
//  Output:  returns code number of climate file's format
//  Purpose: determines what format the climate file is in.
//
{
    char recdType[4] = "";
    char elemType[4] = "";
    char filler[5] = "";
    char staID[80];
    char s[80];
    char line[MAXLINE];

    int  y, m, d, n;

    // --- read first line of file
    if ( fgets(line, MAXLINE, project->Fclimate.file) == NULL ) return UNKNOWN_FORMAT;

    // --- check for TD3200 format
    sstrncpy(recdType, line, 3);
    sstrncpy(filler, &line[23], 4);
    if ( strcmp(recdType, "DLY") == 0 &&
         strcmp(filler, "9999")  == 0 ) return TD3200;

    // --- check for DLY0204 format
    if ( strlen(line) >= 233 )
    {
        sstrncpy(elemType, &line[13], 3);
        n = atoi(elemType);
        if ( n == 1 || n == 2 || n == 151 ) return DLY0204;
    }

    // --- check for USER_PREPARED format
    n = sscanf(line, "%s %d %d %d %s", staID, &y, &m, &d, s);
    if ( n == 5 ) return USER_PREPARED;

    // --- check for GHCND format                                              //(5.1.007)
	if (isGhcndFormat(project, line)) return GHCND;                                   //(5.1.007)

    return UNKNOWN_FORMAT;
}

//=============================================================================

void readFileLine(Project* project, int *y, int *m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from next line of climate file.
//
{
    // --- read next line from climate data file
    while ( strlen(project->FileLine) == 0 )
    {
        if ( fgets(project->FileLine, MAXLINE, project->Fclimate.file) == NULL ) return;
     	if ( project->FileLine[0] == '\n' ) project->FileLine[0] = '\0';
    }

    // --- parse year & month from line
    switch (project->FileFormat)
    {
	case  USER_PREPARED: readUserFileLine(project, y, m);   break;
	case  TD3200:        readTD3200FileLine(project, y, m);  break;
	case  DLY0204:       readDLY0204FileLine(project, y, m); break;
	case  GHCND:         readGhcndFileLine(project, y, m);   break;                      //(5.1.007)
    }
}

//=============================================================================

void readUserFileLine(Project* project, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of User-Prepared climate file.
//
{
    int n;
    char staID[80];
    n = sscanf(project->FileLine, "%s %d %d", staID, y, m);
    if ( n < 3 )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_READ, project->Fclimate.name);
    }
}

//=============================================================================

void readTD3200FileLine(Project* project, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of TD-3200 climate file.
//
{
    char recdType[4] = "";
    char year[5] = "";
    char month[3] = "";
    int  len;

    // --- check for minimum number of characters
    len = strlen(project->FileLine);
    if ( len < 30 )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_READ, project->Fclimate.name);
        return;
    }

    // --- check for proper type of record
    sstrncpy(recdType, project->FileLine, 3);
    if ( strcmp(recdType, "DLY") != 0 )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_READ, project->Fclimate.name);
        return;
    }

    // --- get record's date
    sstrncpy(year,  &project->FileLine[17], 4);
    sstrncpy(month, &project->FileLine[21], 2);
    *y = atoi(year);
    *m = atoi(month);
}

//=============================================================================

void readDLY0204FileLine(Project* project, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of DLY02 or DLY04 climate file.
//
{
    char year[5] = "";
    char month[3] = "";
    int  len;

    // --- check for minimum number of characters
    len = strlen(project->FileLine);
    if ( len < 16 )
    {
        report_writeErrorMsg(project,ERR_CLIMATE_FILE_READ, project->Fclimate.name);
        return;
    }

    // --- get record's date
    sstrncpy(year,  &project->FileLine[7], 4);
    sstrncpy(month, &project->FileLine[11], 2);
    *y = atoi(year);
    *m = atoi(month);
}

//=============================================================================

void readFileValues(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reads next month's worth of data from climate file.
//
{
    int  i, j;
    int  y, m;

    // --- initialize project->FileData array to missing values
    for ( i=0; i<MAXCLIMATEVARS; i++)
    {
        for (j=0; j<MAXDAYSPERMONTH; j++) project->FileData[i][j] = MISSING;
    }

    while ( !project->ErrorCode )
    {
        // --- return when date on line is after current file date
        if ( feof(project->Fclimate.file) ) return;
		readFileLine(project, &y, &m);
        if ( y > project->FileYear || m > project->FileMonth ) return;

        // --- parse climate values from file line
        switch (project->FileFormat)
        {
		case  USER_PREPARED: parseUserFileLine(project );   break;
		case  TD3200:        parseTD3200FileLine(project );  break;
		case  DLY0204:       parseDLY0204FileLine(project ); break;
		case  GHCND:         parseGhcndFileLine(project );   break;                    //(5.1.007)
        }
        strcpy(project->FileLine, "");
    }
}

//=============================================================================

void parseUserFileLine(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: parses climate variable values from a line of a user-prepared
//           climate file.
//
{
    int   n;
    int   y, m, d;
    char  staID[80];
    char  s0[80];
    char  s1[80];
    char  s2[80];
    char  s3[80];
    double x;

    // --- read day, project->Tmax, project->Tmin, Evap, & project->Wind from file line
    n = sscanf(project->FileLine, "%s %d %d %d %s %s %s %s",
        staID, &y, &m, &d, s0, s1, s2, s3);
    if ( n < 4 ) return;
    if ( d < 1 || d > 31 ) return;

    // --- process TMAX
    if ( strlen(s0) > 0 && *s0 != '*' )
    {
        x = atof(s0);
        if ( project->UnitSystem == SI ) x = 9./5.*x + 32.0;
        project->FileData[TMAX][d] =  x;
    }

    // --- process TMIN
    if ( strlen(s1) > 0 && *s1 != '*' )
    {
        x = atof(s1);
        if ( project->UnitSystem == SI ) x = 9./5.*x + 32.0;
        project->FileData[TMIN][d] =  x;
    }

    // --- process EVAP
    if ( strlen(s2) > 0 && *s2 != '*' ) project->FileData[EVAP][d] = atof(s2);

    // --- process WIND
    if ( strlen(s3) > 0 && *s3 != '*' ) project->FileData[WIND][d] = atof(s3);
}

//=============================================================================

void parseTD3200FileLine(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: parses climate variable values from a line of a TD3200 file.
//
{
    int  i;
    char param[5] = "";

    // --- parse parameter name
    sstrncpy(param, &project->FileLine[11], 4);

    // --- see if parameter is temperature, evaporation or wind speed
    for (i=0; i<MAXCLIMATEVARS; i++)
    {
        if (strcmp(param, ClimateVarWords[i]) == 0 ) setTD3200FileValues(project,i);
    }
}

//=============================================================================

void setTD3200FileValues(Project* project,int i)
//
//  Input:   i = climate variable code
//  Output:  none
//  Purpose: reads month worth of values for climate variable from TD-3200 file.
//
{
    char valCount[4] = "";
    char day[3] = "";
    char sign[2] = "";
    char value[6] = "";
    char flag2[2] = "";
    double x;
    int  nValues;
    int  j, k, d;
    int  lineLength;

    // --- parse number of days with data from cols. 27-29 of file line
    sstrncpy(valCount, &project->FileLine[27], 3);
    nValues = atoi(valCount);
    lineLength = strlen(project->FileLine);

    // --- check for enough characters on line
    if ( lineLength >= 12*nValues + 30 )
    {
        // --- for each day's value
        for (j=0; j<nValues; j++)
        {
            // --- parse day, value & flag from file line
            k = 30 + j*12;
            sstrncpy(day,   &project->FileLine[k], 2);
            sstrncpy(sign,  &project->FileLine[k+4], 1);
            sstrncpy(value, &project->FileLine[k+5], 5);
            sstrncpy(flag2, &project->FileLine[k+11], 1);

            // --- if value is valid then store it in project->FileData array
            d = atoi(day);
            if ( strcmp(value, "99999") != 0
                 && ( flag2[0] == '0' || flag2[0] == '1')
                 &&   d > 0
                 &&   d <= 31 )
            {
                // --- convert from string value to numerical value
                x = atof(value);
                if ( sign[0] == '-' ) x = -x;

                // --- convert evaporation from hundreths of inches
                if ( i == EVAP )
                {
                    x /= 100.0;

                    // --- convert to mm if using SI units
                    if ( project->UnitSystem == SI ) x *= MMperINCH;
                }

                // --- convert wind speed from miles/day to miles/hour
                if ( i == WIND ) x /= 24.0;

                // --- store value
                project->FileData[i][d] = x;
            }
        }
    }
}

//=============================================================================

void parseDLY0204FileLine(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: parses a month's worth of climate variable values from a line of
//           a DLY02 or DLY04 climate file.
//
{
    int  j, k, p;
    char param[4] = "";
    char sign[2]  = "";
    char value[6] = "";
    char code[2]  = "";
    double x;

    // --- parse parameter name
    sstrncpy(param, &project->FileLine[13], 3);

    // --- see if parameter is min or max temperature
    p = atoi(param);
    if ( p == 1 ) p = TMAX;
    else if ( p == 2 ) p = TMIN;
    else if ( p == 151 ) p = EVAP;
    else return;

    // --- check for 233 characters on line
    if ( strlen(project->FileLine) < 233 ) return;

    // --- for each of 31 days
    k = 16;
    for (j=1; j<=31; j++)
    {
        // --- parse value & flag from file line
        sstrncpy(sign,  &project->FileLine[k], 1);
        sstrncpy(value, &project->FileLine[k+1], 5);
        sstrncpy(code,  &project->FileLine[k+6], 1);
        k += 7;

        // --- if value is valid then store it in project->FileData array

        if ( strcmp(value, "99999") != 0 && strcmp(value, "     ") != 0 )
        {
            switch (p)
            {
            case TMAX:
            case TMIN:
                // --- convert from integer tenths of a degree C to degrees F
                x = atof(value) / 10.0;
                if ( sign[0] == '-' ) x = -x;
                x = 9./5.*x + 32.0;
                break;
            case EVAP:
                // --- convert from 0.1 mm to inches or mm
                x = atof(value) / 10.0;
                if ( project->UnitSystem == US ) x /= MMperINCH;
                break;
			default: return;
            }
            project->FileData[p][j] = x;
        }
    }
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

int isGhcndFormat(Project* project, char* line)
//
//  Input:   line = first line of text from a climate file
//  Output:  returns TRUE if climate file is in NCDC GHCN Daily format.
//  Purpose: Checks if a climate file is in the NCDC GHCN Daily format
//           and determines the position of each climate variable field.
//
{
    int i;
    char* ptr;

    // --- find starting position of the DATE field
    ptr = strstr(line, "DATE");
    if ( ptr == NULL ) return FALSE;
    project->FileDateFieldPos = ptr - line;

    // --- initialize starting position of each data field
    for ( i = TMIN; i <= WIND; i++) project->FileFieldPos[i] = -1;

    // --- find starting position of each climate variable's data field
    ptr = strstr(line, "TMIN");
    if ( ptr ) project->FileFieldPos[TMIN] = ptr - line;
    ptr = strstr(line, "TMAX");
    if ( ptr ) project->FileFieldPos[TMAX] = ptr - line;
    ptr = strstr(line, "EVAP");
    if ( ptr ) project->FileFieldPos[EVAP] = ptr - line;

    // --- WIND can either be daily movement or average speed
    project->FileWindType = WDMV;
    ptr = strstr(line, "WDMV");
    if ( ptr == NULL )
    {
        project->FileWindType = AWND;
        ptr = strstr(line, "AWND");
    }
    if ( ptr ) project->FileFieldPos[WIND] = ptr - line;

    // --- check if at least one climate variable was found
    for (i = TMIN; i <= WIND; i++) if (project->FileFieldPos[i] >= 0 ) return TRUE;
    return FALSE;
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

void readGhcndFileLine(Project* project, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of a NCDC GHCN Daily climate file.
//
{
    int n = sscanf(&project->FileLine[project->FileDateFieldPos], "%4d%2d", y, m);
    if ( n != 2 )
    {
        *y = -99999;
        *m = -99999;
    }
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

void parseGhcndFileLine(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: parses a line of a NCDC GHCN Daily file for daily
//           values of max/min temperature, pan evaporation and
//           wind speed.
//
{
    int y, m, d, n, v;
    double x;

    // --- parse day of month from date field
    n = sscanf(&project->FileLine[project->FileDateFieldPos], "%4d%2d%2d", &y, &m, &d);
    if ( n < 3 ) return;
    if ( d < 1 || d > 31 ) return;

    // --- parse temperatures (in tenths of deg. C) to deg F
    if ( project->FileFieldPos[TMAX] >= 0 )
    {
        if ( sscanf(&project->FileLine[project->FileFieldPos[TMAX]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
                project->FileData[TMAX][d] = (double)v*0.1*9.0/5.0 + 32.0;
        }
    }
    if ( project->FileFieldPos[TMIN] >= 0 )
    {
        if ( sscanf(&project->FileLine[project->FileFieldPos[TMIN]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
                project->FileData[TMIN][d] = (double)v*0.1*9.0/5.0 + 32.0;
        }
    }

    // -- parse evaporation (in tenths of mm) to user units
    if ( project->FileFieldPos[EVAP] >= 0 )
    {
        if ( sscanf(&project->FileLine[project->FileFieldPos[EVAP]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
            {
                x = (double)v * 0.1;
                if ( project->UnitSystem == US ) x /= MMperINCH;
                project->FileData[EVAP][d] = x;
            }
        }
    }

    // --- parse wind speed (in km/day for WDMV or tenths of m/s for AWND)
    //     to miles/hr
    if ( project->FileFieldPos[WIND] >= 0 )
    {
        if ( sscanf(&project->FileLine[project->FileFieldPos[WIND]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
            {
                if ( project->FileWindType == WDMV ) x = (double)v * 0.62137 / 24.;
                else x = (double)v * 0.1 / 1000. * 0.62137 * 3600.;
                project->FileData[WIND][d] = x;
            }
        }
    }
}

//=============================================================================

////  New function added to release 5.1.010.  ////                             //(5.1.010)

void updateTempMoveAve(Project* project, double tmin, double tmax)
//
//  Input:   tmin = minimum daily temperature (deg F)
//           tmax = maximum daily temperature (deg F)
//  Output:  none
//  Purpose: updates moving averages of average daily temperature
//           and daily temperature range stored in structure project->Tma.
//
{
    double ta,               // new day's average temperature (deg F)
           tr;               // new day's temperature range (deg F)
    int    count = project->Tma.count;

    // --- find ta and tr from new day's min and max temperature
    ta = (tmin + tmax) / 2.0;
    tr = fabs(tmax - tmin);

    // --- if the array used to store previous days' temperatures is full
    if ( count == project->Tma.maxCount )
    {
        // --- update the moving averages with the new day's value
        project->Tma.tAve = (project->Tma.tAve * count + ta - project->Tma.ta[project->Tma.front]) / count;
        project->Tma.tRng = (project->Tma.tRng * count + tr - project->Tma.tr[project->Tma.front]) / count;

        // --- replace the values at the front of the moving average window
        project->Tma.ta[project->Tma.front] = ta;
        project->Tma.tr[project->Tma.front] = tr;

        // --- move the front one position forward
        project->Tma.front++;
        if ( project->Tma.front == count ) project->Tma.front = 0;
    }

    // --- array of previous day's values not full (at start of simulation)
    else
    {
        // --- find new moving averages by adding new values to previous ones
        project->Tma.tAve = (project->Tma.tAve * count + ta) / (count + 1);
        project->Tma.tRng = (project->Tma.tRng * count + tr) / (count + 1);

        // --- save new day's values
        project->Tma.ta[project->Tma.front] = ta;
        project->Tma.tr[project->Tma.front] = tr;

        // --- increment count and front of moving average window
        project->Tma.count++;
        project->Tma.front++;
        if ( project->Tma.count == project->Tma.maxCount ) project->Tma.front = 0;
    }
}
