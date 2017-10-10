//-----------------------------------------------------------------------------
//   funcs.h
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14  (Build 5.1.000)
//             09/15/14  (Build 5.1.007)
//             04/02/15  (Build 5.1.008)
//             08/05/15  (Build 5.1.010)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//
//   Global interfacing functions.
//
//   Build 5.1.007:
//   - climate_readAdjustments() added.
//
//   Build 5.1.008:
//   - Function list was re-ordered and blank lines added for readability.
//   - Pollutant buildup/washoff functions for the new surfqual.c module added.
//   - Several other functions added, re-named or have modified arguments.
//
//   Build 5.1.010:
//   - New roadway_getInflow() function added.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//   Project Manager Methods
//-----------------------------------------------------------------------------

#ifndef FUNCS_H
#define FUNCS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"

void     project_open(Project* project, char *f1, char *f2, char *f3);
void     project_close(Project* project);

void     project_readInput(Project* project);
int      project_readOption(Project* project, char* s1, char* s2);
void     project_validate(Project* project);
int      project_init(Project* project);

int      project_addObject(Project* project, int type, char* id, int n);
int      project_findObject(Project* project, int type, char* id);
char*    project_findID(Project* project, int type, char* id);

double** project_createMatrix(int nrows, int ncols);
void     project_freeMatrix(double** m);

//-----------------------------------------------------------------------------
//   Input Reader Methods
//-----------------------------------------------------------------------------
int     input_countObjects(Project* project);
int     input_readData(Project* project);

//-----------------------------------------------------------------------------
//   Report Writer Methods
//-----------------------------------------------------------------------------
int     report_readOptions(Project* project, char* tok[], int ntoks);

void    report_writeLine(Project* project, char* line);
void    report_writeSysTime(Project* project);
void    report_writeLogo(Project* project);
void    report_writeTitle(Project* project);
void    report_writeOptions(Project* project);
void    report_writeReport(Project* project);

void    report_writeRainStats(Project* project, int gage, TRainStats* rainStats);
void    report_writeRdiiStats(Project* project, double totalRain, double totalRdii);

void    report_writeControlActionsHeading(Project* project);
void    report_writeControlAction(Project* project, DateTime aDate, char* linkID, double value,
        char* ruleID);

void    report_writeRunoffError(Project* project, TRunoffTotals* totals, double area);
void    report_writeLoadingError(Project* project, TLoadingTotals* totals);
void    report_writeGwaterError(Project* project, TGwaterTotals* totals, double area);
void    report_writeFlowError(Project* project, TRoutingTotals* totals);
void    report_writeQualError(Project* project, TRoutingTotals* totals);

void    report_writeMaxStats(Project* project, TMaxStats massBalErrs[], TMaxStats CourantCrit[],
        int nMaxStats);
void    report_writeMaxFlowTurns(Project* project, TMaxStats flowTurns[], int nMaxStats);
void    report_writeSysStats(Project* project, TSysStats* sysStats);

void    report_writeErrorMsg(Project* project, int code, char* msg);
void    report_writeErrorCode(Project* project);
void    report_writeInputErrorMsg(Project* project, int k, int sect, char* line, long lineCount);
void    report_writeWarningMsg(Project* project, char* msg, char* id); 
void    report_writeTseriesErrorMsg(Project* project, int code, TTable *tseries);

void    inputrpt_writeInput(Project* project);
void    statsrpt_writeReport(Project* project);

//-----------------------------------------------------------------------------
//   Temperature/Evaporation Methods
//-----------------------------------------------------------------------------
int      climate_readParams(Project* project, char* tok[], int ntoks);
int      climate_readEvapParams(Project* project, char* tok[], int ntoks);
int      climate_readAdjustments(Project* project, char* tok[], int ntoks);                      //(5.1.007)
void     climate_validate(Project* project);
void     climate_openFile(Project* project);
void     climate_initState(Project* project);
void     climate_setState(Project* project, DateTime aDate);
DateTime climate_getNextEvapDate(Project* project);                                        //(5.1.008)

//-----------------------------------------------------------------------------
//   Rainfall Processing Methods
//-----------------------------------------------------------------------------
void    rain_open(Project* project);
void    rain_close(Project* project);

//-----------------------------------------------------------------------------
//   Snowmelt Processing Methods
//-----------------------------------------------------------------------------
int     snow_readMeltParams(Project* project, char* tok[], int ntoks);
int     snow_createSnowpack(Project* project, int subcacth, int snowIndex);

void    snow_validateSnowmelt(Project* project, int snowIndex);
void    snow_initSnowpack(Project* project, int subcatch);
void    snow_initSnowmelt(Project* project, int snowIndex);

void    snow_getState(Project* project, int subcatch, int subArea, double x[]);
void    snow_setState(Project* project, int subcatch, int subArea, double x[]);

void    snow_setMeltCoeffs(Project* project, int snowIndex, double season);
void    snow_plowSnow(Project* project, int subcatch, double tStep);
double  snow_getSnowMelt(Project* project, int subcatch, double rainfall, double snowfall,
        double tStep, double netPrecip[]);
double  snow_getSnowCover(Project* project, int subcatch);

//-----------------------------------------------------------------------------
//   Runoff Analyzer Methods
//-----------------------------------------------------------------------------
int     runoff_open(Project* project);
void    runoff_execute(Project* project);
void    runoff_close(Project* project);

//-----------------------------------------------------------------------------
//   Conveyance System Routing Methods
//-----------------------------------------------------------------------------
int     routing_open(Project* project);
double  routing_getRoutingStep(Project* project, int routingModel, double fixedStep);
void    routing_execute(Project* project, int routingModel, double routingStep);
void    routing_close(Project* project, int routingModel);

//-----------------------------------------------------------------------------
//   Output Filer Methods
//-----------------------------------------------------------------------------
int     output_open(Project* project);
void    output_end(Project* project);
void    output_close(Project* project);
void    output_checkFileSize(Project* project);
void    output_saveResults(Project* project, double reportTime);
void    output_readDateTime(Project* project, int period, DateTime *aDate);
void    output_readSubcatchResults(Project* project, int period, int area);
void    output_readNodeResults(Project* project, int period, int node);
void    output_readLinkResults(Project* project, int period, int link);

//-----------------------------------------------------------------------------
//   Groundwater Methods
//-----------------------------------------------------------------------------
int     gwater_readAquiferParams(Project* project, int aquifer, char* tok[], int ntoks);
int     gwater_readGroundwaterParams(Project* project, char* tok[], int ntoks);
int     gwater_readFlowExpression(Project* project, char* tok[], int ntoks);
void    gwater_deleteFlowExpression(Project* project, int subcatch);

void    gwater_validateAquifer(Project* project, int aquifer);
void    gwater_validate(Project* project, int subcatch);

void    gwater_initState(Project* project, int subcatch);
void    gwater_getState(Project* project, int subcatch, double x[]);
void    gwater_setState(Project* project, int subcatch, double x[]);

void    gwater_getGroundwater(Project* project, int subcatch, double evap, double infil,
        double tStep);
double  gwater_getVolume(Project* project, int subcatch);

//-----------------------------------------------------------------------------
//   RDII Methods
//-----------------------------------------------------------------------------
int     rdii_readRdiiInflow(Project* project, char* tok[], int ntoks);
void    rdii_deleteRdiiInflow(Project* project, int node);
void    rdii_initUnitHyd(Project* project, int unitHyd);
int     rdii_readUnitHydParams(Project* project, char* tok[], int ntoks);
void    rdii_openRdii(Project* project);
void    rdii_closeRdii(Project* project);
int     rdii_getNumRdiiFlows(Project* project, DateTime aDate);
void    rdii_getRdiiFlow(Project* project, int index, int* node, double* q);

//-----------------------------------------------------------------------------
//   Landuse Methods
//-----------------------------------------------------------------------------
int     landuse_readParams(Project* project, int landuse, char* tok[], int ntoks);
int     landuse_readPollutParams(Project* project, int pollut, char* tok[], int ntoks);
int     landuse_readBuildupParams(Project* project, char* tok[], int ntoks);
int     landuse_readWashoffParams(Project* project, char* tok[], int ntoks);

void    landuse_getInitBuildup(Project* project, TLandFactor* landFactor,  double* initBuildup,
	    double area, double curb);
double  landuse_getBuildup(Project* project, int landuse, int pollut, double area, double curb,
        double buildup, double tStep);

double  landuse_getWashoffLoad(Project* project, int landuse, int p, double area,                //(5.1.008)
        TLandFactor landFactor[], double runoff, double vOutflow);             //(5.1.008)
double  landuse_getAvgBmpEffic(Project* project, int j, int p);
double  landuse_getCoPollutLoad(Project* project, int p, double washoff[]);

//-----------------------------------------------------------------------------
//   Flow/Quality Routing Methods
//-----------------------------------------------------------------------------
void    flowrout_init(Project* project, int routingModel);
void    flowrout_close(Project* project, int routingModel);
double  flowrout_getRoutingStep(Project* project, int routingModel, double fixedStep);
int     flowrout_execute(Project* project, int links[], int routingModel, double tStep);

void    toposort_sortLinks(Project* project, int links[]);
int     kinwave_execute(Project* project, int link, double* qin, double* qout, double tStep);

void    dynwave_validate(Project* project);                                                //(5.1.008)
void    dynwave_init(Project* project);
void    dynwave_close(Project* project);
double  dynwave_getRoutingStep(Project* project, double fixedStep);
int     dynwave_execute(Project* project, double tStep);
void    dwflow_findConduitFlow(Project* project, int j, int steps, double omega, double dt);

void    qualrout_init(Project* project);
void    qualrout_execute(Project* project, double tStep);

//-----------------------------------------------------------------------------
//   Treatment Methods
//-----------------------------------------------------------------------------
int     treatmnt_open(Project* project);
void    treatmnt_close(Project* project);
int     treatmnt_readExpression(Project* project, char* tok[], int ntoks);
void    treatmnt_delete(Project* project, int node);
void    treatmnt_treat(Project* project, int node, double q, double v, double tStep);
void    treatmnt_setInflow(Project* project, double qIn, double wIn[]);

//-----------------------------------------------------------------------------
//   Mass Balance Methods
//-----------------------------------------------------------------------------
int     massbal_open(Project* project);
void    massbal_close(Project* project);
void    massbal_report(Project* project);


void    massbal_updateRunoffTotals(Project* project, int type, double v);                        //(5.1.008)
void    massbal_updateLoadingTotals(Project* project, int type, int pollut, double w);
void    massbal_updateGwaterTotals(Project* project, double vInfil, double vUpperEvap,
        double vLowerEvap, double vLowerPerc, double vGwater);
void    massbal_updateRoutingTotals(Project* project, double tStep);

void    massbal_initTimeStepTotals(Project* project);
void    massbal_addInflowFlow(Project* project, int type, double q);
void    massbal_addInflowQual(Project* project, int type, int pollut, double w);
void    massbal_addOutflowFlow(Project* project, double q, int isFlooded);
void    massbal_addOutflowQual(Project* project, int pollut, double mass, int isFlooded);
void    massbal_addNodeLosses(Project* project, double evapLoss, double infilLoss);
void    massbal_addLinkLosses(Project* project, double evapLoss, double infilLoss);
void    massbal_addReactedMass(Project* project, int pollut, double mass);
void    massbal_addSeepageLoss(Project* project, int pollut, double seepLoss);                   //(5.1.008)
void    massbal_addToFinalStorage(Project* project, int pollut, double mass);                    //(5.1.008)
//   Simulation Statistics Methods
//-----------------------------------------------------------------------------
int     stats_open(Project* project);
void    stats_close(Project* project);
void    stats_report(Project* project);

void    stats_updateCriticalTimeCount(Project* project, int node, int link);
void    stats_updateFlowStats(Project* project, double tStep, DateTime aDate, int stepCount,
        int steadyState);
void    stats_updateSubcatchStats(Project* project, int subcatch, double rainVol, double runonVol,
        double evapVol, double infilVol, double runoffVol, double runoff);
void    stats_updateGwaterStats(Project* project, int j, double infil, double evap,              //(5.1.008)
        double latFlow, double deepFlow, double theta, double waterTable,      //(5.1.008)
        double tStep);                                                         //(5.1.008)
void    stats_updateMaxRunoff(Project* project);
void    stats_updateMaxNodeDepth(Project* project, int node, double depth);                      //(5.1.008)

//-----------------------------------------------------------------------------
//   Raingage Methods
//-----------------------------------------------------------------------------
int      gage_readParams(Project* project, int gage, char* tok[], int ntoks);
void     gage_validate(Project* project, int gage);
void     gage_initState(Project* project, int gage);
void     gage_setState(Project* project, int gage, DateTime aDate);
double   gage_getPrecip(Project* project, int gage, double *rainfall, double *snowfall);
void     gage_setReportRainfall(Project* project, int gage, DateTime aDate);
DateTime gage_getNextRainDate(Project* project, int gage, DateTime aDate);

//-----------------------------------------------------------------------------
//   Subcatchment Methods
//-----------------------------------------------------------------------------
int     subcatch_readParams(Project* project, int subcatch, char* tok[], int ntoks);
int     subcatch_readSubareaParams(Project* project, char* tok[], int ntoks);
int     subcatch_readLanduseParams(Project* project, char* tok[], int ntoks);
int     subcatch_readInitBuildup(Project* project, char* tok[], int ntoks);

void    subcatch_validate(Project* project, int subcatch);
void    subcatch_initState(Project* project, int subcatch);
void    subcatch_setOldState(Project* project, int subcatch);

double  subcatch_getFracPerv(Project* project, int subcatch);
double  subcatch_getStorage(Project* project, int subcatch);
double  subcatch_getDepth(Project* project, int subcatch);

void    subcatch_getRunon(Project* project, int subcatch);
void    subcatch_addRunonFlow(Project* project, int subcatch, double flow);                      //(5.1.008)
double  subcatch_getRunoff(Project* project, int subcatch, double tStep);

double  subcatch_getWtdOutflow(Project* project, int subcatch, double wt);
void    subcatch_getResults(Project* project, int subcatch, double wt, float x[]);

////  New functions added to release 5.1.008.  ////                            //(5.1.008)
//-----------------------------------------------------------------------------
//  Surface Pollutant Buildup/Washoff Methods
//-----------------------------------------------------------------------------
void    surfqual_initState(Project* project, int subcatch);
void    surfqual_getWashoff(Project* project, int subcatch, double runoff, double tStep);
void    surfqual_getBuildup(Project* project, int subcatch, double tStep);
void    surfqual_sweepBuildup(Project* project, int subcatch, DateTime aDate);
double  surfqual_getWtdWashoff(Project* project, int subcatch, int pollut, double wt);

//-----------------------------------------------------------------------------
//   Conveyance System Node Methods
//-----------------------------------------------------------------------------
int     node_readParams(Project* project, int node, int type, int subIndex, char* tok[], int ntoks);
void    node_validate(Project* project, int node);

void    node_initState(Project* project, int node);
void    node_initInflow(Project* project, int node, double tStep);
void    node_setOldHydState(Project* project, int node);
void    node_setOldQualState(Project* project, int node);

void    node_setOutletDepth(Project* project, int node, double yNorm, double yCrit, double z);
void    node_setDividerCutoff(int node, int link);

double  node_getSurfArea(Project* project, int node, double depth);
double  node_getDepth(Project* project, int node, double volume);
double  node_getVolume(Project* project, int node, double depth);
//double  node_getPondedDepth(int node, double volume); removed                //(5.1.008)
double  node_getPondedArea(Project* project, int node, double depth);

double  node_getOutflow(Project* project, int node, int link);
double  node_getLosses(Project* project, int node, double tStep);
double  node_getMaxOutflow(Project* project, int node, double q, double tStep);
double  node_getSystemOutflow(Project* project, int node, int *isFlooded);
void    node_getResults(Project* project, int node, double wt, float x[]);

//-----------------------------------------------------------------------------
//   Conveyance System Inflow Methods
//-----------------------------------------------------------------------------
int     inflow_readExtInflow(Project* project, char* tok[], int ntoks);
int     inflow_readDwfInflow(Project* project, char* tok[], int ntoks);
int     inflow_readDwfPattern(Project* project, char* tok[], int ntoks);

void    inflow_initDwfInflow(Project* project, TDwfInflow* inflow);
void    inflow_initDwfPattern(Project* project, int pattern);

double  inflow_getExtInflow(Project* project, TExtInflow* inflow, DateTime aDate);
double  inflow_getDwfInflow(Project* project, TDwfInflow* inflow, int m, int d, int h);
double  inflow_getPatternFactor(Project* project, int p, int month, int day, int hour);

void    inflow_deleteExtInflows(Project* project, int node);
void    inflow_deleteDwfInflows(Project* project, int node);

//-----------------------------------------------------------------------------
//   Routing Interface File Methods
//-----------------------------------------------------------------------------
int     iface_readFileParams(Project* project, char* tok[], int ntoks);
void    iface_openRoutingFiles(Project* project);
void    iface_closeRoutingFiles(Project* project);
int     iface_getNumIfaceNodes(Project* project, DateTime aDate);
int     iface_getIfaceNode(Project* project, int index);
double  iface_getIfaceFlow(Project* project, int index);
double  iface_getIfaceQual(Project* project, int index, int pollut);
void    iface_saveOutletResults(Project* project, DateTime reportDate, FILE* file);

//-----------------------------------------------------------------------------
//   Hot Start File Methods
//-----------------------------------------------------------------------------
int     hotstart_open(Project* project);
void    hotstart_close(Project* project);

//-----------------------------------------------------------------------------
//   Conveyance System Link Methods
//-----------------------------------------------------------------------------
int     link_readParams(Project* project, int link, int type, int subIndex, char* tok[], int ntoks);
int     link_readXsectParams(Project* project, char* tok[], int ntoks);
int     link_readLossParams(Project* project, char* tok[], int ntoks);

void    link_validate(Project* project, int link);
void    link_initState(Project* project, int link);
void    link_setOldHydState(Project* project, int link);
void    link_setOldQualState(Project* project, int link);

void    link_setTargetSetting(Project* project, int j);
void    link_setSetting(Project* project, int j, double tstep);
int     link_setFlapGate(Project* project, int link, int n1, int n2, double q);

double  link_getInflow(Project* project, int link);
void    link_setOutfallDepth(Project* project, int link);
double  link_getLength(Project* project, int link);
double  link_getYcrit(Project* project, int link, double q);
double  link_getYnorm(Project* project, int link, double q);
double  link_getVelocity(Project* project, int link, double q, double y);
double  link_getFroude(Project* project, int link, double v, double y);
double  link_getPower(Project* project, int link);
double  link_getLossRate(Project* project, int link, double q, double tStep);                    //(5.1.008)
char    link_getFullState(double a1, double a2, double aFull);                 //(5.1.008)

void    link_getResults(Project* project, int link, double wt, float x[]);

//-----------------------------------------------------------------------------
//   Link Cross-Section Methods
//-----------------------------------------------------------------------------
int     xsect_isOpen(int type);
int     xsect_setParams(Project* project, TXsect *xsect, int type, double p[], double ucf);
void    xsect_setIrregXsectParams(Project* project, TXsect *xsect);
void    xsect_setCustomXsectParams(Project* project, TXsect *xsect);
double  xsect_getAmax(TXsect* xsect);

double  xsect_getSofA(Project* project, TXsect* xsect, double area);
double  xsect_getYofA(Project* project, TXsect* xsect, double area);
double  xsect_getRofA(Project* project, TXsect* xsect, double area);
double  xsect_getAofS(Project* project, TXsect* xsect, double sFactor);
double  xsect_getdSdA(Project* project, TXsect* xsect, double area);
double  xsect_getAofY(Project* project, TXsect* xsect, double y);
double  xsect_getRofY(Project* project, TXsect* xsect, double y);
double  xsect_getWofY(Project* project, TXsect* xsect, double y);
double  xsect_getYcrit(Project* project, TXsect* xsect, double q);

//-----------------------------------------------------------------------------
//   Culvert/Roadway Methods                                                   //(5.1.010)
//-----------------------------------------------------------------------------
double  culvert_getInflow(Project* project, int link, double q, double h);
double  roadway_getInflow(Project* project, int link, double dir, double hcrest, double h1,
        double h2);                                                            //(5.1.010)

//-----------------------------------------------------------------------------
//   Force Main Methods
//-----------------------------------------------------------------------------
double  forcemain_getEquivN(Project* project, int j, int k);
double  forcemain_getRoughFactor(Project* project, int j, double lengthFactor);
double  forcemain_getFricSlope(Project* project,int j, double v, double hrad);

//-----------------------------------------------------------------------------
//   Cross-Section Transect Methods
//-----------------------------------------------------------------------------
int     transect_create(Project* project, int n);
void    transect_delete(Project* project);
int     transect_readParams(Project* project, int* count, char* tok[], int ntoks);
void    transect_validate(Project* project, int j);

//-----------------------------------------------------------------------------
//   Custom Shape Cross-Section Methods
//-----------------------------------------------------------------------------
int     shape_validate(TShape *shape, TTable *curve);

//-----------------------------------------------------------------------------
//   Control Rule Methods
//-----------------------------------------------------------------------------
int     controls_create(Project* project, int n);
void    controls_delete(Project* project);
int     controls_addRuleClause(Project* project, int rule, int keyword, char* Tok[], int nTokens);
int     controls_evaluate(Project* project, DateTime currentTime, DateTime elapsedTime, 
        double tStep);

//-----------------------------------------------------------------------------
//   Table & Time Series Methods
//-----------------------------------------------------------------------------
int     table_readCurve(Project* project, char* tok[], int ntoks);
int     table_readTimeseries(Project* project, char* tok[], int ntoks);

int     table_addEntry(TTable* table, double x, double y);
int     table_getFirstEntry(TTable* table, double* x, double* y);
int     table_getNextEntry(TTable* table, double* x, double* y);
void    table_deleteEntries(TTable* table);

void    table_init(TTable* table);
int     table_validate(TTable* table);
//      table_interpolate now defined in table.c                               //(5.1.008)

double  table_lookup(TTable* table, double x);
double  table_lookupEx(TTable* table, double x);
double  table_intervalLookup(TTable* table, double x);
double  table_inverseLookup(TTable* table, double y);

double  table_getSlope(TTable *table, double x);
double  table_getMaxY(TTable *table, double x);
double  table_getArea(TTable* table, double x);
double  table_getInverseArea(TTable* table, double a);

void    table_tseriesInit(TTable *table);
double  table_tseriesLookup(TTable* table, double t, char extend);

//-----------------------------------------------------------------------------
//   Utility Methods
//-----------------------------------------------------------------------------
double   UCF(Project* project, int quantity);                   // units conversion factor
int      getInt(char *s, int *y);             // get integer from string
int      getFloat(char *s, float *y);         // get float from string
int      getDouble(char *s, double *y);       // get double from string
char*    getTempFileName(Project* project, char *s);            // get temporary file name
int      findmatch(char *s, char *keyword[]); // search for matching keyword
int      match(char *str, char *substr);      // true if substr matches part of str
int      strcomp(char *s1, char *s2);         // case insensitive string compare
char*    sstrncpy(char *dest, const char *src,
         size_t maxlen);                      // safe string copy
void     writecon(char *s);                   // writes string to console
DateTime getDateTime(Project* project, double elapsedMsec);     // convert elapsed time to date
void     getElapsedTime(Project* project, DateTime aDate,       // convert elapsed date
         int* days, int* hrs, int* mins);


//-----------------------------------------------------------------------------
//   Coupling Functions
//-----------------------------------------------------------------------------

void setCouplingNodeDepths(Project* project);
//void setCouplingNodeDepth(Project* project, int index);

void setCouplingLateralInflows(Project* project);
//void setCouplingLateralInflow(Project* project, int index);

#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif
