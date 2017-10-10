//-----------------------------------------------------------------------------
//   iface.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//   Author:   L. Rossman
//
//   Routing interface file functions.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
extern double Qcf[];                   // flow units conversion factors
                                       // (see swmm5.c)
////-----------------------------------------------------------------------------                  
////  Shared variables
////-----------------------------------------------------------------------------                  
//static int      project->IfaceFlowUnits;        // flow units for routing interface file
//static int      project->IfaceStep;             // interface file time step (sec)
//static int      project->NumIfacePolluts;       // number of pollutants in interface file
//static int*     project->IfacePolluts;          // indexes of interface file pollutants
//static int      project->NumIfaceNodes;         // number of nodes on interface file
//static int*     project->IfaceNodes;            // indexes of nodes on interface file
//static double** project->OldIfaceValues;        // interface flows & WQ at previous time
//static double** project->NewIfaceValues;        // interface flows & WQ at next time
//static double   project->IfaceFrac;             // fraction of interface file time step
//static DateTime project->OldIfaceDate;          // previous date of interface values
//static DateTime project->NewIfaceDate;          // next date of interface values

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  iface_readFileParams     (called by input_readLine)
//  iface_openRoutingFiles   (called by routing_open)
//  iface_closeRoutingFiles  (called by routing_close)
//  iface_getNumIfaceNodes   (called by addIfaceInflows in routing.c)
//  iface_getIfaceNode       (called by addIfaceInflows in routing.c)
//  iface_getIfaceFlow       (called by addIfaceInflows in routing.c)
//  iface_getIfaceQual       (called by addIfaceInflows in routing.c)
//  iface_saveOutletResults  (called by output_saveResults)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void  openFileForOutput(Project* project);
static void  openFileForInput(Project* project);
static int   getIfaceFilePolluts(Project* project);
static int   getIfaceFileNodes(Project* project);
static void  setOldIfaceValues(Project* project);
static void  readNewIfaceValues(Project* project);
static int   isOutletNode(Project* project, int node);


//=============================================================================

int iface_readFileParams(Project* project, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads interface file information from a line of input data.
//
//  Data format is:
//  USE/SAVE  FileType  FileName
//
{
    char  k;
    int   j;

    // --- determine file disposition and type
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
    k = (char)findmatch(tok[0], FileModeWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);
    j = findmatch(tok[1], FileTypeWords);
    if ( j < 0 ) return error_setInpError(ERR_KEYWORD, tok[1]);
    if ( ntoks < 3 ) return 0;

    // --- process file name
    switch ( j )
    {
      case RAINFALL_FILE:
        project->Frain.mode = k;
        sstrncpy(project->Frain.name, tok[2], MAXFNAME);
        break;

      case RUNOFF_FILE:
        project->Frunoff.mode = k;
        sstrncpy(project->Frunoff.name, tok[2], MAXFNAME);
        break;

      case HOTSTART_FILE:
        if ( k == USE_FILE )
        {
            project->Fhotstart1.mode = k;
            sstrncpy(project->Fhotstart1.name, tok[2], MAXFNAME);
        }
        else if ( k == SAVE_FILE )
        {
            project->Fhotstart2.mode = k;
            sstrncpy(project->Fhotstart2.name, tok[2], MAXFNAME);
        }
        break;

      case RDII_FILE:
        project->Frdii.mode = k;
        sstrncpy(project->Frdii.name, tok[2], MAXFNAME);
        break;

      case INFLOWS_FILE:
        if ( k != USE_FILE ) return error_setInpError(ERR_ITEMS, "");
        project->Finflows.mode = k;
        sstrncpy(project->Finflows.name, tok[2], MAXFNAME);
        break;

      case OUTFLOWS_FILE:
        if ( k != SAVE_FILE ) return error_setInpError(ERR_ITEMS, "");
        project->Foutflows.mode = k;
        sstrncpy(project->Foutflows.name, tok[2], MAXFNAME);
        break;
    }
    return 0;
}

//=============================================================================

void iface_openRoutingFiles(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens routing interface files.
//
{
    // --- initialize shared variables
    project->NumIfacePolluts = 0;
    project->IfacePolluts = NULL;
    project->NumIfaceNodes = 0;
    project->IfaceNodes = NULL;
    project->OldIfaceValues = NULL;
    project->NewIfaceValues = NULL;

    // --- check that inflows & outflows files are not the same
    if ( project->Foutflows.mode != NO_FILE && project->Finflows.mode != NO_FILE )
    {
        if ( strcomp(project->Foutflows.name, project->Finflows.name) )
        {
            report_writeErrorMsg(project,ERR_ROUTING_FILE_NAMES, "");
            return;
        }
    }

    // --- open the file for reading or writing
    if ( project->Foutflows.mode == SAVE_FILE ) openFileForOutput(project);
	if (project->Finflows.mode == USE_FILE) openFileForInput(project);
}

//=============================================================================

void iface_closeRoutingFiles(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: closes routing interface files.
//
{
    FREE(project->IfacePolluts);
    FREE(project->IfaceNodes);
    if ( project->OldIfaceValues != NULL ) project_freeMatrix(project->OldIfaceValues);
    if ( project->NewIfaceValues != NULL ) project_freeMatrix(project->NewIfaceValues);
    if ( project->Finflows.file )  fclose(project->Finflows.file);
    if ( project->Foutflows.file ) fclose(project->Foutflows.file);
}

//=============================================================================

int iface_getNumIfaceNodes(Project* project, DateTime currentDate)
//
//  Input:   currentDate = current date/time
//  Output:  returns number of interface nodes if data exists or
//           0 otherwise
//  Purpose: reads inflow data from interface file at current date.
//
{
    // --- return 0 if file begins after current date
    if ( project->OldIfaceDate > currentDate ) return 0;

    // --- keep updating new interface values until current date bracketed
    while ( project->NewIfaceDate < currentDate && project->NewIfaceDate != NO_DATE )
    {
		setOldIfaceValues(project);
		readNewIfaceValues(project);
    }

    // --- return 0 if no data available
    if ( project->NewIfaceDate == NO_DATE ) return 0;

    // --- find fraction current date is bewteen old & new interface dates
    project->IfaceFrac = (currentDate - project->OldIfaceDate) / (project->NewIfaceDate - project->OldIfaceDate);
    project->IfaceFrac = MAX(0.0, project->IfaceFrac);
    project->IfaceFrac = MIN(project->IfaceFrac, 1.0);

    // --- return number of interface nodes
    return project->NumIfaceNodes;
}

//=============================================================================

int iface_getIfaceNode(Project* project, int index)
//
//  Input:   index = interface file node index
//  Output:  returns project node index
//  Purpose: finds index of project node associated with interface node index
//
{
    if ( index >= 0 && index < project->NumIfaceNodes ) return project->IfaceNodes[index];
    else return -1;
}

//=============================================================================

double iface_getIfaceFlow(Project* project, int index)
//
//  Input:   index = interface file node index
//  Output:  returns inflow to node
//  Purpose: finds interface flow for particular node index.
//
{
    double q1, q2;

    if ( index >= 0 && index < project->NumIfaceNodes )
    {
        // --- interpolate flow between old and new values
        q1 = project->OldIfaceValues[index][0];
        q2 = project->NewIfaceValues[index][0];
        return (1.0 - project->IfaceFrac)*q1 + project->IfaceFrac*q2;
    }
    else return 0.0;
}

//=============================================================================

double iface_getIfaceQual(Project* project, int index, int pollut)
//
//  Input:   index = index of node on interface file
//           pollut = index of pollutant on interface file
//  Output:  returns inflow pollutant concentration
//  Purpose: finds interface concentration for particular node index & pollutant.
//
{
    int    j;
    double c1, c2;

    if ( index >= 0 && index < project->NumIfaceNodes )
    {
        // --- find index of pollut on interface file
        j = project->IfacePolluts[pollut];
        if ( j < 0 ) return 0.0;

        // --- interpolate flow between old and new values
        //     (remember that 1st col. of values matrix is for flow)
        c1 = project->OldIfaceValues[index][j+1];
        c2 = project->NewIfaceValues[index][j+1];
        return (1.0 - project->IfaceFrac)*c1 + project->IfaceFrac*c2;
    }
    else return 0.0;
}

//=============================================================================

void iface_saveOutletResults(Project* project, DateTime reportDate, FILE* file)
//
//  Input:   reportDate = reporting date/time
//           file = ptr. to interface file
//  Output:  none
//  Purpose: saves system outflows to routing interface file.
//
{
    int i, p, yr, mon, day, hr, min, sec;
    char theDate[25];
    datetime_decodeDate(reportDate, &yr, &mon, &day);
    datetime_decodeTime(reportDate, &hr, &min, &sec);
    sprintf(theDate, " %04d %02d  %02d  %02d  %02d  %02d ",
            yr, mon, day, hr, min, sec);
    for (i=0; i<project->Nobjects[NODE]; i++)
    {
        // --- check that node is an outlet node
		if (!isOutletNode(project,i)) continue;

        // --- write node ID, date, flow, and quality to file
        fprintf(file, "\n%-16s", project->Node[i].ID);
        fprintf(file, "%s", theDate);
		fprintf(file, " %-10f", project->Node[i].inflow * UCF(project,FLOW));
        for ( p = 0; p < project->Nobjects[POLLUT]; p++ )
        {
            fprintf(file, " %-10f", project->Node[i].newQual[p]);
        }
    }
}

//=============================================================================

void openFileForOutput(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a routing interface file for writing.
//
{
    int i, n;

    // --- open the routing file for writing text
    project->Foutflows.file = fopen(project->Foutflows.name, "wt");
    if ( project->Foutflows.file == NULL )
    {
        report_writeErrorMsg(project,ERR_ROUTING_FILE_OPEN, project->Foutflows.name);
        return;
    }

    // --- write title & reporting time step to file
    fprintf(project->Foutflows.file, "SWMM5 Interface File");
    fprintf(project->Foutflows.file, "\n%s", project->Title[0]);
    fprintf(project->Foutflows.file, "\n%-4d - reporting time step in sec", project->ReportStep);

    // --- write number & names of each constituent (including flow) to file
    fprintf(project->Foutflows.file, "\n%-4d - number of constituents as listed below:",
            project->Nobjects[POLLUT] + 1);
    fprintf(project->Foutflows.file, "\nFLOW %s", FlowUnitWords[project->FlowUnits]);
    for (i=0; i<project->Nobjects[POLLUT]; i++)
    {
        fprintf(project->Foutflows.file, "\n%s %s", project->Pollut[i].ID,
            QualUnitsWords[project->Pollut[i].units]);
    }

    // --- count number of outlet nodes
    n = 0;
    for (i=0; i<project->Nobjects[NODE]; i++)
    {
		if (isOutletNode(project,i)) n++;
    }

    // --- write number and names of outlet nodes to file
    fprintf(project->Foutflows.file, "\n%-4d - number of nodes as listed below:", n);
    for (i=0; i<project->Nobjects[NODE]; i++)
    {
		if (isOutletNode(project,i))
            fprintf(project->Foutflows.file, "\n%s", project->Node[i].ID);
    }

    // --- write column headings
    fprintf(project->Foutflows.file,
        "\nNode             Year Mon Day Hr  Min Sec FLOW      ");
    for (i=0; i<project->Nobjects[POLLUT]; i++)
    {
        fprintf(project->Foutflows.file, " %-10s", project->Pollut[i].ID);
    }

    // --- if reporting starts immediately, save initial outlet values
    if ( project->ReportStart == project->StartDateTime )
    {
        iface_saveOutletResults(project,project->ReportStart, project->Foutflows.file);
    }
}

//=============================================================================

void openFileForInput(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a routing interface file for reading.
//
{
    int   err;                         // error code
    char  line[MAXLINE+1];             // line from Routing interface file
    char  s[MAXLINE+1];                // general string variable

    // --- open the routing interface file for reading text
    project->Finflows.file = fopen(project->Finflows.name, "rt");
    if ( project->Finflows.file == NULL )
    {
        report_writeErrorMsg(project,ERR_ROUTING_FILE_OPEN, project->Finflows.name);
        return;
    }

    // --- check for correct file type
    fgets(line, MAXLINE, project->Finflows.file);
    sscanf(line, "%s", s);
    if ( !strcomp(s, "SWMM5") )
    {
        report_writeErrorMsg(project,ERR_ROUTING_FILE_FORMAT, project->Finflows.name);
        return;
    }

    // --- skip title line
    fgets(line, MAXLINE, project->Finflows.file);

    // --- read reporting time step (sec)
    project->IfaceStep = 0;
    fgets(line, MAXLINE, project->Finflows.file);
    sscanf(line, "%d", &project->IfaceStep);
    if ( project->IfaceStep <= 0 )
    {
        report_writeErrorMsg(project,ERR_ROUTING_FILE_FORMAT, project->Finflows.name);
        return;
    }

    // --- match constituents in file with those in project
	err = getIfaceFilePolluts(project);
    if ( err > 0 )
    {
        report_writeErrorMsg(project,err, project->Finflows.name);
        return;
    }

    // --- match nodes in file with those in project
	err = getIfaceFileNodes(project);
    if ( err > 0 )
    {
        report_writeErrorMsg(project,err, project->Finflows.name);
        return;
    }

    // --- create matrices for old & new interface flows & WQ values
    project->OldIfaceValues = project_createMatrix(project->NumIfaceNodes,
                                         1+project->NumIfacePolluts);
    project->NewIfaceValues = project_createMatrix(project->NumIfaceNodes,
                                         1+project->NumIfacePolluts);
    if ( project->OldIfaceValues == NULL || project->NewIfaceValues == NULL )
    {
        report_writeErrorMsg(project,ERR_MEMORY, "");
        return;
    }

    // --- read in new interface flows & WQ values
	readNewIfaceValues(project);
    project->OldIfaceDate = project->NewIfaceDate;
}

//=============================================================================

int  getIfaceFilePolluts(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: reads names of pollutants saved on the inflows interface file.
//
{
    int   i, j;
    char  line[MAXLINE+1];             // line from inflows interface file
    char  s1[MAXLINE+1];               // general string variable
    char  s2[MAXLINE+1];         

    // --- read number of pollutants (minus FLOW)
    fgets(line, MAXLINE, project->Finflows.file);
    sscanf(line, "%d", &project->NumIfacePolluts);
    project->NumIfacePolluts--;
    if ( project->NumIfacePolluts < 0 ) return ERR_ROUTING_FILE_FORMAT;

    // --- read flow units
    fgets(line, MAXLINE, project->Finflows.file);
    sscanf(line, "%s %s", s1, s2);
    if ( !strcomp(s1, "FLOW") )  return ERR_ROUTING_FILE_FORMAT;
    project->IfaceFlowUnits = findmatch(s2, FlowUnitWords);
    if ( project->IfaceFlowUnits < 0 ) return ERR_ROUTING_FILE_FORMAT;

    // --- allocate memory for pollutant index array
    if ( project->Nobjects[POLLUT] > 0 )
    {
        project->IfacePolluts = (int *) calloc(project->Nobjects[POLLUT], sizeof(int));
        if ( !project->IfacePolluts ) return ERR_MEMORY;
        for (i=0; i<project->Nobjects[POLLUT]; i++) project->IfacePolluts[i] = -1;
    }

    // --- read pollutant names & units
    if ( project->NumIfacePolluts > 0 )
    {
        // --- check each pollutant name on file with project's pollutants
        for (i=0; i<project->NumIfacePolluts; i++)
        {
            if ( feof(project->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
            fgets(line, MAXLINE, project->Finflows.file);
            sscanf(line, "%s %s", s1, s2);
            if ( project->Nobjects[POLLUT] > 0 )
            {
				j = project_findObject(project, POLLUT, s1);
                if ( j < 0 ) continue;
                if ( !strcomp(s2, QualUnitsWords[project->Pollut[j].units]) )
                    return ERR_ROUTING_FILE_NOMATCH;
                project->IfacePolluts[j] = i;
            }
        }
    }
    return 0;
}

//=============================================================================

int getIfaceFileNodes(Project* project)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: reads names of nodes contained on inflows interface file.
//
{
    int   i;
    char  line[MAXLINE+1];             // line from inflows interface file
    char  s[MAXLINE+1];                // general string variable

    // --- read number of interface nodes
    fgets(line, MAXLINE, project->Finflows.file);
    sscanf(line, "%d", &project->NumIfaceNodes);
    if ( project->NumIfaceNodes <= 0 ) return ERR_ROUTING_FILE_FORMAT;

    // --- allocate memory for interface nodes index array
    project->IfaceNodes = (int *) calloc(project->NumIfaceNodes, sizeof(int));
    if ( !project->IfaceNodes ) return ERR_MEMORY;

    // --- read names of interface nodes from file & save their indexes
    for ( i=0; i<project->NumIfaceNodes; i++ )
    {
        if ( feof(project->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
        fgets(line, MAXLINE, project->Finflows.file);
        sscanf(line, "%s", s);
		project->IfaceNodes[i] = project_findObject(project, NODE, s);
    }

    // --- skip over column headings line
    if ( feof(project->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
    fgets(line, MAXLINE, project->Finflows.file);
    return 0;
}

//=============================================================================

void readNewIfaceValues(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: reads data from inflows interface file for next date.
//
{
    int    i, j;
    char*  s;
    int    yr = 0, mon = 0, day = 0,
		   hr = 0, min = 0, sec = 0;   // year, month, day, hour, minute, second
    char   line[MAXLINE+1];            // line from interface file

    // --- read a line for each interface node
    project->NewIfaceDate = NO_DATE;
    for (i=0; i<project->NumIfaceNodes; i++)
    {
        if ( feof(project->Finflows.file) ) return;
        fgets(line, MAXLINE, project->Finflows.file);

        // --- parse date & time from line
        if ( strtok(line, SEPSTR) == NULL ) return;
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        yr  = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        mon = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        day = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        hr  = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        min = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        sec = atoi(s);

        // --- parse flow value
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        project->NewIfaceValues[i][0] = atof(s) / Qcf[project->IfaceFlowUnits]; 

        // --- parse pollutant values
        for (j=1; j<=project->NumIfacePolluts; j++)
        {
            s = strtok(NULL, SEPSTR);
            if ( s == NULL ) return;
            project->NewIfaceValues[i][j] = atof(s);
        }

    }

    // --- encode date & time values
    project->NewIfaceDate = datetime_encodeDate(yr, mon, day) +
                   datetime_encodeTime(hr, min, sec);
}

//=============================================================================

void setOldIfaceValues(Project* project)
//
//  Input:   none
//  Output:  none
//  Purpose: replaces old values read from routing interface file with new ones. 
//
{
    int i, j;
    project->OldIfaceDate = project->NewIfaceDate;
    for ( i=0; i<project->NumIfaceNodes; i++)
    {
        for ( j=0; j<project->NumIfacePolluts+1; j++ )
        {
            project->OldIfaceValues[i][j] = project->NewIfaceValues[i][j];
        }
    }
}

//=============================================================================

int  isOutletNode(Project* project, int i)
//
//  Input:   i = node index
//  Output:  returns 1 if node is an outlet, 0 if not.
//  Purpose: determines if a node is an outlet point or not.
//
{
    // --- for DW routing only outfalls are outlets
    if ( project->RouteModel == DW )
    {
        return (project->Node[i].type == OUTFALL);
    }

    // --- otherwise outlets are nodes with no outflow links (degree is 0)
    else return (project->Node[i].degree == 0);
}
