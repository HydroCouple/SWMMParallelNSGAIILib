//-----------------------------------------------------------------------------
//  odesolve.h
//
//  Header file for ODE solver contained in odesolve.c
//
//-----------------------------------------------------------------------------
#ifndef ODESOLVE_H
#define ODESOLVE_H

typedef struct Project Project;

// functions that open, close, and use the ODE solver
int  odesolve_open(Project *project, int n);
void odesolve_close(Project *project);
int  odesolve_integrate(Project* project, double ystart[], int n, double x1, double x2,
        double eps, double h1, void(*derivs)(Project*,double, double*, double*));

#endif
