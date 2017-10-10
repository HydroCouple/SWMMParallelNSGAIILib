//-----------------------------------------------------------------------------
//   findroot.h
//
//   Header file for root finding method contained in findroot.c
//
//   Last modified on 11/19/13.
//-----------------------------------------------------------------------------
#ifndef FINDROOT_H
#define FINDROOT_H


typedef struct Project Project;

int findroot_Newton(Project* project, double x1, double x2, double* rts, double xacc,
	void(*func) (Project* , double x, double* f, double* df, void* p),
					void* p);
double findroot_Ridder(Project* project, double x1, double x2, double xacc,
			   double (*func)(Project*, double, void* p), void* p);


#endif
