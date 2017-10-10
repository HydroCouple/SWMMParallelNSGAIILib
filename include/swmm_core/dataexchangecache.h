#ifndef DATAEXCHANGECACHE_H
#define DATAEXCHANGECACHE_H

typedef struct Project Project;

#ifdef __cplusplus
extern "C"
{
#endif

	void addNodeLateralInflow(Project* project, int index, double value);
	int  containsNodeLateralInflow(Project* project, int index, double* value);
	int  removeNodeLateralInflow(Project* project, int index);

	void addNodeDepth(Project* project, int index, double value);
	int containsNodeDepth(Project* project, int index, double* value);
	int removeNodeDepth(Project* project, int index);

	void addSubcatchRain(Project* project, int index, double value);
	int containsSubcatchRain(Project* project, int index, double* value);
	int removeSubcatchRain(Project* project, int index);

	void clearDataCache(Project* project);


#ifdef __cplusplus
}   // matches the linkage specification from above */
#endif

#endif // DATAEXCHANGECACHE_H

