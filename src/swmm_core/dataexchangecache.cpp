//#include "stdafx.h"
#include "dataexchangecache.h"
#include <unordered_map>

using namespace std;

unordered_map<Project*, unordered_map<int, double> > NodeLateralInflows;
unordered_map<Project*, unordered_map<int, double> > NodeDepths;
unordered_map<Project*, unordered_map<int, double> > SubcatchRainfall;

typedef struct OpenMIDataCache OpenMIDataCache;

//node lateral inflow
void addNodeLateralInflow(Project* project, int index, double value)
{
   NodeLateralInflows[project][index] = value;
}

int containsNodeLateralInflow(Project* project, int index, double* const  value)
{
   unordered_map<Project*, unordered_map<int, double> >::iterator it = NodeLateralInflows.find(project);

   if (it != NodeLateralInflows.end())
   {
      unordered_map<int, double > foundProject = (*it).second;

      unordered_map<int, double > ::iterator it1 = foundProject.find(index);

      if (it1 != foundProject.end())
      {
         *value = (*it1).second;
         return 1;
      }
   }

   return 0;
}

//node lateral inflow
int removeNodeLateralInflow(Project* project, int index)
{
   unordered_map<Project*, unordered_map<int, double> >::iterator it = NodeLateralInflows.find(project);

   if (it != NodeLateralInflows.end())
   {
      unordered_map<int, double > foundProject = (*it).second;
      return foundProject.erase(index);
   }

   return 0;
}

//Node Depths
void addNodeDepth(Project* project, int index, double value)
{
   NodeDepths[project][index] = value;
}

int containsNodeDepth(Project* project, int index, double* const value)
{
   unordered_map<Project*, unordered_map<int, double> >::iterator it = NodeDepths.find(project);

   if (it != NodeDepths.end())
   {
      unordered_map<int, double > foundProject = (*it).second;

      unordered_map<int, double > ::iterator it1 = foundProject.find(index);

      if (it1 != foundProject.end())
      {
         *value = (*it1).second;
         return 1;
      }
   }

   return 0;
}

int removeNodeDepth(Project* project, int index)
{
   unordered_map<Project*, unordered_map<int, double> >::iterator it = NodeDepths.find(project);

   if (it != NodeDepths.end())
   {
      unordered_map<int, double > foundProject = (*it).second;
      return foundProject.erase(index);
   }

   return 0;
}

//SubcatchRainfall
void addSubcatchRain(Project* project, int index, double value)
{
   SubcatchRainfall[project][index] = value;
}

int containsSubcatchRain(Project* project, int index, double* const value)
{

   unordered_map<Project*, unordered_map<int, double> >::iterator it = SubcatchRainfall.find(project);

   if (it != SubcatchRainfall.end())
   {
      unordered_map<int, double > foundProject = (*it).second;

      unordered_map<int, double > ::iterator it1 = foundProject.find(index);

      if (it1 != foundProject.end())
      {
         *value = (*it1).second;
         return 1;
      }
   }

   return 0;
}

int removeSubcatchRain(Project* project, int index)
{
   unordered_map<Project*, unordered_map<int, double> >::iterator it = SubcatchRainfall.find(project);

   if (it != SubcatchRainfall.end())
   {
      unordered_map<int, double > foundProject = (*it).second;
      return foundProject.erase(index);
   }

   return 0;
}

void clearDataCache(Project *project)
{
   NodeLateralInflows.erase(project);
   NodeDepths.erase(project);
   SubcatchRainfall.erase(project);
}

