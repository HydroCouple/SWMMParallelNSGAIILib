#ifndef SWMMPARALLELNSGAII_H
#define SWMMPARALLELNSGAII_H

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "swmmparallelnsgaii_global.h"
#include "swmm5_iface.h"

#ifdef __cplusplus
extern "C" {
#endif

std::unordered_map<std::string, std::vector<std::pair<double, double>>> loadedTimeSeriesFiles;
std::unordered_map<std::string, std::unordered_map<std::string, double>> loadedMultiplierTables;
std::string mSWMMProjectFileLines;
std::string swmmExecutablePath;

static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
    return !std::isspace(ch);
  }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

// trim from start (copying)
static inline std::string ltrim_copy(std::string s) {
  ltrim(s);
  return s;
}

// trim from end (copying)
static inline std::string rtrim_copy(std::string s) {
  rtrim(s);
  return s;
}

// trim from both ends (copying)
static inline std::string trim_copy(std::string s) {
  trim(s);
  return s;
}

std::vector<std::string> splitText(const std::string &text, const std::string& delim);

std::string fileExtension(const std::string& filePath);

bool replace(std::string& str, const std::string& from, const std::string& to);

void replaceAll(std::string& str, const std::string& from, const std::string& to);

bool fileExists(const std::string& name);

std::vector<std::pair<double, double>> readTimeSeriesFile(const std::string & timeSeriesFile);

std::vector<std::pair<double, double>> readSWMMTimeSeries(IFaceData *faceData, int iType, int iIndex, int vIndex);

double CalcRMS(const std::vector<std::pair<double, double>> & obs,
               const std::vector<std::pair<double, double>> & sim);

void RMS(const std::vector<std::string> &options, double & value, IFaceData *faceData);

double CalcVolumeError(const std::vector<std::pair<double, double>> & obs,
                       const std::vector<std::pair<double, double>> & sim);

void VolumeError(const std::vector<std::string> &options, double & value, IFaceData *faceData);

void SWMMParallelNSGAII_EXPORT SWMMParallelNSGAII(int gen, int indIndex, int nreal, double *xreal, int nbin, double *xbin, int *nbits,
                                                  int **gene, int nobj, double *obj, int ncon, double *constr, const std::vector<std::string>& optionalArgs);

#ifdef __cplusplus
}
#endif

#endif // SWMMPARALLELNSGAII_H
