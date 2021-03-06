# include <vector>
# include <string>
# include <sys/stat.h>
# include <fstream>
# include <limits>
# include <cmath>
# include <cstring>
# include <unordered_map>
# include <stdlib.h>
# include "swmmparallelnsgaii.h"
# include "datetime.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

std::vector<std::string> splitText(const std::string &text, const std::string& delim)
{
  std::vector<std::string> parsed;

  size_t start = 0;
  size_t end = text.find(delim);

  while (end != std::string::npos)
  {
    parsed.push_back(text.substr(start, end - start));
    start = end + delim.length();
    end = text.find(delim, start);
  }

  parsed.push_back(text.substr(start, end));

  return parsed;
}

std::string fileExtension(const std::string& filePath)
{
  if(filePath.find(".") != std::string::npos)
  {
    return filePath.substr(filePath.find_last_of("."));
  }
  else
  {
    return "";
  }
}

bool replace(std::string& str, const std::string& from, const std::string& to)
{
  size_t start_pos = str.find(from);

  if(start_pos == std::string::npos)
    return false;

  str.replace(start_pos, from.length(), to);

  return true;
}

void replaceAll(std::string& str, const std::string& from, const std::string& to)
{
  if(from.empty())
    return;

  size_t start_pos = 0;

  while((start_pos = str.find(from, start_pos)) != std::string::npos)
  {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

bool fileExists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

std::vector<std::pair<double, double>> readTimeSeriesFile(const std::string & timeSeriesFile)
{
  auto found = loadedTimeSeriesFiles.find(timeSeriesFile);

  if(found != loadedTimeSeriesFiles.end())
  {
    return found->second;
  }
  else
  {

    std::vector<std::pair<double, double>> timeSeries;

#ifdef USE_OPENMP
#pragma omp critical
    {
#endif


      std::fstream readFileStream;
      readFileStream.open(timeSeriesFile,std::fstream::in);

      std::string line;

      if(readFileStream.is_open())
      {
        std::getline(readFileStream,line);

        while (std::getline(readFileStream,line))
        {
          std::vector<std::string> split = splitText(line, ",");
          double value = stof(split[1]);

          replaceAll(split[0],"\"","");

          std::vector<std::string> dateTime = splitText(split[0], " ");
          std::vector<std::string> date = splitText(dateTime[0], "-");
          std::vector<std::string> time = splitText(dateTime[1], ":");

          double dateTimeDbl = datetime_encodeDate(std::stoi(date[0]), std::stoi(date[1]), std::stoi(date[2])) +
              datetime_encodeTime(std::stoi(time[0]), std::stoi(time[1]), std::stoi(time[2]));

          timeSeries.push_back(std::pair<double,double>(dateTimeDbl,value));
        }
      }

      loadedTimeSeriesFiles[timeSeriesFile] = timeSeries;

#ifdef USE_OPENMP
    }
#endif

    return timeSeries;
  }
}

std::vector<std::pair<double, double>> readSWMMTimeSeries(IFaceData *faceData, int iType, int iIndex, int vIndex)
{
  std::vector<std::pair<double, double>> timeSeries;

  for(int i = 0; i < faceData->SWMM_Nperiods; i++)
  {
    double dateTime  = faceData->SWMM_StartDate + (faceData->SWMM_ReportStep * i) /(24.0 * 60 * 60);
    float value = 0.0;
    GetSwmmResult(faceData, iType, iIndex, vIndex,i+1, &value);

    timeSeries.push_back(std::pair<double,double>(dateTime,value));
  }

  return timeSeries;
}

double CalcRMS(const std::vector<std::pair<double, double>> & obs,
               const std::vector<std::pair<double, double>> & sim)
{


  double value = 0.0;
  int currentObs = 0;
  double count = 0;

  for(size_t i = 0; i < sim.size() ; i++)
  {
    double simDT = sim[i].first;

    for(size_t j = currentObs; j < obs.size() - 1; j++)
    {
      double obsDT1 = obs[j].first;
      double obsDT2 = obs[j+1].first;

      if(simDT >= obsDT1 && simDT <= obsDT2)
      {
        double obsv1 = obs[j].second;
        double obsv2 = obs[j+1].second;

        double est = obsv1 + ((simDT - obsDT1)/(obsDT2 - obsDT1))* (obsv2 - obsv1);
        est = sim[i].second - est;
        value += est * est;
        count = count + 1.0;

        currentObs = j;
        break;
      }
    }
  }

  value = count > 0.0 ? std::sqrt(value/count) : std::numeric_limits<double>::max();

  return value;
}

void RMS(const std::vector<std::string> &options, double & value, IFaceData *faceData)
{

  std::vector<std::pair<double,double>> obs = readTimeSeriesFile(options[4]);
  std::vector<std::pair<double,double>> sim = readSWMMTimeSeries(faceData, stoi(options[0]), stoi(options[1]), stoi(options[2]));

  int currentObs = 0;
  double count = 0;

  for(size_t i = 0; i < sim.size() ; i++)
  {
    double simDT = sim[i].first;

    for(size_t j = currentObs; j < obs.size() - 1; j++)
    {
      double obsDT1 = obs[j].first;
      double obsDT2 = obs[j+1].first;

      if(simDT >= obsDT1 && simDT <= obsDT2)
      {
        double obsv1 = obs[j].second;
        double obsv2 = obs[j+1].second;

        double est = obsv1 + ((simDT - obsDT1)/(obsDT2 - obsDT1))* (obsv2 - obsv1);
        est = sim[i].second - est;
        value += est * est;
        count = count + 1.0;

        currentObs = j;
        break;
      }
    }
  }

  value = count > 0.0 ? std::sqrt(value/count) : std::numeric_limits<double>::max();

}

void RMS_EXACT(const std::vector<std::string> &options, double & value, IFaceData *faceData)
{
  std::vector<std::pair<double,double>> obs = readTimeSeriesFile(options[4]);
  std::vector<std::pair<double,double>> sim = readSWMMTimeSeries(faceData, stoi(options[0]), stoi(options[1]), stoi(options[2]));

  double eps = 1.1574074074074073e-5;
  int currentObs = 0;
  double count = 0;

  for(size_t i = 0; i < sim.size() ; i++)
  {
    double simDT = sim[i].first;

    for(size_t j = currentObs; j < obs.size(); j++)
    {
      double obsDT1 = obs[j].first;

      if(fabs(simDT - obsDT1) < eps)
      {
        double est = sim[i].second - obs[j].second;
        value += est * est;
        count = count + 1.0;
        currentObs = j;
        break;
      }
    }
  }

  value = count > 0.0 ? std::sqrt(value/count) : std::numeric_limits<double>::max();
}

void VolumeError(const std::vector<std::string> &options, double & value, IFaceData *faceData)
{
  std::vector<std::pair<double,double>> obs = readTimeSeriesFile(options[4]);
  std::vector<std::pair<double,double>> sim = readSWMMTimeSeries(faceData, stoi(options[0]), stoi(options[1]), stoi(options[2]));

  double minDate = std::min(obs[0].first, sim[0].first);
  double maxDate = std::min(obs[obs.size() -1].first, sim[sim.size() - 1].first);

  double obsVol = 0.0;
  double simVol = 0.0;

  for(size_t i = 1; i < obs.size() ; i++)
  {
    double dt1 = obs[i-1].first;
    double dt2 = obs[i].first;

    if(dt1 >= minDate && dt2 <= maxDate)
    {
      obsVol += fabs(0.5 * (dt2 - dt1) * (obs[i].second - obs[i-1].second)) + fabs( std::min(obs[i].second, obs[i-1].second) * (dt2 - dt1));
    }
  }

  for(size_t i = 1; i < sim.size() ; i++)
  {
    double dt1 = sim[i-1].first;
    double dt2 = sim[i].first;

    if(dt1 >= minDate && dt2 <= maxDate)
    {
      simVol += fabs(0.5 * (dt2 - dt1) * (sim[i].second - sim[i-1].second)) + fabs( std::min(sim[i].second, sim[i-1].second) * (dt2 - dt1));
    }
  }

  value = fabs(simVol - obsVol);
}

void SWMMParallelNSGAIIInitialize(int procRank, const std::vector<std::string>& optionalArgs)
{
  printf("Initializing SWMMParallelNSGAII...\n");

  swmmExecutablePath = optionalArgs[0];
  trim(swmmExecutablePath);
  printf("swmm executable: %s", swmmExecutablePath.c_str());

  swmmInputFile = optionalArgs[1];
  trim(swmmInputFile);

  printf("swmm inputfile: %s", swmmInputFile.c_str());

  if(fileExists(swmmInputFile))
  {
    std::string line = "";
    std::string lines = "";

    std::fstream readFileStream;
    readFileStream.open(swmmInputFile,std::fstream::in);

    if(readFileStream.is_open())
    {
      std::getline(readFileStream,line);
      lines = line;

      while (std::getline(readFileStream,line))
      {
        lines += "\n" + line;
      }
    }

    swmmProjectFileLines = lines;

    readFileStream.close();

    //clean previous sim files
    if(procRank == 0)
    {
      std::string newswmmfile = swmmInputFile;
      std::vector<std::string> fileSplit = splitText(newswmmfile, "/");
      replace(newswmmfile,fileSplit[fileSplit.size() -1], "");

      std::string deleteCommand = "";

      if(newswmmfile.size())
      {
        deleteCommand = "find " + newswmmfile + " -name \"*gen_*_ind*\" -type f -delete";
      }
      else
      {
        deleteCommand = "find . -name \"*gen_*_ind*\" -type f -delete";
      }

      system(deleteCommand.c_str());
    }

    variables = splitText(optionalArgs[2], " ");
    std::vector<std::string> varmultiplierS = splitText(optionalArgs[3], " ");

    if(variables.size() != varmultiplierS.size())
    {
      abort();
    }

    multiplierFound = false;

    for(size_t i = 0 ; i < varmultiplierS.size() ; i++)
    {
      int tempm = stoi(varmultiplierS[i]);
      variableMultiplierIndexes.push_back(tempm);

      if(tempm > 0)
      {
        multiplierFound = true;
      }
    }


    if(multiplierFound == true)
    {
      std::string multiplierFile = optionalArgs[4];
      trim(multiplierFile);

      std::fstream multFileStream;
      multFileStream.open(multiplierFile,std::fstream::in);

      if(multFileStream.is_open())
      {
        std::unordered_map<std::string,double> multiplierTable;

        std::getline(multFileStream, line);

        while (std::getline(multFileStream,line))
        {
          std::vector<std::string> cols = splitText(line, " ");
          trim(cols[0]);
          trim(cols[1]);
          multiplierTable[cols[0]] = stof(cols[1]);
        }

        loadedMultiplierTables[multiplierFile] = multiplierTable;
      }
      else
      {
        printf("Multiplier file does not exist:%s\n", multiplierFile.c_str());
        abort();
      }

      multFileStream.close();

    }
  }
  else
  {
    printf("swmm input file does not exist: %s", swmmInputFile.c_str());
    abort();
  }

  initialized = true;
}

void SWMMParallelNSGAII(int gen, int indIndex, int nreal, double *xreal, int nbin, double *xbin, int *nbits,
                        int **gene, int nobj, double *obj, int ncon, double *constr, const std::vector<std::string>& optionalArgs)
{

  if(!initialized)
  {
    SWMMParallelNSGAIIInitialize(1,optionalArgs);
  }

  std::string allLines = swmmProjectFileLines;

  std::unordered_map<std::string,double> multiplierTable;

  if(multiplierFound == true)
  {
    std::string multiplierFile = optionalArgs[4];
    trim(multiplierFile);

    auto searchLoadedMultiplier = loadedMultiplierTables.find(multiplierFile);

    if(searchLoadedMultiplier != loadedMultiplierTables.end())
    {
      multiplierTable = searchLoadedMultiplier->second;
    }
  }


  std::string extension = fileExtension(swmmInputFile);
  std::string newswmmfile = swmmInputFile;
  replace(newswmmfile,extension, "gen_" + std::to_string(gen) + "_ind" + std::to_string(indIndex) + extension);

  std::string reportfile = newswmmfile;
  replace(reportfile,extension,".rpt");

  std::string outputfile = newswmmfile;
  replace(outputfile,extension,".out");

  std::fstream writeFileStream;
  writeFileStream.open(newswmmfile,std::fstream::out | std::fstream::trunc);\

  for(size_t i = 0; i < variables.size() ; i++)
  {
    std::string variable = variables[i];
    int varmult = variableMultiplierIndexes[i];

    if(varmult > 0)
    {
      for(int f = 0 ; f < varmult; f++)
      {
        std::string newVarName = variable + "_" + std::to_string(f);
        double newValue =  multiplierTable[newVarName] * (1.0 + xreal[i] / 100);
        replaceAll(allLines," " + newVarName + " " , " " + std::to_string(newValue) + " ");
      }
    }
    else
    {
      replaceAll(allLines, " " + variable + " ", " " + std::to_string(xreal[i]) + " ");
    }
  }

  writeFileStream << allLines;
  writeFileStream.flush();
  writeFileStream.close();

  char *inpO = new char[outputfile.length() + 1] ;
  std::strcpy (inpO, outputfile.c_str());

#ifdef USE_EXTERNAL_EXECUTABLE
  std::string execCommand =   swmmExecutablePath + " " + newswmmfile + " " + reportfile + " " + outputfile;

  int status = system(execCommand.c_str());

  if(status)
  {
    printf("The command may not have been issued: %s", execCommand.c_str());
    //    abort();
  }

#else
  char *inpF = new char[newswmmfile.length() + 1] ;
  std::strcpy (inpF, newswmmfile.c_str());


  char *inpR = new char[reportfile.length() + 1] ;
  std::strcpy (inpR, reportfile.c_str());

  RunSwmmDll(inpF,inpR,inpO);

  delete[] inpF;
  delete[] inpR;

#endif

  IFaceData *faceData  =  nullptr;

  if(fileExists(outputfile))
  {
    int error;
    faceData = OpenSwmmOutFile(inpO ,&error);
  }
  else
  {
    printf("Output file was not found: %s", outputfile.c_str());
    abort();
  }

  delete[] inpO;


  for(size_t i = 5; i < optionalArgs.size(); i++)
  {
    std::vector<std::string> args = splitText(optionalArgs[i], " ");

    if(args[3] == "RMS")
    {
      RMS(args, obj[i-5],faceData);
    }
    else if(args[3] == "RMS_EXACT")
    {
      RMS_EXACT(args, obj[i-5],faceData);
    }
    else if(args[3] == "VolError")
    {
      VolumeError(args, obj[i-5],faceData);
    }
  }

  CloseSwmmOutFile(faceData);

}
