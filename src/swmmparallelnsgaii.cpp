
#include <vector>
#include <string>
# include <sys/stat.h>
#include "swmmparallelnsgaii.h"


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




void SWMMParallelNSGAII(int gen, int nreal, double *xreal, int nbin, double *xbin, int *nbits,
                        int **gene, int nobj, double *obj, int ncon, double *constr, const std::list<std::string>& optionalArgs)
{
  int i, j;
  int u[11];
  int v[11];
  double f1, f2, g, h;

  for (i=0; i<11; i++)
  {
    u[i] = 0;
  }

  for (j=0; j<30; j++)
  {
    if (gene[0][j] == 1)
    {
      u[0]++;
    }
  }

  for (i=1; i<11; i++)
  {
    for (j=0; j<4; j++)
    {
      if (gene[i][j] == 1)
      {
        u[i]++;
      }
    }
  }

  f1 = 1.0 + u[0];

  for (i=1; i<11; i++)
  {
    if (u[i] < 5)
    {
      v[i] = 2 + u[i];
    }
    else
    {
      v[i] = 1;
    }
  }

  g = 0;

  for (i=1; i<11; i++)
  {
    g += v[i];
  }

  h = 1.0/f1;
  f2 = g*h;
  obj[0] = f1;
  obj[1] = f2;
  return;
}
