#include <fstream>

#include "nrlib/iotools/fileio.hpp"

#include "src/modelsettings.h"
#include "src/definitions.h"
#include "src/io.h"

IO::IO(void)
{
}

IO::~IO(void)
{
}

void
IO::setOutputPath(const std::string & outputPath) 
{
  outputPath_ = outputPath;
}

void
IO::setFilePrefix(const std::string & filePrefix)               
{
  filePrefix_ = filePrefix;
}

std::string 
IO::makeFullFileName(const std::string & subDir, 
                     const std::string & fileName)
{
  return (outputPath_ + subDir + filePrefix_ + fileName);
}


void   
IO::writeSurfaceToFile(const Surface     & surface,
                       const std::string & baseName,
                       const std::string & path,
                       int                 format)
{
  std::string fileName = IO::makeFullFileName(path, baseName);

  if((format & ASCII) > 0)
    NRLib::WriteIrapClassicAsciiSurf(surface, fileName + SuffixAsciiIrapClassic());
  else  
    NRLib::WriteStormBinarySurf(surface, fileName + SuffixStormBinary());
}


int
IO::findGridType(const std::string & fileName)
{
  if (IsCravaBinaryFile(fileName))
    return(CRAVA);
  else if (IsStormBinaryFile(fileName))
    return(STORM);
  else
    return(SEGY);
  //  else if (NRLib::IsSegyFile(fileName))
  //    return(SEGY);
  //  else
  //    return(UNKNOWN)
}

bool IO::IsCravaBinaryFile(const std::string & fileName)
{
  std::ifstream file;
  NRLib::OpenRead(file, fileName);
  std::string line;
  getline(file,line);         
  file.close();
  std::string label = line.substr(0,20);
  if (label.compare("crava_fftgrid_binary") == 0)
    return true;
  else
    return false;
}

bool IO::IsStormBinaryFile(const std::string & fileName)
{
  std::ifstream file;
  NRLib::OpenRead(file, fileName);
  std::string line;
  getline(file,line);         
  file.close();
  std::string label = line.substr(0,18);  
  if (label.compare("storm_petro_binary") == 0)
    return true;
  else
    return false;
}

std::string IO::outputPath_ = "";
std::string IO::filePrefix_ = "";
