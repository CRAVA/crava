/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

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

  if ((format & ASCII) > 0)
    surface.WriteToFile(fileName + SuffixAsciiIrapClassic(), NRLib::SURF_IRAP_CLASSIC_ASCII);

  // Ensure surfaces in STORM format when STORM grids are requested. Otherwise
  // STORM cubes cannot be imported to RMS.
  if ((format & STORM) > 0 || (format & ASCII) == 0)
    surface.WriteToFile(fileName + SuffixStormBinary(), NRLib::SURF_STORM_BINARY);
}


int
IO::findGridType(const std::string & fileName)
{
  //This conversion may seem strange, but we need our internal formatflags, so we convert
  //from NRLib flags here.
  if (IsCravaBinaryFile(fileName))
    return(CRAVA);
  else {
    int fType = NRLib::FindGridFileType(fileName);
    if (fType == NRLib::STORM_PETRO_BINARY)
      return(STORM);
    else if (fType == NRLib::SEGY)
      return(SEGY);
    else if (fType == NRLib::SGRI)
      return(SGRI);
    else
      return(UNKNOWN);
  }
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
