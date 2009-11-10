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

std::string IO::outputPath_ = "";
std::string IO::filePrefix_ = "CRAVA_";
