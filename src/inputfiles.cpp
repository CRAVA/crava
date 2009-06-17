#include <fstream>

#include "src/inputfiles.h"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/exception/exception.hpp"

InputFiles::InputFiles(void)
  : seedFile_(""),          
    wellFiles_(0),           
    seismicFiles_(0),        
    waveletFiles_(0),        
    waveletEstIntFile_(2),  
    faciesEstIntFile_(2),   
    timeSurfFiles_(0),       
    depthSurfFiles_(2),      
    velocityField_(""),
    backFile_(3),           
    backVelFile_(""),       
    reflMatrFile_(""),
    corrDirFile_(""),
    paramCorrFile_(""),
    tempCorrFile_("")
{
}

InputFiles::~InputFiles(void)
{
}

std::string 
InputFiles::addInputPathAndCheckFiles()
{
  unsigned int i;
  std::string errTxt = addPathAndCheck(seedFile_);
  for(i=0;i<wellFiles_.size();i++)
    errTxt += addPathAndCheck(wellFiles_[i]);
  for(i=0;i<seismicFiles_.size();i++)
    errTxt += addPathAndCheck(seismicFiles_[i]);
  for(i=0;i<waveletFiles_.size();i++)
    errTxt += addPathAndCheck(waveletFiles_[i]);
  for(i=0;i<waveletEstIntFile_.size();i++)
    errTxt += addPathAndCheck(waveletEstIntFile_[i]);
  for(i=0;i<faciesEstIntFile_.size();i++)
    errTxt += addPathAndCheck(faciesEstIntFile_[i]);
  for(i=0;i<timeSurfFiles_.size();i++)
    errTxt += addPathAndCheck(timeSurfFiles_[i], true);
  for(i=0;i<depthSurfFiles_.size();i++)
    errTxt += addPathAndCheck(depthSurfFiles_[i]);
  errTxt += addPathAndCheck(velocityField_);
  for(i=0;i<backFile_.size();i++)
    errTxt += addPathAndCheck(backFile_[i]);
  errTxt += addPathAndCheck(backVelFile_);
  errTxt += addPathAndCheck(corrDirFile_);
  errTxt += addPathAndCheck(reflMatrFile_);
  errTxt += addPathAndCheck(paramCorrFile_);
  errTxt += addPathAndCheck(tempCorrFile_);
  return(errTxt);
}


std::string
InputFiles::addPathAndCheck(std::string & fileName, const bool possiblyNumber)
{
  std::string error = "";
  if(fileName != "" && (possiblyNumber == false || NRLib2::IsNumber(fileName) == false)) {
    fileName = inputDirectory_+fileName;
    std::ifstream infile;
    try {
      NRLib2::OpenRead(infile, fileName);
      infile.close();
    }
    catch (NRLib2::Exception e) {
      error = e.what();
      error += "\n";
    }
  }
  return(error);
}
