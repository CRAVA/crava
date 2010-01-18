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
    wellMoveIntFile_(2),   
    faciesEstIntFile_(2),   
    timeSurfFiles_(0),       
    depthSurfFiles_(2),      
    velocityField_(""),
    backFile_(3),           
    backVelFile_(""),       
    reflMatrFile_(""),
    corrDirFile_(""),
    paramCorrFile_(""),
    tempCorrFile_(""),
    areaSurfaceFile_("")
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
  for(i=0;i<localNoiseFiles_.size();i++)
    errTxt += addPathAndCheck(localNoiseFiles_[i]);
  for(i=0;i<waveletFiles_.size();i++)
    errTxt += addPathAndCheck(waveletFiles_[i]);
  for(i=0;i<waveletShiftFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletShiftFiles_[i]);
  for(i=0;i<waveletScaleFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletScaleFiles_[i]);
  for(i=0;i<waveletFilterFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletFilterFiles_[i]);
  for(i=0;i<waveletCorrFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletCorrFiles_[i]);
  for(i=0;i<waveletEstIntFile_.size();i++)
    errTxt += addPathAndCheck(waveletEstIntFile_[i]);
  for(i=0;i<faciesEstIntFile_.size();i++)
    errTxt += addPathAndCheck(faciesEstIntFile_[i]);
  for(i=0;i<wellMoveIntFile_.size();i++)
    errTxt += addPathAndCheck(wellMoveIntFile_[i], true);
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
  errTxt += addPathAndCheck(refSurfaceFile_,true);
  errTxt += addPathAndCheck(areaSurfaceFile_);

  std::map<std::string, std::string>::iterator j;
  for(j=priorFaciesProb_.begin();j != priorFaciesProb_.end();j++) {
    errTxt += addPathAndCheck(j->second);
  }

  return(errTxt);
}


std::string
InputFiles::addPathAndCheck(std::string & fileName, const bool possiblyNumber)
{
  std::string error = "";
  if(fileName != "" && (possiblyNumber == false || NRLib::IsNumber(fileName) == false)) {
    fileName = inputDirectory_+fileName;
    std::ifstream infile;
    try {
      NRLib::OpenRead(infile, fileName);
      infile.close();
    }
    catch (NRLib::Exception e) {
      error = e.what();
      error += "\n";
    }
  }
  return(error);
}
