/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <fstream>

#include "src/inputfiles.h"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/exception/exception.hpp"

InputFiles::InputFiles(void)
  : seedFile_(""),
    wellFiles_(0),
    wellMoveIntFile_(2),
    faciesEstIntFile_(2),
    timeSurfFiles_(0),
    depthSurfFiles_(2),
    velocityField_(""),
    backFile_(3),
    backVelFile_(""),
    reflMatrFile_(""),
    corrDirFile_(""),
    corrDirTopFile_(""),
    corrDirBaseFile_(""),
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
  int          t;
  int          nTimeLapse           = static_cast<int>(timeLapseSeismicFiles_.size());
  int          nTimeLapseTravelTime = static_cast<int>(rmsVelocities_.size());

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
  for(t=0; t<nTimeLapse; t++){
    for(i=0;i<timeLapseSeismicFiles_[t].size();i++)
      errTxt += addPathAndCheck(timeLapseSeismicFiles_[t][i]);
    for(i=0;i<timeLapseWaveletFiles_[t].size();i++)
      errTxt += addPathAndCheck(timeLapseWaveletFiles_[t][i]);
    for(i=0;i<timeLapseWaveletShiftFiles_[t].size(); i++)
      errTxt += addPathAndCheck(timeLapseWaveletShiftFiles_[t][i]);
    for(i=0;i<timeLapseWaveletScaleFiles_[t].size(); i++)
      errTxt += addPathAndCheck(timeLapseWaveletScaleFiles_[t][i]);
    for(i=0;i<timeLapseLocalNoiseFiles_[t].size();i++)
      errTxt += addPathAndCheck(timeLapseLocalNoiseFiles_[t][i]);
  }
  for(i=0; i<travelTimeHorizons_.size(); i++)
    errTxt += addPathAndCheck(travelTimeHorizons_[i]);
  for(i=0; i<rmsVelocities_.size(); i++)
    errTxt += addPathAndCheck(rmsVelocities_[i]);
  for(t=0; t<nTimeLapseTravelTime; t++) {
    for(i=0; i<timeLapseTravelTimeHorizons_[t].size(); i++)
      errTxt += addPathAndCheck(timeLapseTravelTimeHorizons_[t][i]);
  }
  for(i=0; i<gravimetricLocations_.size(); i++)
    errTxt += addPathAndCheck(gravimetricLocations_[i]);
  for(i=0; i<gravimetricData_.size(); i++)
    errTxt += addPathAndCheck(gravimetricData_[i]);
  for(i=0;i<waveletFilterFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletFilterFiles_[i]);
  for(i=0;i<waveletCorrFiles_.size(); i++)
    errTxt += addPathAndCheck(waveletCorrFiles_[i]);
  for(i=0;i<waveletEstIntFileTop_.size();i++)
    errTxt += addPathAndCheck(waveletEstIntFileTop_[i]);
  for(i=0;i<waveletEstIntFileBase_.size();i++)
    errTxt += addPathAndCheck(waveletEstIntFileBase_[i]);
  for(i=0;i<faciesEstIntFile_.size();i++)
    errTxt += addPathAndCheck(faciesEstIntFile_[i]);
  for(i=0;i<wellMoveIntFile_.size();i++)
    errTxt += addPathAndCheck(wellMoveIntFile_[i], true);
  for(i=0;i<timeSurfFiles_.size();i++)
    errTxt += addPathAndCheck(timeSurfFiles_[i], true);
  for(i=0;i<depthSurfFiles_.size();i++)
    errTxt += addPathAndCheck(depthSurfFiles_[i]);
  for(i=0;i<multizoneSurfaceFiles_.size();i++)
    errTxt += addPathAndCheck(multizoneSurfaceFiles_[i]);
  errTxt += addPathAndCheck(velocityField_);
  for(i=0;i<backFile_.size();i++)
    errTxt += addPathAndCheck(backFile_[i]);

  std::map<std::string, std::string>::iterator j;
  for(j=interval_base_time_surface_.begin();j != interval_base_time_surface_.end();j++) {
    errTxt += addPathAndCheck(j->second, true);
  }
  for(j=interval_base_depth_surface_.begin();j != interval_base_depth_surface_.end();j++) {
    errTxt += addPathAndCheck(j->second, true);
  }
  for(j=interval_corrDirFiles_.begin();j != interval_corrDirFiles_.end();j++) {
    errTxt += addPathAndCheck(j->second);
  }
  for(j=interval_corrDirTopFiles_.begin();j != interval_corrDirTopFiles_.end();j++) {
    errTxt += addPathAndCheck(j->second);
  }
  for(j=interval_corrDirBaseFiles_.begin();j != interval_corrDirBaseFiles_.end();j++) {
    errTxt += addPathAndCheck(j->second);
  }

  errTxt += addPathAndCheck(backVelFile_);
  errTxt += addPathAndCheck(corrDirFile_);
  errTxt += addPathAndCheck(corrDirTopFile_);
  errTxt += addPathAndCheck(corrDirBaseFile_);
  errTxt += addPathAndCheck(reflMatrFile_);
  errTxt += addPathAndCheck(paramCorrFile_);
  errTxt += addPathAndCheck(tempCorrFile_);
  errTxt += addPathAndCheck(refSurfaceFile_,true);
  errTxt += addPathAndCheck(areaSurfaceFile_);

  for(j=priorFaciesProb_.begin();j != priorFaciesProb_.end();j++) {
    errTxt += addPathAndCheck(j->second);
  }
  for(i=0;i<trendCubes_.size();i++)
    errTxt += addPathAndCheck(trendCubes_[i]);

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
    catch (NRLib::Exception & e) {
      error = e.what();
      error += "\n";
    }
  }
  return(error);
}
