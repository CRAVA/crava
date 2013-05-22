/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelgravitystatic.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/analyzelog.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/welldata.h"
#include "src/blockedlogs.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/waveletfilter.h"
#include "src/tasklist.h"


ModelGravityStatic::ModelGravityStatic(ModelSettings        *& modelSettings,
                                       ModelGeneral         *& modelGeneral,
                                       const InputFiles      * inputFiles)
{
  failed_                 = false;
  before_injection_start_ = false; // When do we know what this should be??

  bool failedLoadingModel = false; // Usikker på om trenger to slike failed-variable...
  bool failedReadingFile  = false;
  std::string errText("");

  // Check that there is more than one gravimetric surveys (if any)
  int numberOfVintages = modelSettings->getNumberOfVintages();
  int numberOfGravitySurveys = 0;
  for(int i = 0; i < numberOfVintages; i++)
  {
    if(modelSettings->getGravityTimeLapse(i)){
      numberOfGravitySurveys++;
    }
  }
  if(numberOfGravitySurveys == 1){
    errText+="Not valid with only one gravity survey. Need at least two gravity surveys.";
    failedLoadingModel = true;
  }

  // Read first gravimetric data file
  if(numberOfGravitySurveys > 0 && failedLoadingModel == false){
    // Find first gravity data file
    std::string fileName = inputFiles->getGravimetricData(0);

    int nObs = 30;   // should this be given in input file?
    int nColumns = 5;  // We require data files to have five columns

    observation_location_utmx_.resize(nObs);
    observation_location_utmy_.resize(nObs);
    observation_location_depth_.resize(nObs);
    gravity_response_.resize(nObs);
    gravity_std_dev_.resize(nObs);

    readGravityDataFile(fileName, "gravimetric base survey",
                        nObs, nColumns,
                        observation_location_utmx_,
                        observation_location_utmy_,
                        observation_location_depth_,
                        gravity_response_,
                        gravity_std_dev_,
                        failedReadingFile,
                        errText);
    failedLoadingModel = failedReadingFile;
  }

  if (failedLoadingModel) {
    LogKit::WriteHeader("Error(s) with gravimetric surveys");
    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel || failedReadingFile;
  failed_details_.push_back(failedReadingFile);
 // failed_details_.push_back(failedPriorFacies);
}


ModelGravityStatic::~ModelGravityStatic(void)
{
}

void
ModelGravityStatic::readGravityDataFile(const std::string   & fileName,
                                        const std::string   & readReason,
                                        int                   nObs,
                                        int                   nColumns,
                                        std::vector <float> & obs_loc_utmx,
                                        std::vector <float> & obs_loc_utmy,
                                        std::vector <float> & obs_loc_depth,
                                        std::vector <float> & gravity_response,
                                        std::vector <float> & gravity_std_dev,
                                        bool                  failed,
                                        std::string         & errText)
{
  float * tmpRes = new float[nObs*nColumns+1];
  std::ifstream inFile;
  NRLib::OpenRead(inFile,fileName);
  std::string text = "Reading "+readReason+" from file "+fileName+" ... ";
  LogKit::LogFormatted(LogKit::Low,text);
  std::string storage;
  int index = 0;
  failed = false;

  while(failed == false && inFile >> storage) {
    if(index < nObs*nColumns) {
      try {
        tmpRes[index] = NRLib::ParseType<float>(storage);
      }
      catch (NRLib::Exception & e) {
        errText += "Error in "+fileName+"\n";
        errText += e.what();
        failed = true;
      }
    }
    index++;
  }
  if(failed == false) {
    if(index != nObs*nColumns) {
      failed = true;
      errText += "Found "+NRLib::ToString(index)+" in file "+fileName+", expected "+NRLib::ToString(nObs*nColumns)+".\n";
    }
  }

  if(failed == false) {
    LogKit::LogFormatted(LogKit::Low,"ok.\n");
    index = 0;
    for(int i=0;i<nObs;i++) {
      obs_loc_utmx[i] = tmpRes[index];
      index++;
      obs_loc_utmy[i] = tmpRes[index];
      index++;
      obs_loc_depth[i] = tmpRes[index];
      index++;
      gravity_response[i] = tmpRes[index];
      index++;
      gravity_std_dev[i] = tmpRes[index];
      index++;
    }
  }
  else{
    failed = true;
    LogKit::LogFormatted(LogKit::Low,"failed.\n");
  }
  delete [] tmpRes;
}
