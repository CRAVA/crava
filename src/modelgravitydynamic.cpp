/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelgravitystatic.h"
#include "src/modelgravitydynamic.h"
#include "src/xmlmodelfile.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/io.h"
#include "src/tasklist.h"
#include "src/seismicparametersholder.h"

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/volume/volume.hpp"


ModelGravityDynamic::ModelGravityDynamic(const ModelSettings          * modelSettings,
                                         const ModelGeneral           * modelGeneral,
                                         const InputFiles             * inputFiles,
                                         int                            t,
                                         SeismicParametersHolder      & seismicParameters)

{ //Time lapse constructor
  failed_                 = false;
  thisTimeLapse_          = t;

  bool failedLoadingModel = false;
  bool failedReadingFile = false;
  std::string errText("");

  int nObs = 30;     // should this be given in input file
  int nColumns = 5;  // We require data files to have five columns

  // Check that timeLapse is ok
  if(thisTimeLapse_ < 1 && thisTimeLapse_ >modelSettings->getNumberOfVintages()){
    errText += "Not valid time lapse";
    failedLoadingModel = true;
  }

  if(failedLoadingModel == false){
  // Find first gravity data file
  std::string fileName = inputFiles->getGravimetricData(thisTimeLapse_);

  observation_location_utmx_.resize(nObs);
  observation_location_utmy_.resize(nObs);
  observation_location_depth_.resize(nObs);
  gravity_response_.resize(nObs);
  gravity_std_dev_.resize(nObs);

  ModelGravityStatic::readGravityDataFile(fileName,
                                          "gravimetric survey time lapse "+thisTimeLapse_,
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
    LogKit::WriteHeader("Error(s) with gravimetric surveys.");
    LogKit::LogFormatted(LogKit::Error,"\n"+errText);
    LogKit::LogFormatted(LogKit::Error,"\nAborting\n");
  }

  failed_ = failedLoadingModel || failedReadingFile;
  failed_details_.push_back(failedReadingFile);
  }

ModelGravityDynamic::~ModelGravityDynamic(void)
{
}

