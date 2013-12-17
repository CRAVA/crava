/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#include <cmath>

#include "lib/timekit.hpp"

#include "src/modelgeneral.h"
#include "src/modeltraveltimedynamic.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/seismicparametersholder.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/gridmapping.h"
#include "src/timings.h"

#include "nrlib/iotools/fileio.hpp"
#include "nrlib/surface/regularsurface.hpp"

ModelTravelTimeDynamic::ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                                               const ModelGeneral            * modelGeneral,
                                               const InputFiles              * inputFiles,
                                               const int                     & vintage)
: rms_data_(NULL),
  thisTimeLapse_(vintage)
{
  std::string errTxt = "";

  bool failed_surfaces = false;
  processHorizons(horizons_,
                  inputFiles,
                  errTxt,
                  failed_surfaces);

  bool failed_rms_data = false;
  //processRMSData(rms_data_, //H Commented out since modelGeneral doesn't have TimeCutMapping. Fix when setting up TravelTimeDynamic
  //               modelSettings,
  //               inputFiles,
  //               modelGeneral->getTimeSimbox(),
  //               modelGeneral->getTimeDepthMapping(),
  //               modelGeneral->getTimeCutMapping(),
  //               errTxt,
  //               failed_rms_data);

  bool failed_loading_model = failed_surfaces || failed_rms_data;

  if (failed_loading_model) {
    LogKit::WriteHeader("Error(s) while loading travel time data");
    LogKit::LogMessage(LogKit::Error,"\n"+errTxt);
    LogKit::LogMessage(LogKit::Error,"\nAborting\n");
  }

  failed_ = failed_loading_model;
  failed_details_.push_back(failed_surfaces);
  failed_details_.push_back(failed_rms_data);
}

ModelTravelTimeDynamic::~ModelTravelTimeDynamic()
{
  delete rms_data_;
}

void
ModelTravelTimeDynamic::processHorizons(std::vector<Surface>   & horizons,
                                        const InputFiles       * inputFiles,
                                        std::string            & errTxt,
                                        bool                   & failed)
{

  const std::vector<std::string> & travel_time_horizons = inputFiles->getTravelTimeHorizons(thisTimeLapse_);

  int n_horizons = static_cast<int>(travel_time_horizons.size());

  if(n_horizons == 1) {
    if(travel_time_horizons[0] != "") {
      errTxt += "Only one surface is given for inversion of the horizons in the travel time data. At least two surfaces should be given\n";
      failed = true;
    }
  }

  else {
    horizons.resize(n_horizons);
    for(int i=0; i<n_horizons; i++)
      horizons[i] = Surface(travel_time_horizons[i]);
  }

}

void
ModelTravelTimeDynamic::processRMSData(FFTGrid                 *& rms_data,
                                       const ModelSettings      * modelSettings,
                                       const InputFiles         * inputFiles,
                                       const Simbox             * timeSimbox,
                                       GridMapping              * timeDepthMapping,
                                       const GridMapping        * timeCutMapping,
                                       std::string              & errTxt,
                                       bool                     & failed)

{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  LogKit::WriteHeader("Reading RMS travel time data");

  const SegyGeometry * geometry = new const SegyGeometry;
  rms_data                      = new FFTGrid;

  float offset = modelSettings->getTravelTimeSegyOffset(thisTimeLapse_);
  if(offset < 0)
    offset = modelSettings->getSegyOffset(thisTimeLapse_);

  const Simbox * timeCutSimbox = NULL;
  if (timeCutMapping != NULL)
    timeCutSimbox = timeCutMapping->getSimbox(); // For the got-enough-data test
  else
    timeCutSimbox = timeSimbox;

  const std::string & file_name = inputFiles->getRmsVelocities(thisTimeLapse_);
  std::string         data_name = "RMS data";
  std::string         tmpErrText = "";

  //ModelGeneral::readGridFromFile(file_name,
  //                               data_name,
  //                               offset,
  //                               rms_data,
  //                               geometry,
  //                               modelSettings->getTravelTimeTraceHeaderFormat(thisTimeLapse_),
  //                               FFTGrid::DATA,
  //                               timeSimbox,
  //                               timeCutSimbox,
  //                               modelSettings,
  //                               tmpErrText);

  if(tmpErrText != "") {
    tmpErrText += "\nReading of file \'"+file_name+"\' for "+data_name+" failed.\n";
    errTxt += tmpErrText;
    failed = true;
  }


  LogKit::LogFormatted(LogKit::Low,"\n");

  if(failed == false) {
    bool segyVolumesRead = false;

    if (geometry != NULL)
      segyVolumesRead = true;

    if (segyVolumesRead == true) {
      LogKit::LogFormatted(LogKit::Low,"\nArea/resolution           x0           y0            lx         ly     azimuth         dx      dy\n");
      LogKit::LogFormatted(LogKit::Low,"-------------------------------------------------------------------------------------------------\n");

      if (geometry != NULL) {
        double geoAngle = (-1)*timeSimbox->getAngle()*(180/M_PI);
        if (geoAngle < 0)
          geoAngle += 360.0;
        LogKit::LogFormatted(LogKit::Low,"RMS travel time data   %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                             geometry->GetX0(), geometry->GetY0(),
                             geometry->Getlx(), geometry->Getly(), geoAngle,
                             geometry->GetDx(), geometry->GetDy());
      }
    }

    if((modelSettings->getOutputGridsOther() & IO::RMS_VELOCITIES) > 0) {
      std::string baseName = IO::PrefixTravelTimeData();
      std::string sgriLabel = std::string("RMS travel time data");

      rms_data->writeFile(baseName,
                          IO::PathToTravelTimeData(),
                          timeSimbox,
                          sgriLabel,
                          offset,
                          timeDepthMapping,
                          timeCutMapping,
                          *modelSettings->getTraceHeaderFormatOutput());
    }

    if (geometry != NULL)
      delete geometry;
  }

  Timings::setTimeSeismic(wall,cpu);
}
