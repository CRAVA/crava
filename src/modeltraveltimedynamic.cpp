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
#include "src/rmstrace.h"

#include "nrlib/iotools/fileio.hpp"
#include "nrlib/surface/regularsurface.hpp"

ModelTravelTimeDynamic::ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                                               const ModelGeneral            * modelGeneral,
                                               const InputFiles              * inputFiles,
                                               const int                     & vintage)
: thisTimeLapse_(vintage)
{
  std::string errTxt = "";

  bool failed_surfaces = false;
  processHorizons(horizons_,
                  inputFiles,
                  errTxt,
                  failed_surfaces);

  bool failed_rms_data = false;
  processRMSData(modelSettings,
                 inputFiles,
                 modelGeneral->getTimeSimbox(),
                 errTxt,
                 failed_rms_data);

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
ModelTravelTimeDynamic::processRMSData(const ModelSettings      * modelSettings,
                                       const InputFiles         * inputFiles,
                                       const Simbox             * timeSimbox,
                                       std::string              & errTxt,
                                       bool                     & failed)

{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  LogKit::WriteHeader("Reading RMS travel time data");

  const std::string & file_name  = inputFiles->getRmsVelocities(thisTimeLapse_);
  std::string         tmpErrText = "";

  readRMSData(file_name, tmpErrText);

  int n_rms_traces     = static_cast<int>(rms_traces_.size());
  int n_layers_above   = modelSettings->getRMSnLayersAbove();
  int n_layers_below   = modelSettings->getRMSnLayersBelow();
  int n_layers_simbox  = timeSimbox->getnz();
  int n_layers_padding = modelSettings->getNZpad();
  int n_layers         = n_layers_above + n_layers_padding + n_layers_below;

  double t_top;
  double t_bot;
  double dt_simbox;

  for(int i=0; i<n_rms_traces; i++) {

    getCoordinates(timeSimbox,
                   rms_traces_[i],
                   t_top,
                   t_bot,
                   dt_simbox);

    NRLib::Grid2D<double> G = calculateG(rms_traces_[i].getTime(),
                                         t_top,
                                         t_bot,
                                         dt_simbox,
                                         n_layers,
                                         n_layers_above,
                                         n_layers_below,
                                         n_layers_simbox,
                                         n_layers_padding);
  }

  if(tmpErrText != "") {
    errTxt += tmpErrText;
    failed = true;
  }


  Timings::setTimeSeismic(wall,cpu);
}

//----------------------------------------------------------------------------
void
ModelTravelTimeDynamic::readRMSData(const std::string & fileName,
                                    std::string       & errTxt)
{

  int error = 0;

  std::ifstream file;
  NRLib::OpenRead(file, fileName);

  if(file == 0) {
    error = 1;
    errTxt += "Could not open RMS data file "+fileName+" for reading.\n";
  }


  int line = 0;

  while(line < 32) {
    NRLib::DiscardRestOfLine(file,line,false);
  }

  std::string endWord = "::Goodbye::";
  std::string token   = "";

  int j       = 0;
  int IL      = IMISSING;
  int XL      = IMISSING;
  double utmx = RMISSING;
  double utmy = RMISSING;

  std::vector<double> time;
  std::vector<double> velocity;

  while(NRLib::CheckEndOfFile(file) == false) {

    try {

      NRLib::ReadNextToken(file,token,line);

      if(token != endWord) {
        int    this_il       = NRLib::ParseType<int>(token);
        int    this_xl       = NRLib::ReadNext<int>(file,line);
        double this_utmx     = NRLib::ReadNext<double>(file,line);
        double this_utmy     = NRLib::ReadNext<double>(file,line);
        double this_time     = NRLib::ReadNext<double>(file,line);
        double this_velocity = NRLib::ReadNext<double>(file,line);


        if(time.size() > 1) {

          if( (IL != this_il) || (XL != this_xl) ) {
            RMSTrace trace(IL,
                           XL,
                           utmx,
                           utmy,
                           time,
                           velocity);

            rms_traces_.push_back(trace);

            time.clear();
            velocity.clear();

            j = -1;
          }
        }

        IL = this_il;
        XL = this_xl;
        utmx = this_utmx;
        utmy = this_utmy;
        time.push_back(this_time);
        velocity.push_back(this_velocity);

        j++;
      }

      else {
        RMSTrace trace(IL,
                       XL,
                       utmx,
                       utmy,
                       time,
                       velocity);

        rms_traces_.push_back(trace);

        break;
      }

    }
    catch (NRLib::IOError e) {
      std::string text;
      text += std::string("\nERROR: Reading of RMS data \'") + fileName + "\' failed.\n";
      text += std::string("\nERROR message is \'") + e.what() + "\'";
      LogKit::LogMessage(LogKit::Error,text);
      errTxt += text;
      exit(1);
    }
  }

  file.close();
  file.clear();

  int n_traces = static_cast<int>(rms_traces_.size());

  LogKit::LogFormatted(LogKit::Low, "\nRead "+NRLib::ToString(n_traces)+" RMS traces in file "+fileName+"\n");

}

double
ModelTravelTimeDynamic::findMaxTime() const
{

  int n_rms_traces = static_cast<int>(rms_traces_.size());

  double max_time = 0;

  for(int i=0; i<n_rms_traces; i++) {
    const std::vector<double> & rms_time = rms_traces_[i].getTime();
    double max = rms_time[rms_time.size()-1];

    if(max > max_time)
      max_time = max;
  }

  return max_time;
}

NRLib::Grid2D<double>
ModelTravelTimeDynamic::calculateG(const std::vector<double> & rms_time,
                                   const double              & t_top,
                                   const double              & t_bot,
                                   const double              & dt_simbox,
                                   const int                 & n_layers,
                                   const int                 & n_layers_above,
                                   const int                 & n_layers_below,
                                   const int                 & n_layers_simbox,
                                   const int                 & n_layers_padding) const
{
  std::vector<double> t(n_layers + 1, 0);
  std::vector<double> dt(n_layers + 1, 0);

  double max_time = findMaxTime();

  double dt_above  = static_cast<double>(  t_top             / n_layers_above);
  double dt_below  = static_cast<double>( (max_time - t_bot) / n_layers_below);

  for(int j=0; j<n_layers + 1; j++) {
    if(j < n_layers_above) {
      t[j]  = j * dt_above;
      dt[j] = dt_above;
    }
    else if(j >= n_layers_above && j < n_layers_above + n_layers_simbox) {
      t[j]  = t_top + (j - n_layers_above) * dt_simbox;
      dt[j] = dt_simbox;
    }
    else if(j >= n_layers_above + n_layers_padding) {
      t[j]  = t_bot + (j - n_layers_above - n_layers_padding) * dt_below;
      dt[j] = dt_below;
    }
  }

  int n_rms_data = static_cast<int>(rms_time.size());

  NRLib::Grid2D<double> G(n_rms_data, n_layers, 0);

  for(int j=0; j<n_rms_data; j++) {
    int k=0;
    while(rms_time[j] >= t[k] && k < n_layers) {
      G(j,k) = dt[k] / rms_time[j];
      k++;
    }
    if(k < n_layers)
      G(j,k) = (rms_time[j] - t[k-1]) / rms_time[j];
  }

  return G;
}

void
ModelTravelTimeDynamic::getCoordinates(const Simbox   * timeSimbox,
                                       const RMSTrace & rms_trace,
                                       double         & t_top,
                                       double         & t_bot,
                                       double         & dt_simbox) const
{

  const double x = rms_trace.getUtmx();
  const double y = rms_trace.getUtmy();

  int i_ind;
  int j_ind;

  timeSimbox->getIndexes(x, y, i_ind, j_ind);

  t_top     = timeSimbox->getTop(i_ind, j_ind);
  t_bot     = timeSimbox->getBot(i_ind, j_ind);
  dt_simbox = timeSimbox->getdz(i_ind, j_ind);
}
