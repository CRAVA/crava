/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "lib/timekit.hpp"
#include "src/modeltraveltimedynamic.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/rmstrace.h"
#include "src/simbox.h"
#include "src/definitions.h"

ModelTravelTimeDynamic::ModelTravelTimeDynamic(const ModelSettings * modelSettings,
                                               const InputFiles    * inputFiles,
                                               const Simbox        * timeSimbox,
                                               const int           & vintage)
: push_down_horizons_(0),
  initial_horizons_(0),
  horizon_names_(0, ""),
  horizon_standard_deviation_(0),
  rms_traces_(0, NULL),
  rms_standard_deviation_(RMISSING),
  n_layers_above_(0),
  n_layers_below_(0),
  mean_vp_top_(RMISSING),
  mean_vp_base_(RMISSING),
  var_vp_above_(RMISSING),
  var_vp_below_(RMISSING),
  range_above_(RMISSING),
  range_below_(RMISSING),
  failed_(false),
  failed_details_(0),
  this_time_lapse_(vintage),
  lz_limit_(RMISSING),
  simbox_above_(NULL),
  simbox_below_(NULL),
  rms_data_given_(false),
  horizon_data_given_(false)
{
  std::string errTxt = "";

  bool failed_surfaces = false;
  if (vintage > 0)
    processHorizons(modelSettings,
                    inputFiles,
                    errTxt,
                    failed_surfaces);

  bool failed_rms_data = false;
  processRMSData(modelSettings,
                 inputFiles,
                 timeSimbox,
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

 //-------------------------------------------------------------------------------------------//

ModelTravelTimeDynamic::~ModelTravelTimeDynamic()
{
  for (size_t i = 0; i < rms_traces_.size(); i++)
    delete rms_traces_[i];

  delete simbox_above_;
  delete simbox_below_;
}

//-------------------------------------------------------------------------------------------//

void
ModelTravelTimeDynamic::processHorizons(const ModelSettings  * modelSettings,
                                        const InputFiles     * inputFiles,
                                        std::string          & errTxt,
                                        bool                 & failed)
{
  std::string tmpErrText = "";

  const std::vector<std::string> & push_down_horizons = inputFiles->getTravelTimeHorizons(this_time_lapse_);
  const std::vector<std::string> & initial_horizons   = inputFiles->getTravelTimeHorizons(0);

  int n_horizons = static_cast<int>(push_down_horizons.size());

  if (n_horizons > 0) {
    LogKit::WriteHeader("Reading horizon travel time data");

    horizon_data_given_ = true;

    horizon_names_              = modelSettings->getTimeLapseTravelTimeHorizons(this_time_lapse_);
    horizon_standard_deviation_ = modelSettings->getTimeLapseTravelTimeHorizonSD(this_time_lapse_);

    initial_horizons_.resize(n_horizons);
    push_down_horizons_.resize(n_horizons);

    for (int i = 0; i < n_horizons; i++) {
      LogKit::LogFormatted(LogKit::Low, "\nHorizon \""+horizon_names_[i]+"\":\n");

      LogKit::LogFormatted(LogKit::Low, "  Reading horizon file file "+initial_horizons[i]+"\n");
      initial_horizons_[i] = Surface(initial_horizons[i]); //Finn riktig fil, trenger ikke være gitt i samme rekkefølge

      LogKit::LogFormatted(LogKit::Low, "  Reading push down file "+push_down_horizons[i]+"\n");
      push_down_horizons_[i] = Surface(push_down_horizons[i]);
    }
  }

  if (tmpErrText != "") {
    errTxt += tmpErrText;
    failed = true;
  }

}

//-------------------------------------------------------------------------------------------//

void
ModelTravelTimeDynamic::processRMSData(const ModelSettings * modelSettings,
                                       const InputFiles    * inputFiles,
                                       const Simbox        * timeSimbox,
                                       std::string         & errTxt,
                                       bool                & failed)

{
  std::string tmpErrText = "";

  const std::string & file_name  = inputFiles->getRmsVelocities(this_time_lapse_);

  if (file_name != "") {
    LogKit::WriteHeader("Reading RMS travel time data");

    rms_data_given_ = true;

    readRMSData(file_name, timeSimbox, tmpErrText);

    rms_standard_deviation_ = modelSettings->getRMSStandardDeviation();

    n_layers_above_         = modelSettings->getRMSnLayersAbove();
    n_layers_below_         = modelSettings->getRMSnLayersBelow();

    mean_vp_top_            = modelSettings->getRMSMeanVpTop();
    mean_vp_base_           = modelSettings->getRMSMeanVpBase();

    var_vp_above_           = modelSettings->getRMSVarianceVpAbove();
    var_vp_below_           = modelSettings->getRMSVarianceVpBelow();

    range_above_            = static_cast<float>(modelSettings->getRMSTemporalCorrelationRangeAbove());
    range_below_            = static_cast<float>(modelSettings->getRMSTemporalCorrelationRangeBelow());

    lz_limit_               = modelSettings->getLzLimit();

    setupSimboxAbove(timeSimbox,
                     modelSettings->getOutputGridFormat(),
                     modelSettings->getOutputGridDomain(),
                     modelSettings->getOtherOutputFlag(),
                     lz_limit_,
                     tmpErrText);

    setupSimboxBelow(timeSimbox,
                     modelSettings->getOutputGridFormat(),
                     modelSettings->getOutputGridDomain(),
                     modelSettings->getOtherOutputFlag(),
                     tmpErrText);

  }

  if (tmpErrText != "") {
    errTxt += tmpErrText;
    failed = true;
  }
}

//----------------------------------------------------------------------------
void
ModelTravelTimeDynamic::readRMSData(const std::string & fileName,
                                    const Simbox      * timeSimbox,
                                    std::string       & errTxt)
{

  int error = 0;

  std::ifstream file;
  NRLib::OpenRead(file, fileName);

  if (file == 0) {
    error = 1;
    errTxt += "Could not open RMS data file "+fileName+" for reading.\n";
  }


  int line = 0;

  while (line < 32) {
    NRLib::DiscardRestOfLine(file, line, false);
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

  while (NRLib::CheckEndOfFile(file) == false) {

    try {

      NRLib::ReadNextToken(file, token, line);

      if (token != endWord) {
        int    this_il       = NRLib::ParseType<int>(token);
        int    this_xl       = NRLib::ReadNext<int>(file, line);
        double this_utmx     = NRLib::ReadNext<double>(file, line);
        double this_utmy     = NRLib::ReadNext<double>(file, line);
        double this_time     = NRLib::ReadNext<double>(file, line);
        double this_velocity = NRLib::ReadNext<double>(file, line);


        if (time.size() > 1) {

          if ((IL != this_il) || (XL != this_xl)) {

            int i_ind;
            int j_ind;
            timeSimbox->getIndexes(this_utmx, this_utmy, i_ind, j_ind);

            if (i_ind != IMISSING && j_ind != IMISSING) {
              RMSTrace * trace = new RMSTrace(IL,
                                              XL,
                                              i_ind,
                                              j_ind,
                                              utmx,
                                              utmy,
                                              time,
                                              velocity);

              rms_traces_.push_back(trace);
            }

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

        int i_ind;
        int j_ind;
        timeSimbox->getIndexes(utmx, utmy, i_ind, j_ind);

        if (i_ind != IMISSING && j_ind != IMISSING) {
          RMSTrace * trace = new RMSTrace(IL,
                                          XL,
                                          i_ind,
                                          j_ind,
                                          utmx,
                                          utmy,
                                          time,
                                          velocity);

          rms_traces_.push_back(trace);

          break;
        }
      }

    }
    catch (NRLib::IOError e) {
      std::string text;
      text += std::string("\nERROR: Reading of RMS data \'") + fileName + "\' failed.\n";
      text += std::string("\nERROR message is \'") + e.what() + "\'";
      LogKit::LogMessage(LogKit::Error, text);
      errTxt += text;
      exit(1);
    }
  }

  file.close();
  file.clear();

  int n_traces = static_cast<int>(rms_traces_.size());

  LogKit::LogFormatted(LogKit::Low, "\nRead "+NRLib::ToString(n_traces)+" RMS traces in file "+fileName+"\n");

}

//----------------------------------------------------------------------------

void
ModelTravelTimeDynamic::setupSimboxAbove(const Simbox  * timeSimbox,
                                         const int     & outputFormat,
                                         const int     & outputDomain,
                                         const int     & otherOutput,
                                         const double  & lz_limit,
                                         std::string   & errTxt)
{

  std::string tmpErrTxt = "";

  Surface top_simbox_surface(dynamic_cast<const Surface &> (timeSimbox->GetTopSurface()));

  double xmin, xmax, ymin, ymax;

  timeSimbox->getMinAndMaxXY(xmin, xmax, ymin, ymax);

  double lx = xmax-xmin;
  double ly = ymax-ymin;

  int nx = timeSimbox->getnx();
  int ny = timeSimbox->getny();

  NRLib::Grid2D<double> z_grid_top(nx, ny, 0);

  Surface top_surface = Surface(xmin, ymin, lx, ly, z_grid_top);

  //
  // Make new simbox
  //
  simbox_above_ = new Simbox(timeSimbox);
  simbox_above_->setDepth(top_surface, top_simbox_surface, n_layers_above_);

  int status = simbox_above_->calculateDz(lz_limit, tmpErrTxt);

  if (status == Simbox::INTERNALERROR) {
    errTxt += "A problem was encountered for the simbox above the reservoir in the RMS inversion\n";
    errTxt += tmpErrTxt;
  }


  if ((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Above_Reservoir";
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Above_Reservoir";
    simbox_above_->writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToBackground(),
                                    outputFormat);
    simbox_above_->setTopBotName(topSurf, baseSurf, outputFormat);
  }

}

//----------------------------------------------------------------------------

void
ModelTravelTimeDynamic::setupSimboxBelow(const Simbox  * timeSimbox,
                                         const int     & outputFormat,
                                         const int     & outputDomain,
                                         const int     & otherOutput,
                                         std::string   & errTxt)
{
  int n_rms_traces = rms_traces_.size();

  double max_time = 0;

  for (int i = 0; i < n_rms_traces; i++) {
    const std::vector<double> & rms_time = rms_traces_[i]->getTime();
    double max = rms_time[rms_time.size() - 1];

    if (max > max_time)
      max_time = max;
  }

  Surface bot_simbox_surface(dynamic_cast<const Surface &> (timeSimbox->GetBotSurface()));

  double z_max = timeSimbox->getBotZMax();
  double lz    = max_time - z_max;

  //
  // Make new simbox
  //
  simbox_below_ = new Simbox(timeSimbox);

  double dz      = static_cast<double>(lz / n_layers_below_);
  double z_shift = 0;

  simbox_below_->setDepth(bot_simbox_surface, z_shift, lz, dz);

  if ((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Below_Reservoir";
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Below_Reservoir";
    simbox_below_->writeTopBotGrids(topSurf,
                                    baseSurf,
                                    IO::PathToBackground(),
                                    outputFormat);
    simbox_below_->setTopBotName(topSurf, baseSurf, outputFormat);

  }

  if (lz < 0) {
    errTxt += "A problem was encountered for the simbox below the reservoir in the RMS inversion\n";
    errTxt += "There are no RMS data below the reservoir\n";
  }
}
