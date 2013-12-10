/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/modeltraveltimestatic.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/simbox.h"
#include "src/definitions.h"

ModelTravelTimeStatic::ModelTravelTimeStatic(const ModelSettings * modelSettings,
                                             const InputFiles    * inputFiles,
                                             const Simbox        * timeSimbox)
: initial_horizons_(0),
  horizon_names_(0, ""),
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
  lz_limit_(RMISSING),
  format_flag_(modelSettings->getOutputGridFormat()),
  simbox_above_(NULL)
{
  std::string errTxt = "";

  bool failed_surfaces = false;
  bool failed_rms_data = false;

  bool doTravelTimeInversion = false;
  if (modelSettings->getTravelTimeTimeLapse(0) == true)
    doTravelTimeInversion = true;


  if (doTravelTimeInversion == true) {

    LogKit::WriteHeader("Reading Static Travel Time Data");

    processHorizons(modelSettings,
                    inputFiles,
                    errTxt,
                    failed_surfaces);

    processRMSData(modelSettings,
                   timeSimbox,
                   errTxt,
                   failed_rms_data);
  }

  bool failed_loading_model = failed_surfaces || failed_rms_data;

  if (failed_loading_model) {
    LogKit::WriteHeader("Error(s) while loading static travel time data");
    LogKit::LogMessage(LogKit::Error,"\n"+errTxt);
    LogKit::LogMessage(LogKit::Error,"\nAborting\n");
  }

  failed_ = failed_loading_model;
  failed_details_.push_back(failed_surfaces);
  failed_details_.push_back(failed_rms_data);
}

 //-------------------------------------------------------------------------------------------//

ModelTravelTimeStatic::~ModelTravelTimeStatic()
{
  delete simbox_above_;
}

//-------------------------------------------------------------------------------------------//

void
ModelTravelTimeStatic::processHorizons(const ModelSettings  * modelSettings,
                                       const InputFiles     * inputFiles,
                                       std::string          & errTxt,
                                       bool                 & failed)
{
  std::string tmpErrText = "";

  const std::vector<std::string> & initial_horizons = inputFiles->getTravelTimeHorizons(0);

  int n_horizons = static_cast<int>(initial_horizons.size());

  if (n_horizons > 0) {

    LogKit::LogFormatted(LogKit::Low, "\nReading horizon travel time data:\n");

    horizon_names_ = modelSettings->getTimeLapseTravelTimeHorizons(0);

    int n_horizons = static_cast<int>(horizon_names_.size());

    initial_horizons_.resize(n_horizons);

    for (int i = 0; i < n_horizons; i++) {
      LogKit::LogFormatted(LogKit::Low, "\nHorizon \""+horizon_names_[i]+"\":\n");

      LogKit::LogFormatted(LogKit::Low, "  Reading horizon file "+initial_horizons[i]+"\n");
      initial_horizons_[i] = Surface(initial_horizons[i]);

    }
  }

  if (tmpErrText != "") {
    errTxt += tmpErrText;
    failed = true;
  }
}

//-------------------------------------------------------------------------------------------//

void
ModelTravelTimeStatic::processRMSData(const ModelSettings * modelSettings,
                                      const Simbox        * timeSimbox,
                                      std::string         & errTxt,
                                      bool                & failed)

{
  std::string tmpErrText = "";

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

  if (tmpErrText != "") {
    errTxt += tmpErrText;
    failed = true;
  }
}

//----------------------------------------------------------------------------

void
ModelTravelTimeStatic::setupSimboxAbove(const Simbox  * timeSimbox,
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
