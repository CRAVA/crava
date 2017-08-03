/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/modeltraveltimestatic.h"
#include "src/modeltraveltimedynamic.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/rmstrace.h"
#include "src/simbox.h"
#include "src/definitions.h"

ModelTravelTimeDynamic::ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                                               const InputFiles              * inputFiles,
                                               const ModelTravelTimeStatic   * modelTravelTimeStatic,
                                               const Simbox                  * timeSimbox,
                                               const int                     & vintage)
: push_down_horizons_(0),
  horizon_names_(0, ""),
  horizon_standard_deviation_(0),
  rms_traces_(0, NULL),
  rms_standard_deviation_(RMISSING),
  failed_(false),
  failed_details_(0),
  this_time_lapse_(vintage),
  simbox_below_(NULL),
  rms_data_given_(false),
  horizon_data_given_(false)
{
  LogKit::WriteHeader("Reading Travel Time Data");

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
                 modelTravelTimeStatic,
                 timeSimbox,
                 errTxt,
                 failed_rms_data);

  errorCorrXY_ = setErrorCorrXYGrid(timeSimbox,  modelSettings);

  bool failed_loading_model = failed_surfaces || failed_rms_data;

  if (failed_loading_model) {
    LogKit::WriteHeader("Error(s) while loading dynamic travel time data");
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

  if (errorCorrXY_ != NULL)
    delete errorCorrXY_;

  delete simbox_below_;
}

//-------------------------------------------------------------------------------------------//

void
ModelTravelTimeDynamic::processHorizons(const ModelSettings         * modelSettings,
                                        const InputFiles            * inputFiles,
                                        std::string                 & errTxt,
                                        bool                        & failed)
{
  std::string tmpErrText = "";

  const std::vector<std::string> & push_down_horizons = inputFiles->getTravelTimeHorizons(this_time_lapse_);

  int n_horizons = static_cast<int>(push_down_horizons.size());

  if (n_horizons > 0) {
    LogKit::LogFormatted(LogKit::Low, "\nReading horizon travel time data:\n");

    horizon_data_given_ = true;

    horizon_names_              = modelSettings->getTimeLapseTravelTimeHorizons(this_time_lapse_);
    horizon_standard_deviation_ = modelSettings->getTimeLapseTravelTimeHorizonSD(this_time_lapse_);

    push_down_horizons_.resize(n_horizons);

    for (int i = 0; i < n_horizons; i++) {
      LogKit::LogFormatted(LogKit::Low, "\nHorizon \""+horizon_names_[i]+"\":\n");

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
ModelTravelTimeDynamic::processRMSData(const ModelSettings         * modelSettings,
                                       const InputFiles            * inputFiles,
                                       const ModelTravelTimeStatic * modelTravelTimeStatic,
                                       const Simbox                * timeSimbox,
                                       std::string                 & errTxt,
                                       bool                        & failed)

{
  std::string tmpErrText = "";

  const std::string & file_name  = inputFiles->getRmsVelocities(this_time_lapse_);

  if (file_name != "") {
    LogKit::LogFormatted(LogKit::Low, "\nReading RMS travel time data:\n");

    rms_data_given_ = true;

    readRMSData(file_name, timeSimbox, tmpErrText);

    rms_standard_deviation_ = modelSettings->getRMSStandardDeviation(this_time_lapse_);

    setupSimboxBelow(timeSimbox,
                     modelSettings->getOutputGridFormat(),
                     modelSettings->getOutputGridDomain(),
                     modelSettings->getOtherOutputFlag(),
                     modelTravelTimeStatic->getNumberOfLayersBelow(),
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
  std::ifstream file;
  NRLib::OpenRead(file, fileName);

  if (!file)
    errTxt += "Could not open RMS data file "+fileName+" for reading.\n";


  int line = 0;

  while (line < 32) {  // NBNB can be made more robust
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
    catch (NRLib::IOError & e) {
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
ModelTravelTimeDynamic::setupSimboxBelow(const Simbox  * timeSimbox,
                                         const int     & outputFormat,
                                         const int     & outputDomain,
                                         const int     & otherOutput,
                                         const int     & n_layers,
                                         std::string   & errTxt)
{
  // Dynamic as timeSimbox changes in different time lapses

  int n_rms_traces = static_cast<int>(rms_traces_.size());

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
  simbox_below_ = new Simbox(*timeSimbox);

  double dz      = static_cast<double>(lz / n_layers);
  double z_shift = 0;

  simbox_below_->setDepth(bot_simbox_surface, z_shift, lz, dz);

  if ((otherOutput & IO::EXTRA_SURFACES) > 0 && (outputDomain & IO::TIMEDOMAIN) > 0) {
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixTime() + "_Below_Reservoir";
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixTime() + "_Below_Reservoir";
//    simbox_below_->writeTopBotGrids(topSurf,
//                                    baseSurf,
//                                    IO::PathToBackground(),
//                                    outputFormat);
    simbox_below_->setTopBotName(topSurf, baseSurf, outputFormat);

  }

  if (lz < 0) {
    errTxt += "A problem was encountered for the simbox below the reservoir in the RMS inversion\n";
    errTxt += "There are no RMS data below the reservoir\n";
  }
}

//----------------------------------------------------------------------------

Surface *
ModelTravelTimeDynamic::setErrorCorrXYGrid(const Simbox        * timeSimbox,
                                           const ModelSettings * modelSettings) const
{
  float dx  = static_cast<float>(timeSimbox->getdx());
  float dy  = static_cast<float>(timeSimbox->getdy());

  //NBNB Next two lines should be active, but gave errors, so leave with dummies until activation.
  //int   nx  = modelSettings->getNXpad();
  //int   ny  = modelSettings->getNYpad();
  int nx=-1;
  int ny=-1;


  Surface * grid = new Surface(0, 0, dx * nx, dy * ny, nx, ny, 0.0, RMISSING);

  Vario * vario = modelSettings->getLateralTravelTimeErrorCorr(this_time_lapse_);

  if (vario != NULL) {
    int refi,refj;

    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        if (i < (nx / 2 + 1))
          refi = i;
        else
          refi = i - nx;

        if ( j < (ny / 2 + 1))
          refj = j;
        else
          refj = j - ny;

        (*grid)(j * nx + i) = vario->corr(refi * dx, refj * dy);
      }
    }
  }
  return(grid);
}
