/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "lib/timekit.hpp"
#include "src/modeltraveltimedynamic.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/timings.h"
#include "src/rmstrace.h"

ModelTravelTimeDynamic::ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                                               const InputFiles              * inputFiles,
                                               const int                     & vintage)
: rms_traces_(0),
  thisTimeLapse_(vintage)
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
  for(size_t i=0; i<rms_traces_.size(); i++)
    delete rms_traces_[i];
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
                                       std::string              & errTxt,
                                       bool                     & failed)

{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  LogKit::WriteHeader("Reading RMS travel time data");

  const std::string & file_name  = inputFiles->getRmsVelocities(thisTimeLapse_);
  std::string         tmpErrText = "";

  readRMSData(file_name, tmpErrText);

  n_layers_above_   = modelSettings->getRMSnLayersAbove();
  n_layers_below_   = modelSettings->getRMSnLayersBelow();

  mean_vp_top_      = modelSettings->getRMSMeanVpTop();
  mean_vp_base_     = modelSettings->getRMSMeanVpBase();

  var_vp_above_     = modelSettings->getRMSVarianceVpAbove();
  var_vp_below_     = modelSettings->getRMSVarianceVpBelow();

  range_above_ = static_cast<float>(modelSettings->getRMSTemporalCorrelationRangeAbove());
  range_below_ = static_cast<float>(modelSettings->getRMSTemporalCorrelationRangeBelow());

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
            RMSTrace * trace = new RMSTrace(IL,
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
        RMSTrace * trace = new RMSTrace(IL,
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
