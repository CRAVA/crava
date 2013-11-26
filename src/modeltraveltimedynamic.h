/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODEL_TRAVEL_TIME_DYNAMIC_H
#define MODEL_TRAVEL_TIME_DYNAMIC_H

#include "src/definitions.h"

class InputFiles;
class ModelSettings;
class ModelTravelTimeStatic;
class RMSTrace;
class Simbox;

class ModelTravelTimeDynamic
{
public:
  ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                         const InputFiles              * inputFiles,
                         const ModelTravelTimeStatic   * modelTravelTimeStatic,
                         const Simbox                  * timeSimbox,
                         const int                     & vintage);

  ~ModelTravelTimeDynamic();


  bool                             getFailed()                         const { return failed_                          ;}
  std::vector<bool>                getFailedDetails()                  const { return failed_details_                  ;}

  const std::vector<Surface>       getPushDownHorizons()               const { return push_down_horizons_              ;}
  const std::vector<double>        getHorizonStandardDeviation()       const { return horizon_standard_deviation_      ;}
  const bool                       getHorizonDataGiven()               const { return horizon_data_given_              ;}

  const std::vector<RMSTrace *>    getRMSTraces()                      const { return rms_traces_                      ;}
  const double                     getRMSStandardDeviation()           const { return rms_standard_deviation_          ;}
  const Simbox *                   getSimboxBelow()                    const { return simbox_below_                    ;}
  const int                        getThisTimeLapse()                  const { return this_time_lapse_                 ;}
  const bool                       getRMSDataGiven()                   const { return rms_data_given_                  ;}
  const Surface *                  getErrorCorrXY()                    const { return errorCorrXY_                     ;}
  const std::vector<std::string>   getHorizonNames()                   const { return horizon_names_                   ;}

private:

  void                          processHorizons(const ModelSettings         * modelSettings,
                                                const InputFiles            * inputFiles,
                                                std::string                 & errTxt,
                                                bool                        & failed);

  void                          processRMSData(const ModelSettings         * modelSettings,
                                               const InputFiles            * inputFiles,
                                               const ModelTravelTimeStatic * modelTravelTimeStatic,
                                               const Simbox                * timeSimbox,
                                               std::string                 & errTxt,
                                               bool                        & failed);

  void                          readRMSData(const std::string & fileName,
                                            const Simbox      * timeSimbox,
                                            std::string       & errTxt);

  void                          setupSimboxBelow(const Simbox  * timeSimbox,
                                                 const int     & outputFormat,
                                                 const int     & outputDomain,
                                                 const int     & otherOutput,
                                                 const int     & n_layers,
                                                 std::string   & errTxt);

  Surface *                     setErrorCorrXYGrid(const Simbox        * timeSimbox,
                                                   const ModelSettings * modelSettings) const;

  std::vector<Surface>      push_down_horizons_;                    ///< Push down horizons used for horizon inversion
  std::vector<std::string>  horizon_names_;                         ///< Names corresponding to the horizons

  std::vector<double>       horizon_standard_deviation_;            ///< Observation error for the horizon data

  Surface *                 errorCorrXY_ ;                          ///< Lateral error coorrelation for traveltime data

  std::vector<RMSTrace *>   rms_traces_;

  double                    rms_standard_deviation_;                ///< Observation error for the RMS data

  bool                      failed_;                                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;                        ///< Detailed failed information.

  int                       this_time_lapse_;                       ///< Time lapse of the current travel time data set

  Simbox *                  simbox_below_;                          ///< Simbox to be used below the reservoir

  bool                      rms_data_given_;                        ///< True if RMS data are given
  bool                      horizon_data_given_;                    ///< True if horizon data are given

};

#endif
