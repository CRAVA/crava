/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELTRAVELTIMEDYNAMIC_H
#define MODELTRAVELTIMEDYNAMIC_H

#include "src/definitions.h"

class InputFiles;
class ModelSettings;
class RMSTrace;
class Simbox;

class ModelTravelTimeDynamic
{
public:
  ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                         const InputFiles              * inputFiles,
                         const Simbox                  * timeSimbox,
                         const int                     & vintage);

  ~ModelTravelTimeDynamic();


  bool                          getFailed()                         const { return failed_                          ;}
  std::vector<bool>             getFailedDetails()                  const { return failed_details_                  ;}

  const std::vector<Surface>    getInitialHorizons()                const { return initial_horizons_                ;}
  const std::vector<Surface>    getPushDownHorizons()               const { return push_down_horizons_              ;}
  const std::vector<double>     getHorizonStandardDeviation()       const { return horizon_standard_deviation_      ;}
  const bool                    getHorizonDataGiven()               const { return horizon_data_given_              ;}

  const std::vector<RMSTrace *> getRMSTraces()                      const { return rms_traces_                      ;}
  const double                  getMeanVpTop()                      const { return mean_vp_top_                     ;}
  const double                  getMeanVpBase()                     const { return mean_vp_base_                    ;}
  const double                  getVarVpAbove()                     const { return var_vp_above_                    ;}
  const double                  getVarVpBelow()                     const { return var_vp_below_                    ;}
  const double                  getRangeAbove()                     const { return range_above_                     ;}
  const double                  getRangeBelow()                     const { return range_below_                     ;}
  const double                  getRMSStandardDeviation()           const { return rms_standard_deviation_          ;}
  const Simbox *                getSimboxAbove()                    const { return simbox_above_                    ;}
  const Simbox *                getSimboxBelow()                    const { return simbox_below_                    ;}
  const int                     getThisTimeLapse()                  const { return this_time_lapse_                 ;}
  const double                  getLzLimit()                        const { return lz_limit_                        ;}
  const bool                    getRMSDataGiven()                   const { return rms_data_given_                  ;}

private:

  void                          processHorizons(const ModelSettings   * modelSettings,
                                                const InputFiles      * inputFiles,
                                                std::string           & errTxt,
                                                bool                  & failed);

  void                          processRMSData(const ModelSettings * modelSettings,
                                               const InputFiles    * inputFiles,
                                               const Simbox        * timeSimbox,
                                               std::string         & errTxt,
                                               bool                & failed);

  void                          readRMSData(const std::string & fileName,
                                            const Simbox      * timeSimbox,
                                            std::string       & errTxt);

  void                          setupSimboxAbove(const Simbox  * timeSimbox,
                                                 const int     & outputFormat,
                                                 const int     & outputDomain,
                                                 const int     & otherOutput,
                                                 const double  & lz_limit,
                                                 std::string   & errTxt);

  void                          setupSimboxBelow(const Simbox  * timeSimbox,
                                                 const int     & outputFormat,
                                                 const int     & outputDomain,
                                                 const int     & otherOutput,
                                                 std::string   & errTxt);

  std::vector<Surface>      push_down_horizons_;                    ///< Push down horizons used for horizon inversion
  std::vector<Surface>      initial_horizons_;                      ///< TP0 being the initial horizon before push down
  std::vector<std::string>  horizon_names_;                         ///< Names corresponding to the horizons

  std::vector<double>       horizon_standard_deviation_;            ///< Observation error for the horizon data

  std::vector<RMSTrace *>   rms_traces_;

  double                    rms_standard_deviation_;                ///< Observation error for the RMS data

  int                       n_layers_above_;                        ///< Number of layers to be used in the RMS inversion above the reservoir
  int                       n_layers_below_;                        ///< Number of layers to be used in the RMS inversion below the reservoir

  double                    mean_vp_top_;                           ///< E(Vp) at the top of the zone above the reservoir, that is, at sea level
  double                    mean_vp_base_;                          ///< E(Vp) at the base of the zone below the reservoir

  double                    var_vp_above_;                          ///< Var(Vp) above the reservoir
  double                    var_vp_below_;                          ///< Var(Vp) below the reservoir

  double                    range_above_;                           ///< Range of the temporal corralation function used above the reservoir
  double                    range_below_;                           ///< Range of the temporal corralation function used below the reservoir

  bool                      failed_;                                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;                        ///< Detailed failed information.

  int                       this_time_lapse_;                       ///< Time lapse of the current travel time data set

  double                    lz_limit_;                              ///< Minimum allowed value for (min interval thickness)/(max interval thickness)
                                                                    ///< Also stored in modelSettings. Needed in rmsinversion.cpp

  Simbox *                  simbox_above_;                          ///< Simbox to be used above the reservoir
  Simbox *                  simbox_below_;                          ///< Simbox to be used below the reservoir

  bool                      rms_data_given_;                        ///< True if RMS data are given
  bool                      horizon_data_given_;                    ///< True if horizon data are given

};

#endif
