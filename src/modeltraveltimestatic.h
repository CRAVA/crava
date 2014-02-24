/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODEL_TRAVEL_TIME_STATIC_H
#define MODEL_TRAVEL_TIME_STATIC_H

#include "src/definitions.h"

class InputFiles;
class ModelSettings;
class Simbox;

class ModelTravelTimeStatic
{
public:
  ModelTravelTimeStatic(const ModelSettings           * modelSettings,
                        const InputFiles              * inputFiles,
                        const Simbox                  * timeSimbox);

  ~ModelTravelTimeStatic();


  bool                             getFailed()                         const { return failed_                          ;}
  std::vector<bool>                getFailedDetails()                  const { return failed_details_                  ;}

  const double                     getMeanVpTop()                      const { return mean_vp_top_                     ;}
  const double                     getMeanVpBase()                     const { return mean_vp_base_                    ;}
  const double                     getVarVpAbove()                     const { return var_vp_above_                    ;}
  const double                     getVarVpBelow()                     const { return var_vp_below_                    ;}
  const double                     getRangeAbove()                     const { return range_above_                     ;}
  const double                     getRangeBelow()                     const { return range_below_                     ;}
  const Simbox *                   getSimboxAbove()                    const { return simbox_above_                    ;}
  const double                     getLzLimit()                        const { return lz_limit_                        ;}
  const int                        getOutputGridFormat()               const { return format_flag_                     ;}
  const std::vector<Surface>       getInitialHorizons()                const { return initial_horizons_                ;}
  const std::vector<std::string>   getInitialHorizonNames()            const { return horizon_names_                   ;}
  const int                        getNumberOfLayersBelow()            const { return n_layers_below_                  ;}

private:

  void                          processHorizons(const ModelSettings   * modelSettings,
                                                const InputFiles      * inputFiles,
                                                std::string           & errTxt,
                                                bool                  & failed);

  void                          processRMSData(const ModelSettings * modelSettings,
                                               const Simbox        * timeSimbox,
                                               std::string         & errTxt,
                                               bool                & failed);


  void                          setupSimboxAbove(const Simbox  * timeSimbox,
                                                 const int     & outputFormat,
                                                 const int     & outputDomain,
                                                 const int     & otherOutput,
                                                 const double  & lz_limit,
                                                 std::string   & errTxt);

  std::vector<Surface>      initial_horizons_;                      ///< TP0 being the initial horizon before push down
  std::vector<std::string>  horizon_names_;                         ///< Names corresponding to the horizons

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

  double                    lz_limit_;                              ///< Minimum allowed value for (min interval thickness)/(max interval thickness)
                                                                    ///< Also stored in modelSettings. Needed in rmsinversion.cpp

  int                       format_flag_;                           ///< Decides output format, see modelSettings where it is also used

  Simbox *                  simbox_above_;                          ///< Simbox to be used above the reservoir

};

#endif
