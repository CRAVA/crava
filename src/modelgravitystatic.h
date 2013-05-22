#ifndef MODELGRAVITYSTATIC_H
#define MODELGRAVITYSTATIC_H

#include <stdio.h>

#include "src/modelsettings.h"
#include "src/inputfiles.h"


class ModelGravityStatic
{
public:
  ModelGravityStatic(ModelSettings        *& modelSettings,
                     ModelGeneral         *& modelGeneral,
                     const InputFiles      * inputFiles);
  ~ModelGravityStatic();

  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  // To be used by ModelGravityDynamic as well
  static void readGravityDataFile(const std::string   & fileName,
                                  const std::string   & readReason,
                                  int                   nObs,
                                  int                   nColumns,
                                  std::vector <float> & obs_loc_utmx,
                                  std::vector <float> & obs_loc_utmy,
                                  std::vector <float> & obs_loc_depth,
                                  std::vector <float> & gravity_response,
                                  std::vector <float> & gravity_std_dev,
                                  bool                  failed,
                                  std::string         & errText);

private:
  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  bool before_injection_start_;   // Flagg marking if this first vintage is befor start of injection

  // Vectors to store data in
  std::vector<float> observation_location_utmx_;
  std::vector<float> observation_location_utmy_;
  std::vector<float> observation_location_depth_;
  std::vector<float> gravity_response_;
  std::vector<float> gravity_std_dev_;
};

#endif
