/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELGRAVITYDYNAMIC_H
#define MODELGRAVITYDYNAMIC_H

#include <stdio.h>

//#include "src/definitions.h"
//#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"

class Simbox;
class FFTGrid;
class GridMapping;
class InputFiles;
class ModelGravityStatic;
class ModelGeneral;
class Simbox;
class SeismicParametersHolder;

class ModelGravityDynamic
{
public:
  ModelGravityDynamic(const ModelSettings          * modelSettings,
                      const ModelGeneral           * modelGeneral,
                      ModelGravityStatic           * modelGravityStatic,
                      const InputFiles             * inputFiles,
                      int                            t,
                      SeismicParametersHolder      & seismicParameters);

  ModelGravityDynamic(const ModelSettings          * modelSettings,
                      const ModelGeneral           * modelGeneral,
                      ModelGravityStatic           * modelGravityStatic,
                      CommonData                   * commonData,
                      const InputFiles             * inputFiles,
                      int                            t,
                      SeismicParametersHolder      & seismicParameters);

  //ModelGravityDynamic(ModelGravityStatic       * modelGravityStatic,
  //                    Simbox                   * simbox,
  //                    SeismicParametersHolder  & seismicParameters,
  //                    const std::vector<float> & observation_location_utmx,
  //                    const std::vector<float> & observation_location_utmy,
  //                    const std::vector<float> & observation_location_depth,
  //                    const std::vector<float> & gravity_response,
  //                    const std::vector<float> & gravity_std_dev);

  // Possible other parameters?
  // GridMapping timeDepthMapping
  // Simbox timeSimbox
  // Simbox timeBGSimbox

  ~ModelGravityDynamic();

  bool                          GetFailed()                const { return failed_                 ;}
  std::vector<bool>             GetFailedDetails()         const { return failed_details_         ;}

  std::vector<float>            GetGravityResponse()       const { return gravity_response_       ;}
  std::vector<float>            GetGravityStdDev()         const { return gravity_std_dev_        ;}
  NRLib::Matrix                 GetGMatrix()               const { return G_                      ;}
  NRLib::Matrix                 GetGMatrixFullSize()       const { return G_fullsize_             ;}

  void                          SetSyntheticData(NRLib::Vector data) { synthetic_data_ = data     ;}

private:
  bool                      debug_;

  bool                      failed_;              ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;      ///< Detailed failed information.

  int                       thisTimeLapse_;

  std::vector<float> observation_location_utmx_;  ///< Vectors to store observation location coordinates
  std::vector<float> observation_location_utmy_;
  std::vector<float> observation_location_depth_;
  std::vector<float> gravity_response_;           ///< Vector to store gravimetric response
  std::vector<float> gravity_std_dev_;            ///< Vector to store standard deviation

  NRLib::Matrix G_;                               ///< Forward model gravity matrix to be used in inversion. Upscaled.
  NRLib::Matrix G_fullsize_;                      ///< Forward model gravity matrix to be used for calculating synthetic response. Full size.

  NRLib::Vector synthetic_data_;                  ///< Synthetic data response after inversion at this time vintage

  const ModelGeneral * modelGeneral_;

  void BuildGMatrix(ModelGravityStatic      * modelGravityStatic,
                    SeismicParametersHolder & seismicParameters,
                    Simbox                  * simbox = NULL);
};

#endif
