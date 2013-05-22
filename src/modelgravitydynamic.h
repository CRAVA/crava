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
class SeismicParametersHolder;

class ModelGravityDynamic
{
public:
  ModelGravityDynamic(const ModelSettings          * modelSettings,
                      const ModelGeneral           * modelGeneral,
                      const InputFiles             * inputFiles,
                      int                            t,
                      SeismicParametersHolder      & seismicParameters);

  // GridMapping timeDepthMapping
  // Simbox timeSimbox
  // Simbox timeBGSimbox

  ~ModelGravityDynamic();

  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

private:
  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  int                       thisTimeLapse_;

  // Vectors to store data in
  std::vector<float> observation_location_utmx_;
  std::vector<float> observation_location_utmy_;
  std::vector<float> observation_location_depth_;
  std::vector<float> gravity_response_;
  std::vector<float> gravity_std_dev_;

};

#endif
