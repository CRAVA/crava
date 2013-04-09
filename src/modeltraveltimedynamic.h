/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELTRAVELTIMEDYNAMIC_H
#define MODELTRAVELTIMEDYNAMIC_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"

class Simbox;
class WellData;
class FFTGrid;
class InputFiles;
class ModelGeneral;
class SeismicParametersHolder;

class ModelTravelTimeDynamic
{
public:
  ModelTravelTimeDynamic(const ModelSettings           * modelSettings,
                         const ModelGeneral            * modelGeneral,
                         const InputFiles              * inputFiles,
                         const int                     & vintage);

  ~ModelTravelTimeDynamic();

  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}


private:

  void                          processHorizons(std::vector<Surface>   & horizons,
                                                const InputFiles       * inputFiles,
                                                std::string            & errTxt,
                                                bool                   & failed);

  void                          processRMSData(FFTGrid                 *& rms_data,
                                               const ModelSettings      * modelSettings,
                                               const InputFiles         * inputFiles,
                                               const Simbox             * timeSimbox,
                                               GridMapping              * timeDepthMapping,
                                               const GridMapping        * timeCutMapping,
                                               std::string              & errTxt,
                                               bool                     & failed);

  std::vector<Surface>      horizons_;              ///< Horizons used for horizon inversion
  FFTGrid                 * rms_data_;              ///< RMS data U^2

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.

  int                       thisTimeLapse_;         ///< Time lapse of the current travel time data set

};

#endif
