/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELAVOSTATIC_H
#define MODELAVOSTATIC_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/background.h" //or move getAlpha & co to cpp-file.
#include "src/modelsettings.h"
#include "src/inputfiles.h"

struct irapgrid;
class Wavelet;
class Vario;
class Simbox;
class FFTGrid;
class GridMapping;
class InputFiles;
class ModelAVOStatic
{
public:
  ModelAVOStatic(ModelSettings        *& modelSettings,
                 ModelGeneral         *& modelGeneral,
                 const InputFiles      * inputFiles,
                 GridMapping           * timeCutMapping,
                 Simbox                * timeSimbox,
                 Simbox               *& timeBGSimbox,
                 Simbox                * timeSimboxConstThick,
                 std::vector<WellData *> wells);
  ~ModelAVOStatic();

  const std::vector<Surface*> & getFaciesEstimInterval()   const { return faciesEstimInterval_    ;}

  /*const*/ std::vector<Surface *> & getWaveletEstimInterval()  /*const*/ { return waveletEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getFaciesEstimInterval()   /*const*/ { return faciesEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getWellMoveInterval()      /*const*/ { return wellMoveInterval_   ;}


  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  void                          writeWells(       std::vector<WellData *> wells, ModelSettings * modelSettings) const;
  void                          writeBlockedWells(std::vector<WellData *> wells, ModelSettings * modelSettings, std::vector<std::string> facies_name, std::vector<int> facies_label) const;

  void             addSeismicLogs(std::vector<WellData *>     wells,
                                  FFTGrid                  ** seisCube,
                                  const ModelSettings             * modelSettings,
                                  int                         nAngles);                              // Changes wells

  void             generateSyntheticSeismic(Wavelet              ** wavelet,
                                            std::vector<WellData *> wells,
                                            const float   * const * reflectionMatrix,
                                            const Simbox          * timeSimbox,
                                            const ModelSettings   * modelSettings,
                                            int                     nAngles);                    // Changes wells

  void             deleteDynamicWells(std::vector<WellData *> wells,
                                      int         nWells);

private:
  void             blockLogs(std::vector<WellData *> & wells,
                                Simbox               * timeSimbox,
                                Simbox               * timeBGSimbox,
                                Simbox               * timeSimboxConstThick,
                                ModelSettings       *& modelSettings);

  void             loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
                                     std::vector<Surface *> & faciesEstimInterval,
                                     std::vector<Surface *> & wellMoveInterval,
                                     Simbox                 * timeSimbox,
                                     const InputFiles       * inputFiles,
                                     std::string            & errText,
                                     bool                   & failed);

  void             checkAvailableMemory(Simbox              * timeSimbox,
                                        ModelSettings       * modelSettings,
                                        const InputFiles    * inputFiles);

  bool                      forwardModeling_;

  std::vector<Surface *>    waveletEstimInterval_;  ///< Grids giving the wavelet estimation interval.
  std::vector<Surface *>    faciesEstimInterval_;   ///< Grids giving the facies estimation intervals.
  std::vector<Surface *>    wellMoveInterval_;      ///< Grids giving the facies estimation intervals.

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.
};

#endif
