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
class Corr;
class Wavelet;
class Vario;
class Simbox;
class WellData;
class FFTGrid;
class GridMapping;
class InputFiles;

class ModelAVOStatic
{
public:
  ModelAVOStatic(ModelSettings      *& modelSettings,
                 const InputFiles    * inputFiles,
                 std::vector<bool>     failedGeneralDetails,
                 GridMapping         * timeCutMapping,
                 Simbox              * timeSimbox,
                 Simbox             *& timeBGSimbox,
                 Simbox              * timeSimboxConstThick);
  ~ModelAVOStatic();

  WellData                  **& getWells()                 /*const*/ { return wells_                  ;}
  const float                 * getPriorFacies()           const { return priorFacies_            ;}
  FFTGrid                    ** getPriorFaciesCubes()      const { return priorFaciesProbCubes_   ;}
  const std::vector<Surface*> & getFaciesEstimInterval()   const { return faciesEstimInterval_    ;}

  /*const*/ std::vector<Surface *> & getWaveletEstimInterval()  /*const*/ { return waveletEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getFaciesEstimInterval()   /*const*/ { return faciesEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getWellMoveInterval()      /*const*/ { return wellMoveInterval_   ;}


  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  void                          writeWells(       WellData ** wells, ModelSettings * modelSettings) const;
  void                          writeBlockedWells(WellData ** wells, ModelSettings * modelSettings) const;

  void             addSeismicLogs(WellData ** wells, FFTGrid ** seisCube,
                                  ModelSettings * modelSettings);                              // Changes wells
  void             generateSyntheticSeismic(Wavelet      ** wavelet,
                                            WellData     ** wells,
                                            float        ** reflectionMatrix,
                                            Simbox        * timeSimbox,
                                            ModelSettings * modelSettings);                    // Changes wells
  void             processWellLocation(FFTGrid                     ** seisCube,
                                       WellData                    ** wells,
                                       float                       ** reflectionMatrix,
                                       Simbox                       * timeSimbox,
                                       ModelSettings                * modelSettings,
                                       const std::vector<Surface *> & interval);              // Changes wells

private:
  void             processWells(WellData       **& wells,
                                Simbox           * timeSimbox,
                                Simbox           * timeBGSimbox,
                                Simbox           * timeSimboxConstThick,
                                ModelSettings   *& modelSettings,
                                const InputFiles * inputFiles,
                                std::string      & errText,
                                bool             & failed);

  void             processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
                                          float                       *& priorFacies,
                                          WellData                    ** wells,
                                          Simbox                       * timeSimbox,
                                          Simbox                       * timeCutSimbox,
                                          ModelSettings                * modelSettings,
                                          bool                         & failed,
                                          std::string                  & errTxt,
                                          const InputFiles             * inputFiles);
  void             readPriorFaciesProbCubes(const InputFiles  * inputFiles,
                                            ModelSettings     * modelSettings,
                                            FFTGrid         **& priorFaciesProbCubes,
                                            Simbox            * timeSimbox,
                                            Simbox            * timeCutSimbox,
                                            std::string       & errTxt,
                                            bool              & failed);
  void             loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
                                     std::vector<Surface *> & faciesEstimInterval,
                                     std::vector<Surface *> & wellMoveInterval,
                                     Simbox                 * timeSimbox,
                                     const InputFiles       * inputFiles,
                                     std::string            & errText,
                                     bool                   & failed);
  void             checkFaciesNames(WellData        ** wells,
                                    ModelSettings   *& modelSettings,
                                    const InputFiles * inputFiles,
                                    std::string      & tmpErrText,
                                    int              & error);

  void             checkAvailableMemory(Simbox              * timeSimbox,
                                        ModelSettings       * modelSettings,
                                        const InputFiles    * inputFiles);

  bool                      forwardModeling_;
  int                       numberOfWells_;

  WellData               ** wells_;                 ///< Well data

  std::vector<Surface *>    waveletEstimInterval_;  ///< Grids giving the wavelet estimation interval.
  std::vector<Surface *>    faciesEstimInterval_;   ///< Grids giving the facies estimation intervals.
  std::vector<Surface *>    wellMoveInterval_;      ///< Grids giving the facies estimation intervals.

  float                   * priorFacies_;           ///<
  FFTGrid                ** priorFaciesProbCubes_;  ///< Cubes for prior facies probabilities

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.
};

#endif
