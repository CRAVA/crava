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
class FFTGrid;
class GridMapping;
class InputFiles;
class ModelAVOStatic
{
public:
  ModelAVOStatic(ModelSettings        *& modelSettings,
                 const InputFiles      * inputFiles,
                 std::vector<bool>       failedGeneralDetails,
                 GridMapping           * timeCutMapping,
                 Simbox                * timeSimbox,
                 Simbox               *& timeBGSimbox,
                 Simbox                * timeSimboxConstThick,
                 std::vector<WellData *> wells);
  ~ModelAVOStatic();

  const float                 * getPriorFacies()           const { return priorFacies_            ;}
  FFTGrid                    ** getPriorFaciesCubes()      const { return priorFaciesProbCubes_   ;}
  const std::vector<Surface*> & getFaciesEstimInterval()   const { return faciesEstimInterval_    ;}

  /*const*/ std::vector<Surface *> & getWaveletEstimInterval()  /*const*/ { return waveletEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getFaciesEstimInterval()   /*const*/ { return faciesEstimInterval_   ;}
  /*const*/ std::vector<Surface *> & getWellMoveInterval()      /*const*/ { return wellMoveInterval_   ;}


  bool                          getFailed()                const { return failed_                 ;}
  std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  void                          writeWells(       std::vector<WellData *> wells, ModelSettings * modelSettings) const;
  void                          writeBlockedWells(std::vector<WellData *> wells, ModelSettings * modelSettings) const;

  void             addSeismicLogs(std::vector<WellData *>     wells,
                                  FFTGrid                  ** seisCube,
                                  ModelSettings             * modelSettings,
                                  int                         nAngles);                              // Changes wells

  void             generateSyntheticSeismic(Wavelet              ** wavelet,
                                            std::vector<WellData *> wells,
                                            float                ** reflectionMatrix,
                                            Simbox                * timeSimbox,
                                            ModelSettings         * modelSettings,
                                            int                     nAngles);                    // Changes wells

  void             deleteDynamicWells(std::vector<WellData *> wells,
                                      int         nWells);

private:
  void             blockLogs(std::vector<WellData *> & wells,
                                Simbox               * timeSimbox,
                                Simbox               * timeBGSimbox,
                                Simbox               * timeSimboxConstThick,
                                ModelSettings       *& modelSettings);

  void             processPriorFaciesProb(const std::vector<Surface*>  & faciesEstimInterval,
                                          float                       *& priorFacies,
                                          std::vector<WellData *>        wells,
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

  void             processFaciesInformation(ModelSettings     *& modelSettings,
                                            const InputFiles   * inputFiles,
                                            std::string        & tmpErrText,
                                            int                & error) const;

  void             checkAvailableMemory(Simbox              * timeSimbox,
                                        ModelSettings       * modelSettings,
                                        const InputFiles    * inputFiles);

  bool                      forwardModeling_;

  std::vector<Surface *>    waveletEstimInterval_;  ///< Grids giving the wavelet estimation interval.
  std::vector<Surface *>    faciesEstimInterval_;   ///< Grids giving the facies estimation intervals.
  std::vector<Surface *>    wellMoveInterval_;      ///< Grids giving the facies estimation intervals.

  float                   * priorFacies_;           ///<
  FFTGrid                ** priorFaciesProbCubes_;  ///< Cubes for prior facies probabilities

  bool                      failed_;                ///< Indicates whether errors occured during construction.
  std::vector<bool>         failed_details_;        ///< Detailed failed information.
};

#endif
