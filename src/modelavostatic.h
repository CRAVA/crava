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
  //ModelAVOStatic(ModelSettings        *& modelSettings,
  //               ModelGeneral         *& modelGeneral,
  //               const InputFiles      * inputFiles,
  //               GridMapping           * timeCutMapping,
  //               Simbox                * timeSimbox,
  //               Simbox               *& timeBGSimbox,
  //               Simbox                * timeSimboxConstThick,
  //               std::vector<WellData *> wells);

  ModelAVOStatic(ModelSettings        *& model_settings,
                 //ModelGeneral         *& model_general,
                 const InputFiles      * input_files,
                 CommonData            * common_data,
                 Simbox                * simbox,
                 int                     i_intervals);

  ~ModelAVOStatic();

  const std::vector<Surface*> & getFaciesEstimInterval()   const { return facies_estim_interval_    ;}

  /*const*/ std::vector<Surface *> & getWaveletEstimInterval()  /*const*/ { return wavelet_estim_interval_ ;}
  /*const*/ std::vector<Surface *> & getFaciesEstimInterval()   /*const*/ { return facies_estim_interval_  ;}
  /*const*/ std::vector<Surface *> & getWellMoveInterval()      /*const*/ { return well_move_interval_     ;}

  FFTGrid                          * getErrCorr()                         { return err_corr_              ;}

  //bool                          getFailed()                const { return failed_                 ;}
  //std::vector<bool>             getFailedDetails()         const { return failed_details_         ;}

  //void                          writeWells(       std::vector<WellData *> wells, ModelSettings * modelSettings) const;
  //void                          writeBlockedWells(std::vector<WellData *> wells, ModelSettings * modelSettings, std::vector<std::string> facies_name, std::vector<int> facies_label) const;
  void             writeBlockedWells(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                     ModelSettings                            * model_settings,
                                     std::vector<std::string>                   facies_name,
                                     std::vector<int>                           facies_label) const;

  //void             addSeismicLogs(std::vector<WellData *>     wells,
  //                                FFTGrid                  ** seisCube,
  //                                const ModelSettings       * modelSettings,
  //                                int                         nAngles);                              // Changes wells

  void             addSeismicLogs(std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                  std::vector<FFTGrid *>                     seis_cube,
                                  //FFTGrid                                 ** seis_cube,
                                  //const ModelSettings                      * model_settings,
                                  int                                        n_angles);                              // Changes wells

  //void             generateSyntheticSeismic(Wavelet              ** wavelet,
  //                                          std::vector<WellData *> wells,
  //                                          const float   * const * reflectionMatrix,
  //                                          const Simbox          * timeSimbox,
  //                                          const ModelSettings   * modelSettings,
  //                                          int                     nAngles);                    // Changes wells

  void             generateSyntheticSeismic(std::vector<Wavelet *>                   & wavelet,
                                            std::map<std::string, BlockedLogsCommon *> blocked_wells,
                                            const float *                      const * reflection_matrix,
                                            const Simbox                             * time_simbox,
                                            const ModelSettings                      * model_settings,
                                            int                                        n_angles);

  static FFTGrid *        createFFTGrid(int nx,  int ny,  int nz,
                                        int nxp, int nyp, int nzp,
                                        bool file_grid);

  //void             deleteDynamicWells(std::vector<WellData *> wells,
  //                                    int         nWells);

private:
  //void             blockLogs(std::vector<WellData *> & wells,
  //                              Simbox               * timeSimbox,
  //                              Simbox               * timeBGSimbox,
  //                              Simbox               * timeSimboxConstThick,
  //                              ModelSettings       *& modelSettings);

  //void             loadExtraSurfaces(std::vector<Surface *> & waveletEstimInterval,
  //                                   std::vector<Surface *> & faciesEstimInterval,
  //                                   std::vector<Surface *> & wellMoveInterval,
  //                                   Simbox                 * timeSimbox,
  //                                   const InputFiles       * inputFiles,
  //                                   std::string            & errText,
  //                                   bool                   & failed);

  void             checkAvailableMemory(Simbox              * time_simbox,
                                        ModelSettings       * model_settings,
                                        const InputFiles    * input_files);

  bool                      forward_modeling_;

  std::vector<Surface *>    wavelet_estim_interval_;  ///< Grids giving the wavelet estimation interval.
  std::vector<Surface *>    facies_estim_interval_;   ///< Grids giving the facies estimation intervals.
  std::vector<Surface *>    well_move_interval_;      ///< Grids giving the facies estimation intervals.

  FFTGrid                 * err_corr_;

  //bool                      failed_;                ///< Indicates whether errors occured during construction.
  //std::vector<bool>         failed_details_;        ///< Detailed failed information.
};

#endif
