/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELAVOSTATIC_H
#define MODELAVOSTATIC_H

#include <stdio.h>

#include "nrlib/surface/regularsurface.hpp"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/inputfiles.h"
#include "src/commondata.h"

#include "nrlib/well/well.hpp"


struct irapgrid;
class Wavelet;
class Vario;
class Simbox;
class FFTGrid;
class InputFiles;
class ModelAVOStatic
{
public:

  ModelAVOStatic(ModelSettings        *& model_settings,
                 const InputFiles      * input_files,
                 CommonData            * common_data,
                 const Simbox          * simbox,
                 int                     i_intervals);

  ~ModelAVOStatic();

  const std::vector<Surface*>      & GetFaciesEstimInterval()       const { return facies_estim_interval_  ;}

  /*const*/ std::vector<Surface *> & GetWaveletEstimInterval()  /*const*/ { return wavelet_estim_interval_ ;}
  /*const*/ std::vector<Surface *> & GetFaciesEstimInterval()   /*const*/ { return facies_estim_interval_  ;}
  /*const*/ std::vector<Surface *> & GetWellMoveInterval()      /*const*/ { return well_move_interval_     ;}

  FFTGrid                          * GetErrCorr()                         { return err_corr_               ;}

  bool                               GetForwardModeling()                 { return forward_modeling_       ;}

  //void             AddSeismicLogs(std::map<std::string, BlockedLogsCommon *> blocked_wells,
  //                                std::vector<FFTGrid *>                     seis_cube,
  //                                int                                        n_angles);                              // Changes wells

  static FFTGrid *        CreateFFTGrid(int nx,  int ny,  int nz,
                                        int nxp, int nyp, int nzp,
                                        bool file_grid);

private:

  void             CheckAvailableMemory(const Simbox              * time_simbox,
                                        ModelSettings       * model_settings,
                                        const InputFiles    * input_files);

  bool                      forward_modeling_;

  std::vector<Surface *>    wavelet_estim_interval_;  ///< Grids giving the wavelet estimation interval.
  std::vector<Surface *>    facies_estim_interval_;   ///< Grids giving the facies estimation intervals.
  std::vector<Surface *>    well_move_interval_;      ///< Grids giving the facies estimation intervals.

  FFTGrid                 * err_corr_;

};

#endif
