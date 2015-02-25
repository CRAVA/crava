/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "nrlib/random/beta.hpp"
#include "src/blockedlogscommon.h"

class Vario;
class Simbox;
class FFTGrid;
class CovGrid2D;
class GridMapping;
class KrigingData3D;
class KrigingData2D;
class ModelSettings;
class BlockedLogsCommon;
class MultiIntervalGrid;

//Special note on the use of Background:
//All pointers used here are also used externally, so no deletion happens.

class Background
{
public:

  static
  void SetupBackground(std::vector<NRLib::Grid<float> *>                & parameters,
                       std::vector<std::vector<double> >                & vertical_trends,
                       NRLib::Grid<float>                               * velocity,
                       const Simbox                                     * time_simbox,
                       const Simbox                                     * time_bg_simbox,
                       const std::map<std::string, BlockedLogsCommon *> & bl,
                       const std::map<std::string, BlockedLogsCommon *> & bg_bl,
                       const ModelSettings                              * modelSettings,
                       const std::string                                & interval_name,
                       std::string                                      & err_text);

  static
  void SetupBackground(FFTGrid ** grids);

  static
  void         FillInVerticalTrend(FFTGrid                   * bgTrend,
                                   const std::vector<double> & trend);


private:

  static
  void         GenerateBackgroundModel(NRLib::Grid<float>                               * bg_vp,
                                       NRLib::Grid<float>                               * bg_vs,
                                       NRLib::Grid<float>                               * bg_rho,
                                       std::vector<std::vector<double> >                & vertical_trends,
                                       NRLib::Grid<float>                               * velociy,
                                       const Simbox                                     * simbox,
                                       const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                       const ModelSettings                              * model_settings,
                                       const std::string                                & interval_name,
                                       std::string                                      & err_text);

  static
  void         ResampleBackgroundModel(NRLib::Grid<float>  * & bg_vp,
                                       NRLib::Grid<float>  * & bg_vs,
                                       NRLib::Grid<float>  * & bg_rho,
                                       const Simbox        *   bg_simbox,
                                       const Simbox        *   simbox);

  static
  void         CalculateBackgroundTrend(std::vector<double>               & trend,
                                        std::vector<double>               & avgDev,
                                        const int                           nz,
                                        const float                         dz,
                                        float                               logMin,
                                        float                               logMax,
                                        float                               maxHz,
                                        std::vector<std::vector<double> > & wellTrend,
                                        std::vector<std::vector<double> > & highCutWellTrend,
                                        const std::string                 & name);

  static
  void         GetKrigingWellTrends(std::vector<std::vector<double> >                & bl_vp,
                                    std::vector<std::vector<double> >                & bl_vs,
                                    std::vector<std::vector<double> >                & bl_rho,
                                    std::vector<std::vector<double> >                & vt_vp,
                                    std::vector<std::vector<double> >                & vt_vs,
                                    std::vector<std::vector<double> >                & vt_rho,
                                    std::vector<const std::vector<int> *>            & ipos,
                                    std::vector<const std::vector<int> *>            & jpos,
                                    std::vector<const std::vector<int> *>            & kpos,
                                    std::vector<int>                                 & n_blocks,
                                    int                                              & tot_blocks,
                                    const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                    const int                                        & n_wells);

  static
  void         GetWellTrends(std::vector<std::vector<double> >                & well_trend,
                             std::vector<std::vector<double> >                & high_cut_Well_trend,
                             const std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                             const int                                        & nz,
                             const std::string                                & name,
                             std::string                                      & err_text);


  static
  void         WriteTrendsToFile(std::vector<double> & trend,
                                 const Simbox        * simbox,
                                 bool                  write1D,
                                 bool                  write3D,
                                 bool                  hasVelocityTrend,
                                 const std::string   & name,
                                 bool                  isFile,
                                 const std::string   & interval_name);


  static
  const CovGrid2D & MakeCovGrid2D(const Simbox      * simbox,
                                  Vario             * vario,
                                  int                 debugFlag,
                                  const std::string & interval_name);


  static
  void         SetupKrigingData2D(std::vector<KrigingData2D>                 & kriging_data_vp,
                                  std::vector<KrigingData2D>                 & kriging_data_vs,
                                  std::vector<KrigingData2D>                 & kriging_data_rho,
                                  std::vector<double>                        & trend_vp,
                                  std::vector<double>                        & trend_vs,
                                  std::vector<double>                        & trend_rho,
                                  const int                                    output_flag,
                                  const int                                  & nz,
                                  const float                                & dz,
                                  const int                                  & tot_blocks,
                                  const std::vector<int>                     & n_blocks,
                                  const std::vector<std::vector<double> >    & bl_vp,
                                  const std::vector<std::vector<double> >    & bl_vs,
                                  const std::vector<std::vector<double> >    & bl_rho,
                                  const std::vector<std::vector<double> >    & vt_vp,
                                  const std::vector<std::vector<double> >    & vt_vs,
                                  const std::vector<std::vector<double> >    & vt_rho,
                                  const std::vector<const std::vector<int> *>  ipos,
                                  const std::vector<const std::vector<int> *>  jpos,
                                  const std::vector<const std::vector<int> *>  kpos,
                                  const std::string                          & interval_name);


  static
  void         MakeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                    NRLib::Grid<float>               * bg_grid,
                                    std::vector<double>              & trend,
                                    const Simbox                     * simbox,
                                    const CovGrid2D                  & cov_grid_2D,
                                    const std::string                & type,
                                    int                                n_threads);


  static
  void         CalculateVelocityDeviations(NRLib::Grid<float>                               * velocity,
                                           const Simbox                                     * simbox,
                                           const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                           std::vector<double>                              & trend_vel,
                                           std::vector<double>                              & avg_dev_vel,
                                           std::vector<double>                              & avg_dev_vp,
                                           //int                                                outputFlag,
                                           int                                                n_wells);


  static
  void         ResampleParameter(NRLib::Grid<float> *& p_new, // Resample to
                                 NRLib::Grid<float> *  p_old, // Resample from
                                 const Simbox       *  simbox_new,
                                 const Simbox       *  simbox_old);


  static
  void         CalculateVerticalTrend(std::vector<std::vector<double> > & wellTrend,
                                      std::vector<double>               & trend,
                                      float                               logMin,
                                      float                               logMax,
                                      float                               maxHz,
                                      int                                 nz,
                                      float                               dz,
                                      const std::string                 & name);


  static
  void         WriteVerticalTrend(std::vector<double> & trend,
                                  float                 dz,
                                  int                   nz,
                                  std::string           name,
                                  const std::string   & interval_name);


  static
  void         CalculateDeviationFromVerticalTrend(std::vector<std::vector<double> > & wellTrend,
                                                   const std::vector<double>         & trend,
                                                   std::vector<double>               & avg_dev,
                                                   const int                           nd);


  static
  void         WriteDeviationsFromVerticalTrend(const std::vector<double>                        & avg_dev_vp,
                                                const std::vector<double>                        & avg_dev_vs,
                                                const std::vector<double>                        & avg_dev_rho,
                                                const std::vector<double>                        & trend_vp,
                                                const std::vector<double>                        & trend_vs,
                                                const std::vector<double>                        & trend_rho,
                                                const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                                const int                                          n_wells,
                                                const int                                          nz);


  static
  void         SmoothTrendWithLocalLinearRegression(std::vector<double> & trend,
                                                    int                 * count,
                                                    int                   nWells,
                                                    int                   nz,
                                                    float                 dz,
                                                    float                 min_value,
                                                    float                 max_value,
                                                    std::string           parName);



  static
  FFTGrid    * CopyFFTGrid(FFTGrid   * origGrid,
                           const bool  expTrans,
                           const bool  fileGrid);

};

#endif
