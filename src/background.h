/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "nrlib/random/beta.hpp"

#include "src/blockedlogscommon.h"
#include "src/multiintervalgrid.h"

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

  Background(std::vector<NRLib::Grid<float> *>                & parameters,
             NRLib::Grid<float>                               * velocity,
             const Simbox                                     * time_simbox,
             const Simbox                                     * time_bg_simbox,
             const std::map<std::string, BlockedLogsCommon *> & bl,
             const std::map<std::string, BlockedLogsCommon *> & bg_bl,
             const ModelSettings                              * modelSettings,
             std::string                                      & err_text);

  //Background(std::vector<NRLib::Grid<float> *> & parameters,
  //           const std::vector<NRLib::Well>    & wells,
  //           const Simbox                      * timeSimbox,
  //           const ModelSettings               * modelSettings,
  //           const std::vector<std::string>    & surface_files,
  //           std::string                       & err_text);

  //Background(std::vector<std::vector<NRLib::Grid<float> *> > & parameters,
  //           const std::vector<NRLib::Well>                  & wells,
  //           MultiIntervalGrid                               * multiple_interval_grid,
  //           const ModelSettings                             * model_settings,
  //           std::string                                     & err_text);

  Background(FFTGrid ** grids);
  ~Background(void);

  //FFTGrid    * getAlpha(void) { return back_model_[0]; }
  //FFTGrid    * getBeta(void)  { return back_model_[1]; }

  //FFTGrid    * getVp(void)  { return back_model_[0]; }
  //FFTGrid    * getVs(void)  { return back_model_[1]; }
  //FFTGrid    * getRho(void) { return back_model_[2]; }

  //double       getMeanVsVp() const { return vsvp_;}

  //void         setClassicVsVp(); //For debugging purposes.

  //void         writeBackgrounds(const Simbox            * simbox,
  //                              GridMapping             * depthMapping,
  //                              const GridMapping       * timeMapping,
  //                              const bool                isFile,
  //                              const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)) const;

  //void         releaseGrids(); //backModel grids are now taken care of by other classes.

private:

  void         GenerateBackgroundModel(NRLib::Grid<float>                               * bg_vp,
                                       NRLib::Grid<float>                               * bg_vs,
                                       NRLib::Grid<float>                               * bg_rho,
                                       NRLib::Grid<float>                               * velociy,
                                       const Simbox                                     * simbox,
                                       const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                       const ModelSettings                              * model_settings,
                                       std::string                                      & err_text);

  //void         GenerateMultizoneBackgroundModel(NRLib::Grid<float>             * bg_vp,
  //                                              NRLib::Grid<float>             * bg_vs,
  //                                              NRLib::Grid<float>             * bg_rho,
  //                                              const std::vector<NRLib::Well> & wells,
  //                                              const Simbox                   * simbox,
  //                                              const ModelSettings            * modelSettings,
  //                                              const std::vector<std::string> & surface_files,
  //                                              std::string                    & err_text);

  //void         GenerateMultiIntervalBackgroundModel(std::vector<std::vector<NRLib::Grid<float> *> > & bg_vp,
  //                                                  const std::vector<NRLib::Well>                  & wells,
  //                                                  MultiIntervalGrid                               * multiple_interval_grid,
  //                                                  const ModelSettings                             * model_settings,
  //                                                  std::string                                     & err_text);

  void         ResampleBackgroundModel(NRLib::Grid<float>  * bg_vp,
                                       NRLib::Grid<float>  * bg_vs,
                                       NRLib::Grid<float>  * bg_rho,
                                       const Simbox        * bg_simbox,
                                       const Simbox        * simbox);

  //void         padAndSetBackgroundModel(FFTGrid * bgAlpha,
  //                                      FFTGrid * bgBeta,
  //                                      FFTGrid * bgRho);

  //void         createPaddedParameter(FFTGrid *& pNew,
  //                                   FFTGrid  * pOld);

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
                                    const int                                        & n_wells) const;

  //void         getKrigingWellTrendsZone(std::vector<BlockedLogsCommon *>     & blocked_logs,
  //                                      std::vector<std::vector<double> >    & bl_vp,
  //                                      std::vector<std::vector<double> >    & bl_vs,
  //                                      std::vector<std::vector<double> >    & bl_rho,
  //                                      std::vector<std::vector<double> >    & vt_vp,
  //                                      std::vector<std::vector<double> >    & vt_vs,
  //                                      std::vector<std::vector<double> >    & vt_rho,
  //                                      std::vector<const std::vector<int> > & ipos,
  //                                      std::vector<const std::vector<int> > & jpos,
  //                                      std::vector<const std::vector<int> > & kpos,
  //                                      std::vector<int>                     & n_blocks,
  //                                      int                                  & tot_blocks,
  //                                      const int                            & nz) const;

  void         GetWellTrends(std::vector<std::vector<double> >                & well_trend,
                             std::vector<std::vector<double> >                & high_cut_Well_trend,
                             const std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                             const int                                        & nz,
                             const std::string                                & name,
                             std::string                                      & err_text) const;

  //void         getWellTrendsZone(std::vector<BlockedLogsCommon *>  & bl,
  //                               std::vector<std::vector<double> > & wellTrend,
  //                               std::vector<std::vector<double> > & highCutWellTrend,
  //                               const std::vector<NRLib::Well>    & wells,
  //                               StormContGrid                     & eroded_zone,
  //                               const std::vector<bool>           & hitZone,
  //                               const int                         & nz,
  //                               const std::string                 & name,
  //                               const int                         & i,
  //                               std::string                       & err_text) const;

  //void        checkWellHitsZone(std::vector<bool>              & hitZone,
  //                              const std::vector<NRLib::Well> & wells,
  //                              StormContGrid                  & eroded_zone,
  //                              const int                      & nWells) const;

  void         WriteTrendsToFile(std::vector<double> & trend,
                                 const Simbox        * simbox,
                                 bool                  write1D,
                                 bool                  write3D,
                                 bool                  hasVelocityTrend,
                                 const std::string &   name,
                                 bool                  isFile);

  //void         writeMultizoneTrendsToFile(const std::vector<std::vector<double> > & vp_zones,
  //                                        const std::vector<std::vector<double> > & vs_zones,
  //                                        const std::vector<std::vector<double> > & rho_zones,
  //                                        std::vector<StormContGrid>              & alpha_trend_zone,
  //                                        std::vector<StormContGrid>              & beta_trend_zone,
  //                                        std::vector<StormContGrid>              & rho_trend_zone,
  //                                        const Simbox                            * simbox,
  //                                        const std::vector<int>                  & erosion_priority,
  //                                        const std::vector<Surface>              & surface,
  //                                        const std::vector<double>               & surface_uncertainty,
  //                                        const bool                                isFile) const;

  //void         writeMultiIntervalTrendsToFile(const std::vector<std::vector<double> >   & vp_zones,
  //                                            const std::vector<std::vector<double> >   & vs_zones,
  //                                            const std::vector<std::vector<double> >   & rho_zones,
  //                                            std::vector<StormContGrid>                & vp_trend_zone,
  //                                            std::vector<StormContGrid>                & vs_trend_zone,
  //                                            std::vector<StormContGrid>                & rho_trend_zone,
  //                                            MultiIntervalGrid                         * multiple_interval_grid,
  //                                            std::vector<const NRLib::Surface<double> *> surfaces,
  //                                            //const std::vector<double>                 & surface_uncertainty,
  //                                            const bool                                  isFile) const;

  const CovGrid2D & MakeCovGrid2D(const Simbox * simbox,
                                  Vario        * vario,
                                  int            debugFlag);

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
                                  const std::vector<const std::vector<int> *>  kpos) const;

  void         MakeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                    NRLib::Grid<float>               * bg_grid,
                                    std::vector<double>              & trend,
                                    const Simbox                     * simbox,
                                    const CovGrid2D                  & cov_grid_2D,
                                    const std::string                & type) const;

  //void         makeTrendZone(const std::vector<double> & trend,
  //                           StormContGrid             & trend_zone) const;

  //void         makeKrigedZone(const std::vector<KrigingData2D> & krigingData,
  //                            const std::vector<double>        & trend,
  //                            StormContGrid                    & kriged_zone,
  //                            const CovGrid2D                  & covGrid2D) const;

  //void         MakeMultizoneBackground(NRLib::Grid<float>               * bg_vp,
  //                                     NRLib::Grid<float>               * bg_vs,
  //                                     NRLib::Grid<float>               * bg_rho,
  //                                     const std::vector<StormContGrid> & vp_zones,
  //                                     const std::vector<StormContGrid> & vs_zones,
  //                                     const std::vector<StormContGrid> & rho_zones,
  //                                     const Simbox                     * simbox,
  //                                     const std::vector<int>           & erosion_priority,
  //                                     const std::vector<Surface>       & surface,
  //                                     const std::vector<double>        & surface_uncertainty,
  //                                     const std::string                & type) const;

  //void         MakeMultiIntervalBackground(std::vector<std::vector<NRLib::Grid<float> *> > & parameters,
  //                                         const std::vector<StormContGrid>                & vp_zones,
  //                                         const std::vector<StormContGrid>                & vs_zones,
  //                                         const std::vector<StormContGrid>                & rho_zones,
  //                                         MultiIntervalGrid                               * multiple_interval_grid,
  //                                         std::vector<const NRLib::Surface<double> *>       surfaces,
  //                                         //const std::vector<double>                       & surface_uncertainty,
  //                                         const std::string                               & type) const;

  void         CalculateVelocityDeviations(NRLib::Grid<float>                               * velocity,
                                           const Simbox                                     * simbox,
                                           const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                           std::vector<double>                              & trend_vel,
                                           std::vector<double>                              & avg_dev_vel,
                                           std::vector<double>                              & avg_dev_vp,
                                           //int                                                outputFlag,
                                           int                                                n_wells);

  void         ResampleParameter(NRLib::Grid<float> *& p_new, // Resample to
                                 NRLib::Grid<float> *  p_old, // Resample from
                                 const Simbox       *  simbox_new,
                                 const Simbox       *  simbox_old);

  void         CalculateVerticalTrend(std::vector<std::vector<double> > & wellTrend,
                                      std::vector<double>               & trend,
                                      float                               logMin,
                                      float                               logMax,
                                      float                               maxHz,
                                      int                                 nz,
                                      float                               dz,
                                      const std::string                 & name);

  void         WriteVerticalTrend(std::vector<double> & trend,
                                  float                 dz,
                                  int                   nz,
                                  std::string           name);

  void         CalculateDeviationFromVerticalTrend(std::vector<std::vector<double> > & wellTrend,
                                                   const std::vector<double>         & trend,
                                                   std::vector<double>               & avg_dev,
                                                   const int                           nd);

  void         WriteDeviationsFromVerticalTrend(const std::vector<double>                        & avg_dev_vp,
                                                const std::vector<double>                        & avg_dev_vs,
                                                const std::vector<double>                        & avg_dev_rho,
                                                const std::vector<double>                        & trend_vp,
                                                const std::vector<double>                        & trend_vs,
                                                const std::vector<double>                        & trend_rho,
                                                const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                                const int                                          n_wells,
                                                const int                                          nz);

  void         SmoothTrendWithLocalLinearRegression(std::vector<double> & trend,
                                                    int                 * count,
                                                    int                   nWells,
                                                    int                   nz,
                                                    float                 dz,
                                                    float                 min_value,
                                                    float                 max_value,
                                                    std::string           parName);

  void         FillInVerticalTrend(FFTGrid                   * bgTrend,
                                   const std::vector<double> & trend);

  //void         findMeanVsVp(FFTGrid * Vp,
  //                          FFTGrid * Vs);

  FFTGrid    * CopyFFTGrid(FFTGrid   * origGrid,
                           const bool  expTrans,
                           const bool  fileGrid) const;

  //void        ComputeZoneProbability(const std::vector<double>      & z,
  //                                   const std::vector<NRLib::Beta> & horizon_distributions,
  //                                   const std::vector<int>         & erosion_priority,
  //                                   std::vector<double>            & zone_probability) const;

  //void         ErodeSurface(Surface       *& surface,
  //                          const Surface *  priority_surface,
  //                          const Simbox  *  simbox,
  //                          const bool    &  compare_upward) const;

  //void         ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
  //                              const std::vector<int>     & erosion_priority,
  //                              const std::vector<Surface> & surface,
  //                              const Simbox               * simbox) const;

  //void         BuildErodedZones(StormContGrid                & eroded_zone,
  //                              const std::vector<Surface *> & eroded_surfaces,
  //                              const int                    & nz,
  //                              const Simbox                 * simbox,
  //                              const int                    & i) const;

  //void         BuildErodedIntervals(StormContGrid                & eroded_zone,
  //                                  //const std::vector<Surface>   & eroded_surfaces,
  //                                  const int                    & nz,
  //                                  const Simbox                 * simbox) const;

  //void         BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
  //                                       std::vector<StormContGrid> & vs_zones,
  //                                       std::vector<StormContGrid> & rho_zones,
  //                                       const std::vector<Surface> & surf,
  //                                       const std::vector<int>     & correlation_structure,
  //                                       const Simbox               * simbox,
  //                                       const float                & dz,
  //                                       std::string                & err_text) const;

  //void         BuildSeismicPropertyIntervals(std::vector<StormContGrid> & vp_zones,
  //                                           std::vector<StormContGrid> & vs_zones,
  //                                           std::vector<StormContGrid> & rho_zones,
  //                                           MultiIntervalGrid          * multiple_interval_grid,
  //                                           std::string                & err_text) const;

  //FFTGrid    * back_model_[3];       // Cubes for background model files.

  int          DataTarget_;         // Number of data requested in a kriging block
  //double       vsvp_;               // Average ratio between vs and vp.
};

#endif
