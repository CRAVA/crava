/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "nrlib/random/beta.hpp"
#include "nrlib/well/well.hpp"

#include "src/blockedlogscommon.h"
#include "src/multiintervalgrid.h"

class Vario;
class Simbox;
class FFTGrid;
class CovGrid2D;
class WellData;
//class Well;
class GridMapping;
class KrigingData3D;
class KrigingData2D;
class ModelSettings;
class BlockedLogsForZone;
class BlockedLogsCommon;
class BlockedLogs;
class MultiIntervalGrid;


//Special note on the use of Background:
//All pointers used here are also used externally, so no deletion happens.

class Background
{
public:

  Background(std::vector<NRLib::Grid<double> >          & parameters,
             const std::vector<NRLib::Well>             & wells,
             NRLib::Grid<double>                        & velocity,
             const Simbox                               * time_simbox,
             const Simbox                               * time_bg_simbox,
             std::map<std::string, BlockedLogsCommon *> & bl,
             std::map<std::string, BlockedLogsCommon *> & bg_bl,
             const ModelSettings                        * modelSettings);

  Background(std::vector<NRLib::Grid<double> > & parameters,
             const std::vector<NRLib::Well>    & wells,
             const Simbox                      * timeSimbox,
             const ModelSettings               * modelSettings,
             const std::vector<std::string>    & surface_files,
             std::string                       & err_text);

  Background(std::vector<std::vector<NRLib::Grid<double> > > & parameters,
             const std::vector<NRLib::Well>                  & wells,
             MultiIntervalGrid                               * multiple_interval_grid,
             const ModelSettings                             * model_settings,
             std::string                                     & err_text);

  Background(FFTGrid ** grids);
  ~Background(void);

  //NRLib::Grid<double> & GetVp(void) { return back_model_[0]; }
  //NRLib::Grid<double> & GetVs(void) { return back_model_[1]; }
  //NRLib::Grid<double> & GetRho(void) { return back_model_[2]; }

  FFTGrid    * getAlpha(void) { return back_model_[0]; }
  FFTGrid    * getBeta(void)  { return back_model_[1]; }

  FFTGrid    * getVp(void)  { return back_model_[0]; }
  FFTGrid    * getVs(void)  { return back_model_[1]; }
  FFTGrid    * getRho(void) { return back_model_[2]; }

  double       getMeanVsVp() const { return vsvp_;}

  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(const Simbox            * simbox,
                                GridMapping             * depthMapping,
                                const GridMapping       * timeMapping,
                                const bool                isFile,
                                const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)) const;

  void         releaseGrids(); //backModel grids are now taken care of by other classes.

  //NRLib::Grid<double> FFTGridRealToGrid(const FFTGrid * fft_grid);

private:
  //void         generateBackgroundModel(FFTGrid                      *& bgAlpha,
  //                                     FFTGrid                      *& bgBeta,
  //                                     FFTGrid                      *& bgRho,
  //                                     FFTGrid                      *& velocity,
  //                                     const std::vector<WellData *> & wells,
  //                                     const Simbox                  * simbox,
  //                                     const ModelSettings           * modelSettings);

  void         generateBackgroundModel(NRLib::Grid<double>                        & bg_vp,
                                       NRLib::Grid<double>                        & bg_vs,
                                       NRLib::Grid<double>                        & bg_rho,
                                       NRLib::Grid<double>                        & velociy,
                                       const std::vector<NRLib::Well>             & wells,
                                       const Simbox                               * simbox,
                                       std::map<std::string, BlockedLogsCommon *> & bl,
                                       std::map<std::string, BlockedLogsCommon *> & bg_bl,
                                       const ModelSettings                        * model_settings);

  //void         generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
  //                                              FFTGrid                       *& bgBeta,
  //                                              FFTGrid                       *& bgRho,
  //                                              const std::vector<WellData *>  & wells,
  //                                              const Simbox                   * simbox,
  //                                              const ModelSettings            * modelSettings,
  //                                              const std::vector<std::string> & surface_files);

  void         GenerateMultizoneBackgroundModel(NRLib::Grid<double>            & bg_vp,//FFTGrid                       *& bgAlpha,
                                                NRLib::Grid<double>            & bg_vs,
                                                NRLib::Grid<double>            & bg_rho,
                                                const std::vector<NRLib::Well> & wells,
                                                const Simbox                   * simbox,
                                                const ModelSettings            * modelSettings,
                                                const std::vector<std::string> & surface_files,
                                                std::string                    & err_text);

  void         GenerateMultiIntervalBackgroundModel(std::vector<NRLib::Grid<double> > & bg_vp,
                                                    std::vector<NRLib::Grid<double> > & bg_vs,
                                                    std::vector<NRLib::Grid<double> > & bg_rho,
                                                    const std::vector<NRLib::Well>    & wells,
                                                    MultiIntervalGrid                 * multiple_interval_grid,
                                                    const ModelSettings               * model_settings,
                                                    std::string                       & err_text);

  void         resampleBackgroundModel(NRLib::Grid<double> & bg_vp, //FFTGrid      *& bgAlpha,
                                       NRLib::Grid<double> & bg_vs,
                                       NRLib::Grid<double> & bg_rho,
                                       const Simbox  * timeBGsimbox,
                                       const Simbox  * timeSimbox,
                                       const ModelSettings * modelSettings);
  void         padAndSetBackgroundModel(FFTGrid * bgAlpha,
                                        FFTGrid * bgBeta,
                                        FFTGrid * bgRho);
  //void         padAndSetBackgroundModelInterval(std::vector<FFTGrid *> bg_vp,
  //                                              std::vector<FFTGrid *> bg_vs,
  //                                              std::vector<FFTGrid *> bg_rho);
  void         createPaddedParameter(FFTGrid *& pNew,
                                     FFTGrid  * pOld);

  void         calculateBackgroundTrend(float              * trend,
                                        float              * avgDev,
                                        const int            nz,
                                        const float          dz,
                                        float                logMin,
                                        float                logMax,
                                        float                maxHz,
                                        std::vector<float *> wellTrend,
                                        std::vector<float *> highCutWellTrend,
                                        const std::string  & name);

  //void         getKrigingWellTrends(std::vector<float *>          & blAlpha,
  //                                  std::vector<float *>          & blBeta,
  //                                  std::vector<float *>          & blRho,
  //                                  std::vector<float *>          & vtAlpha,
  //                                  std::vector<float *>          & vtBeta,
  //                                  std::vector<float *>          & vtRho,
  //                                  std::vector<const int *>      & ipos,
  //                                  std::vector<const int *>      & jpos,
  //                                  std::vector<const int *>      & kpos,
  //                                  std::vector<int>              & nBlocks,
  //                                  int                           & totBlocks,
  //                                  const std::vector<WellData *> & wells,
  //                                  const int                     & nWells) const;

  void         getKrigingWellTrends(std::vector<std::vector<double> >          & blAlpha,
                                    std::vector<std::vector<double> >          & blBeta,
                                    std::vector<std::vector<double> >          & blRho,
                                    std::vector<std::vector<double> >          & vtAlpha,
                                    std::vector<std::vector<double> >          & vtBeta,
                                    std::vector<std::vector<double> >          & vtRho,
                                    std::vector<const std::vector<int> >       & ipos,
                                    std::vector<const std::vector<int> >       & jpos,
                                    std::vector<const std::vector<int> >       & kpos,
                                    std::vector<int>                           & nBlocks,
                                    int                                        & totBlocks,
                                    const std::vector<NRLib::Well>             & wells,
                                    std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                                    const int                                  & nWells) const;

  //void         getKrigingWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
  //                                      std::vector<float *>              & blAlpha,
  //                                      std::vector<float *>              & blBeta,
  //                                      std::vector<float *>              & blRho,
  //                                      std::vector<float *>              & vtAlpha,
  //                                      std::vector<float *>              & vtBeta,
  //                                      std::vector<float *>              & vtRho,
  //                                      std::vector<const int *>          & ipos,
  //                                      std::vector<const int *>          & jpos,
  //                                      std::vector<const int *>          & kpos,
  //                                      std::vector<int>                  & nBlocks,
  //                                      int                               & totBlocks,
  //                                      const int                         & nz) const;

  void         getKrigingWellTrendsZone(std::vector<BlockedLogsCommon *>     & blocked_logs,
                                        std::vector<std::vector<double> >    & bl_vp,
                                        std::vector<std::vector<double> >    & bl_vs,
                                        std::vector<std::vector<double> >    & bl_rho,
                                        std::vector<std::vector<double> >    & vt_vp,
                                        std::vector<std::vector<double> >    & vt_vs,
                                        std::vector<std::vector<double> >    & vt_rho,
                                        std::vector<const std::vector<int> > & ipos,
                                        std::vector<const std::vector<int> > & jpos,
                                        std::vector<const std::vector<int> > & kpos,
                                        std::vector<int>                     & n_blocks,
                                        int                                  & tot_blocks,
                                        const int                            & nz) const;

  //void         getWellTrends(std::vector<float *>          & wellTrend,
  //                           std::vector<float *>          & highCutWellTrend,
  //                           const std::vector<WellData *> & wells,
  //                           const int                     & nz,
  //                           const std::string             & name) const;

  void         getWellTrends(std::vector<float *>           & wellTrend,
                             std::vector<float *>           & highCutWellTrend,
                             const std::vector<NRLib::Well> & wells,
                             std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                             const int                      & nz,
                             const std::string              & name) const;

  //void         getWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
  //                               std::vector<float *>              & wellTrend,
  //                               std::vector<float *>              & highCutWellTrend,
  //                               const std::vector<WellData *>     & wells,
  //                               StormContGrid                     & eroded_zone,
  //                               const std::vector<bool>           & hitZone,
  //                               const int                         & nz,
  //                               const std::string                 & name,
  //                               const int                         & i) const;

  void         getWellTrendsZone(const ModelSettings               * modelSettings,
                                 std::vector<BlockedLogsCommon *>  & bl,
                                 std::vector<float *>              & wellTrend,
                                 std::vector<float *>              & highCutWellTrend,
                                 const std::vector<NRLib::Well>    & wells,
                                 StormContGrid                     & eroded_zone,
                                 //const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                                 const std::vector<bool>           & hitZone,
                                 const int                         & nz,
                                 const std::string                 & name,
                                 const int                         & i) const;

  //void        checkWellHitsZone(std::vector<bool>             & hitZone,
  //                              const std::vector<WellData *> & wells,
  //                              StormContGrid                 & eroded_zone,
  //                              const int                     & nWells) const;

  void        checkWellHitsZone(std::vector<bool>              & hitZone,
                                const std::vector<NRLib::Well> & wells,
                                StormContGrid                  & eroded_zone,
                                const int                      & nWells) const;

  void         writeTrendsToFile(float             * trend,
                                 const Simbox      * simbox,
                                 bool                write1D,
                                 bool                write3D,
                                 bool                hasVelocityTrend,
                                 const std::string & name,
                                 bool                isFile);

  void         writeMultizoneTrendsToFile(const std::vector<float *>   alpha_zones,
                                          const std::vector<float *>   beta_zones,
                                          const std::vector<float *>   rho_zones,
                                          std::vector<StormContGrid> & alpha_trend_zone,
                                          std::vector<StormContGrid> & beta_trend_zone,
                                          std::vector<StormContGrid> & rho_trend_zone,
                                          const Simbox               * simbox,
                                          const std::vector<int>     & erosion_priority,
                                          const std::vector<Surface> & surface,
                                          const std::vector<double>  & surface_uncertainty,
                                          const bool                   isFile) const;

  void         writeMultiIntervalTrendsToFile(const std::vector<float *>   vp_zones,
                                              const std::vector<float *>   vs_zones,
                                              const std::vector<float *>   rho_zones,
                                              std::vector<StormContGrid> & vp_trend_zone,
                                              std::vector<StormContGrid> & vs_trend_zone,
                                              std::vector<StormContGrid> & rho_trend_zone,
                                              //const Simbox             * simbox,
                                              //const std::vector<int>   & erosion_priority,
                                              MultiIntervalGrid          * multiple_interval_grid,
                                              //const std::vector<int>     & erosion_priority,
                                              //const std::vector<NRLib::Surface<double> > & surfaces,
                                              //std::vector<const NRLib::Surface<double>& > surface,
                                              std::vector<const NRLib::Surface<double> *> surfaces,
                                              const std::vector<double>  & surface_uncertainty,
                                              const bool                   isFile) const;

  const CovGrid2D & makeCovGrid2D(const Simbox * simbox,
                                  Vario        * vario,
                                  int            debugFlag);

  void         setupKrigingData2D(std::vector<KrigingData2D>     & kriging_data_vp,
                                  std::vector<KrigingData2D>     & kriging_data_vs,
                                  std::vector<KrigingData2D>     & kriging_data_rho,
                                  float                          * trend_vp,
                                  float                          * trend_vs,
                                  float                          * trend_rho,
                                  const int                        output_flag,
                                  const int                      & nz,
                                  const float                    & dz,
                                  const int                      & tot_blocks,
                                  const std::vector<int>         & n_blocks,
                                  const std::vector<std::vector<double> > & bl_vp,
                               //const std::vector<float *>     & blAlpha,
                                  const std::vector<std::vector<double> > & bl_vs,
                                  const std::vector<std::vector<double> > & bl_rho,
                                  const std::vector<std::vector<double> > & vt_vp,
                                  const std::vector<std::vector<double> > & vt_vs,
                                  const std::vector<std::vector<double> > & vt_rho,
                                  const std::vector<const std::vector<int> >   ipos,
                                  const std::vector<const std::vector<int> >   jpos,
                                  const std::vector<const std::vector<int> >   kpos) const;

  void         makeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                    NRLib::Grid<double>              & bg_grid,
                                    const float                      * trend,
                                    const Simbox                     * simbox,
                                    const CovGrid2D                  & cov_grid_2D,
                                    const std::string                & type) const;
                                    //bool                               isFile) const;

  void         makeTrendZone(const float   * trend,
                             StormContGrid & trend_zone) const;

  void         makeKrigedZone(const std::vector<KrigingData2D> & krigingData,
                              const float                      * trend,
                              StormContGrid                    & kriged_zone,
                              const CovGrid2D                  & covGrid2D) const;

  void         MakeMultizoneBackground(NRLib::Grid<double>              & bg_vp, //FFTGrid                                   *& bgAlpha,
                                       NRLib::Grid<double>              & bg_vs,
                                       NRLib::Grid<double>              & bg_rho,
                                       const std::vector<StormContGrid> & vp_zones,
                                       const std::vector<StormContGrid> & vs_zones,
                                       const std::vector<StormContGrid> & rho_zones,
                                       const Simbox                     * simbox,
                                       const std::vector<int>           & erosion_priority,
                                       const std::vector<Surface>       & surface,
                                       const std::vector<double>        & surface_uncertainty,
                                       const bool                         isFile,
                                       const std::string                & type) const;

  void         MakeMultiIntervalBackground(std::vector<NRLib::Grid<double> > & bg_vp,
                                           std::vector<NRLib::Grid<double> > & bg_vs,
                                           std::vector<NRLib::Grid<double> > & bg_rho,
                                           const std::vector<StormContGrid>  & vp_zones,
                                           const std::vector<StormContGrid>  & vs_zones,
                                           const std::vector<StormContGrid>  & rho_zones,
                                           MultiIntervalGrid                 * multiple_interval_grid,
                                           //const Simbox                     * simbox,
                                           //const std::vector<int>           & erosion_priority,
                                           //const std::vector<NRLib::Surface<double> > & surface,
                                           //std::vector<const NRLib::Surface<double>& > surface_test,
                                           std::vector<const NRLib::Surface<double> *> surfaces,
                                           const std::vector<double>        & surface_uncertainty,
                                           const bool                         is_file,
                                           const std::string                & type) const;

  //void         calculateVelocityDeviations(FFTGrid                       * velocity,
  //                                         const std::vector<WellData *> & wells,
  //                                         const Simbox                  * simbox,
  //                                         float                        *& trendVel,
  //                                         float                        *& avgDevVel,
  //                                         float                         * avgDevAlpha,
  //                                         int                             outputFlag,
  //                                         int                             nWells);

  void         calculateVelocityDeviations(NRLib::Grid<double>            & velocity, //FFTGrid                        * velocity,
                                           const std::vector<NRLib::Well> & wells,
                                           const Simbox                   * simbox,
                                           std::map<std::string, BlockedLogsCommon *> & bl,
                                           std::map<std::string, BlockedLogsCommon *> & bg_bl,
                                           float                         *& trendVel,
                                           float                         *& avgDevVel,
                                           float                          * avgDevAlpha,
                                           int                              outputFlag,
                                           int                              nWells);

  //void         resampleParameter(FFTGrid *& parameterNew,
  //                               FFTGrid  * parameterOld,
  //                               const Simbox   * simboxNew,
  //                               const Simbox   * simboxOld,
  //                               bool       isFile);

  void         resampleParameter(NRLib::Grid<double> & p_new, //FFTGrid *& pNew,        // Resample to
                                 NRLib::Grid<double> & p_old, //FFTGrid  * pOld,        // Resample from
                                 const Simbox        * simbox_new,
                                 const Simbox        * simbox_old);

  void         calculateVerticalTrend(std::vector<float *> wellTrend,
                                      float              * trend,
                                      float                logMin,
                                      float                logMax,
                                      float                maxHz,
                                      int                  nz,
                                      float                dz,
                                      const std::string  & name);

  void         writeVerticalTrend(float      * trend,
                                  float        dz,
                                  int          nz,
                                  std::string  name);
  void         calculateDeviationFromVerticalTrend(std::vector<float *>  wellTrend,
                                                   const float         * trend,
                                                   float               * avg_dev,
                                                   const int             nd);

  void         writeDeviationsFromVerticalTrend(const float                   * avg_dev_alpha,
                                                const float                   * avg_dev_beta,
                                                const float                   * avg_dev_rho,
                                                const float                   * trend_alpha,
                                                const float                   * trend_beta,
                                                const float                   * trend_rho,
                                                const std::vector<WellData *> & wells,
                                                const int                       nWells,
                                                const int                       nz);

  void         writeDeviationsFromVerticalTrend(const float                    * avg_dev_alpha,
                                                const float                    * avg_dev_beta,
                                                const float                    * avg_dev_rho,
                                                const float                    * trend_alpha,
                                                const float                    * trend_beta,
                                                const float                    * trend_rho,
                                                const std::vector<NRLib::Well> & wells,
                                                const int                        nWells,
                                                const int                        nz);

  void         smoothTrendWithLocalLinearRegression(float      * trend,
                                                    int        * count,
                                                    int          nWells,
                                                    int          nz,
                                                    float        dz,
                                                    float        min_value,
                                                    float        max_value,
                                                    std::string  parName);
  void         fillInVerticalTrend(FFTGrid     * bgTrend,
                                   const float * trend);
  void         findMeanVsVp(FFTGrid * Vp,
                            FFTGrid * Vs);
  FFTGrid    * copyFFTGrid(FFTGrid   * origGrid,
                           const bool  expTrans,
                           const bool  fileGrid) const;

  void        ComputeZoneProbability(const std::vector<double>      & z,
                                     const std::vector<NRLib::Beta> & horizon_distributions,
                                     const std::vector<int>         & erosion_priority,
                                     std::vector<double>            & zone_probability) const;

  void         ErodeSurface(Surface       *& surface,
                            const Surface *  priority_surface,
                            const Simbox  *  simbox,
                            const bool    &  compare_upward) const;

  void         ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
                                const std::vector<int>     & erosion_priority,
                                const std::vector<Surface> & surface,
                                const Simbox               * simbox) const;

  void         BuildErodedZones(StormContGrid                & eroded_zone,
                                const std::vector<Surface *> & eroded_surfaces,
                                const int                    & nz,
                                const Simbox                 * simbox,
                                const int                    & i) const;

  void         BuildErodedIntervals(StormContGrid                & eroded_zone,
                                    //const std::vector<Surface>   & eroded_surfaces,
                                    const int                    & nz,
                                    const Simbox                 * simbox) const;

  void         BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
                                         std::vector<StormContGrid> & vs_zones,
                                         std::vector<StormContGrid> & rho_zones,
                                         const std::vector<Surface> & surf,
                                         const std::vector<int>     & correlation_structure,
                                         const Simbox               * simbox,
                                         const float                & dz,
                                         std::string                & err_text) const;

  void         BuildSeismicPropertyIntervals(std::vector<StormContGrid> & vp_zones,
                                             std::vector<StormContGrid> & vs_zones,
                                             std::vector<StormContGrid> & rho_zones,
                                             MultiIntervalGrid          * multiple_interval_grid,
                                             std::string                & err_text) const;

  //std::vector<FFTGrid *> back_model_interval_[3];
  FFTGrid    * back_model_[3];       // Cubes for background model files.
  //std::vector<NRLib::Grid<double> > back_model_;

  int          DataTarget_;         // Number of data requested in a kriging block
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
