/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "nrlib/random/beta.hpp"
#include "nrlib/well/well.hpp"

#include "src/blockedlogscommon.h"

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

class BlockedLogs;


//Special note on the use of Background:
//All pointers used here are also used externally, so no deletion happens.

class Background
{
public:
  //Background(FFTGrid                       ** grids,
  //           const std::vector<WellData *>  & wells,
  //           FFTGrid                       *& velocity,
  //           const Simbox                   * timeSimbox,
  //           const Simbox                   * timeBGSimbox,
  //           const ModelSettings            * modelSettings);

  Background(FFTGrid                          ** grids, //std::vector<NRLib::Grid<double> > & grids,
             const std::vector<NRLib::Well>    & wells,
             FFTGrid                          *& velocity,
             const Simbox                      * time_simbox,
             const Simbox                      * time_bg_simbox,
             std::map<std::string, BlockedLogsCommon *> & bl,
             std::map<std::string, BlockedLogsCommon *> & bg_bl,
             const ModelSettings               * modelSettings);

  //Background(FFTGrid                       ** grids,
  //           const std::vector<WellData *>  & wells,
  //           const Simbox                   * timeSimbox,
  //           const ModelSettings            * modelSettings,
  //           const std::vector<std::string> & correlation_files);

  Background(FFTGrid                          ** grids, //std::vector<NRLib::Grid<double> > & grids, //FFTGrid                       ** grids,
             const std::vector<NRLib::Well>    & wells,
             const Simbox                      * timeSimbox,
             const ModelSettings               * modelSettings,
             const std::vector<std::string>    & surface_files);
             //std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs);

  Background(FFTGrid ** grids);
  ~Background(void);

  //NRLib::Grid<double> & GetVp(void) { return back_model_[0]; }
  //NRLib::Grid<double> & GetVs(void) { return back_model_[1]; }
  //NRLib::Grid<double> & GetRho(void) { return back_model_[2]; }

  FFTGrid    * getVp(void) { return back_model_[0]; }
  FFTGrid    * getVs(void)  { return back_model_[1]; }
  FFTGrid    * getRho(void)   { return back_model_[2]; }

  double       getMeanVsVp() const { return vsvp_;}

  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(const Simbox            * simbox,
                                GridMapping             * depthMapping,
                                const GridMapping       * timeMapping,
                                const bool                isFile,
                                const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)) const;

  void         releaseGrids(); //backModel grids are now taken care of by other classes.

private:
  //void         generateBackgroundModel(FFTGrid                      *& bgAlpha,
  //                                     FFTGrid                      *& bgBeta,
  //                                     FFTGrid                      *& bgRho,
  //                                     FFTGrid                      *& velocity,
  //                                     const std::vector<WellData *> & wells,
  //                                     const Simbox                  * simbox,
  //                                     const ModelSettings           * modelSettings);

  void         generateBackgroundModel(FFTGrid                      *& bg_vp,
                                       FFTGrid                      *& bg_vs,
                                       FFTGrid                      *& bg_rho,
                                       FFTGrid                       *& velocity,
                                       const std::vector<NRLib::Well> & wells,
                                       const Simbox                   * simbox,
                                       std::map<std::string, BlockedLogsCommon *> & bl,
                                       std::map<std::string, BlockedLogsCommon *> & bg_bl,
                                       const ModelSettings            * modelSettings);

  //void         generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
  //                                              FFTGrid                       *& bgBeta,
  //                                              FFTGrid                       *& bgRho,
  //                                              const std::vector<WellData *>  & wells,
  //                                              const Simbox                   * simbox,
  //                                              const ModelSettings            * modelSettings,
  //                                              const std::vector<std::string> & surface_files);

  void         GenerateMultizoneBackgroundModel(FFTGrid                                   *& bgAlpha,
                                                FFTGrid                                   *& bgBeta,
                                                FFTGrid                                   *& bgRho,
                                                const std::vector<NRLib::Well> & wells,
                                                const Simbox                   * simbox,
                                                const ModelSettings            * modelSettings,
                                                const std::vector<std::string> & surface_files);

  void         resampleBackgroundModel(FFTGrid      *& bgAlpha,
                                       FFTGrid      *& bgBeta,
                                       FFTGrid      *& bgRho,
                                       const Simbox  * timeBGsimbox,
                                       const Simbox  * timeSimbox,
                                       const ModelSettings * modelSettings);
  void         padAndSetBackgroundModel(FFTGrid * bgAlpha,
                                        FFTGrid * bgBeta,
                                        FFTGrid * bgRho);
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
                                    //std::vector<const int *>       & ipos,
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

  void         getKrigingWellTrendsZone(std::vector<BlockedLogsCommon *>  & blocked_logs,
                                     std::vector<std::vector<double> >           & bl_vp,
                                     std::vector<std::vector<double> >           & bl_vs,
                                     std::vector<std::vector<double> >           & bl_rho,
                                     std::vector<std::vector<double> >           & vt_vp,
                                     std::vector<std::vector<double> >           & vt_vs,
                                     std::vector<std::vector<double> >           & vt_rho,
                                     //std::vector<float *>              & blAlpha,
                                     //std::vector<float *>              & blBeta,
                                     //std::vector<float *>              & blRho,
                                     //std::vector<float *>              & vtAlpha,
                                     //std::vector<float *>              & vtBeta,
                                     //std::vector<float *>              & vtRho,
                                     std::vector<const std::vector<int> >       & ipos,
                                     std::vector<const std::vector<int> >       & jpos,
                                     std::vector<const std::vector<int> >       & kpos,
                                     //std::vector<const int *>          & ipos,
                                     //std::vector<const int *>          & jpos,
                                     //std::vector<const int *>          & kpos,
                                     std::vector<int>                  & n_blocks,
                                     int                               & tot_blocks,
                                     const int                         & nz) const;

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

  void         makeKrigedBackground(const std::vector<KrigingData2D> & krigingData,
                                    FFTGrid                         *& bgGrid,
                                    const float                      * trend,
                                    const Simbox                     * simbox,
                                    const CovGrid2D                  & covGrid2D,
                                    const std::string                & type,
                                    bool                               isFile) const;

  void         makeTrendZone(const float   * trend,
                             StormContGrid & trend_zone) const;

  void         makeKrigedZone(const std::vector<KrigingData2D> & krigingData,
                              const float                      * trend,
                              StormContGrid                    & kriged_zone,
                              const CovGrid2D                  & covGrid2D) const;

  void         MakeMultizoneBackground(FFTGrid                                   *& bgAlpha,
                                       FFTGrid                                   *& bgBeta,
                                       FFTGrid                                   *& bgRho,
                                       const std::vector<StormContGrid> & vp_zones,
                                       const std::vector<StormContGrid> & vs_zones,
                                       const std::vector<StormContGrid> & rho_zones,
                                       const Simbox                     * simbox,
                                       const std::vector<int>           & erosion_priority,
                                       const std::vector<Surface>       & surface,
                                       const std::vector<double>        & surface_uncertainty,
                                       const bool                         isFile,
                                       const std::string                & type) const;

  //void         calculateVelocityDeviations(FFTGrid                       * velocity,
  //                                         const std::vector<WellData *> & wells,
  //                                         const Simbox                  * simbox,
  //                                         float                        *& trendVel,
  //                                         float                        *& avgDevVel,
  //                                         float                         * avgDevAlpha,
  //                                         int                             outputFlag,
  //                                         int                             nWells);

  void         calculateVelocityDeviations(FFTGrid                        * velocity,
                                           const std::vector<NRLib::Well> & wells,
                                           const Simbox                   * simbox,
                                           std::map<std::string, BlockedLogsCommon *> & bl,
                                           std::map<std::string, BlockedLogsCommon *> & bg_bl,
                                           float                         *& trendVel,
                                           float                         *& avgDevVel,
                                           float                          * avgDevAlpha,
                                           int                              outputFlag,
                                           int                              nWells);

  void         resampleParameter(FFTGrid *& parameterNew,
                                 FFTGrid  * parameterOld,
                                 const Simbox   * simboxNew,
                                 const Simbox   * simboxOld,
                                 bool       isFile);
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

  void         BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
                                         std::vector<StormContGrid> & vs_zones,
                                         std::vector<StormContGrid> & rho_zones,
                                         const std::vector<Surface> & surf,
                                         const std::vector<int>     & correlation_structure,
                                         const Simbox               * simbox,
                                         const float                & dz) const;

  FFTGrid    * back_model_[3];       // Cubes for background model files.
  //std::vector<NRLib::Grid<double> > back_model_;

  int          DataTarget_;         // Number of data requested in a kriging block
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
