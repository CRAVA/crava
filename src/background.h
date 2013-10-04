/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "nrlib/random/beta.hpp"

class Vario;
class Simbox;
class FFTGrid;
class CovGrid2D;
class WellData;
class GridMapping;
class KrigingData3D;
class KrigingData2D;
class ModelSettings;
class BlockedLogsForZone;

//Special note on the use of Background:
//All pointers used here are also used externally, so no deletion happens.

class Background
{
public:
  Background(FFTGrid                       ** grids,
             const std::vector<WellData *>  & wells,
             FFTGrid                       *& velocity,
             const Simbox                   * timeSimbox,
             const Simbox                   * timeBGSimbox,
             const ModelSettings            * modelSettings);

  Background(FFTGrid                       ** grids,
             const std::vector<WellData *>  & wells,
             const Simbox                   * timeSimbox,
             const ModelSettings            * modelSettings,
             const std::vector<std::string> & correlation_files);
  Background(FFTGrid ** grids);
  ~Background(void);

  FFTGrid    * getAlpha(void) { return backModel_[0]; }
  FFTGrid    * getBeta(void)  { return backModel_[1]; }
  FFTGrid    * getRho(void)   { return backModel_[2]; }
  double       getMeanVsVp() const { return vsvp_;}

  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(const Simbox            * simbox,
                                GridMapping             * depthMapping,
                                const GridMapping       * timeMapping,
                                const bool                isFile,
                                const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)) const;

  void         releaseGrids(); //backModel grids are now taken care of by other classes.

  static void  createPaddedParameter(FFTGrid *& pNew,
                                     FFTGrid  * pOld);
private:
  void         generateBackgroundModel(FFTGrid              *& bgAlpha,
                                       FFTGrid              *& bgBeta,
                                       FFTGrid              *& bgRho,
                                       FFTGrid              *& velocity,
                                       const std::vector<WellData *> & wells,
                                       const Simbox          * simbox,
                                       const ModelSettings   * modelSettings);

  void         generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
                                                FFTGrid                       *& bgBeta,
                                                FFTGrid                       *& bgRho,
                                                const std::vector<WellData *>  & wells,
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

  void         getKrigingWellTrends(std::vector<float *>          & blAlpha,
                                    std::vector<float *>          & blBeta,
                                    std::vector<float *>          & blRho,
                                    std::vector<float *>          & vtAlpha,
                                    std::vector<float *>          & vtBeta,
                                    std::vector<float *>          & vtRho,
                                    std::vector<const int *>      & ipos,
                                    std::vector<const int *>      & jpos,
                                    std::vector<const int *>      & kpos,
                                    std::vector<int>              & nBlocks,
                                    int                           & totBlocks,
                                    const std::vector<WellData *> & wells,
                                    const int                     & nWells) const;

  void         getKrigingWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
                                        std::vector<float *>              & blAlpha,
                                        std::vector<float *>              & blBeta,
                                        std::vector<float *>              & blRho,
                                        std::vector<float *>              & vtAlpha,
                                        std::vector<float *>              & vtBeta,
                                        std::vector<float *>              & vtRho,
                                        std::vector<const int *>          & ipos,
                                        std::vector<const int *>          & jpos,
                                        std::vector<const int *>          & kpos,
                                        std::vector<int>                  & nBlocks,
                                        int                               & totBlocks,
                                        const int                         & nz) const;

  void         getWellTrends(std::vector<float *>          & wellTrend,
                             std::vector<float *>          & highCutWellTrend,
                             const std::vector<WellData *> & wells,
                             const int                     & nz,
                             const std::string             & name) const;

  void         getWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
                                 std::vector<float *>              & wellTrend,
                                 std::vector<float *>              & highCutWellTrend,
                                 const std::vector<WellData *>     & wells,
                                 StormContGrid                     & background_zone,
                                 NRLib::Volume                     & eroded_zone,
                                 const std::vector<bool>           & hitZone,
                                 const int                         & nz,
                                 const std::string                 & name,
                                 const int                         & i) const;

  void        checkWellHitsZone(std::vector<bool>             & hitZone,
                                const std::vector<WellData *> & wells,
                                NRLib::Volume                 & eroded_zone,
                                const int                     & nWells) const;

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

  void         setupKrigingData2D(std::vector<KrigingData2D>     & krigingDataAlpha,
                                  std::vector<KrigingData2D>     & krigingDataBeta,
                                  std::vector<KrigingData2D>     & krigingDataRho,
                                  float                          * trendAlpha,
                                  float                          * trendBeta,
                                  float                          * trendRho,
                                  const int                        outputFlag,
                                  const int                      & nz,
                                  const float                    & dz,
                                  const int                      & totBlocks,
                                  const std::vector<int>         & nBlocks,
                                  const std::vector<float *>     & blAlpha,
                                  const std::vector<float *>     & blBeta,
                                  const std::vector<float *>     & blRho,
                                  const std::vector<float *>     & vtAlpha,
                                  const std::vector<float *>     & vtBeta,
                                  const std::vector<float *>     & vtRho,
                                  const std::vector<const int *>   ipos,
                                  const std::vector<const int *>   jpos,
                                  const std::vector<const int *>   kpos) const;

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

  void         MakeMultizoneBackground(FFTGrid                         *& bgAlpha,
                                       FFTGrid                         *& bgBeta,
                                       FFTGrid                         *& bgRho,
                                       const std::vector<StormContGrid> & alpha_zones,
                                       const std::vector<StormContGrid> & beta_zones,
                                       const std::vector<StormContGrid> & rho_zones,
                                       const Simbox                     * simbox,
                                       const std::vector<int>           & erosion_priority,
                                       const std::vector<Surface>       & surface,
                                       const std::vector<double>        & surface_uncertainty,
                                       const bool                         isFile,
                                       const std::string                & type) const;

  void         calculateVelocityDeviations(FFTGrid               * velocity,
                                           const std::vector<WellData *> & wells,
                                           const Simbox          * simbox,
                                           float                *& trendVel,
                                           float                *& avgDevVel,
                                           float                 * avgDevAlpha,
                                           int                     outputFlag,
                                           int                     nWells);
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

  void         RegularizeZoneSurfaces(const std::vector<Surface> & surface,
                                      const Simbox * simbox,
                                      std::vector<Surface>       & regular_surface);

  void         ErodeSurface(Surface       *& surface,
                            const Surface *  priority_surface,
                            const bool    &  compare_upward) const;

  void         ErodeAllSurfaces(std::vector<Surface *>     & eroded_surfaces,
                                const std::vector<int>     & erosion_priority,
                                const std::vector<Surface> & surface) const;

  void         BuildErodedZones(NRLib::Volume                & eroded_zone,
                                const std::vector<Surface *> & eroded_surfaces,
                                const Simbox                 * simbox,
                                const int                    & i) const;

  void         BuildSeismicPropertyZones(std::vector<StormContGrid> & vp_zones,
                                         std::vector<StormContGrid> & vs_zones,
                                         std::vector<StormContGrid> & rho_zones,
                                         const std::vector<Surface> & surf,
                                         const std::vector<int>     & correlation_structure,
                                         const Simbox               * simbox,
                                         const float                & dz) const;

  FFTGrid    * backModel_[3];       // Cubes for background model files.
  int          DataTarget_;         // Number of data requested in a kriging block
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
