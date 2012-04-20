#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

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
  Background(FFTGrid       ** grids,
             WellData      ** wells,
             FFTGrid       *& velocity,
             Simbox         * timeSimbox,
             Simbox         * timeBGSimbox,
             ModelSettings  * modelSettings);
  Background(FFTGrid                       ** grids,
             WellData                      ** wells,
             Simbox                         * timeSimbox,
             ModelSettings                  * modelSettings,
             const std::vector<std::string> & correlation_files);
  Background(FFTGrid ** grids);
  ~Background(void);

  FFTGrid    * getAlpha(void) { return backModel_[0]; }
  FFTGrid    * getBeta(void)  { return backModel_[1]; }
  FFTGrid    * getRho(void)   { return backModel_[2]; }
  double       getMeanVsVp()  { return vsvp_;}

  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(Simbox                  * simbox,
                                GridMapping             * depthMapping,
                                GridMapping             * timeMapping,
                                const bool                isFile,
                                const TraceHeaderFormat & thf = TraceHeaderFormat(TraceHeaderFormat::SEISWORKS)) const;

  void         releaseGrids(); //backModel grids are now taken care of by other classes.

private:
  void         generateBackgroundModel(FFTGrid      *& bgAlpha,
                                       FFTGrid      *& bgBeta,
                                       FFTGrid      *& bgRho,
                                       FFTGrid      *& velocity,
                                       WellData     ** wells,
                                       Simbox        * simbox,
                                       ModelSettings * modelSettings);

  void         generateMultizoneBackgroundModel(FFTGrid                       *& bgAlpha,
                                                FFTGrid                       *& bgBeta,
                                                FFTGrid                       *& bgRho,
                                                WellData                      ** wells,
                                                Simbox                         * simbox,
                                                ModelSettings                  * modelSettings,
                                                const std::vector<std::string> & surface_files);

  void         resampleBackgroundModel(FFTGrid      *& bgAlpha,
                                       FFTGrid      *& bgBeta,
                                       FFTGrid      *& bgRho,
                                       Simbox        * timeBGsimbox,
                                       Simbox        * timeSimbox,
                                       ModelSettings * modelSettings);
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

  void         getKrigingWellTrends(std::vector<float *>     & blAlpha,
                                    std::vector<float *>     & blBeta,
                                    std::vector<float *>     & blRho,
                                    std::vector<float *>     & vtAlpha,
                                    std::vector<float *>     & vtBeta,
                                    std::vector<float *>     & vtRho,
                                    std::vector<const int *> & ipos,
                                    std::vector<const int *> & jpos,
                                    std::vector<const int *> & kpos,
                                    std::vector<int>         & nBlocks,
                                    int                      & totBlocks,
                                    WellData                ** wells,
                                    const int                & nWells) const;

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
                                        int                               & totBlocks) const;

  void         getWellTrends(std::vector<float *> wellTrend,
                             std::vector<float *> highCutWellTrend,
                             WellData          ** wells,
                             const int          & nWells,
                             const std::string  & name) const;

  void         getWellTrendsZone(std::vector<BlockedLogsForZone *> & bl,
                                 std::vector<float *>              & wellTrend,
                                 std::vector<float *>              & highCutWellTrend,
                                 WellData                         ** wells,
                                 StormContGrid                     * eroded_zone,
                                 const std::string                  & name,
                                 const int                          & i) const;

  void         writeTrendsToFile(float             * trend,
                                 Simbox            * simbox,
                                 bool                write1D,
                                 bool                write3D,
                                 bool                hasVelocityTrend,
                                 const std::string & name,
                                 bool                isFile);
  const CovGrid2D & makeCovGrid2D(Simbox * simbox,
                                  Vario  * vario,
                                  int      debugFlag);
  void         setupKrigingData2D(std::vector<KrigingData2D> & krigingDataAlpha,
                                  std::vector<KrigingData2D> & krigingDataBeta,
                                  std::vector<KrigingData2D> & krigingDataRho,
                                  float                      * trendAlpha,
                                  float                      * trendBeta,
                                  float                      * trendRho,
                                  int                          outputFlag,
                                  const int                  & nz,
                                  const float                & dz,
                                  const int                  & totBlocks,
                                  std::vector<int>           & nBlocks,
                                  std::vector<float *>       & blAlpha,
                                  std::vector<float *>       & blBeta,
                                  std::vector<float *>       & blRho,
                                  std::vector<float *>       & vtAlpha,
                                  std::vector<float *>       & vtBeta,
                                  std::vector<float *>       & vtRho,
                                  std::vector<const int *>    ipos,
                                  std::vector<const int *>    jpos,
                                  std::vector<const int *>    kpos) const;

  void         makeKrigedBackground(const std::vector<KrigingData2D> & krigingData,
                                    FFTGrid                         *& bgGrid,
                                    float                            * trend,
                                    Simbox                           * simbox,
                                    const CovGrid2D                  & covGrid2D,
                                    const std::string                & type,
                                    bool                               isFile) const;

  void         makeKrigedBackgroundZone(const std::vector<KrigingData2D> & krigingData,
                                        float                            * trend,
                                        StormContGrid                    * correlation_zone,
                                        const CovGrid2D                  & covGrid2D) const;

  void         calculateVelocityDeviations(FFTGrid   * velocity,
                                           WellData ** wells,
                                           Simbox    * simbox,
                                           float    *& trendVel,
                                           float    *& avgDevVel,
                                           float     * avgDevAlpha,
                                           int         outputFlag,
                                           int         nWells);
  void         resampleParameter(FFTGrid *& parameterNew,
                                 FFTGrid  * parameterOld,
                                 Simbox   * simboxNew,
                                 Simbox   * simboxOld,
                                 bool       isFile);
  void         calculateVerticalTrend(std::vector<float *> wellTrend,
                                      float              * trend,
                                      float                logMin,
                                      float                logMax,
                                      float                maxHz,
                                      int                  nz,
                                      float                dz,
                                      const std::string  & name);

  void         calculateVerticalTrendZone(WellData         ** wells,
                                          float             * trend,
                                          float               logMin,
                                          float               logMax,
                                          float               maxHz,
                                          int                 nWells,
                                          int                 nz,
                                          float               dz,
                                          const std::string & name) const;

  void         writeVerticalTrend(float      * trend,
                                  float        dz,
                                  int          nz,
                                  std::string  name);
  void         calculateDeviationFromVerticalTrend(std::vector<float *>  wellTrend,
                                                   const float         * trend,
                                                   float               * avg_dev,
                                                   const int             nd);

  void         writeDeviationsFromVerticalTrend(const float *  avg_dev_alpha,
                                                const float *  avg_dev_beta,
                                                const float *  avg_dev_rho,
                                                const float *  trend_alpha,
                                                const float *  trend_beta,
                                                const float *  trend_rho,
                                                WellData    ** wells,
                                                const int      nWells,
                                                const int      nz);
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

  void         ErodeSurface(Surface       *& surface,
                            const Surface  * priority_surface,
                            const Simbox   * simbox,
                            const bool     & compare_upward) const;

  void         BuildErodedZones(std::vector<StormContGrid *> & eroded_zones,
                                const std::vector<int>       & erosion_priority,
                                const std::vector<Surface>   & surf,
                                const std::vector<int>       & nz_zone,
                                const Simbox                 * simbox) const;

  void         BuildCorrelationZones(std::vector<StormContGrid *> & correlation_zones,
                                     const std::vector<Surface>   & surf,
                                     const std::vector<int>       & correlation_structure,
                                     const Simbox                 * simbox) const;

  FFTGrid    * backModel_[3];       // Cubes for background model files.
  int          DataTarget_;         // Number of data requested in a kriging block
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
