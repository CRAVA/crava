#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "lib/global_def.h"

class Vario;
class Simbox;
class FFTGrid;
class WellData;
class GridMapping;
class KrigingData3D;
class ModelSettings;

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
  Background(FFTGrid ** grids);
  ~Background(void);

  FFTGrid    * getAlpha(void) { return backModel_[0]; }
  FFTGrid    * getBeta(void)  { return backModel_[1]; }
  FFTGrid    * getRho(void)   { return backModel_[2]; }
  double       getMeanVsVp()  { return vsvp_;}
  
  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(Simbox * simbox, GridMapping * depthMapping, GridMapping * timeMapping) const;

private:
  void         generateBackgroundModel(FFTGrid      *& bgAlpha,
                                       FFTGrid      *& bgBeta,
                                       FFTGrid      *& bgRho,
                                       FFTGrid      *& velocity,
                                       WellData     ** wells,
                                       Simbox        * simbox,
                                       ModelSettings * modelSettings);
  void         calculateBackgroundTrend(FFTGrid          *& trendGrid,
                                        WellData         ** wells,
                                        Simbox            * simbox,
                                        float             * trend,
                                        float             * avgDev,
                                        float               logMin, 
                                        float               logMax,
                                        float               maxHz, 
                                        int                 outputFlag,
                                        int                 nWells, 
                                        bool                hasVelocityTrend,
                                        const std::string & name);
  void         setupKrigingData(KrigingData3D *& krigingData,
                                WellData      ** wells,
                                Simbox         * simbox,
                                float          * trendAlpha,
                                float          * trendBeta, 
                                float          * trendRho , 
                                const int        nWells);
  void         interpolateBackgroundTrend(KrigingData3D * krigingData,
                                          FFTGrid       * bgAlpha,
                                          FFTGrid       * bgBeta,
                                          FFTGrid       * bgRho,
                                          Simbox        * simbox,
                                          Vario         * vario,
                                          int             debugFlag);
  void         calculateVelocityDeviations(FFTGrid   * velocity,
                                           WellData ** wells,
                                           Simbox    * simbox,                               
                                           float     * trendVel,
                                           float     * avgDevVel,
                                           float     * avgDevAlpha,
                                           float       logMin, 
                                           float       logMax,
                                           int         outputFlag,
                                           int         nWells);
  void         resampleParameter(FFTGrid *& parameterNew,
                                 FFTGrid  * parameterOld,
                                 Simbox   * simboxNew,
                                 Simbox   * simboxOld);
  void         calculateVerticalTrend(WellData         ** wells,
                                      float             * trend, 
                                      float               logMin,
                                      float               logMax,
                                      float               maxHz,
                                      int                 nWells,
                                      int                 nz,
                                      float               dz,
                                      const std::string & name);
  void         writeVerticalTrend(float      * trend, 
                                  float        dz,
                                  int          nz,
                                  std::string  name);
  void         calculateDeviationFromVerticalTrend(WellData    ** wells,
                                                   const float  * trend, 
                                                   float        * avg_dev,
                                                   int            nWells,
                                                   int            nd,
                                                   std::string    name);
  void         writeDeviationsFromVerticalTrend(const float *  avg_dev_alpha,
                                                const float *  avg_dev_beta,
                                                const float *  avg_dev_rho,
                                                const float *  trend_alpha, 
                                                const float *  trend_beta, 
                                                const float *  trend_rho,
                                                WellData    ** wells,
                                                const int      nWells,
                                                const int      nz);
  void         smoothTrendWithMovingAverage(float * trend, 
                                            int   * count,
                                            int     nWells,
                                            int     nz);
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
  void         extrapolateTrend(std::string  pName, 
                                float      * log,
                                int          nz);
  void         findMeanVsVp(FFTGrid * Vp,
                            FFTGrid * Vs);

  FFTGrid    * backModel_[3];       // Cubes for background model files.
  int          DataTarget_;         // Number of data requested in a kriging block 
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
