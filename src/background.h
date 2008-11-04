#ifndef BACKGROUND_H
#define BACKGROUND_H

#include <stdio.h>

#include "lib/global_def.h"

class WellData;
class Simbox;
class FFTGrid;
class ModelSettings;

//Special note on the use of Background:
//All pointers used here are also used externally, so no deletion happens. 

class Background
{
public:
  Background(FFTGrid       ** grids,
             WellData      ** wells,
             Simbox         * simbox,
             ModelSettings  * modelSettings);
  Background(FFTGrid ** grids);
  ~Background(void);

  FFTGrid    * getAlpha(void) { return backModel_[0]; }
  FFTGrid    * getBeta(void)  { return backModel_[1]; }
  FFTGrid    * getRho(void)   { return backModel_[2]; }
  double       getMeanVsVp()  { return vsvp_;}
  
  void         setClassicVsVp(); //For debugging purposes.

  void         writeBackgrounds(Simbox * simbox) const;

private:
  void         generateBackgroundModel(WellData      ** wells,
                                       Simbox         * simbox,
                                       ModelSettings  * modelSettings);
  void         calculateVerticalTrend(WellData   ** wells,
                                      float       * trend, 
                                      float         logMin,
                                      float         logMax,
                                      float         maxHz,
                                      int           nWells,
                                      int           nz,
                                      float         dz,
                                      const char  * name);
  void         writeVerticalTrend(float      * trend, 
                                  float        dz,
                                  int          nz,
                                  std::string  name);
  void         createTrendCube(Simbox      *  simbox,
                               FFTGrid     *& pFFTGrid, 
                               const float *  trend);
  void         calculateDeviationFromVerticalTrend(WellData    ** wells,
                                                   const float  * trend, 
                                                   float        * avg_dev,
                                                   int            nWells,
                                                   int            nd,
                                                   const char   * name);
  void         writeDeviationsFromVerticalTrend(WellData    ** wells,
                                                const float *  avg_dev_alpha,
                                                const float *  avg_dev_beta,
                                                const float *  avg_dev_rho,
                                                const float *  trend_alpha, 
                                                const float *  trend_beta, 
                                                const float *  trend_rho,
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
                                                    const char * parName);
  void         extrapolateTrend(const char * pName, 
                                float      * log,
                                int          nz);
  void         findMeanVsVp();

  FFTGrid    * backModel_[3];       // Cubes for background model files.
  int          DataTarget_;         // Number of data requested in a kriging block 
  double       vsvp_;               // Average ratio between vs and vp.
};

#endif
