#ifndef BLOCKEDLOGS_H
#define BLOCKEDLOGS_H

#include "lib/log.h"
#include <stdlib.h>


class RandomGen;
class FFTGrid;
class WellData;
class Simbox;

class BlockedLogs
{
public:
  BlockedLogs(WellData  * well, 
              Simbox    * simbox,
              RandomGen * random);
  ~BlockedLogs(void);

  const char   * getWellname(void)                  const { return wellname_ ;} 
  const int      getNumberOfBlocks(void)            const { return nBlocks_ ;}  
  const double * getXpos(void)                      const { return xpos_ ;}
  const double * getYpos(void)                      const { return ypos_ ;}
  const double * getZpos(void)                      const { return zpos_ ;}
  const int    * getIpos(void)                      const { return ipos_ ;}
  const int    * getJpos(void)                      const { return jpos_ ;}
  const int    * getKpos(void)                      const { return kpos_ ;}
  const float  * getAlpha(void)                     const { return alpha_ ;}
  const float  * getBeta(void)                      const { return beta_ ;}
  const float  * getRho(void)                       const { return rho_ ;}
  const int    * getFacies(void)                    const { return facies_ ;}
  const float  * getAlphaBackgroundResolution(void) const { return alpha_background_resolution_ ;} 
  const float  * getBetaBackgroundResolution(void)  const { return beta_background_resolution_ ;}
  const float  * getRhoBackgroundResolution(void)   const { return rho_background_resolution_ ;}
  const float  * getAlphaSeismicResolution(void)    const { return alpha_seismic_resolution_ ;}
  const float  * getBetaSeismicResolution(void)     const { return beta_seismic_resolution_ ;} 
  const float  * getRhoSeismicResolution(void)      const { return rho_seismic_resolution_ ;}  
  void           getVerticalTrend(const float * blockedLog, float * trend);
  void           getVerticalTrend(const int * blockedLog,int * trend, RandomGen * random);
  void           getBlockedGrid(FFTGrid * grid, float * blockedLog);
  void           writeToFile(float dz, int type, bool exptrans) const;

private:
  void           blockWell(WellData  * well, 
                           Simbox    * simbox,
                           RandomGen * random);
  void           blockContinuousLog(const int   *  bInd,
                                    const float *  wellLog,
                                    float       *& blockedLog);
  void           blockDiscreteLog(const int * bInd,
                                  const int *  wellLog,
                                  const int *  faciesNumbers,
                                  int          nFacies,
                                  int       *& blockedLog,
                                  RandomGen * random);
  int            findMostProbable(const int * count,
                                  int         nFacies,
                                  RandomGen * random); 
  void           findSizeAndBlockPointers(WellData * well,
                                          Simbox   * simbox,
                                          int      * bInd);
  void           findBlockIJK(WellData  * well,
                              Simbox    * simbox,
                              const int * bInd);
  void           findBlockXYZ(Simbox * simbox);

  char        * wellname_;                    //| Name of well   

  double      * xpos_;                        //|
  double      * ypos_;                        //| Simbox XYZ value for block
  double      * zpos_;                        //|

  int         * ipos_;                        //|
  int         * jpos_;                        //| Simbox IJK value for block
  int         * kpos_;                        //|

  float       * alpha_;                       //|
  float       * beta_;                        //| Raw logs (log-domain)
  float       * rho_;                         //|
  int 	      * facies_;                      //| Facies numbers using *internal* numbering

  float       * alpha_background_resolution_; //| 
  float       * beta_background_resolution_;  //| Logs filtered to background resolution (log-domain)
  float       * rho_background_resolution_;   //| 

  float       * alpha_seismic_resolution_;    //| 
  float       * beta_seismic_resolution_;     //| Logs filtered to seismic resolution (log-domain)
  float       * rho_seismic_resolution_;      //| 

  int           firstM_;                      //| First well log entry contributing to blocked well
  int           lastM_;                       //| Last well log entry contributing to blocked well

  int           firstB_;                      //| First block with contribution from well log
  int           lastB_;                       //| Last block with contribution from well log

  int           nBlocks_;                     //| Number of blocks
  int           nLayers_;                     //| Number of vertical layers (nz)

  const int   * faciesNumbers_;
  int           nFacies_;
};

#endif
