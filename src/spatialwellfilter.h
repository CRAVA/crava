#ifndef SPATIALWELLFILTER_H
#define SPATIALWELLFILTER_H

#include "src/definitions.h"

class Corr;
class FFTGrid;
class Crava;
class WellData;
class BlockedLogs;

class SpatialWellFilter
{
public:
  SpatialWellFilter(int nwells);
  ~SpatialWellFilter();
  
  void                     setPriorSpatialCorr(FFTGrid *parSpatialCorr, WellData *well, int wellnr);
  void                     doFiltering(Corr                        * corr, 
                                       WellData                   ** wells, 
                                       int                           nWells, 
                                       bool                          useVpRhoFilter, 
                                       const std::vector<int>      & filterWell,
                                       int                           nAngles,
                                       const Crava                 * cravaResult, 
                                       const std::vector<Grid2D *> & noiseScale);
  std::vector<double **> & getSigmae(void) {return sigmae_;}
  
private:
  std::vector<double **> sigmae_;
  std::vector<double **> sigmaeVpRho_;
  int                    nData_;   ///< sum no blocks in all wells
  int                    nWells_;
  int                  * n_;
  double             *** priorSpatialCorr_;

  void doVpRhoFiltering(const double ** sigmapri, const double ** sigmapost, int n, 
                        BlockedLogs * blockedlogs);
  void updateSigmaE(double ** filter, double ** postCov, int n);
  void completeSigmaE(int lastn, const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale);
  
  void updateSigmaEVpRho(double ** filter, double ** postCov,  int nDim, int n);
  void completeSigmaEVpRho(int lastn, const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale);

  void computeSigmaEAdjusted(double** sigmae ,double** sigmaE0,double** sigmaETmp,int n,double** sigmaEAdj);

  void adjustDiagSigma(double ** sigmae, int n);

  void calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n, bool useVs);
  void MakeInterpolatedResiduals(const float * bwLog, const float * bwLogBG, const int n, const int offset, double ** residuals);
  void fillValuesInSigmapost(double **sigmapost, const int *ipos, const int *jpos, const int *kpos, FFTGrid *covgrid, int n, int ni, int nj);

};
#endif

