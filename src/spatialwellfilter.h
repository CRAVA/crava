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
  void                     doFiltering(Corr *corr, WellData **wells, int nWells, bool useVpRhoFilter, int nAngles,
                                       const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale);
  int                      getNdata(void) const {return nData_;}
  std::vector<double **> & getSigmae(void) {return sigmae_;}
  
private:
  std::vector<double **> sigmae_;
  std::vector<double **> sigmaeVpRho_;
  int                    nData_;   ///< sum no blocks in all wells
  int                    nWells_;
  int                  * n_;
  double             *** priorSpatialCorr_;

  void doVpRhoFiltering(const double ** sigmapri, const double ** sigmapost, int n, 
                        BlockedLogs * blockedlogs, const Crava * cravaResult,
                        const std::vector<Grid2D *> & noiseScale);
  void updateSigmaE(double ** filter, double ** priCov, double ** postCov, int n, 
                    const Crava * cravaResult, const std::vector<Grid2D *> & noiseScale);
  void completeSigmaE(int lastn);
  void updateSigmaEVpRho(double ** filter, 
                         double ** priCov, 
                         double ** postCov, 
                         int nDim, 
                         int n,
                         const Crava * cravaResult, 
                         const std::vector<Grid2D *> & noiseScale);
  void completeSigmaEVpRho(int lastn);

  void adjustDiagSigma(double ** sigmae, int n);
  void calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n, bool useVs);
  void MakeInterpolatedResiduals(const float * bwLog, const float * bwLogBG, const int n, const int offset, double ** residuals);


};
#endif

