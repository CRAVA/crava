#ifndef SPATIALWELLFILTER_H
#define SPATIALWELLFILTER_H

class Corr;
class FFTGrid;

class SpatialWellFilter
{
public:
  SpatialWellFilter(int nwells);
  ~SpatialWellFilter();
  
  void             setPriorSpatialCorr(FFTGrid *parSpatialCorr, WellData *well, int wellnr);
  void             doFiltering(Corr *corr, WellData **wells, int nWells);
  const float    * getAlphaFiltered(void) const {return alphaFiltered_;}
  const float    * getBetaFiltered(void) const {return betaFiltered_;}
  const float    * getRhoFiltered(void) const {return rhoFiltered_;}
  const int        getNdata(void) const {return nData_;}
  double        ** getSigmae(void) {return sigmae_;}
  
private:
  void             adjustDiagSigma(double ** sigmae, int n);
  void             calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n);
  void             MakeInterpolatedResiduals(const float * bwLog, const float * bwLogBG, const int n, const int offset, double ** residuals);
  double       *** priorSpatialCorr_;
  float          * alphaFiltered_;
  float          * betaFiltered_; 
  float          * rhoFiltered_;

  //float       ** alpha_;
  //float       ** beta_;
  //float       ** rho_; 
  double        ** sigmae_;
  int              nData_;   ///< sum no blocks in all wells
  int              nWells_;
  int            * n_;
};
#endif

