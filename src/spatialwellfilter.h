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
  void             doFiltering(Corr *corr, WellData **wells, int nWells, bool useVpRhoFilter);
  const int        getNdata(void) const {return nData_;}
  double        ** getSigmae(void) {return sigmae_;}
  
private:
  void             doVpRhoFiltering(const double ** sigmapri, const double ** sigmapost, int n, 
                                    BlockedLogs * blockedlogs);

  void             adjustDiagSigma(double ** sigmae, int n);
  void             calculateFilteredLogs(double **Aw, BlockedLogs *blockedlogs, int n, bool useVs);
  void             MakeInterpolatedResiduals(const float * bwLog, const float * bwLogBG, const int n, const int offset, double ** residuals);
  double       *** priorSpatialCorr_;

  double        ** sigmae_;
  double        ** sigmaeVpRho_;
  int              nData_;   ///< sum no blocks in all wells
  int              nWells_;
  int            * n_;
};
#endif

