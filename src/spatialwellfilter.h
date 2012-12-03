/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

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

  void                     setPriorSpatialCorr(FFTGrid  * parSpatialCorr,
                                               WellData * well,
                                               int        wellnr);

  //-----------------------------------------------------------------------------
  void                     doFiltering(Corr                        * corr,
                                       WellData                   ** wells,
                                       int                           nWells,
                                       bool                          useVpRhoFilter,
                                       int                           nAngles,
                                       const Crava                 * cravaResult,
                                       const std::vector<Grid2D *> & noiseScale);
  //-----------------------------------------------------------------------------

  std::vector<double **> & getSigmae(void) { return sigmae_ ;}

private:
  //-------------------------------------------------------
  void doVpRhoFiltering(const double ** sigmapri,
                        const double ** sigmapost,
                        int             n,
                        BlockedLogs  *  blockedLogs);
  //-------------------------------------------------------
  void updateSigmaE(const NRLib::Matrix &  filter,
                    const NRLib::Matrix &  PostCov,
                    double              ** sigmae,
                    int                    n);
  //-------------------------------------------------------


  void completeSigmaE(int                           lastn,
                      const Crava                 * cravaResult,
                      const std::vector<Grid2D *> & noiseScale);

  void updateSigmaEVpRho(double ** filter,
                         double ** postCov,
                         int       nDim,
                         int       n);

  void completeSigmaEVpRho(int                           lastn,
                           const Crava                 * cravaResult,
                           const std::vector<Grid2D *> & noiseScale);

  void computeSigmaEAdjusted(double ** sigmae ,
                             double ** sigmaE0,
                             double ** sigmaETmp,
                             int       n,
                             double ** sigmaEAdj);

  void adjustDiagSigma(double ** sigmae,
                       int       n);



  //----------------------------------------------------------------
  void calculateFilteredLogs(const NRLib::Matrix & Aw,
                             BlockedLogs         * blockedlogs,
                             int                   n,
                             bool                  useVs);
  //----------------------------------------------------------------
  void MakeInterpolatedResiduals(const float   * bwLog,
                                 const float   * bwLogBG,
                                 const int       n,
                                 const int       offset,
                                 NRLib::Vector & residuals);
  //----------------------------------------------------------------

  void fillValuesInSigmapost(double    ** sigmapost,
                             const int *  ipos,
                             const int *  jpos,
                             const int *  kpos,
                             FFTGrid   *  covgrid,
                             int          n,
                             int          ni,
                             int          nj);

  std::vector<double **>     sigmae_;
  std::vector<double **>     sigmaeVpRho_;
  int                        nData_;   ///< sum no blocks in all wells
  int                        nWells_;
  int                    *   n_;
  double                 *** priorSpatialCorr_;

};
#endif

