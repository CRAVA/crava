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

  std::vector<NRLib::Matrix> & getSigmae(void) { return sigmae_ ;}

private:

  //----------------------------------------------------------------
  void doVpRhoFiltering(const NRLib::SymmetricMatrix & sigmapri,
                        const NRLib::SymmetricMatrix & sigmapost,
                        const int                      n,
                        BlockedLogs                  * blockedLogs);
  //----------------------------------------------------------------
  void updateSigmaE(const NRLib::Matrix & filter,
                    const NRLib::Matrix & Spost,
                    NRLib::Matrix       & sigmae,
                    int                   n);
  //-------------------------------------------------------


  void completeSigmaE(int                           lastn,
                      const Crava                 * cravaResult,
                      const std::vector<Grid2D *> & noiseScale);

  void updateSigmaEVpRho(const NRLib::Matrix & filter,
                         const NRLib::Matrix & postCov,
                         int                   nDim,
                         int                   n);

  void completeSigmaEVpRho(int                           lastn,
                           const Crava                 * cravaResult,
                           const std::vector<Grid2D *> & noiseScale);

  void computeSigmaEAdjusted(NRLib::Matrix &  sigmaeX,
                             double        ** sigmaE0,
                             double        ** sigmaETmp,
                             int              n,
                             double        ** sigmaEAdj);

  void adjustDiagSigma(NRLib::Matrix & sigmae,
                       int             n);



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

  std::vector<NRLib::Matrix> sigmae_;
  std::vector<NRLib::Matrix> sigmaeVpRho_;


  int                        nData_;   ///< sum no blocks in all wells
  int                        nWells_;
  int                    *   n_;
  double                 *** priorSpatialCorr_;

};
#endif

