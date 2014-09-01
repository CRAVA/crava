/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SPATIALWELLFILTER_H
#define SPATIALWELLFILTER_H

#include "src/definitions.h"
#include "rplib/syntwelldata.h"
#include "libs/nrlib/flens/nrlib_flens.hpp"

class FFTGrid;
class AVOInversion;
class BlockedLogsCommon;
class SeismicParametersHolder;

class SpatialWellFilter
{
public:
  SpatialWellFilter();

  SpatialWellFilter(int       nwells,
                    bool      cov_estimated);

  ~SpatialWellFilter();


  std::vector<NRLib::Matrix> &  getSigmae()           { return sigmae_                           ;}


protected:

  void doVpRhoFiltering(std::vector<NRLib::Matrix>      & sigmaeVpRho,
                        double                         ** sigmapri,
                        double                         ** sigmapost,
                        const int                         n,
                        BlockedLogsCommon               * blockedLogs);

  void completeSigmaE(std::vector<NRLib::Matrix>        & sigmae,
                      int                                 lastn,
                      const AVOInversion                * avoInversionResult,
                      const std::vector<Grid2D *>       & noiseScale);

  void updateSigmaE(NRLib::Matrix                       & sigmae,
                    const NRLib::Matrix                 & Filter,
                    const NRLib::Matrix                 & PostCov,
                    int                                   n);

  void updateSigmaEVpRho(std::vector<NRLib::Matrix>     & sigmaeVpRho,
                         const NRLib::Matrix            & Aw,
                         const NRLib::Matrix            & Spost,
                         int                              nDim,
                         int                              n);

  void completeSigmaEVpRho(std::vector<NRLib::Matrix>   & sigmaeVpRho,
                           int                            lastn,
                           const AVOInversion           * avoInversionResult,
                           const std::vector<Grid2D *>  & noiseScale);

  void computeSigmaEAdjusted(NRLib::Matrix              & sigmae,
                             NRLib::Matrix              & sigmaE0,
                             NRLib::Matrix              & sigmaETmp,
                             NRLib::Matrix              & sigmaEAdj);

  void adjustDiagSigma(NRLib::Matrix                    & sigmae);

  void calculateFilteredLogs(const NRLib::Matrix        & Aw,
                             BlockedLogsCommon          * blockedlogs,
                             int                          n,
                             bool                         useVs);

  void MakeInterpolatedResiduals(const std::vector<double>  & bwLog,
                                 const std::vector<double>  & bwLogBG,
                                 const int                    n,
                                 const int                    offset,
                                 NRLib::Vector              & residuals);

  // ****************************************************************

  std::vector<NRLib::Matrix>                            sigmae_;

  int                                                   nData_;             ///< sum of blocks in all wells
  int                                                   nWells_;
  int                                                 * n_;
  double                                            *** priorSpatialCorr_;  // vector (one entry per blocked log) of spatial 2D corr matrices
  std::vector<NRLib::Grid2D<double> >                    prior_cov_vp_;      // vector (one entry per blocked log) of 2D covariance matrices
  std::vector<NRLib::Grid2D<double> >                    prior_cov_vs_;
  std::vector<NRLib::Grid2D<double> >                    prior_cov_rho_;
  std::vector<NRLib::Grid2D<double> >                    prior_cov_vpvs_;
  std::vector<NRLib::Grid2D<double> >                    prior_cov_vprho_;
  std::vector<NRLib::Grid2D<double> >                    prior_cov_vsrho_;

  std::vector<double **>                                sigmaeVpRho_;

};
#endif
