/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef SPATIALREALWELLFILTER_H
#define SPATIALREALWELLFILTER_H

#include "src/definitions.h"
#include "src/spatialwellfilter.h"
#include "src/blockedlogscommon.h"
#include "src/timings.h"
#include "src/blockedlogscommon.h"
#include "src/seismicparametersholder.h"
#include "lib/timekit.hpp"

class AVOInversion;
class FFTGrid;


class SpatialRealWellFilter: public SpatialWellFilter
{
public:
  SpatialRealWellFilter();

  SpatialRealWellFilter(int       nwells,
                        bool      cov_estimated);

  ~SpatialRealWellFilter();

  // when prior covariance is estimated
  void                    SetPriorSpatialCovariance(const BlockedLogsCommon   * blocked_log,
                                                    int                         wellnr,
                                                    const FFTGrid             * cov_vp,
                                                    const FFTGrid             * cov_vs,
                                                    const FFTGrid             * cov_rho,
                                                    const FFTGrid             * cov_vpvs,
                                                    const FFTGrid             * cov_vprho,
                                                    const FFTGrid             * cov_vsrho);

  // when prior correlation is given as input
  void                     setPriorSpatialCorr(FFTGrid           * parSpatialCorr,
                                               BlockedLogsCommon * blocked_log,
                                               int                 wellnr);

  void                    doFiltering(std::map<std::string, BlockedLogsCommon *> blocked_logs,
                                      bool                                       useVpRhoFilter,
                                      int                                        nAngles,
                                      const AVOInversion                       * avoInversionResult,
                                      const std::vector<Grid2D *>              & noiseScale,
                                      SeismicParametersHolder                  & seismicParameters);


private:

  void adjustDiagSigma(NRLib::Matrix & sigmae);

  /*
  void calculateFilteredLogs(const NRLib::Matrix & Aw,
                             BlockedLogsCommon   * blockedlogs,
                             int                   n,
                             bool                  useVs);
  */

  void MakeInterpolatedResiduals(const std::vector<double> & bwLog,
                                 const std::vector<double> & bwLogBG,
                                 const int                   n,
                                 const int                   offset,
                                 NRLib::Vector &             residuals);

  
  void fillValuesInSigmapost(double    ** sigmapost,
                             const int *  ipos,
                             const int *  jpos,
                             const int *  kpos,
                             FFTGrid   *  covgrid,
                             int          n,
                             int          ni,
                             int          nj);


};
#endif
