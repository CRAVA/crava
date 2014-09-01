/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef SPATIALSYNTWELLFILTER_H
#define SPATIALSYNTWELLFILTER_H

#include "src/definitions.h"
#include "src/spatialwellfilter.h"
#include "rplib/syntwelldata.h"
#include "rplib/distributionsrock.h"
#include "src/avoinversion.h"
#include "src/timings.h"
#include "src/seismicparametersholder.h"
#include "lib/timekit.hpp"

class SpatialSyntWellFilter: public SpatialWellFilter
{
public:
  SpatialSyntWellFilter();

  SpatialSyntWellFilter(int         n_synt_wells,
                        bool        cov_estimated);

  ~SpatialSyntWellFilter();

  void                    SetPriorSpatialCovarianceSyntWell(const FFTGrid               * cov_vp,
                                                            const FFTGrid               * cov_vs,
                                                            const FFTGrid               * cov_rho,
                                                            const FFTGrid               * cov_vpvs,
                                                            const FFTGrid               * cov_vprho,
                                                            const FFTGrid               * cov_vsrho,
                                                            int                           wellnr);

  void                     SetPriorSpatialCorrSyntWell(FFTGrid                          * parSpatialCorr,
                                                       int                                wellnr);

  void                     DoFilteringSyntWells(SeismicParametersHolder                  & seismicParameters,
                                                const NRLib::Matrix                      & priorVar0);

  void                    GenerateSyntWellData (const std::map<std::string, DistributionsRock *>       & rock_distributions,
                                                const std::vector<std::string>                         & facies_names,
                                                const std::vector<double>                              & trend_min,
                                                const std::vector<double>                              & trend_max,
                                                double                                                   dz,
                                                int                                                      nWellsPerCombinationOfTrendParams,
                                                int                                                      n_bins_trend,
                                                int                                                      syntWellLength);

  void    SetNumberOfSyntWellsToBeFiltered(int n_synt_wells)                      { nWellsToBeFiltered_ = n_synt_wells                ;}
  void    SetNumberOfSyntWellsPerCombinationOfTrendParams(int n_synt_wells)       { nWellsPerCombinationOfTrendParams_ = n_synt_wells ;}

  const std::vector<SyntWellData *> & GetSyntWellData()                                                 const { return syntWellData_                     ;}
  const std::vector<double>         & GetTrend1()                                                       const { return trend_1_                          ;}
  const std::vector<double>         & GetTrend2()                                                       const { return trend_2_                          ;}
  int                                 GetNumberOfSyntWellsToBeFiltered()                                const { return nWellsToBeFiltered_               ;}
  int                                 GetNumberOfSyntWellsPerCombinationOfTrendParams()                 const { return nWellsPerCombinationOfTrendParams_;}



private:


  void    FillValuesInSigmapostSyntWell(double     ** sigmapost,
                                         const int  *  ipos,
                                         const int  *  jpos,
                                         const int  *  kpos,
                                         FFTGrid    *  covgrid,
                                         int           n,
                                         int           ni,
                                         int           nj);

  //In the case of synthetic wells
  std::vector<SyntWellData *>             syntWellData_;
  int                                     n_bins_trend_;
  std::vector<double>                     trend_1_;
  std::vector<double>                     trend_2_;
  double                                  trend_1_bin_size_;
  double                                  trend_2_bin_size_;

  int                                     nWellsToBeFiltered_;
  int                                     nWellsPerCombinationOfTrendParams_;


};
#endif