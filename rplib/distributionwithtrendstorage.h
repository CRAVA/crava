#ifndef RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H
#define RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H

#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/trend/trendstorage.hpp"


class DistributionWithTrendStorage
{
public:
  DistributionWithTrendStorage();

  virtual ~DistributionWithTrendStorage();

  virtual DistributionWithTrend        * GenerateDistributionWithTrend(const std::string                       & /*path*/,
                                                                       const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                                       const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                                       std::string                             & /*errTxt*/)                    = 0;

  virtual NRLib::TrendStorage          * CloneMean() const = 0;

  virtual bool                           GetIsShared() const = 0;

};

//--------------------------------------------------------------//

class DeltaDistributionWithTrendStorage : public DistributionWithTrendStorage
{
public:

  DeltaDistributionWithTrendStorage();

  DeltaDistributionWithTrendStorage(double mean,
                                    bool   is_shared);

  DeltaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                    bool                        is_shared);

  virtual ~DeltaDistributionWithTrendStorage();

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual bool                             GetIsShared() const   { return(is_shared_) ;}

private:
  const NRLib::TrendStorage            * mean_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
};

//--------------------------------------------------------------//

class NormalDistributionWithTrendStorage : public DistributionWithTrendStorage
{
public:

  NormalDistributionWithTrendStorage();

  NormalDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                     const NRLib::TrendStorage * variance,
                                     bool                        is_shared);

  virtual ~NormalDistributionWithTrendStorage();

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual bool                             GetIsShared() const   { return(is_shared_) ;}

private:
  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
};

//--------------------------------------------------------------//

class BetaDistributionWithTrendStorage : public DistributionWithTrendStorage
{
public:

  BetaDistributionWithTrendStorage();

  BetaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                   const NRLib::TrendStorage * variance,
                                   bool                        is_shared);

  virtual ~BetaDistributionWithTrendStorage();

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual bool                             GetIsShared() const   { return(is_shared_) ;}

private:
  void                                     CheckBetaConsistency(NRLib::Trend * mean,
                                                                NRLib::Trend * variance,
                                                                std::string  & errTxt) const;

  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
};
#endif
