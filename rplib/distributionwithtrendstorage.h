#ifndef RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H
#define RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H

#include <string>
#include <vector>

#include "rplib/distributionwithtrend.h"


class DistributionWithTrend;

namespace NRLib {
  class TrendStorage;
  class Trend;
}


class DistributionWithTrendStorage
{
public:
  DistributionWithTrendStorage();

  virtual ~DistributionWithTrendStorage();

  virtual DistributionWithTrendStorage * Clone() const = 0;

  virtual DistributionWithTrend        * GenerateDistributionWithTrend(const std::string                       & /*path*/,
                                                                       const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                                       const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                                       std::string                             & /*errTxt*/)                    = 0;

  virtual NRLib::TrendStorage          * CloneMean() const = 0;

  virtual bool                           GetIsShared() const = 0;

  virtual void                           SetVintageYear(int year) = 0;

  virtual void                           SetOneYearCorrelation(double correlation) = 0;

  virtual int                            GetVintageYear() = 0;

  virtual double                         GetOneYearCorrelation() = 0;

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

  DeltaDistributionWithTrendStorage(const DeltaDistributionWithTrendStorage & dist);

  virtual ~DeltaDistributionWithTrendStorage();

  virtual DistributionWithTrendStorage   * Clone() const                             { return new DeltaDistributionWithTrendStorage(*this) ;}

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual const NRLib::TrendStorage      * GetMean()   const                         { return mean_                        ;}

  virtual void                             SetVintageYear(int year)                  { vintage_year_ = year                ;}

  virtual void                             SetOneYearCorrelation(double correlation) { one_year_correlation_ = correlation ;}

  virtual bool                             GetIsShared() const                       { return is_shared_                   ;}

  virtual int                              GetVintageYear()                          { return vintage_year_                ;}

  virtual double                           GetOneYearCorrelation()                   { return one_year_correlation_        ;}

private:
  const NRLib::TrendStorage            * mean_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
  int                                    vintage_year_;
  double                                 one_year_correlation_;

};

//--------------------------------------------------------------//

class NormalDistributionWithTrendStorage : public DistributionWithTrendStorage
{
public:

  NormalDistributionWithTrendStorage();

  NormalDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                     const NRLib::TrendStorage * variance,
                                     bool                        is_shared);

  NormalDistributionWithTrendStorage(const NormalDistributionWithTrendStorage & dist);

  virtual ~NormalDistributionWithTrendStorage();

  virtual DistributionWithTrendStorage   * Clone() const                             { return new NormalDistributionWithTrendStorage(*this) ;}

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual void                             SetVintageYear(int year)                  { vintage_year_ = year                ;}

  virtual void                             SetOneYearCorrelation(double correlation) { one_year_correlation_ = correlation ;}

  virtual bool                             GetIsShared() const                       { return is_shared_                   ;}

  virtual int                              GetVintageYear()                          { return vintage_year_                ;}

  virtual double                           GetOneYearCorrelation()                   { return one_year_correlation_        ;}

private:
  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
  int                                    vintage_year_;
  double                                 one_year_correlation_;
};

//--------------------------------------------------------------//

class BetaDistributionWithTrendStorage : public DistributionWithTrendStorage
{
public:

  BetaDistributionWithTrendStorage();

  BetaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                   const NRLib::TrendStorage * variance,
                                   const double              & lower_limit,
                                   const double              & upper_limit,
                                   bool                        is_shared);

  BetaDistributionWithTrendStorage(const BetaDistributionWithTrendStorage & dist);

  virtual ~BetaDistributionWithTrendStorage();

  virtual DistributionWithTrendStorage   * Clone() const                             { return new BetaDistributionWithTrendStorage(*this) ;}

  virtual DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                         std::string                             & errTxt);

  virtual NRLib::TrendStorage            * CloneMean() const;

  virtual void                             SetVintageYear(int year)                  { vintage_year_ = year                ;}

  virtual void                             SetOneYearCorrelation(double correlation) { one_year_correlation_ = correlation ;}

  virtual bool                             GetIsShared() const                       { return is_shared_                   ;}

  virtual int                              GetVintageYear()                          { return vintage_year_                ;}

  virtual double                           GetOneYearCorrelation()                   { return one_year_correlation_        ;}

  virtual const double                     GetLowerLimit() const                     { return lower_limit_                 ;}
  virtual const double                     GetUpperLimit() const                     { return upper_limit_                 ;}

private:
  void                                     CheckBetaConsistency(NRLib::Trend * mean,
                                                                NRLib::Trend * variance,
                                                                std::string  & errTxt) const;

  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
  double                                 lower_limit_;
  double                                 upper_limit_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_shared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
  int                                    vintage_year_;
  double                                 one_year_correlation_;
};
#endif
