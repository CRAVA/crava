#ifndef RPLIB_BETADISTRIBUTIONWITHTREND_H
#define RPLIB_BETADISTRIBUTIONWITHTREND_H

#include "nrlib/trend/trend.hpp"
#include "nrlib/random/distribution.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"

class BetaDistributionWithTrend : public DistributionWithTrend {
 public:

   BetaDistributionWithTrend();

   BetaDistributionWithTrend(const NRLib::Trend * mean,
                             const NRLib::Trend * var,
                             bool                 shared);

   BetaDistributionWithTrend(const BetaDistributionWithTrend & dist);

   virtual ~BetaDistributionWithTrend();

   virtual DistributionWithTrend    * Clone() const                           { return new BetaDistributionWithTrend(*this) ;}

   virtual bool                       GetIsShared() const                     { return(is_shared_)                          ;}
   virtual bool                       GetIsDistribution() const               { return(true)                                ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                     ;}

   virtual double                     ReSample(double s1, double s2) const;
   virtual double                     GetQuantileValue(double u, double s1, double s2) const;

 private:
   void CalculateAlpha(double mean, double var, double & alpha) const;
   void CalculateBeta(double mean, double var, double & beta) const;

  NRLib::Grid2D<NRLib::Distribution<double> *>              * beta_distribution_;
  const NRLib::Trend                                        * mean_;
  const NRLib::Trend                                        * var_;
  const bool                                                  is_shared_;       // Use is_shared_ like in DistributionWithTrendStorage to know if we have a reservoir variable.
  std::vector<bool>                                           use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used
  int                                                         ni_;
  int                                                         nj_;
  std::vector<double>                                         mean_sampling_;
  std::vector<double>                                         var_sampling_;
  int                                                         n_samples_mean_;
  int                                                         n_samples_var_;

};
#endif
