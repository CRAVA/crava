#ifndef RPLIB_BETADISTRIBUTIONWITHTREND_H
#define RPLIB_BETADISTRIBUTIONWITHTREND_H

#include "rplib/distributionwithtrend.h"

namespace NRLib {
  class Trend;
  template <typename T>
  class Grid2D;
}

class BetaDistributionWithTrend : public DistributionWithTrend {
 public:

   BetaDistributionWithTrend();

   BetaDistributionWithTrend(const NRLib::Trend * mean,
                             const NRLib::Trend * var,
                             const double       & lower_limit,
                             const double       & upper_limit,
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
   void CalculateAlpha(const double & mean,
                       const double & var,
                       const double & lower_limit,
                       const double & upper_limit,
                       double       & alpha) const;

   void CalculateBeta(const double & mean,
                      const double & var,
                      const double & lower_limit,
                      const double & upper_limit,
                      double       & beta) const;

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
