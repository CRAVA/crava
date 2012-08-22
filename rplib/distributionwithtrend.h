#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include "nrlib/trend/trend.hpp"
#include "nrlib/random/distribution.hpp"

class DistributionWithTrend {
 public:
   DistributionWithTrend(const NRLib::Distribution<double> * distr,
                         const NRLib::Trend                * mean,
                         const NRLib::Trend                * sd,
                         bool                                sheared,
                         bool                                is_distribution);
   virtual ~DistributionWithTrend();


   double ReSample(double s1, double s2) const {
     double dummy = 0;
     return sd_->GetValue(s1, s2, dummy)*distribution_->Draw() + mean_->GetValue(s1, s2, dummy);
   }

   bool                       GetIsDistribution() const               { return(is_distribution_)                    ;}
   std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                     ;}

 private:

  void  FindUseTrendCube(int dim, int reference);

  const NRLib::Distribution<double> * distribution_; // in standard form
  const NRLib::Trend                * mean_;
  const NRLib::Trend                * sd_;
  const bool                          is_sheared_;       // Use is_shared_ like in DistributionWithTrendStorage to know if we have a reservoir variable.
  const bool                          is_distribution_;  // True if distribution_ is other distribution than dirac.
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used


};
#endif
