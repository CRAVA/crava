#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include "nrlib/trend/trend.hpp"
#include "nrlib/random/distribution.hpp"

class DistributionWithTrend {
 public:
   DistributionWithTrend(const NRLib::Distribution<double> * distr,
                         const NRLib::Trend                * mean,
                         const NRLib::Trend                * sd);
   virtual ~DistributionWithTrend();


   double ReSample(double s1, double s2) const {
     double dummy = 0;
     return sd_->GetValue(s1, s2, dummy)*distribution_->Draw() + mean_->GetValue(s1, s2, dummy);
   }


 private:
  const NRLib::Distribution<double> * distribution_; // in standard form
  const NRLib::Trend                * mean_;
  const NRLib::Trend                * sd_;

};
#endif
