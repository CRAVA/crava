#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include "rplib/trend.h"

#include "nrlib/random/distribution.hpp"

class DistributionWithTrend {
 public:
   DistributionWithTrend(const NRLib::Distribution<double>&  distr,
                         const Trend&                        mean,
                         const Trend&                        sd);
   virtual ~DistributionWithTrend();


   double ReSample(double s1, double s2) const {
     return sd_.GetValue(s1, s2)*distr_.Draw() + mean_.GetValue(s1, s2);
   }


 private:
  const NRLib::Distribution<double>&  distr_; // in standard form
  const Trend&                        mean_;
  const Trend&                        sd_;

};
#endif
