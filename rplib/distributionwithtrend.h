#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include "rplib/trend.h"

#include "nrlib/random/distribution.hpp"

class DistributionWithTrend {
 public:
   DistributionWithTrend(NRLib::Distribution<double>*  distr,
                         Trend*                        mean,
                         Trend*                        sd);
   virtual ~DistributionWithTrend();
   

   double ReSample(int s1, int s2) { 
     return sd_->GetValue(s1, s2)*distr_->Draw() + mean_->GetValue(s1, s2); 
   }
  

 private:
  NRLib::Distribution<double>*  distr_; // in standard form
  Trend*                        mean_;
  Trend*                        sd_;
   
};
#endif
