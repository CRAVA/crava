#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend(const NRLib::Distribution<double> * distr,
                                             const NRLib::Trend                * mean,
                                             const NRLib::Trend                * sd) 
: distribution_(distr), mean_(mean), sd_(sd) 
{
//Marit: Har tatt inn variance, ikke sd
}

DistributionWithTrend::~DistributionWithTrend()
{
  delete distribution_;
  delete mean_;
  delete sd_;
}
