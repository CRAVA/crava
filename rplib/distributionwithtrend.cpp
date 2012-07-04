#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend(const NRLib::Distribution<double> * distr,
                                             const NRLib::Trend                * mean,
                                             const NRLib::Trend                * sd,
                                             bool                                sheared)
: distribution_(distr), mean_(mean), sd_(sd), is_sheared_(sheared)
{
}

DistributionWithTrend::~DistributionWithTrend()
{
  delete distribution_;
  delete mean_;
  delete sd_;
}
