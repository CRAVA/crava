#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend(const NRLib::Distribution<double>&  distr,
                                             const NRLib::Trend&                 mean,
                                             const NRLib::Trend&                 sd) :
  distr_(distr), mean_(mean), sd_(sd) {


}

DistributionWithTrend::~DistributionWithTrend()
{

}
