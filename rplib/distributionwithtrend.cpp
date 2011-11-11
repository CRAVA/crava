#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend(NRLib::Distribution<double>*  distr,
                                             Trend*                        mean,
                                             Trend*                        sd) :
  distr_(distr), mean_(mean), sd_(sd) {


}

DistributionWithTrend::~DistributionWithTrend()
{
  
}
