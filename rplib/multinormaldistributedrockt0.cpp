#include "rplib/multinormaldistributedrockt0.h"

#include "rplib/multinormalrock.h"
#include "rplib/multinormalwithtrend.h"

MultiNormalDistributedRockT0::MultiNormalDistributedRockT0(const MultiNormalWithTrend& mult_normal_distr) :
  mult_normal_distr_(mult_normal_distr)
{
}


MultiNormalDistributedRockT0::~MultiNormalDistributedRockT0()
{
}

Rock* 
MultiNormalDistributedRockT0::GenerateSample(const std::vector<double> & trend_params) const {
  std::vector<double>  param(3, 0); 
  std::vector<double>  saturation;
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  mult_normal_distr_.ReSample(s1, s2, param[0], param[1], param[2]);
  return new MultiNormalRock(param, saturation);
}