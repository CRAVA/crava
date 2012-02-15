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
  
  NRLib::Normal vp0, vs0, density0;
  mult_normal_distr_.ReSample(vp0, vs0, density0, s1, s2, param[0], param[1], param[2]);
  return new MultiNormalRock(param, saturation);
}

std::vector<double>
MultiNormalDistributedRockT0::GetExpectation(const std::vector<double> & trend_params) const
{
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  std::vector<double> expectation;
  mult_normal_distr_.CalculateExpectation(s1, s2, expectation);
  return(expectation);
}

NRLib::Matrix
MultiNormalDistributedRockT0::GetCovariance(const std::vector<double> & trend_params) const
{
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  NRLib::Matrix covariance;
  mult_normal_distr_.CalculateCovariance(s1, s2, covariance);
  return(covariance);
}
