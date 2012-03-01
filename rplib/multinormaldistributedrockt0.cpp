#include "rplib/multinormaldistributedrockt0.h"
#include "rplib/multinormalrock.h"
#include "rplib/trinormalwith2dtrend.h"
#include "rplib/pdf3dgaussian.h"
#include "rplib/pdf3d.h"

#include "nrlib/grid/grid2d.hpp"

MultiNormalDistributedRockT0::MultiNormalDistributedRockT0(const TriNormalWith2DTrend& mult_normal_distr) :
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
  std::vector<double> expectation(3);
  mult_normal_distr_.GetExpectation(s1, s2, expectation);
  return(expectation);
}

NRLib::Grid2D<double>
MultiNormalDistributedRockT0::GetCovariance(const std::vector<double> & trend_params) const
{
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  NRLib::Grid2D<double> covariance(3,3,0);
  mult_normal_distr_.GetCovariance(s1, s2, covariance);
  return(covariance);
}

Pdf3D* 
MultiNormalDistributedRockT0::GeneratePdf(void) const
{

  return new Pdf3DGaussian(mult_normal_distr_);

}
