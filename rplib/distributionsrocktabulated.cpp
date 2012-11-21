#include "rplib/distributionsrocktabulated.h"
#include "rplib/rocktabulatedvelocity.h"
#include "rplib/demmodelling.h"

#include "nrlib/grid/grid2d.hpp"

DistributionsRockTabulated::DistributionsRockTabulated(const DistributionWithTrend * elastic1,
                                                       const DistributionWithTrend * elastic2,
                                                       const DistributionWithTrend * density,
                                                       double                        corr_elastic1_elastic2,
                                                       double                        corr_elastic1_density,
                                                       double                        corr_elastic2_density,
                                                       DEMTools::TabulatedMethod     method,
                                                       const std::vector<double>   & alpha,
                                                       const std::vector<double>   & s_min,
                                                       const std::vector<double>   & s_max)
  : DistributionsRock(),
    elastic1_(elastic1),
    elastic2_(elastic2),
    density_(density),
    corr_elastic1_elastic2_(corr_elastic1_elastic2),
    corr_elastic1_density_(corr_elastic1_density),
    corr_elastic2_density_(corr_elastic2_density),
    tabulated_method_(method)
{
  alpha_ = alpha;
  s_min_ = s_min;
  s_max_ = s_max;

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = elastic1_;
  elastic_variables[1] = elastic2_;
  elastic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_elastic1_elastic2_;
  corr_matrix(1,0) = corr_elastic1_elastic2_;
  corr_matrix(0,2) = corr_elastic1_density_;
  corr_matrix(2,0) = corr_elastic1_density_;
  corr_matrix(1,2) = corr_elastic2_density_;
  corr_matrix(2,1) = corr_elastic2_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

  SetupExpectationAndCovariances(expectation_,
                                 covariance_,
                                 tabulated_s0_,
                                 tabulated_s1_,
                                 s_min_,
                                 s_max_);
}

DistributionsRockTabulated::DistributionsRockTabulated(const DistributionsRockTabulated & dist)
  : DistributionsRock(dist),
    corr_elastic1_elastic2_(dist.corr_elastic1_elastic2_),
    corr_elastic1_density_(dist.corr_elastic1_density_),
    corr_elastic2_density_(dist.corr_elastic2_density_),
    tabulated_method_(dist.tabulated_method_)
{
  elastic1_ = dist.elastic1_->Clone();
  elastic2_ = dist.elastic1_->Clone();
  density_  = dist.density_ ->Clone();

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = elastic1_;
  elastic_variables[1] = elastic2_;
  elastic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_elastic1_elastic2_;
  corr_matrix(1,0) = corr_elastic1_elastic2_;
  corr_matrix(0,2) = corr_elastic1_density_;
  corr_matrix(2,0) = corr_elastic1_density_;
  corr_matrix(1,2) = corr_elastic2_density_;
  corr_matrix(2,1) = corr_elastic2_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

  alpha_       = dist.alpha_;
  s_min_       = dist.s_min_;
  s_max_       = dist.s_max_;
  expectation_ = dist.expectation_;
  covariance_  = dist.covariance_;
}

DistributionsRockTabulated::~DistributionsRockTabulated()
{
  if(elastic1_->GetIsShared() == false)
    delete elastic1_;

  if(elastic2_->GetIsShared() == false)
    delete elastic2_;

  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

DistributionsRock *
DistributionsRockTabulated::Clone() const
{
  return new DistributionsRockTabulated(*this);
}

Rock *
DistributionsRockTabulated::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(3);

  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Rock * new_rock = GetSample(u, trend_params);

  return new_rock;
}

Rock *
DistributionsRockTabulated::GetSample(const std::vector<double> & u,
                                      const std::vector<double> & trend_params) const
{
  std::vector<double> sample;

  sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double sample_elastic1 = sample[0];
  double sample_elastic2 = sample[1];
  double sample_density  = sample[2];

  double sample_vp;
  double sample_vs;

  if(tabulated_method_ == DEMTools::Modulus)
    DEMTools::CalcSeismicParamsFromElasticParams(sample_elastic1, sample_elastic2, sample_density, sample_vp, sample_vs);
  else {
    sample_vp = sample_elastic1;
    sample_vs = sample_elastic2;
  }

  Rock * new_rock = new RockTabulatedVelocity(sample_vp, sample_vs, sample_density, u);

  return new_rock;
}

bool
DistributionsRockTabulated::HasDistribution() const
{
  bool has_distribution = false;

  if(elastic1_->GetIsDistribution() == true || elastic2_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;

}

std::vector<bool>
DistributionsRockTabulated::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> elastic1_trend = elastic1_->GetUseTrendCube();
  std::vector<bool> elastic2_trend = elastic2_->GetUseTrendCube();
  std::vector<bool> density_trend  = density_ ->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(elastic1_trend[i] == true || elastic2_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return(has_trend);
}

bool
DistributionsRockTabulated::GetIsOkForBounding() const
{
  bool is_ok_for_bounding = false;

  bool              has_distribution = HasDistribution();
  std::vector<bool> has_trend        = HasTrend();

  if(has_distribution == false && has_trend[0] == false && has_trend[1] == false)
    is_ok_for_bounding = true;

  return(is_ok_for_bounding);
}

Rock *
DistributionsRockTabulated::UpdateSample(double                      corr_param,
                                         bool                        param_is_time,
                                         const std::vector<double> & trend,
                                         const Rock                * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);
  Rock * updated_sample = GetSample(u, trend);

  return updated_sample;
}
