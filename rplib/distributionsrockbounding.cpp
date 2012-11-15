#include "rplib/distributionsrockbounding.h"
#include "rplib/rockbounding.h"

#include "rplib/rock.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"

#include <cassert>

DistributionsRockBounding::DistributionsRockBounding(const DistributionsRock     * upper_rock,
                                                     const DistributionsRock     * lower_rock,
                                                     const DistributionWithTrend * porosity,
                                                     const DistributionWithTrend * bulk_weight,
                                                     const DistributionWithTrend * p_wave_weight,
                                                     double                        correlation_weights,
                                                     const std::vector<double>   & alpha,
                                                     const std::vector<double>   & s_min,
                                                     const std::vector<double>   & s_max)
: upper_rock_(upper_rock),
  lower_rock_(lower_rock),
  porosity_(porosity),
  K_weight_(bulk_weight),
  M_weight_(p_wave_weight),
  correlation_weights_(correlation_weights)
{
  alpha_ = alpha;
  s_min_ = s_min;
  s_max_ = s_max;

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> variables(3);
  variables[0] = porosity_;
  variables[1] = K_weight_;
  variables[2] = M_weight_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      if(i == j)
        corr_matrix(i,j) = 1;
      else
        corr_matrix(i,j) = 0;
    }
  }

  corr_matrix(1,2) = correlation_weights_;
  corr_matrix(2,1) = correlation_weights_;

  tabulated_ = new Tabulated(variables, corr_matrix);

}

DistributionsRockBounding::DistributionsRockBounding(const DistributionsRockBounding & dist)
: DistributionsRock(dist),
  correlation_weights_(dist.correlation_weights_),
  tabulated_(dist.tabulated_)
{
  upper_rock_ = dist.upper_rock_->Clone();
  lower_rock_ = dist.lower_rock_->Clone();
  porosity_   = dist.porosity_  ->Clone();
  K_weight_   = dist.K_weight_  ->Clone();
  M_weight_   = dist.M_weight_  ->Clone();
  alpha_      = dist.alpha_;
  s_min_      = dist.s_min_;
  s_max_      = dist.s_max_;
}

DistributionsRockBounding::~DistributionsRockBounding()
{
  delete upper_rock_;
  delete lower_rock_;

  if(porosity_->GetIsShared() == false)
    delete porosity_;

  if(K_weight_->GetIsShared() == false)
    delete K_weight_;

  if(M_weight_->GetIsShared() == false)
    delete M_weight_;

  delete tabulated_;
}

DistributionsRock *
DistributionsRockBounding::Clone() const
{
  return new DistributionsRockBounding(*this);
}

Rock *
DistributionsRockBounding::GenerateSample(const std::vector<double> & trend_params) const
{
  Rock * sample_upper_rock = upper_rock_->GenerateSample(trend_params);
  Rock * sample_lower_rock = lower_rock_->GenerateSample(trend_params);

  std::vector<double> u(3);
  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Rock * sample_rock = GetSample(u, trend_params, sample_upper_rock, sample_lower_rock);

  delete sample_upper_rock;
  delete sample_lower_rock;

  return sample_rock;
}

Rock *
DistributionsRockBounding::GetSample(const std::vector<double> & u,
                                     const std::vector<double> & trend_params,
                                     const Rock                * sample_upper_rock,
                                     const Rock                * sample_lower_rock) const
{
  std::vector<double> sample;

  sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double sample_porosity = sample[0];
  double sample_K_weight = sample[1];
  double sample_M_weight = sample[2];

  Rock * new_rock = new RockBounding(sample_upper_rock, sample_lower_rock, sample_porosity, sample_K_weight, sample_M_weight, u);

  return new_rock;
}

std::vector<double>
DistributionsRockBounding::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  std::vector<double> dummy;
  return(dummy);
}

NRLib::Grid2D<double>
DistributionsRockBounding::GetCovariance(const std::vector<double> & /*trend_params*/) const
{
  NRLib::Grid2D<double> dummy;
  return(dummy);
}

bool
DistributionsRockBounding::HasDistribution() const
{
  bool has_distribution = false;

  if(upper_rock_->HasDistribution() == true || lower_rock_->HasDistribution() == true ||
    porosity_->GetIsDistribution() == true || K_weight_->GetIsDistribution() == true || M_weight_->GetIsDistribution() == true) {
      has_distribution = true;
  }

  return has_distribution;

}

std::vector<bool>
DistributionsRockBounding::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> upper_rock_trend = upper_rock_->HasTrend();
  std::vector<bool> lower_rock_trend = lower_rock_->HasTrend();
  std::vector<bool> poro_trend       = porosity_->GetUseTrendCube();
  std::vector<bool> K_trend          = K_weight_->GetUseTrendCube();
  std::vector<bool> M_trend          = M_weight_->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(upper_rock_trend[i] == true || lower_rock_trend[i] == true ||
      poro_trend[i] == true || K_trend[i] == true || M_trend[i] == true) {
        has_trend[i] = true;
    }
  }

  return has_trend;
}

bool
DistributionsRockBounding::GetIsOkForBounding() const
{
  bool is_ok_for_bounding = false;

  if(upper_rock_->GetIsOkForBounding() == true && lower_rock_->GetIsOkForBounding() == true)
    is_ok_for_bounding = true;

  return is_ok_for_bounding;
}

Rock *
DistributionsRockBounding::UpdateSample(double                      corr_param,
                                        bool                        param_is_time,
                                        const std::vector<double> & trend,
                                        const Rock                * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

  assert(typeid(sample) == typeid(RockBounding));
  const RockBounding * core_sample = dynamic_cast<const RockBounding *>(sample);

  Rock * updated_upper_rock = upper_rock_->UpdateSample(corr_param,
                                                        param_is_time,
                                                        trend,
                                                        core_sample->GetUpperRock());
  Rock * updated_lower_rock = lower_rock_->UpdateSample(corr_param,
                                                        param_is_time,
                                                        trend,
                                                        core_sample->GetLowerRock());

  Rock * updated_sample     = GetSample(u, trend, updated_upper_rock, updated_lower_rock);

  return updated_sample;
}
