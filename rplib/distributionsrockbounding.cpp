#include "rplib/distributionsrockbounding.h"
#include "rplib/rockbounding.h"

#include "rplib/rock.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"


DistributionsRockBounding::DistributionsRockBounding(const DistributionsRock     * upper_rock,
                                                     const DistributionsRock     * lower_rock,
                                                     const DistributionWithTrend * porosity,
                                                     const DistributionWithTrend * bulk_weight,
                                                     const DistributionWithTrend * p_wave_weight,
                                                     double                        correlation_weights)
: upper_rock_(upper_rock),
  lower_rock_(lower_rock),
  porosity_(porosity),
  K_weight_(bulk_weight),
  M_weight_(p_wave_weight),
  correlation_weights_(correlation_weights)
{
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

  // Find has_distribution_
  if(upper_rock_->HasDistribution() == true || lower_rock_->HasDistribution() == true ||
    porosity_->GetIsDistribution() == true || K_weight_->GetIsDistribution() == true || M_weight_->GetIsDistribution() == true) {
      has_distribution_ = true;
  }
  else
    has_distribution_ = false;

  // Find has_trend_
  std::vector<bool> upper_rock_trend = upper_rock_->HasTrend();
  std::vector<bool> lower_rock_trend = lower_rock_->HasTrend();
  std::vector<bool> poro_trend       = porosity_->GetUseTrendCube();
  std::vector<bool> K_trend          = K_weight_->GetUseTrendCube();
  std::vector<bool> M_trend          = M_weight_->GetUseTrendCube();

  has_trend_.resize(2);
  for(int i=0; i<2; i++) {
    if(upper_rock_trend[i] == true || lower_rock_trend[i] == true ||
      poro_trend[i] == true || K_trend[i] == true || M_trend[i] == true) {
        has_trend_[i] = true;
    }
    else
      has_trend_[i] = false;
  }


}

DistributionsRockBounding::~DistributionsRockBounding()
{
  if(porosity_->GetIsShared() == false)
    delete porosity_;

  if(K_weight_->GetIsShared() == false)
    delete K_weight_;

  if(M_weight_->GetIsShared() == false)
    delete M_weight_;

  delete tabulated_;
}

Rock *
DistributionsRockBounding::GenerateSample(const std::vector<double> & trend_params) const
{

  std::vector<double> u;
  std::vector<double> sample;

  sample = tabulated_->GenerateSample(u, trend_params[0], trend_params[1]);

  double sample_porosity = sample[0];
  double sample_K_weight = sample[1];
  double sample_M_weight = sample[2];

  Rock * sample_upper_rock = upper_rock_->GenerateSample(trend_params);
  Rock * sample_lower_rock = lower_rock_->GenerateSample(trend_params);

  Rock * new_rock = new RockBounding(sample_upper_rock, sample_lower_rock, sample_porosity, sample_K_weight, sample_M_weight, u);

  delete sample_upper_rock;
  delete sample_lower_rock;

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

Pdf3D *
DistributionsRockBounding::GeneratePdf() const
{
  Pdf3D * dummy = NULL;
  return(dummy);
}

bool
DistributionsRockBounding::HasDistribution() const
{
  return(has_distribution_);
}

std::vector<bool>
DistributionsRockBounding::HasTrend() const
{
  return(has_trend_);
}
