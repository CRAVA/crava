#include "rplib/distributionssolidtabulated.h"
#include "rplib/solidtabulatedmodulus.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"

#include "nrlib/random/distribution.hpp"

DistributionsSolidTabulated::DistributionsSolidTabulated(const DistributionWithTrend * elastic1,
                                                         const DistributionWithTrend * elastic2,
                                                         const DistributionWithTrend * density,
                                                         double                        corr_elastic1_elastic2,
                                                         double                        corr_elastic1_density,
                                                         double                        corr_elastic2_density,
                                                         DEMTools::TabulatedMethod     method,
                                                         std::vector<double>         & alpha)
: DistributionsSolid()
{
  elastic1_               = elastic1;
  elastic2_               = elastic2;
  density_                = density;
  corr_elastic1_elastic2_ = corr_elastic1_elastic2;
  corr_elastic1_density_  = corr_elastic1_density;
  corr_elastic2_density_  = corr_elastic2_density;
  tabulated_method_       = method;
  alpha_                  = alpha;

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

}

DistributionsSolidTabulated::DistributionsSolidTabulated(const DistributionsSolidTabulated & dist)
: DistributionsSolid(dist),
  corr_elastic1_elastic2_(dist.corr_elastic1_elastic2_),
  corr_elastic1_density_(dist.corr_elastic1_density_),
  corr_elastic2_density_(dist.corr_elastic2_density_),
  tabulated_method_(dist.tabulated_method_)
{
  if(dist.elastic1_->GetIsShared() == false)
    elastic1_ = dist.elastic1_->Clone();
  else
    elastic1_ = dist.elastic1_;

  if(dist.elastic2_->GetIsShared() == false)
    elastic2_ = dist.elastic2_->Clone();
  else
    elastic2_ = dist.elastic2_;

  if(dist.density_->GetIsShared() == false)
    density_  = dist.density_ ->Clone();
  else
    density_ = dist.density_;

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

  alpha_ = dist.alpha_;
}

DistributionsSolidTabulated::~DistributionsSolidTabulated()
{
  if(elastic1_->GetIsShared() == false)
    delete elastic1_;
  if(elastic2_->GetIsShared() == false)
    delete elastic2_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

DistributionsSolid *
DistributionsSolidTabulated::Clone() const
{
  return new DistributionsSolidTabulated(*this);
}

Solid *
DistributionsSolidTabulated::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(3);

  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Solid * solid = GetSample(u, trend_params);

  return solid;
}

Solid *
DistributionsSolidTabulated::GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const
{

  std::vector<double> seismic_sample;

  seismic_sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double elastic1_sample = seismic_sample[0];
  double elastic2_sample = seismic_sample[1];
  double density_sample  = seismic_sample[2];

  double bulk_sample;
  double shear_sample;

  if(tabulated_method_ == DEMTools::Velocity)
    DEMTools::CalcElasticParamsFromSeismicParams(elastic1_sample, elastic2_sample, density_sample, bulk_sample, shear_sample);
  else {
    bulk_sample  = elastic1_sample;
    shear_sample = elastic2_sample;
  }

  Solid * solid = new SolidTabulatedModulus(bulk_sample, shear_sample, density_sample, u);

  return solid;
}

bool
DistributionsSolidTabulated::HasDistribution() const
{
  bool has_distribution = false;

  if(elastic1_->GetIsDistribution() == true || elastic2_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}


std::vector<bool>
DistributionsSolidTabulated::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> elastic1_trend = elastic1_->GetUseTrendCube();
  std::vector<bool> elastic2_trend = elastic2_->GetUseTrendCube();
  std::vector<bool> density_trend  = density_ ->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(elastic1_trend[i] == true || elastic2_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

Solid *
DistributionsSolidTabulated::UpdateSample(double                      corr_param,
                                          bool                        param_is_time,
                                          const std::vector<double> & trend,
                                          const Solid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Solid * updated_sample = GetSample(u, trend);

  return updated_sample;
}
