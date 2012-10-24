#include "rplib/distributionssolidtabulatedvelocity.h"
#include "rplib/solidtabulatedmodulus.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"

DistributionsSolidTabulatedVelocity::DistributionsSolidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * vs,
                                                                         const DistributionWithTrend * density,
                                                                         double                        corr_vp_vs,
                                                                         double                        corr_vp_density,
                                                                         double                        corr_vs_density,
                                                                         std::vector<double>         & alpha)
: DistributionsSolid()
{
  vp_               = vp;
  vs_               = vs;
  density_          = density;
  corr_vp_vs_       = corr_vp_vs;
  corr_vp_density_  = corr_vp_density;
  corr_vs_density_  = corr_vs_density;
  alpha_            = alpha;

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = vp_;
  elastic_variables[1] = vs_;
  elastic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_vp_vs_;
  corr_matrix(1,0) = corr_vp_vs_;
  corr_matrix(0,2) = corr_vp_density_;
  corr_matrix(2,0) = corr_vp_density_;
  corr_matrix(1,2) = corr_vs_density_;
  corr_matrix(2,1) = corr_vs_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

}

DistributionsSolidTabulatedVelocity::~DistributionsSolidTabulatedVelocity()
{
  if(vp_->GetIsShared() == false)
    delete vp_;
  if(vs_->GetIsShared() == false)
    delete vs_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

Solid *
DistributionsSolidTabulatedVelocity::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(3);

  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Solid * solid = GetSample(u, trend_params);

  return solid;
}

Solid *
DistributionsSolidTabulatedVelocity::GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const
{

  std::vector<double> seismic_sample;

  seismic_sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double vp_sample      = seismic_sample[0];
  double vs_sample      = seismic_sample[1];
  double density_sample = seismic_sample[2];

  double bulk_sample;
  double shear_sample;

  DEMTools::CalcElasticParamsFromSeismicParams(vp_sample, vs_sample, density_sample, bulk_sample, shear_sample);

  Solid * solid = new SolidTabulatedModulus(bulk_sample, shear_sample, density_sample, u);

  return solid;
}

bool
DistributionsSolidTabulatedVelocity::HasDistribution() const
{
  bool has_distribution = false;

  if(vp_->GetIsDistribution() == true || vs_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}


std::vector<bool>
DistributionsSolidTabulatedVelocity::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> vp_trend      = vp_     ->GetUseTrendCube();
  std::vector<bool> vs_trend      = vs_     ->GetUseTrendCube();
  std::vector<bool> density_trend = density_->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(vp_trend[i] == true || vs_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

Solid *
DistributionsSolidTabulatedVelocity::UpdateSample(double                      corr_param,
                                                  bool                        param_is_time,
                                                  const std::vector<double> & trend,
                                                  const Solid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Solid * updated_sample = GetSample(u, trend);

  return updated_sample;
}
