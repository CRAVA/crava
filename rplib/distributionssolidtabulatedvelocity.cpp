#include "rplib/distributionssolidtabulatedvelocity.h"
#include "rplib/mineral.h"
#include "rplib/tabulated.h"

DistributionsSolidTabulatedVelocity::DistributionsSolidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * vs,
                                                                         const DistributionWithTrend * density,
                                                                         double                        corr_vp_vs,
                                                                         double                        corr_vp_density,
                                                                         double                        corr_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  corr_vp_vs_(corr_vp_vs),
  corr_vp_density_(corr_vp_density),
  corr_vs_density_(corr_vs_density)
{
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

   // Find has_distribution_
  if(vp_->GetIsDistribution() == true || vs_->GetIsDistribution() == true || density_->GetIsDistribution() == true) {
    has_distribution_ = true;
  }
  else
    has_distribution_ = false;

  // Find has_trend_
  std::vector<bool> vp_trend      = vp_     ->GetUseTrendCube();
  std::vector<bool> vs_trend      = vs_     ->GetUseTrendCube();
  std::vector<bool> density_trend = density_->GetUseTrendCube();

  has_trend_.resize(2);
  for(int i=0; i<2; i++) {
    if(vp_trend[i] == true || vs_trend[i] == true || density_trend[i] == true)
      has_trend_[i] = true;
    else
      has_trend_[i] = false;
  }

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
  std::vector<double> u;
  std::vector<double> seismic_sample;

  seismic_sample = tabulated_->GenerateSample(u, trend_params[0], trend_params[1]);

  double vp_sample      = seismic_sample[0];
  double vs_sample      = seismic_sample[1];
  double density_sample = seismic_sample[2];

  double bulk_sample  = density_sample * (std::pow(vp_sample,2) - 4/3 * std::pow(vs_sample,2));
  double shear_sample = density_sample *  std::pow(vs_sample,2);

  Solid * solid = new Mineral(bulk_sample, shear_sample, density_sample, u, NULL);

  return solid;
}

bool
DistributionsSolidTabulatedVelocity::HasDistribution() const
{
  return(has_distribution_);
}

std::vector<bool>
DistributionsSolidTabulatedVelocity::HasTrend() const
{
  return(has_trend_);
}

Solid *
DistributionsSolidTabulatedVelocity::UpdateSample(const std::vector< double > &/*corr*/,
                                                  const Solid                 & /*solid*/) const {

  return NULL;
}
