#include "rplib/distributionsdryrocktabulated.h"
#include "rplib/dryrocktabulatedmodulus.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"
#include "rplib/distributionssolid.h"

#include "nrlib/random/distribution.hpp"

DistributionsDryRockTabulated::DistributionsDryRockTabulated(DistributionWithTrend     * elastic1,
                                                             DistributionWithTrend     * elastic2,
                                                             DistributionWithTrend     * density,
                                                             DistributionsSolid        * mineral,
                                                             DistributionWithTrend     * total_porosity,
                                                             double                      corr_elastic1_elastic2,
                                                             double                      corr_elastic1_density,
                                                             double                      corr_elastic2_density,
                                                             DEMTools::TabulatedMethod   method,
                                                             std::vector<double>       & alpha)
: DistributionsDryRock()
{
  if (elastic1->GetIsShared() == true)
    elastic1_ = elastic1;
  else
    elastic1_ = elastic1->Clone();

  if (elastic2->GetIsShared() == true)
    elastic2_ = elastic2;
  else
    elastic2_ = elastic2->Clone();

  if (density->GetIsShared() == true)
    density_ = density;
  else
    density_ = density->Clone();

  mineral_ = mineral->Clone();

  if (total_porosity->GetIsShared() == true)
    total_porosity_ = total_porosity;
  else
    total_porosity_ = total_porosity->Clone();

  corr_elastic1_elastic2_ = corr_elastic1_elastic2;
  corr_elastic1_density_  = corr_elastic1_density;
  corr_elastic2_density_  = corr_elastic2_density;
  tabulated_method_       = method;
  alpha_                  = alpha;

  // Generate tabulated_
  std::vector<DistributionWithTrend *> elastic_variables(3);
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

DistributionsDryRockTabulated::DistributionsDryRockTabulated(const DistributionsDryRockTabulated & dist)
: DistributionsDryRock(dist),
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

  mineral_          = dist.mineral_->Clone();

  if(dist.total_porosity_->GetIsShared() == false)
    total_porosity_ = dist.total_porosity_->Clone();
  else
    total_porosity_ = dist.total_porosity_;
  // Generate tabulated_
  std::vector<DistributionWithTrend *> elastic_variables(3);
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

DistributionsDryRockTabulated::~DistributionsDryRockTabulated()
{
  if(elastic1_->GetIsShared() == false)
    delete elastic1_;
  if(elastic2_->GetIsShared() == false)
    delete elastic2_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete mineral_;

  if(total_porosity_->GetIsShared() == false)
    delete total_porosity_;

  delete tabulated_;
}

DistributionsDryRock *
DistributionsDryRockTabulated::Clone() const
{
  return new DistributionsDryRockTabulated(*this);
}

DryRock *
DistributionsDryRockTabulated::GenerateSample(const std::vector<double> & trend_params)
{
  std::vector<double> u(3);

  for(unsigned int i=0; i<u.size(); i++)
    u[i] = NRLib::Random::Unif01();

  std::vector<double> u2(1);

  for(unsigned int i=0; i<u2.size(); i++)
    u2[i] = NRLib::Random::Unif01();

  DryRock * dryrock = GetSample(u, u2, trend_params);

  return dryrock;
}

DryRock *
DistributionsDryRockTabulated::GetSample(const std::vector<double> & u,
                                         const std::vector<double> & u2,
                                         const std::vector<double> & trend_params)
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
  Solid * mineral_sample = mineral_->GenerateSample(trend_params);
  double mineral_moduli_k_sample, dummy1, dummy2;
  mineral_sample->GetElasticParams(mineral_moduli_k_sample, dummy1, dummy2);
  delete mineral_sample;

  double total_porosity_sample    = total_porosity_->GetQuantileValue(u2[0], trend_params[0], trend_params[1]);

  std::vector<double> u_final(4);

  u_final[0] = u[0]; u_final[1] = u[1]; u_final[2] = u[2];
  u_final[3] = u2[0];

  DryRock * dryrock = new DryRockTabulatedModulus(bulk_sample, shear_sample, density_sample,
                                                  total_porosity_sample, mineral_moduli_k_sample, u_final);

  return dryrock;
}

bool
DistributionsDryRockTabulated::HasDistribution() const
{
  if (mineral_->HasDistribution() == true || elastic1_->GetIsDistribution() == true || elastic2_->GetIsDistribution() == true ||
     density_->GetIsDistribution() == true  || total_porosity_->GetIsDistribution() == true)
    return true;

  return false;
}


std::vector<bool>
DistributionsDryRockTabulated::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> elastic1_trend          = elastic1_->GetUseTrendCube();
  std::vector<bool> elastic2_trend          = elastic2_->GetUseTrendCube();
  std::vector<bool> density_trend           = density_ ->GetUseTrendCube();
  std::vector<bool> mineral_moduli_k_trend  = mineral_ ->HasTrend();
  std::vector<bool> total_porosity_trend    = total_porosity_ ->GetUseTrendCube();


  for(int i=0; i<2; i++) {
    if(elastic1_trend[i] == true || elastic2_trend[i] == true || density_trend[i] == true ||
       mineral_moduli_k_trend[i] == true || total_porosity_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

DryRock *
DistributionsDryRockTabulated::UpdateSample(double                      corr_param,
                                            bool                        param_is_time,
                                            const std::vector<double> & trend,
                                            const DryRock             * sample)
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);

  std::vector<double> u1(3);
  std::vector<double> u2(2);

  u1[0] = u[0]; u1[1] = u[1]; u1[2] = u[2];
  u2[0] = u[3];
  DryRock * updated_sample = GetSample(u1, u2, trend);

  return updated_sample;
}
