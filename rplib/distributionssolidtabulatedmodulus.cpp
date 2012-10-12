#include "rplib/distributionssolidtabulatedmodulus.h"

#include "rplib/distributionwithtrend.h"

DistributionsSolidTabulatedModulus::DistributionsSolidTabulatedModulus(const DistributionWithTrend         * distr_k,
                                                                       const DistributionWithTrend         * distr_mu,
                                                                       const DistributionWithTrend         * distr_rho,
                                                                       const double                          corr_k_mu,
                                                                       const double                          corr_k_rho,
                                                                       const double                          corr_mu_rho)
: DistributionsSolid(),
  distr_k_(distr_k),
  distr_mu_(distr_mu),
  distr_rho_(distr_rho),
  corr_k_mu_(corr_k_mu),
  corr_k_rho_(corr_k_rho),
  corr_mu_rho_(corr_mu_rho)
{
  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = distr_k_;
  elastic_variables[1] = distr_mu_;
  elastic_variables[2] = distr_rho_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_k_mu_;
  corr_matrix(1,0) = corr_k_mu_;
  corr_matrix(0,2) = corr_k_rho_;
  corr_matrix(2,0) = corr_k_rho_;
  corr_matrix(1,2) = corr_mu_rho_;
  corr_matrix(2,1) = corr_mu_rho_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

}

DistributionsSolidTabulatedModulus::~DistributionsSolidTabulatedModulus()
{
  if(distr_k_->GetIsShared() == false)
    delete distr_k_;
  if(distr_mu_->GetIsShared() == false)
    delete distr_mu_;
  if(distr_rho_->GetIsShared() == false)
    delete distr_rho_;

  delete tabulated_;
}

Solid *
DistributionsSolidTabulatedModulus::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(3);

  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Solid * solid = GetSample(u, trend_params);

  return solid;
}

Solid *
DistributionsSolidTabulatedModulus::GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const
{

  std::vector<double> sample;

  sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double k   = sample[0];
  double mu  = sample[1];
  double rho = sample[2];

  Solid * solid = new SolidTabulatedModulus(k, mu, rho, u);

  return solid;
}

bool
DistributionsSolidTabulatedModulus::HasDistribution() const
{
  bool has_distribution = false;

  if(distr_k_->GetIsDistribution() == true || distr_mu_->GetIsDistribution() == true || distr_rho_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}

std::vector<bool>
DistributionsSolidTabulatedModulus::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> k_trend       = distr_k_  ->GetUseTrendCube();
  std::vector<bool> mu_trend      = distr_mu_ ->GetUseTrendCube();
  std::vector<bool> density_trend = distr_rho_->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(k_trend[i] == true || mu_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

Solid *
DistributionsSolidTabulatedModulus::UpdateSample(const std::vector< double > &/*corr*/,
                                                 const Solid                 & /*solid*/) const {

  return NULL;
}
