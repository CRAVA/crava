#include "rplib/distributionsmineral.h"

#include "rplib/distributionwithtrend.h"

DistributionsMineral::DistributionsMineral(const DistributionWithTrend         * distr_k,
                                           const DistributionWithTrend         * distr_mu,
                                           const DistributionWithTrend         * distr_rho,
                                           const double                          corr_k_mu,
                                           const double                          corr_k_rho,
                                           const double                          corr_mu_rho,
                                           DistributionsMineralEvolution       * distr_evolution)
: DistributionsSolid(),
  distr_k_(distr_k),
  distr_mu_(distr_mu),
  distr_rho_(distr_rho),
  corr_k_mu_(corr_k_mu),
  corr_k_rho_(corr_k_rho),
  corr_mu_rho_(corr_mu_rho),
  distr_evolution_(distr_evolution)
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

  tabulated_ = Tabulated(elastic_variables, corr_matrix);

  // Find has_distribution_
  if(distr_k_->GetIsDistribution() == true || distr_mu_->GetIsDistribution() == true || distr_rho_->GetIsDistribution() == true) {
    has_distribution_ = true;
  }
  else
    has_distribution_ = false;

  // Find has_trend_
  std::vector<bool> k_trend       = distr_k_  ->GetUseTrendCube();
  std::vector<bool> mu_trend      = distr_mu_ ->GetUseTrendCube();
  std::vector<bool> density_trend = distr_rho_->GetUseTrendCube();

  has_trend_.resize(2);
  for(int i=0; i<2; i++) {
    if(k_trend[i] == true || mu_trend[i] == true || density_trend[i] == true)
      has_trend_[i] = true;
    else
      has_trend_[i] = false;
  }

}

DistributionsMineral::~DistributionsMineral()
{
  if(distr_k_->GetIsShared() == false)
    delete distr_k_;
  if(distr_mu_->GetIsShared() == false)
    delete distr_mu_;
  if(distr_rho_->GetIsShared() == false)
    delete distr_rho_;
}

Solid *
DistributionsMineral::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u;

  std::vector<double> sample = tabulated_.GenerateSample(u, trend_params[0], trend_params[1]);

  double k   = sample[0];
  double mu  = sample[1];
  double rho = sample[2];

  Solid * solid = new Mineral(k, mu, rho, u, distr_evolution_);

  return solid;
}

bool
DistributionsMineral::HasDistribution() const
{
  return(has_distribution_);
}

std::vector<bool>
DistributionsMineral::HasTrend() const
{
  return(has_trend_);
}

Solid *
DistributionsMineral::UpdateSample(const std::vector< double > &/*corr*/,
                                   const Solid                 & /*solid*/) const {

  return NULL;
}
