#include "rplib/distributionsmineral.h"

#include "rplib/distributionwithtrend.h"

DistributionsMineral::DistributionsMineral(const DistributionWithTrend         * distr_k,
                                           const DistributionWithTrend         * distr_mu,
                                           const DistributionWithTrend         * distr_rho,
                                           const DistributionWithTrend         * corr_k_mu,
                                           const DistributionWithTrend         * corr_k_rho,
                                           const DistributionWithTrend         * corr_mu_rho,
                                                 DistributionsMineralEvolution * distr_evolution)
: DistributionsSolid(),
  distr_k_(distr_k),
  distr_mu_(distr_mu),
  distr_rho_(distr_rho),
  corr_k_mu_(corr_k_mu),
  corr_k_rho_(corr_k_rho),
  corr_mu_rho_(corr_mu_rho),
  distr_evolution_(distr_evolution)
{

}

DistributionsMineral::~DistributionsMineral(){}

Solid *
DistributionsMineral::GenerateSample(const std::vector<double> & trend_params) const
{
  double k   = distr_k_->ReSample(trend_params[0], trend_params[1]);
  double mu  = distr_mu_->ReSample(trend_params[0], trend_params[1]);
  double rho = distr_rho_->ReSample(trend_params[0], trend_params[1]);
  Solid * solid = new Mineral(k, mu, rho, distr_evolution_);
  return solid;
}

bool
DistributionsMineral::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsMineral::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Solid *
DistributionsMineral::UpdateSample(const std::vector< double > &/*corr*/,
                                   const Solid                 & /*solid*/) const {

  return NULL;
}
