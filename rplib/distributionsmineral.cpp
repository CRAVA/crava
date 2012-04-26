#include "rplib/distributionsmineral.h"

DistributionsMineral::DistributionsMineral(NRLib::Distribution<double>   * distr_k,
                                           NRLib::Distribution<double>   * distr_mu,
                                           NRLib::Distribution<double>   * distr_rho,
                                           DistributionsMineralEvolution * distr_evolution)
: DistributionsSolid()
{
  distr_k_          = distr_k;
  distr_mu_         = distr_mu;
  distr_rho_        = distr_rho;
  distr_evolution_  = distr_evolution;;
}

DistributionsMineral::~DistributionsMineral(){}

Solid *
DistributionsMineral::GenerateSample() const
{
  double k   = distr_k_->Draw();
  double mu  = distr_mu_->Draw();
  double rho = distr_rho_->Draw();
  Solid * solid = new Mineral(k, mu, rho, distr_evolution_);
  return solid;
}
