#ifndef DISTRIBUTIONSMINERAL_H
#define DISTRIBUTIONSMINERAL_H

#include "rplib/distributionssolid.h"
#include "rplib/mineral.h"
#include "rplib/distributionsmineralevolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsMineral : public DistributionsSolid {
public:

  DistributionsMineral(NRLib::Distribution<double>   * distr_k,
                       NRLib::Distribution<double>   * distr_mu,
                       NRLib::Distribution<double>   * distr_rho,
                       DistributionsMineralEvolution * distr_evolution = NULL)
  : DistributionsSolid()
  {
    distr_k_          = distr_k;
    distr_mu_         = distr_mu;
    distr_rho_        = distr_rho;
    distr_evolution_  = distr_evolution;;
  }

  virtual ~DistributionsMineral(){}

  virtual Solid * GenerateSample() const {
    double k   = distr_k_->Draw();
    double mu  = distr_mu_->Draw();
    double rho = distr_rho_->Draw();
    Solid * solid = new Mineral(k, mu, rho, distr_evolution_);
    return solid;
  }

private:
  NRLib::Distribution<double>   * distr_k_;
  NRLib::Distribution<double>   * distr_mu_;
  NRLib::Distribution<double>   * distr_rho_;
  DistributionsMineralEvolution * distr_evolution_;
};

#endif
