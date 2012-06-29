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
                       DistributionsMineralEvolution * distr_evolution = NULL);

  virtual ~DistributionsMineral();

  virtual Solid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

private:
  NRLib::Distribution<double>   * distr_k_;         // Pointer to external object.
  NRLib::Distribution<double>   * distr_mu_;        // Pointer to external object.
  NRLib::Distribution<double>   * distr_rho_;       // Pointer to external object.
  DistributionsMineralEvolution * distr_evolution_; // Pointer to external object.
};

#endif
