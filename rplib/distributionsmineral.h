#ifndef DISTRIBUTIONSMINERAL_H
#define DISTRIBUTIONSMINERAL_H

#include "rplib/distributionssolid.h"
#include "rplib/mineral.h"

class DistributionWithTrend;
class DistributionsBrineEvolution;

class DistributionsMineral : public DistributionsSolid {
public:

  DistributionsMineral(DistributionWithTrend         * distr_k,
                       DistributionWithTrend         * distr_mu,
                       DistributionWithTrend         * distr_rho,
                       DistributionsMineralEvolution * distr_evolution = NULL);

  virtual ~DistributionsMineral();

  virtual Solid * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

private:
  DistributionWithTrend         * distr_k_;         // Pointer to external object.
  DistributionWithTrend         * distr_mu_;        // Pointer to external object.
  DistributionWithTrend         * distr_rho_;       // Pointer to external object.
  DistributionsMineralEvolution * distr_evolution_; // Pointer to external object.
};

#endif
