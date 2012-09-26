#ifndef RPLIB_DISTRIBUTIONSSOLIDMIX_H
#define RPLIB_DISTRIBUTIONSSOLIDMIX_H

#include "rplib/distributionssolid.h"

#include "nrlib/random/distribution.hpp"
#include "rplib/demmodelling.h"

class DistributionWithTrend;
class DistributionsSolidMixEvolution;
class DistributionsSolid;

class DistributionsSolidMix : public DistributionsSolid {
public:

                                                  DistributionsSolidMix(std::vector< DistributionsSolid * >          & distr_solid,
                                                                        std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                                                        DEMTools::MixMethod                            mix_method,
                                                                        DistributionsSolidMixEvolution               * distr_evolution = NULL);

  virtual                                         ~DistributionsSolidMix();

  virtual                                         Solid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

  virtual bool                                    HasDistribution() const;

  virtual std::vector<bool>                       HasTrend() const;

protected:
  virtual Solid *                                 UpdateSample(const std::vector< double > & /*corr*/,
                                                               const Solid                 & /*solid*/) const;

private:
  std::vector< DistributionsSolid * >             distr_solid_;     // Pointers to external objects.
  std::vector< DistributionWithTrend * >          distr_vol_frac_;  // Pointers to external objects.
  DEMTools::MixMethod                             mix_method_;
  DistributionsSolidMixEvolution               *  distr_evolution_; // Pointer to external object.
};

#endif
