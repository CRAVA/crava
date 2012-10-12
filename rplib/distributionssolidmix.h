#ifndef RPLIB_DISTRIBUTIONSSOLIDMIX_H
#define RPLIB_DISTRIBUTIONSSOLIDMIX_H

#include "rplib/distributionssolid.h"

#include "nrlib/random/distribution.hpp"
#include "rplib/demmodelling.h"

class DistributionWithTrend;

class DistributionsSolidMix : public DistributionsSolid {
public:

  DistributionsSolidMix(std::vector< DistributionsSolid * >          & distr_solid,
                        std::vector< DistributionWithTrend * >       & distr_vol_frac,
                        DEMTools::MixMethod                            mix_method);

  virtual ~DistributionsSolidMix();

  virtual                                         Solid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

  virtual bool                                    HasDistribution() const;

  virtual std::vector<bool>                       HasTrend() const;

protected:
  virtual Solid *                                 UpdateSample(const std::vector< double > & /*corr*/,
                                                               const Solid                 & /*solid*/) const;

private:

  Solid                                         * GetSample(const std::vector<double>  & u,
                                                            const std::vector<double>  & trend_params,
                                                            const std::vector<Solid *> & solid_samples) const;

  std::vector< DistributionsSolid * >             distr_solid_;     // Pointers to external objects.
  std::vector< DistributionWithTrend * >          distr_vol_frac_;  // Pointers to external objects.
  DEMTools::MixMethod                             mix_method_;
};

#endif
