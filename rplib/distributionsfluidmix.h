#ifndef RPLIB_DISTRIBUTIONSFLUIDMIX_H
#define RPLIB_DISTRIBUTIONSFLUIDMIX_H

#include "rplib/distributionsfluid.h"

#include "rplib/demmodelling.h"

#include "nrlib/random/distribution.hpp"

class DistributionWithTrend;
class FluidMix;

class DistributionsFluidMix : public DistributionsFluid {
public:

  DistributionsFluidMix(std::vector< DistributionsFluid * >          & distr_fluid,
                        std::vector< DistributionWithTrend * >       & distr_vol_frac,
                        DEMTools::MixMethod                            mix_method);

  virtual ~DistributionsFluidMix();

  virtual Fluid                                 * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                                    HasDistribution() const;

  virtual std::vector<bool>                       HasTrend() const;

protected:
  virtual Fluid                                 * UpdateSample(const std::vector< double > & /*corr*/,
                                                               const Fluid                 & /*fluid*/) const;

private:

  Fluid                                         * GetSample(const std::vector<double>  & u,
                                                            const std::vector<double>  & trend_params,
                                                            const std::vector<Fluid *> & fluid_samples) const;

  std::vector< DistributionsFluid * >             distr_fluid_;     // Pointers to external objects.
  std::vector< DistributionWithTrend * >          distr_vol_frac_;  // Pointers to external objects.
  DEMTools::MixMethod                             mix_method_;
};

#endif
