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
                        DEMTools::MixMethod                            mix_method,
                        std::vector<double>                          & alpha);

  DistributionsSolidMix(const DistributionsSolidMix & dist);

  virtual ~DistributionsSolidMix();

  virtual DistributionsSolid                    * Clone() const;

  virtual Solid                                 * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                                    HasDistribution() const;

  virtual std::vector<bool>                       HasTrend() const;

  virtual Solid *                                 UpdateSample(double                      corr_param,
                                                               bool                        param_is_time,
                                                               const std::vector<double> & trend,
                                                               const Solid               * sample) const;
protected:

private:

  Solid                                         * GetSample(const std::vector<double>  & u,
                                                            const std::vector<double>  & trend_params,
                                                            const std::vector<Solid *> & solid_samples) const;

  std::vector< DistributionsSolid * >             distr_solid_;     // Pointers to external objects.
  std::vector< DistributionWithTrend * >          distr_vol_frac_;  // Pointers to external objects.
  DEMTools::MixMethod                             mix_method_;
};

#endif
