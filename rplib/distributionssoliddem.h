#ifndef RPLIB_DISTRIBUTIONSSOLIDDEM_H
#define RPLIB_DISTRIBUTIONSSOLIDDEM_H

#include "rplib/distributionssolid.h"

#include "nrlib/random/distribution.hpp"

class DistributionWithTrend;

class DistributionsSolidDEM : public DistributionsSolid {
public:

  DistributionsSolidDEM(DistributionsSolid                           * distr_solid,
                        DistributionsSolid                           * distr_solid_inc,
                        std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                        std::vector< DistributionWithTrend * >       & distr_incl_concentration);

  virtual                       ~DistributionsSolidDEM();

  virtual Solid               * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const; //dummy function that needs to be implemented

  virtual std::vector<bool>     HasTrend() const; //dummy function that needs to be implemented

protected:
  virtual Solid               * UpdateSample(const std::vector< double > & /*corr*/,
                                             const Solid                 & /*solid*/)            const;

private:

  DistributionsSolid                           * distr_solid_;              // Pointer to external object.
  DistributionsSolid                           * distr_solid_inc_;          // Pointer to external object.
  std::vector< DistributionWithTrend * >         distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< DistributionWithTrend * >         distr_incl_concentration_; // Pointers to external objects.
};

#endif
