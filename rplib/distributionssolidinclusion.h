#ifndef RPLIB_DISTRIBUTIONS_SOLID_INCLUSION_H
#define RPLIB_DISTRIBUTIONS_SOLID_INCLUSION_H

#include "rplib/distributionssolid.h"
#include "rplib/distributionssolid.h"
#include "rplib/solid.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/grid/grid2d.hpp"

class DistributionsSolidInclusion : public DistributionsSolid {
public:

  DistributionsSolidInclusion(std::vector<DistributionsSolid *>              distr_solid,
                              std::vector<double>                            distr_incl_spectrum,
                              std::vector<double>                            distr_incl_concentration);

  virtual ~DistributionsSolidInclusion();

  virtual Solid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

protected:
  virtual Solid * UpdateSample(const std::vector< double > & /*corr*/,
                               const Solid                 & /*solid*/)            const;


private:
  std::vector<DistributionsSolid *>              distr_solid_;              // Pointer to external object.
  std::vector<double>                            distr_incl_spectrum_;      // Pointers to external objects.
  std::vector<double>                            distr_incl_concentration_; // Pointers to external objects.

  std::vector<double>                            expectation_;
  NRLib::Grid2D<double>                          covariance_;

};

#endif
