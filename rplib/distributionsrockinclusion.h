#ifndef RPLIB_DISTRIBUTIONSROCKINCLUSION_H
#define RPLIB_DISTRIBUTIONSROCKINCLUSION_H

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/rockinclusion.h"
#include "rplib/distributionsrockinclusionevolution.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/grid/grid2d.hpp"

class DistributionsRockInclusion : public DistributionsRock {
public:

  DistributionsRockInclusion(DistributionsSolid                           * distr_solid,
                             DistributionsFluid                           * distr_fluid,
                             std::vector< NRLib::Distribution<double> * >   distr_incl_spectrum,
                             std::vector< NRLib::Distribution<double> * >   distr_incl_concentration,
                             NRLib::Distribution<double>                  * distr_porosity,
                             DistributionsRockInclusionEvolution          * distr_evolution = NULL);

  virtual ~DistributionsRockInclusion();

  virtual Rock * GenerateSample(const std::vector<double> & /*trend_params*/) const;

  virtual std::vector<double>   GetExpectation(const std::vector<double> & /*trend_params*/) const;

  virtual NRLib::Grid2D<double> GetCovariance(const std::vector<double> & /*trend_params*/)  const;

  virtual Pdf3D * GeneratePdf(void) const; // Returns NULL.

private:
  void SampleVpVsRhoExpectationAndCovariance(std::vector<double>   & expectation,
                                             NRLib::Grid2D<double> & covariance);

  DistributionsSolid                           * distr_solid_;              // Pointer to external object.
  DistributionsFluid                           * distr_fluid_;              // Pointer to external object.
  std::vector< NRLib::Distribution<double> * >   distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< NRLib::Distribution<double> * >   distr_incl_concentration_; // Pointers to external objects.
  NRLib::Distribution<double>                  * distr_porosity_;           // Pointer to external object.
  DistributionsRockInclusionEvolution          * distr_evolution_;          // Pointer to external object.

  std::vector<double>   expectation_;
  NRLib::Grid2D<double> covariance_;

};

#endif
