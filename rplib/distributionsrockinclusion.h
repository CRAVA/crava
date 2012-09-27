#ifndef RPLIB_DISTRIBUTIONSROCKINCLUSION_H
#define RPLIB_DISTRIBUTIONSROCKINCLUSION_H

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/rockinclusion.h"
#include "rplib/distributionsrockinclusionevolution.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/grid/grid2d.hpp"

class DistributionWithTrend;

//NBNB fjellvoll not finished yet
//TODO:covariance and expectations functions and variables.

class DistributionsRockInclusion : public DistributionsRock {
public:

                                                 DistributionsRockInclusion(DistributionsSolid                           * distr_solid,
                                                                            DistributionsFluid                           * distr_fluid,
                                                                            std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                                                            std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                                                                            DistributionWithTrend                        * distr_porosity,
                                                                            DistributionsRockInclusionEvolution          * distr_evolution = NULL);

  virtual                                        ~DistributionsRockInclusion();

  virtual Rock                                 * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                                 * UpdateSample(const std::vector<double> & /*corr*/,
                                                              const Rock                & /*rock*/)                    const { return NULL; }

  virtual Pdf3D                                * GeneratePdf(void) const; // Returns NULL.

  virtual bool                                   HasDistribution() const; //dummy function that needs to be implemented

  virtual std::vector<bool>                      HasTrend() const; //dummy function that needs to be implemented

  virtual std::vector<double>                    GetExpectation(const std::vector<double> & /*trend_params*/) const { return expectation_old_; }

  virtual NRLib::Grid2D<double>                  GetCovariance(const std::vector<double> & /*trend_params*/)  const { return covariance_old_; }

private:
  void                                           SampleVpVsRhoExpectationAndCovariance(std::vector<double>   & expectation,
                                                                                       NRLib::Grid2D<double> & covariance);

  DistributionsSolid                           * distr_solid_;              // Pointer to external object.
  DistributionsFluid                           * distr_fluid_;              // Pointer to external object.
  std::vector< DistributionWithTrend * >         distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< DistributionWithTrend * >         distr_incl_concentration_; // Pointers to external objects.
  DistributionWithTrend                        * distr_porosity_;           // Pointer to external object.
  DistributionsRockInclusionEvolution          * distr_evolution_;          // Pointer to external object.

  //NBNB fjellvoll dummy to be removed soon
  std::vector<double>                            expectation_old_;
  NRLib::Grid2D<double>                          covariance_old_;

};

#endif
