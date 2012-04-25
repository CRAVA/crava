#ifndef RPLIB_DISTRIBUTIONSROCKINCLUSION_H
#define RPLIB_DISTRIBUTIONSROCKINCLUSION_H

#include "rplib/distributionsrock.h"
#include "rplib/rockinclusion.h"
#include "rplib/distributionsrockinclusionevolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsRockInclusion : public DistributionsRock {
public:

  DistributionsRockInclusion(DistributionsSolid                           * distr_solid,
                             DistributionsFluid                           * distr_fluid,
                             std::vector< NRLib::Distribution<double> * >   distr_incl_spectrum,
                             std::vector< NRLib::Distribution<double> * >   distr_incl_concentration,
                             NRLib::Distribution<double>                  * distr_porosity,
                             DistributionsRockInclusionEvolution          * distr_evolution = NULL)
  : DistributionsRock()
  {
    assert( distr_incl_spectrum.size() == distr_incl_concentration.size() );

    distr_solid_              = distr_solid;
    distr_fluid_              = distr_fluid;
    distr_incl_spectrum_      = distr_incl_spectrum;
    distr_incl_concentration_ = distr_incl_concentration;
    distr_porosity_           = distr_porosity;
    distr_evolution_          = distr_evolution;
  }

  virtual ~DistributionsRockInclusion(){}

  virtual Rock * GenerateSample(const std::vector<double> & /*trend_params*/) const
  {
    Solid * solid = distr_solid_->GenerateSample();
    Fluid * fluid = distr_fluid_->GenerateSample();
    size_t n_incl = distr_incl_spectrum_.size();
    std::vector<double> inclusion_spectrum(n_incl);
    std::vector<double> inclusion_concentration(n_incl);

    for (size_t i = 0; i < n_incl; ++i) {
      inclusion_spectrum[i] = distr_incl_spectrum_[i]->Draw();
      inclusion_concentration[i] = distr_incl_concentration_[i]->Draw();
    }
    double porosity = distr_porosity_->Draw();
    Rock * new_rock = new RockInclusion(solid, fluid, inclusion_spectrum, inclusion_concentration, porosity, distr_evolution_);

    // Deep copy taken by constructor of RockInclusion, hence delete
    // solid and fluid here:
    delete solid;
    delete fluid;

    return new_rock;
  }

  virtual std::vector<double>   GetExpectation(const std::vector<double> & /*trend_params*/) const
  {
    // The NRLib class Distributions reports on quantiles, not expectations.
    // Until that is included, this implementation of GetExpectation return the median.
    // For Delta and Normal distributions this equals the expectation.
    size_t n_incl = distr_incl_spectrum_.size();
    std::vector<double> expectation(2*n_incl,0.0);
    for (size_t i = 0; i < n_incl; ++i) {
      expectation[i]          = distr_incl_spectrum_[i]->Quantile(0.5);
      expectation[i + n_incl] = distr_incl_concentration_[i]->Quantile(0.5);
    }
    return(expectation);
  }

  virtual NRLib::Grid2D<double> GetCovariance(const std::vector<double> & /*trend_params*/)  const
  {
    size_t n_incl = distr_incl_spectrum_.size();
    NRLib::Grid2D<double> covariance(2*n_incl, 2*n_incl, 0);
    return(covariance);
  }

  virtual Pdf3D * GeneratePdf(void) const
  {
    Pdf3D * pdf3D = NULL; //FAKE
    return pdf3D;
  }

private:
  DistributionsSolid                           * distr_solid_;
  DistributionsFluid                           * distr_fluid_;
  std::vector< NRLib::Distribution<double> * >   distr_incl_spectrum_;
  std::vector< NRLib::Distribution<double> * >   distr_incl_concentration_;
  NRLib::Distribution<double>                  * distr_porosity_;
  DistributionsRockInclusionEvolution          * distr_evolution_;
};

#endif
