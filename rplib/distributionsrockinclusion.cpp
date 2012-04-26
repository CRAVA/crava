#include "rplib/distributionsrockinclusion.h"

DistributionsRockInclusion::DistributionsRockInclusion(DistributionsSolid                           * distr_solid,
                                                       DistributionsFluid                           * distr_fluid,
                                                       std::vector< NRLib::Distribution<double> * >   distr_incl_spectrum,
                                                       std::vector< NRLib::Distribution<double> * >   distr_incl_concentration,
                                                       NRLib::Distribution<double>                  * distr_porosity,
                                                       DistributionsRockInclusionEvolution          * distr_evolution)
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

DistributionsRockInclusion::~DistributionsRockInclusion(){}

Rock *
DistributionsRockInclusion::GenerateSample(const std::vector<double> & /*trend_params*/) const
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

std::vector<double>
DistributionsRockInclusion::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  size_t n_incl = distr_incl_spectrum_.size();
  std::vector<double> expectation(2*n_incl,0.0);
  for (size_t i = 0; i < n_incl; ++i) {
    expectation[i]          = distr_incl_spectrum_[i]->Quantile(0.5);
    expectation[i + n_incl] = distr_incl_concentration_[i]->Quantile(0.5);
  }
  return(expectation);
}

NRLib::Grid2D<double>
DistributionsRockInclusion::GetCovariance(const std::vector<double> & /*trend_params*/)  const
{
  size_t n_incl = distr_incl_spectrum_.size();
  NRLib::Grid2D<double> covariance(2*n_incl, 2*n_incl, 0.0);
  return(covariance);
}

Pdf3D *
DistributionsRockInclusion::GeneratePdf(void) const
{
  Pdf3D * pdf3D = NULL;
  return pdf3D;
}
