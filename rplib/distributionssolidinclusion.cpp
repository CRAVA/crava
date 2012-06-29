
#include "rplib/distributionssolid.h"
#include "rplib/distributionssolidinclusion.h"
#include "rplib/solid.h"

//NBNB Marit: Uferdig klasse

DistributionsSolidInclusion::DistributionsSolidInclusion(std::vector<DistributionsSolid  *>   distr_solid, //Marit: Blande mange solider. Hvor mange, og hvordan?
                                                         std::vector<double>                  distr_incl_spectrum,
                                                         std::vector<double>                  distr_incl_concentration)
                                                         //DistributionssolidInclusionEvolution          * distr_evolution)
: DistributionsSolid()
{
  assert( distr_incl_spectrum.size() == distr_incl_concentration.size() );

  distr_solid_              = distr_solid;
  distr_incl_spectrum_      = distr_incl_spectrum;
  distr_incl_concentration_ = distr_incl_concentration;
  //distr_porosity_           = distr_porosity;
  //distr_evolution_          = distr_evolution;
}

DistributionsSolidInclusion::~DistributionsSolidInclusion(){}

Solid *
DistributionsSolidInclusion::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Solid * solid1 = distr_solid_[0]->GenerateSample(trend_params);//Marit: Hardkodet, blir feil
  //Solid * solid2 = distr_solid_[1]->GenerateSample(trend_params);
  size_t n_incl  = distr_incl_spectrum_.size();

  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl);

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i] = distr_incl_spectrum_[i];
    inclusion_concentration[i] = distr_incl_concentration_[i];
  }

  Solid * new_solid = NULL;
  //Rock * new_rock = new RockInclusion(solid, solid, inclusion_spectrum, inclusion_concentration, porosity, distr_evolution_);

  return new_solid;
}
