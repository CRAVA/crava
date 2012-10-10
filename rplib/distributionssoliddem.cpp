#include "rplib/distributionssoliddem.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/soliddem.h"

DistributionsSolidDEM::DistributionsSolidDEM(DistributionsSolid                           * distr_solid,
                                             DistributionsSolid                           * distr_solid_inc,
                                             std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                             std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                                             DistributionWithTrend                        * distr_porosity)
: DistributionsSolid()
{
  assert( distr_incl_spectrum.size() == distr_incl_concentration.size() );

  distr_solid_              = distr_solid;
  distr_solid_inc_          = distr_solid_inc;
  distr_incl_spectrum_      = distr_incl_spectrum;
  distr_incl_concentration_ = distr_incl_concentration;
  distr_porosity_           = distr_porosity;
}

DistributionsSolidDEM::~DistributionsSolidDEM(){}

Solid *
DistributionsSolidDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  Solid * solid     = distr_solid_->GenerateSample(trend_params);
  Solid * solid_inc = distr_solid_inc_->GenerateSample(trend_params);
  size_t  n_incl    = distr_incl_spectrum_.size();

  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl);

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i]      = distr_incl_spectrum_[i]->ReSample(trend_params[0], trend_params[1]);
    inclusion_concentration[i] = distr_incl_concentration_[i]->ReSample(trend_params[0], trend_params[1]);
  }
  double  porosity  = distr_porosity_->ReSample(trend_params[0], trend_params[1]);
  Solid * new_solid = new SolidDEM(solid, solid_inc, inclusion_spectrum, inclusion_concentration, porosity);

  // Deep copy taken by constructor of SolidDEM, hence delete
  // solid and solid_inc here:
  delete solid;
  delete solid_inc;

  return new_solid;
}

bool
DistributionsSolidDEM::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsSolidDEM::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Solid *
DistributionsSolidDEM::UpdateSample(const std::vector< double > & /*corr*/,
                                    const Solid                 & /*solid*/) const {

  return NULL;
}
