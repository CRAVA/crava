#include "rplib/distributionsrockmix.h"
#include "rplib/rockmix.h"

#include "rplib/rock.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrend.h"


DistributionsRockMixOfRock::DistributionsRockMixOfRock(std::vector< DistributionsRock * >           & distr_rock,
                                                       std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                                       DEMTools::MixMethod                            mix_method)
: distr_rock_(distr_rock),
  distr_vol_frac_(distr_vol_frac),
  mix_method_(mix_method)
{

}

DistributionsRockMixOfRock::~DistributionsRockMixOfRock()
{

}

Rock *
DistributionsRockMixOfRock::GenerateSample(const std::vector<double> & trend_params) const
{

  size_t n_rocks      =    distr_rock_.size();
  std::vector<Rock*>       rock(n_rocks);
  std::vector<double>      volume_fraction(n_rocks, 0.0);

  size_t missing_index = n_rocks;
  for(size_t i = 0; i < n_rocks; ++i) {
    rock[i] = distr_rock_[i]->GenerateSample(trend_params);
    if (distr_vol_frac_[i])
      volume_fraction[i] = distr_vol_frac_[i]->ReSample(trend_params[0], trend_params[1]);
    else
      missing_index    = i;
  }

  if (missing_index != n_rocks) {
    double sum = 0.0;
    for (size_t i = 0; i < volume_fraction.size(); ++i)
      sum += volume_fraction[i];

    volume_fraction[missing_index] = 1.0 - sum;
  }

  Rock * rock_mixed = new RockMixOfRock(rock, volume_fraction, mix_method_);

  // Deep copy taken by constructor of RockMixOfRock, hence delete rock here:
  for(size_t i = 0; i < n_rocks; ++i)
    delete rock[i];

  return rock_mixed;
}

std::vector<double>
DistributionsRockMixOfRock::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  std::vector<double> dummy;
  return(dummy);
}

NRLib::Grid2D<double>
DistributionsRockMixOfRock::GetCovariance(const std::vector<double> & /*trend_params*/) const
{
  NRLib::Grid2D<double> dummy;
  return(dummy);
}

Pdf3D *
DistributionsRockMixOfRock::GeneratePdf() const
{
  Pdf3D * dummy = NULL;
  return(dummy);
}

bool
DistributionsRockMixOfRock::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsRockMixOfRock::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}
