#include "rplib/distributionsrockmix.h"
#include "rplib/rockmix.h"

#include "rplib/rock.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrend.h"

//This file contains two classes DistributionsRockMixOfRock and DistributionsRockMixOfSolidAndFluid.

//-------------------------------------- DistributionsRockMixOfRock ---------------------------------------------------------


DistributionsRockMixOfRock::DistributionsRockMixOfRock(const std::vector< DistributionsRock * >           & distr_rock,
                                                       const std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                                       DEMTools::MixMethod                                  mix_method)
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

//-------------------------------------- DistributionsRockMixOfSolidAndFluid ---------------------------------------------------------


DistributionsRockMixOfSolidAndFluid::DistributionsRockMixOfSolidAndFluid(const std::vector< DistributionsSolid * >           & distr_solid,
                                                                         const std::vector< DistributionsFluid * >           & distr_fluid,
                                                                         const std::vector< DistributionWithTrend * >        & distr_vol_frac_solid,
                                                                         const std::vector< DistributionWithTrend * >        & distr_vol_frac_fluid,
                                                                         DEMTools::MixMethod                                   mix_method)
: distr_solid_(distr_solid),
  distr_fluid_(distr_fluid),
  distr_vol_frac_solid_(distr_vol_frac_solid),
  distr_vol_frac_fluid_(distr_vol_frac_fluid),
  mix_method_(mix_method)
{
}

DistributionsRockMixOfSolidAndFluid::~DistributionsRockMixOfSolidAndFluid()
{

}

Rock *
DistributionsRockMixOfSolidAndFluid::GenerateSample(const std::vector<double> & trend_params) const
{
  size_t n_fluids      =   distr_fluid_.size();
  size_t n_solids      =   distr_solid_.size();

  std::vector<Solid*>      solid(n_solids);
  std::vector<double>      volume_fraction_solid(n_solids, 0.0);

  size_t missing_index = n_solids + n_fluids;

  for(size_t i = 0; i < n_solids; ++i) {
    solid[i] = distr_solid_[i]->GenerateSample(trend_params);

    if (distr_vol_frac_solid_[i])
      volume_fraction_solid[i] = distr_vol_frac_solid_[i]->ReSample(trend_params[0], trend_params[1]);
    else
      missing_index    = i;
  }

  std::vector<Fluid*>      fluid(n_fluids);
  std::vector<double>      volume_fraction_fluid(n_fluids, 0.0);

  for(size_t i = 0; i < n_fluids; ++i) {
    fluid[i] = distr_fluid_[i]->GenerateSample(trend_params);

    if (distr_vol_frac_fluid_[i])
      volume_fraction_fluid[i] = distr_vol_frac_fluid_[i]->ReSample(trend_params[0], trend_params[1]);
    else
      missing_index    = i + n_solids;
  }

  if (missing_index != n_fluids + n_solids) {
    double sum = 0.0;

    for (size_t i = 0; i < volume_fraction_solid.size(); ++i)
      sum += volume_fraction_solid[i];

    for (size_t i = 0; i < volume_fraction_fluid.size(); ++i)
      sum += volume_fraction_fluid[i];

   if(missing_index < n_solids)
     volume_fraction_solid[missing_index] = 1.0 - sum;
   else
    volume_fraction_fluid[missing_index - n_solids] = 1.0 - sum;
  }

  Rock * rock_mixed = new RockMixOfSolidAndFluid(solid, fluid, volume_fraction_solid, volume_fraction_fluid, mix_method_);

  // Deep copy taken by constructor of RockMixOfSolidAndFluid, hence delete solid here:
  for(size_t i = 0; i < n_solids; ++i)
    delete solid[i];

  for(size_t i = 0; i < n_fluids; ++i)
    delete fluid[i];

  return rock_mixed;
}

std::vector<double>
DistributionsRockMixOfSolidAndFluid::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  std::vector<double> dummy;
  return(dummy);
}

NRLib::Grid2D<double>
DistributionsRockMixOfSolidAndFluid::GetCovariance(const std::vector<double> & /*trend_params*/) const
{
  NRLib::Grid2D<double> dummy;
  return(dummy);
}

Pdf3D *
DistributionsRockMixOfSolidAndFluid::GeneratePdf() const
{
  Pdf3D * dummy = NULL;
  return(dummy);
}

bool
DistributionsRockMixOfSolidAndFluid::HasDistribution() const
{
  bool has_distribution = false;

  size_t n_fluids = distr_fluid_.size();
  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_fluids; i++) {
    if(distr_fluid_[i]->HasDistribution() == true)
      has_distribution = true;
    else if(distr_vol_frac_fluid_[i] != NULL && distr_vol_frac_fluid_[i]->GetIsDistribution() == true)
      has_distribution = true;
  }
  for(size_t i=0; i<n_solids; i++) {
    if(distr_solid_[i]->HasDistribution() == true)
      has_distribution = true;
    else if(distr_vol_frac_solid_[i] != NULL && distr_vol_frac_solid_[i]->GetIsDistribution() == true)
      has_distribution = true;
  }

  return has_distribution;

}

std::vector<bool>
DistributionsRockMixOfSolidAndFluid::HasTrend() const
{
  std::vector<bool> has_trend(2);

  for(int i=0; i<2; i++)
    has_trend[i] = false;

  size_t n_fluids = distr_fluid_.size();

  for(size_t i=0; i<n_fluids; i++) {
    std::vector<bool> fluid_trend  = distr_fluid_[i]->HasTrend();

    std::vector<bool> volume_trend(2,false);
    if(distr_vol_frac_fluid_[i] != NULL)
       volume_trend = distr_vol_frac_fluid_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {
      if(fluid_trend[j] == true)
        has_trend[j] = true;
      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_solids; i++) {
    std::vector<bool> solid_trend  = distr_solid_[i]->HasTrend();

    std::vector<bool> volume_trend(2,false);
    if(distr_vol_frac_solid_[i] != NULL)
       volume_trend = distr_vol_frac_solid_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {
      if(solid_trend[j] == true)
        has_trend[j] = true;
      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  return has_trend;
}

bool
DistributionsRockMixOfSolidAndFluid::GetIsOkForBounding() const
{
  bool is_ok_for_bounding = false;

  if(mix_method_ == DEMTools::Reuss || mix_method_ == DEMTools::Voigt) {
    if(distr_fluid_.size() == 1 && distr_solid_.size() == 1) {
      is_ok_for_bounding = true;

      bool fluid_distr = distr_fluid_[0]->HasDistribution();
      bool solid_distr = distr_solid_[0]->HasDistribution();
      if(fluid_distr == true || solid_distr == true)
        is_ok_for_bounding = false;

      std::vector<bool> fluid_trend = distr_fluid_[0]->HasTrend();
      std::vector<bool> solid_trend = distr_solid_[0]->HasTrend();

      for(int i=0; i<2; i++) {
        if(fluid_trend[i] == true || solid_trend[i] == true)
          is_ok_for_bounding = false;
      }
    }
  }

  return(is_ok_for_bounding);
}

