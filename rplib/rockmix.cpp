#include "rplib/rockmix.h"

#include "nrlib/exception/exception.hpp"

#include "rplib/solid.h"
#include "rplib/fluid.h"

#include <numeric>

//This file contains two classes RockMixOfRock and RockMixOfSolidAndFluid.

//-------------------------------------- RockMixOfRock ---------------------------------------------------------

RockMixOfRock::RockMixOfRock(const std::vector<Rock*>      & rock,
                             const std::vector<double>     & volume_fraction,
                             DEMTools::MixMethod             mix_method)
: Rock()
{
  // Deep copy of rock:
  rock_.resize(rock.size());
  for (size_t i = 0; i < rock.size(); ++i) {
    rock_[i] = rock[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;

  vp_   = vs_ = rho_ = porosity_ = 0;

  if (rock.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> rho(rock.size()), porosity(rock.size()), k(rock.size()), mu(rock.size());

    for (size_t i = 0; i < rock.size(); i++) {
      double vp, vs;
      rock[i]->GetSeismicParams(vp, vs, rho[i]);
      DEMTools::CalcElasticParamsFromSeismicParams(vs, vp, rho[i], k[i], mu[i]);
      porosity[i] = rock[i]->GetPorosity();
    }

    double k_eff = 0, mu_eff = 0;

    switch (mix_method_) {
      case DEMTools::Hill :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingHill(mu, volume_fraction);
        break;
      case DEMTools::Reuss :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);     // homogeneous
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingReuss(mu, volume_fraction);
        break;
      case DEMTools::Voigt :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mu, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid rock mixing algorithm specified.");
    }

    rho_        = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
    porosity_   = DEMTools::CalcEffectivePorosity(porosity, volume_fraction);

    DEMTools::CalcSeismicParamsFromElasticParams(k_eff, mu_eff, rho_, vp_, vs_);
  }
}


RockMixOfRock::~RockMixOfRock()
{
  for (size_t i = 0; i < rock_.size(); ++i)
    delete rock_[i];
}

Rock *
RockMixOfRock::Clone() const
{
  // Provide base class variables.
  RockMixOfRock * s = new RockMixOfRock(*this);

  // Provide variables specific to RockMixOfRock.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;
  s->porosity_        = this->porosity_;

  size_t n_rocks = this->rock_.size();
  s->rock_.resize(n_rocks);
  for (size_t i = 0; i < n_rocks; ++i)
    s->rock_[i] = this->rock_[i]->Clone();

  return s;
}

RockMixOfRock& RockMixOfRock::operator=(const RockMixOfRock& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;
    porosity_        = rhs.porosity_;

    size_t n_rocks_old = rock_.size();
    for (size_t i = 0; i < n_rocks_old; ++i)
      delete rock_[i];

    size_t n_rocks = rhs.rock_.size();
    rock_.resize(n_rocks);
    for (size_t i = 0; i < n_rocks; ++i)
      rock_[i] = rhs.rock_[i]->Clone();
  }
  return *this;
}

Rock *
RockMixOfRock::Evolve(const std::vector<int>          & delta_time,
                      const std::vector< Rock * >     & rock) const
{
  size_t n_rocks = rock_.size();
  std::vector<Rock*> rock_new(n_rocks);
  for (size_t i = 0; i < n_rocks; ++i)
    rock_new[i] = rock_[i]->Evolve(delta_time, rock);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Rock * rock_mixed_new = new RockMixOfRock(rock_new, volume_fraction, mix_method_);

  // Deep copy taken by constructor of RockMixOfRock, hence delete rock_new here:
  for (size_t i = 0; i < n_rocks; ++i)
    delete rock_new[i];

  return rock_mixed_new;
}

double
RockMixOfRock::GetPorosity() const
{
  return porosity_;
}

void
RockMixOfRock::SetPorosity(double porosity)
{
  porosity_ = porosity;
}

//----------------------------------- RockMixOfSolidAndFluid -------------------------------------------------

RockMixOfSolidAndFluid::RockMixOfSolidAndFluid(const std::vector<Solid*> &  solid,
                                               const std::vector<Fluid*> &  fluid,
                                               const std::vector<double> &  volume_fraction,
                                               double                       porosity,
                                               DEMTools::MixMethod          mix_method)
: Rock()
{

  // Deep copy of solid:
  solid_.resize(solid.size());
  for (size_t i = 0; i < solid.size(); ++i) {
    solid_[i] = solid[i]->Clone();
  }

  // Deep copy of fluid:
  fluid_.resize(fluid.size());
  for (size_t i = 0; i < fluid.size(); ++i) {
    fluid_[i] = fluid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;

  vp_   = vs_ = rho_ = 0;

  porosity_ = porosity;

  size_t size = solid.size() + fluid.size();

  if (size != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> k(size), mu(size), rho(size);
    size_t solid_size = solid.size();
    for (size_t i = 0; i < solid_size; i++)
      solid[i]->GetElasticParams(k[i], mu[i], rho[i]);

    for (size_t i = solid.size(); i < size; i++) {
      fluid[i - solid_size]->GetElasticParams(k[i], rho[i]);
      mu[i] = 0;
    }

    double k_eff = 0, mu_eff = 0;

    switch (mix_method_) {
      case DEMTools::Hill :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingHill(mu, volume_fraction);
        break;
      case DEMTools::Reuss :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);     // homogeneous
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingReuss(mu, volume_fraction);
        break;
      case DEMTools::Voigt :
        k_eff     = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        mu_eff    = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mu, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid rock mixing algorithm specified.");
    }

    rho_        = DEMTools::CalcEffectiveDensity(rho, volume_fraction);

    DEMTools::CalcSeismicParamsFromElasticParams(k_eff, mu_eff, rho_, vp_, vs_);
  }
}


RockMixOfSolidAndFluid::~RockMixOfSolidAndFluid()
{
  for (size_t i = 0; i < solid_.size(); ++i)
    delete solid_[i];

  for (size_t i = 0; i < fluid_.size(); ++i)
    delete fluid_[i];
}

Rock *
RockMixOfSolidAndFluid::Clone() const
{
  // Provide base class variables.
  RockMixOfSolidAndFluid * s = new RockMixOfSolidAndFluid(*this);

  // Provide variables specific to RockMixOfSolidAndFluid.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;
  s->porosity_        = this->porosity_;

  size_t n_solids = this->solid_.size();
  s->solid_.resize(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    s->solid_[i] = this->solid_[i]->Clone();

  size_t n_fluids = this->fluid_.size();
  s->fluid_.resize(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    s->fluid_[i] = this->fluid_[i]->Clone();

  return s;
}

RockMixOfSolidAndFluid& RockMixOfSolidAndFluid::operator=(const RockMixOfSolidAndFluid& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;
    porosity_        = rhs.porosity_;

    size_t n_solids_old = solid_.size();
    for (size_t i = 0; i < n_solids_old; ++i)
      delete solid_[i];

    size_t n_solids = rhs.solid_.size();
    solid_.resize(n_solids);
    for (size_t i = 0; i < n_solids; ++i)
      solid_[i] = rhs.solid_[i]->Clone();

    size_t n_fluids_old = fluid_.size();
    for (size_t i = 0; i < n_fluids_old; ++i)
      delete fluid_[i];

    size_t n_fluids = rhs.fluid_.size();
    fluid_.resize(n_fluids);
    for (size_t i = 0; i < n_fluids; ++i)
      fluid_[i] = rhs.fluid_[i]->Clone();
  }
  return *this;
}

Rock *
RockMixOfSolidAndFluid::Evolve(const std::vector<int>          & delta_time,
                               const std::vector< Rock * >     & /*rock*/) const
{
  size_t n_solids = solid_.size();
  std::vector<Solid*> solid_new(n_solids);

  //NBNB fjellvoll convert the rock to solid and fluid
  std::vector< const Solid* > solid;
  std::vector< const Fluid* > fluid;

  for (size_t i = 0; i < n_solids; ++i)
    solid_new[i] = solid_[i]->Evolve(delta_time, solid);

  size_t n_fluids = fluid_.size();
  std::vector<Fluid*> fluid_new(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    fluid_new[i] = fluid_[i]->Evolve(delta_time, fluid);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Rock * rock_mixed_new = new RockMixOfSolidAndFluid(solid_new, fluid_new, volume_fraction, porosity_, mix_method_);

  // Deep copy taken by constructor of RockMixOfSolidAndFluid, hence delete solid_new, fluid_new here:
  for (size_t i = 0; i < n_solids; ++i)
    delete solid_new[i];

  for (size_t i = 0; i < n_fluids; ++i)
    delete fluid_new[i];

  return rock_mixed_new;
}

double
RockMixOfSolidAndFluid::GetPorosity() const
{
  return porosity_;
}

void
RockMixOfSolidAndFluid::SetPorosity(double porosity)
{
  porosity_ = porosity;
}
