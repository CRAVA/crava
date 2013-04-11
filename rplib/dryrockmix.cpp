#include "rplib/dryrockmix.h"

#include "rplib/solid.h"

#include "nrlib/exception/exception.hpp"

#include <numeric>
#include <algorithm>

//This file contains two classes DryRockMix and DryRockMixOfDryRockAndSolid.

//-------------------------------------- DryRockMix ---------------------------------------------------------

DryRockMix::DryRockMix(const std::vector<DryRock*>    & dryrock,
                       const std::vector<double>      & volume_fraction,
                       const std::vector<double>      & u,
                       DEMTools::MixMethod              mix_method)
: DryRock()
{
  u_ = u;  // u contains independent samples used in quantiles of volume_fraction. Use RMISSING for fraction to be calculated from other fractions

  // Deep copy of dryrock:
  dryrock_.resize(dryrock.size());
  for (size_t i = 0; i < dryrock.size(); ++i) {
    dryrock_[i] = dryrock[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;

  k_   = mu_ = rho_ = 0;

  if (dryrock.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> k(dryrock.size()), mu(dryrock.size()), rho(dryrock.size()), total_porosity(dryrock.size()), mineral_moduli_k(dryrock.size());
    for (size_t i = 0; i < dryrock.size(); i++) {
      dryrock[i]->GetElasticParams(k[i], mu[i], rho[i]);
      total_porosity[i]   = dryrock[i]->GetTotalPorosity();
      mineral_moduli_k[i] = dryrock[i]->GetMineralModuliK();
     }

    switch (mix_method_) {
      case DEMTools::Hill :
        k_                = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        mu_               = DEMTools::CalcEffectiveElasticModuliUsingHill(mu, volume_fraction);
        mineral_moduli_k_ = DEMTools::CalcEffectiveElasticModuliUsingHill(mineral_moduli_k, volume_fraction);
        break;
      case DEMTools::Reuss :
        k_                = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);     // homogeneous
        mu_               = DEMTools::CalcEffectiveElasticModuliUsingReuss(mu, volume_fraction);
        mineral_moduli_k_ = DEMTools::CalcEffectiveElasticModuliUsingReuss(mineral_moduli_k, volume_fraction);
        break;
      case DEMTools::Voigt :
        k_                = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        mu_               = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mu, volume_fraction);
        mineral_moduli_k_ = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mineral_moduli_k, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid dryrock mixing algorithm specified.");
    }
    rho_            = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
    total_porosity_ = DEMTools::CalcEffectivePorosity(total_porosity, volume_fraction);

  }
}


DryRockMix::~DryRockMix()
{
  for (size_t i = 0; i < dryrock_.size(); ++i)
    delete dryrock_[i];
}

DryRock *
DryRockMix::Clone() const {
  // Provide base class variables.
  DryRockMix * s = new DryRockMix(*this);

  // Provide variables specific to DryRockMix.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;

  size_t n_dryrocks = this->dryrock_.size();
  s->dryrock_.resize(n_dryrocks);
  for (size_t i = 0; i < n_dryrocks; ++i)
    s->dryrock_[i] = this->dryrock_[i]->Clone();

  return s;
}

DryRockMix& DryRockMix::operator=(const DryRockMix& rhs)
{
  if (this != &rhs) {
    DryRock::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;

    size_t n_dryrocks_old = dryrock_.size();
    for (size_t i = 0; i < n_dryrocks_old; ++i)
      delete dryrock_[i];

    size_t n_dryrocks = rhs.dryrock_.size();
    dryrock_.resize(n_dryrocks);
    for (size_t i = 0; i < n_dryrocks; ++i)
      dryrock_[i] = rhs.dryrock_[i]->Clone();
  }
  return *this;
}

//-------------------------------------- DryRockMixOfDryRockAndSolid ---------------------------------------------------------


DryRockMixOfDryRockAndSolid::DryRockMixOfDryRockAndSolid(const std::vector<DryRock*> &  dryrock,
                                                         const std::vector<Solid*>   &  solid,
                                                         const std::vector<double>   &  volume_fraction_dryrock,
                                                         const std::vector<double>   &  volume_fraction_solid,                                                         
                                                         const std::vector<double>   &  u,
                                                         DEMTools::MixMethod            mix_method)
: DryRock()
{
  // Deep copy of dryrock:
  dryrock_.resize(dryrock.size());
  for (size_t i = 0; i < dryrock.size(); ++i) {
    dryrock_[i] = dryrock[i]->Clone();
  }

  // Deep copy of solid:
  solid_.resize(solid.size());
  for (size_t i = 0; i < solid.size(); ++i) {
    solid_[i] = solid[i]->Clone();
  }

  volume_fraction_solid_ = volume_fraction_solid;
  volume_fraction_dryrock_ = volume_fraction_dryrock;
  u_                     = u; // u contains independent samples used in quantiles of (volume_fraction_dryrock, volume_fraction_solid). Use RMISSING for fraction to be calculated from other fractions
  mix_method_            = mix_method;

  if (solid_.size() != volume_fraction_solid_.size())
    throw NRLib::Exception("Invalid arguments: Number of solid properties are different from number of volume fractions.");

  else if (dryrock_.size() != volume_fraction_dryrock_.size())
    throw NRLib::Exception("Invalid arguments: Number of dryrock properties are different from number of volume fractions.");


  ComputeSeismicVariables();
}


DryRockMixOfDryRockAndSolid::~DryRockMixOfDryRockAndSolid()
{
  for (size_t i = 0; i < solid_.size(); ++i)
    delete solid_[i];

  for (size_t i = 0; i < dryrock_.size(); ++i)
    delete dryrock_[i];
}

DryRock *
DryRockMixOfDryRockAndSolid::Clone() const
{
  // Provide base class variables.
  DryRockMixOfDryRockAndSolid * s = new DryRockMixOfDryRockAndSolid(*this);

  // Provide variables specific to DryRockMixOfDryRockAndSolid.
  s->volume_fraction_solid_ = this->volume_fraction_solid_;
  s->volume_fraction_dryrock_ = this->volume_fraction_dryrock_;
  s->mix_method_            = this->mix_method_;

  size_t n_dryrocks = this->dryrock_.size();
  s->dryrock_.resize(n_dryrocks);
  for (size_t i = 0; i < n_dryrocks; ++i)
    s->dryrock_[i] = this->dryrock_[i]->Clone();

  size_t n_solids = this->solid_.size();
  s->solid_.resize(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    s->solid_[i] = this->solid_[i]->Clone();

  return s;
}

DryRockMixOfDryRockAndSolid& DryRockMixOfDryRockAndSolid::operator=(const DryRockMixOfDryRockAndSolid& rhs)
{
  if (this != &rhs) {
    DryRock::operator=(rhs);

    volume_fraction_solid_ = rhs.volume_fraction_solid_;
    volume_fraction_dryrock_ = rhs.volume_fraction_dryrock_;
    mix_method_      = rhs.mix_method_;

    size_t n_solids_old = solid_.size();
    for (size_t i = 0; i < n_solids_old; ++i)
      delete solid_[i];

    size_t n_solids = rhs.solid_.size();
    solid_.resize(n_solids);
    for (size_t i = 0; i < n_solids; ++i)
      solid_[i] = rhs.solid_[i]->Clone();

    size_t n_dryrocks_old = dryrock_.size();
    for (size_t i = 0; i < n_dryrocks_old; ++i)
      delete dryrock_[i];

    size_t n_dryrocks = rhs.dryrock_.size();
    dryrock_.resize(n_dryrocks);
    for (size_t i = 0; i < n_dryrocks; ++i)
      dryrock_[i] = rhs.dryrock_[i]->Clone();
  }
  return *this;
}

void
DryRockMixOfDryRockAndSolid::ComputeSeismicVariables()
{
  size_t size = solid_.size() + dryrock_.size();

  std::vector<double> volume_fraction(size);
  std::vector<double> k(size);
  std::vector<double> mu(size);
  std::vector<double> rho(size);

  std::vector<double> mineral_moduli_k_all;
  mineral_moduli_k_all.reserve(size);

  size_t solid_size = solid_.size();
  for (size_t i = 0; i < solid_size; i++) {
    solid_[i]->GetElasticParams(k[i], mu[i], rho[i]);
    volume_fraction[i] = volume_fraction_solid_[i];
    mineral_moduli_k_all.push_back(k[i]);
  }

  std::vector<double> total_porosity(size, 0.0);

  for (size_t i = solid_.size(); i < size; i++) {
    dryrock_[i - solid_size]->GetElasticParams(k[i], mu[i], rho[i]);
    total_porosity[i] = dryrock_[i - solid_size]->GetTotalPorosity();
    volume_fraction[i] = volume_fraction_dryrock_[i - solid_size];
    mineral_moduli_k_all.push_back(dryrock_[i - solid_size]->GetMineralModuliK());
  }

  switch (mix_method_) {
      case DEMTools::Hill :
        k_        = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        mu_       = DEMTools::CalcEffectiveElasticModuliUsingHill(mu, volume_fraction);
        break;

      case DEMTools::Reuss :
        k_        = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);
        mu_       = DEMTools::CalcEffectiveElasticModuliUsingReuss(mu, volume_fraction);
        break;

      case DEMTools::Voigt :
        k_        = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        mu_       = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mu, volume_fraction);
        break;

      default :
        throw NRLib::Exception("Invalid rock mixing algorithm specified.");
  }

  rho_              = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
  total_porosity_   = DEMTools::CalcEffectivePorosity(total_porosity, volume_fraction);
  mineral_moduli_k_ = DEMTools::CalcEffectiveElasticModuliUsingHill(mineral_moduli_k_all, volume_fraction); //NBNB fjellvoll check with Erling
}


