#include "rplib/solidmix.h"

#include "nrlib/exception/exception.hpp"

#include <numeric>

SolidMix::SolidMix(const std::vector<Solid*>      & solid,
                   const std::vector<double>      & volume_fraction,
                   const std::vector<double>      & u,
                   DEMTools::MixMethod              mix_method)
: Solid()
{
  u_ = u;  // u contains independent samples used in quantiles of volume_fraction. Use RMISSING for fraction to be calculated from other fractions

  // Deep copy of solid:
  solid_.resize(solid.size());
  for (size_t i = 0; i < solid.size(); ++i) {
    solid_[i] = solid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;

  k_   = mu_ = rho_ = 0;

  if (solid.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> k(solid.size()), mu(solid.size()), rho(solid.size());
    for (size_t i = 0; i < solid.size(); i++)
      solid[i]->GetElasticParams(k[i], mu[i], rho[i]);

    switch (mix_method_) {
      case DEMTools::Hill :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        mu_   = DEMTools::CalcEffectiveElasticModuliUsingHill(mu, volume_fraction);
        break;
      case DEMTools::Reuss :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);     // homogeneous
        mu_   = DEMTools::CalcEffectiveElasticModuliUsingReuss(mu, volume_fraction);
        break;
      case DEMTools::Voigt :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        mu_   = DEMTools::CalcEffectiveElasticModuliUsingVoigt(mu, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid solid mixing algorithm specified.");
    }
    rho_    = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
  }
}


SolidMix::~SolidMix()
{
  for (size_t i = 0; i < solid_.size(); ++i)
    delete solid_[i];
}

Solid *
SolidMix::Clone() const {
  // Provide base class variables.
  SolidMix * s = new SolidMix(*this);

  // Provide variables specific to SolidMix.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;

  size_t n_solids = this->solid_.size();
  s->solid_.resize(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    s->solid_[i] = this->solid_[i]->Clone();

  return s;
}

SolidMix& SolidMix::operator=(const SolidMix& rhs)
{
  if (this != &rhs) {
    Solid::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;

    size_t n_solids_old = solid_.size();
    for (size_t i = 0; i < n_solids_old; ++i)
      delete solid_[i];

    size_t n_solids = rhs.solid_.size();
    solid_.resize(n_solids);
    for (size_t i = 0; i < n_solids; ++i)
      solid_[i] = rhs.solid_[i]->Clone();
  }
  return *this;
}

