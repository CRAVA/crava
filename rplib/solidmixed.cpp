#include "rplib/solidmixed.h"

#include "nrlib/exception/exception.hpp"

SolidMixed::SolidMixed(const std::vector<Solid*>      & solid,
                       const std::vector<double>      & volume_fraction,
                       DEMTools::MixMethod              mix_method,
                       DistributionsSolidMixEvolution * distr_evolution)
: Solid()
{
  // Deep copy of solid:
  solid_.resize(solid.size());
  for (size_t i = 0; i < solid.size(); ++i) {
    solid_[i] = solid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;
  distr_evolution_ = distr_evolution;

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


SolidMixed::~SolidMixed()
{
  for (size_t i = 0; i < solid_.size(); ++i)
    delete solid_[i];
}

Solid *
SolidMixed::Clone() const {
  // Provide base class variables.
  SolidMixed * s = new SolidMixed(*this);

  // Provide variables specific to SolidMixed.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;
  s->distr_evolution_ = this->distr_evolution_;

  size_t n_solids = this->solid_.size();
  s->solid_.resize(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    s->solid_[i] = this->solid_[i]->Clone();

  return s;
}

SolidMixed& SolidMixed::operator=(const SolidMixed& rhs)
{
  if (this != &rhs) {
    Solid::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;
    distr_evolution_ = rhs.distr_evolution_;

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

Solid *
SolidMixed::Evolve(const std::vector<int>             & delta_time,
                   const std::vector< const Solid * > & solid) const {
  size_t n_solids = solid_.size();
  std::vector<Solid*> solid_new(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    solid_new[i] = solid_[i]->Evolve(delta_time, solid);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Solid * solid_mixed_new = new SolidMixed(solid_new, volume_fraction, mix_method_, distr_evolution_);

  // Deep copy taken by constructor of SolidMixed, hence delete solid_new here:
  for (size_t i = 0; i < n_solids; ++i)
    delete solid_new[i];

  return solid_mixed_new;
}

