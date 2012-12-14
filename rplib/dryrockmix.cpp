#include "rplib/dryrockmix.h"

#include "nrlib/exception/exception.hpp"

#include <numeric>

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
    std::vector<double> k(dryrock.size()), mu(dryrock.size()), rho(dryrock.size());
    for (size_t i = 0; i < dryrock.size(); i++)
      dryrock[i]->GetElasticParams(k[i], mu[i], rho[i]);

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
        throw NRLib::Exception("Invalid dryrock mixing algorithm specified.");
    }
    rho_    = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
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

DryRock *
DryRockMix::Evolve(const std::vector<int>               & delta_time,
                   const std::vector< const DryRock * > & dryrock) const {
  size_t n_dryrocks = dryrock_.size();
  std::vector<DryRock*> dryrock_new(n_dryrocks);
  for (size_t i = 0; i < n_dryrocks; ++i)
    dryrock_new[i] = dryrock_[i]->Evolve(delta_time, dryrock);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  DryRock * dryrock_mixed_new = new DryRockMix(dryrock_new, volume_fraction, u_, mix_method_);

  // Deep copy taken by constructor of DryRockMix, hence delete dryrock_new here:
  for (size_t i = 0; i < n_dryrocks; ++i)
    delete dryrock_new[i];

  return dryrock_mixed_new;
}

