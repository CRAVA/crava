#include "rplib/solidmixed.h"

SolidMixed::SolidMixed(const std::vector<Solid*>      & solid,
                       const std::vector<double>      & volume_fraction,
                       DistributionsSolidMixEvolution * distr_evolution)
: Solid()
{
  // Deep copy of solid:
  solid_.resize(solid.size());
  for (size_t i = 0; i < solid.size(); ++i) {
    solid_[i] = solid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  distr_evolution_ = distr_evolution;

  k_   = mu_ = rho_ = 0;

  if (solid.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else if (solid.size() != 2) //NBNB fjellvoll this is only tmp.
    throw NRLib::Exception("Invalid arguments:Number of volume fractions is different than 2.");
  else {
    std::vector<double> k(solid.size()), mu(solid.size()), rho(solid.size());
    for (size_t i = 0; i < solid.size(); i++)
      solid[i]->ComputeElasticParams(k[i], mu[i], rho[i]);

    //NBNB fjellvoll only pair-mix allowed at the moment
    k_      = DEMTools::CalcEffectiveElasticModuliUsingHill(k[0], volume_fraction[0], k[1]);
    mu_     = DEMTools::CalcEffectiveElasticModuliUsingHill(mu[0], volume_fraction[0], mu[1]);
    rho_    = DEMTools::CalcEffectiveDensity(rho[0], volume_fraction[0], rho[1]);
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
  s->k_               = this->k_;
  s->mu_              = this->mu_;
  s->rho_             = this->rho_;
  s->volume_fraction_ = this->volume_fraction_;
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

    k_               = rhs.k_;
    mu_              = rhs.mu_;
    rho_             = rhs.rho_;
    volume_fraction_ = rhs.volume_fraction_;
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

void
SolidMixed::ComputeElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

Solid *
SolidMixed::Evolve(const std::vector<int>             & delta_time,
                   const std::vector< const Solid * > & solid) const {
  size_t n_solids = solid_.size();
  std::vector<Solid*> solid_new(n_solids);
  for (size_t i = 0; i < n_solids; ++i)
    solid_new[i] = solid_[i]->Evolve(delta_time, solid);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Solid * solid_mixed_new = new SolidMixed(solid_new, volume_fraction, distr_evolution_);

  // Deep copy taken by constructor of SolidMixed, hence delete solid_new here:
  for (size_t i = 0; i < n_solids; ++i)
    delete solid_new[i];

  return solid_mixed_new;
}

