#include "rplib/fluidmixed.h"

FluidMixed::FluidMixed(const std::vector<Fluid*>      & fluid,
                       const std::vector<double>      & volume_fraction,
                       DEMTools::MixMethod              mix_method,
                       DistributionsFluidMixEvolution * distr_evolution)
: Fluid()
{

  // Deep copy of fluid:
  fluid_.resize(fluid.size());
  for (size_t i = 0; i < fluid.size(); ++i) {
    fluid_[i] = fluid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  mix_method_      = mix_method;
  distr_evolution_ = distr_evolution;

  k_ = rho_ = 0;

  if (fluid.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else {
    std::vector<double> k(fluid.size()), rho(fluid.size());
    for (size_t i = 0; i < fluid.size(); i++)
      fluid[i]->GetElasticParams(k[i], rho[i]);

    switch (mix_method_) {
      case DEMTools::Hill :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingHill(k, volume_fraction);
        break;
      case DEMTools::Reuss :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingReuss(k, volume_fraction);     // homogeneous
        break;
      case DEMTools::Voigt :
        k_    = DEMTools::CalcEffectiveElasticModuliUsingVoigt(k, volume_fraction);
        break;
      default :
        throw NRLib::Exception("Invalid fluid mixing algorithm specified.");
    }
    rho_  = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
  }
}

FluidMixed::~FluidMixed()
{
  for (size_t i = 0; i < fluid_.size(); ++i)
    delete fluid_[i];
}

Fluid *
FluidMixed::Clone() const
{
  // Provide base class variables.
  FluidMixed * s = new FluidMixed(*this);

  // Provide variables specific to FluidMixed.
  s->volume_fraction_ = this->volume_fraction_;
  s->mix_method_      = this->mix_method_;
  s->distr_evolution_ = this->distr_evolution_;

  size_t n_fluids = this->fluid_.size();
  s->fluid_.resize(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    s->fluid_[i] = this->fluid_[i]->Clone();

  return s;
}

FluidMixed& FluidMixed::operator=(const FluidMixed& rhs)
{
  if (this != &rhs) {
    Fluid::operator=(rhs);

    volume_fraction_ = rhs.volume_fraction_;
    mix_method_      = rhs.mix_method_;
    distr_evolution_ = rhs.distr_evolution_;

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


Fluid *
FluidMixed::Evolve(const std::vector<int>             & delta_time,
                   const std::vector< const Fluid * > & fluid) const
{
  size_t n_fluids = fluid_.size();
  std::vector<Fluid*> fluid_new(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    fluid_new[i] = fluid_[i]->Evolve(delta_time, fluid);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Fluid * fluid_mixed_new = new FluidMixed(fluid_new, volume_fraction, mix_method_, distr_evolution_);

  // Deep copy taken by constructor of FluidMixed, hence delete fluid_new here:
  for (size_t i = 0; i < n_fluids; ++i)
    delete fluid_new[i];

  return fluid_mixed_new;
}


