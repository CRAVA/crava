#include "rplib/fluidmixed.h"

FluidMixed::FluidMixed(const std::vector<Fluid*>      & fluid,
                       const std::vector<double>      & volume_fraction,
                       DistributionsFluidMixEvolution * distr_evolution)
: Fluid()
{

  // Deep copy of fluid:
  fluid_.resize(fluid.size());
  for (size_t i = 0; i < fluid.size(); ++i) {
    fluid_[i] = fluid[i]->Clone();
  }

  volume_fraction_ = volume_fraction;
  distr_evolution_ = distr_evolution;

  k_ = rho_ = 0;

  if (fluid.size() != volume_fraction.size())
    throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
  else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
    throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
  else if (fluid.size() != 2) //NBNB fjellvoll this is only tmp.
    throw NRLib::Exception("Invalid arguments:Number of volume fractions is different than 2.");
  else {
    //NBNB fjellvoll only pair-mix allowed at the moment
    std::vector<double> k(fluid.size()), rho(fluid.size());
    for (size_t i = 0; i < fluid.size(); i++)
      fluid[i]->GetElasticParams(k[i], rho[i]);
    k_    = DEMTools::CalcEffectiveElasticModuliUsingReuss(k[0], volume_fraction[0], k[1]);     // homogeneous
    rho_  = DEMTools::CalcEffectiveDensity(rho[0], volume_fraction[0], rho[1]);
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
  s->k_ = this->k_;
  s->rho_ = this->rho_;
  s->volume_fraction_ = this->volume_fraction_;
  s->distr_evolution_ = this->distr_evolution_;

  size_t n_fluids = this->fluid_.size();
  s->fluid_.resize(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    s->fluid_[i] = this->fluid_[i]->Clone();

  return s;
}

void
FluidMixed::ComputeElasticParams(double /*temp*/, double /*pore_pressure*/)
{
  // Done in constructor.
}

void
FluidMixed::GetElasticParams(double& k, double& rho) const
{
  k     = k_;
  rho   = rho_;
}

Fluid *
FluidMixed::Evolve(const std::vector<int>       & delta_time,
                   const std::vector< Fluid * > & fluid) const
{
  size_t n_fluids = fluid_.size();
  std::vector<Fluid*> fluid_new(n_fluids);
  for (size_t i = 0; i < n_fluids; ++i)
    fluid_new[i] = fluid_[i]->Evolve(delta_time, fluid);

  std::vector<double> volume_fraction = volume_fraction_; // Evolve when model is defined.
  Fluid * fluid_mixed_new = new FluidMixed(fluid_new, volume_fraction, distr_evolution_);

  // Deep copy taken by constructor of FluidMixed, hence delete fluid_new here:
  for (size_t i = 0; i < n_fluids; ++i)
    delete fluid_new[i];

  return fluid_mixed_new;
}


