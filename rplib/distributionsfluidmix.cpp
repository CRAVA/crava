
#include "rplib/fluid.h"
#include "rplib/distributionsfluidmix.h"

DistributionsFluidMix::DistributionsFluidMix(std::vector< DistributionsFluid * >            distr_fluid,
                                             std::vector< NRLib::Distribution<double> * >   distr_vol_frac,
                                             DEMTools::MixMethod                            mix_method,
                                             DistributionsFluidMixEvolution               * distr_evolution)
: DistributionsFluid()
{
  assert(distr_fluid.size() == distr_vol_frac.size());
  distr_fluid_      = distr_fluid;
  distr_vol_frac_   = distr_vol_frac;
  mix_method_       = mix_method;
  distr_evolution_  = distr_evolution;
}

DistributionsFluidMix::~DistributionsFluidMix(){}

Fluid *
DistributionsFluidMix::GenerateSample(const std::vector<double> & trend_params) const
{

  size_t n_fluids = distr_fluid_.size();
  std::vector<Fluid*> fluid(n_fluids);
  std::vector<double> volume_fraction(n_fluids);

  for(size_t i = 0; i < n_fluids; ++i) {
    fluid[i] = distr_fluid_[i]->GenerateSample(trend_params);
    volume_fraction[i] = distr_vol_frac_[i]->Draw();
  }
  Fluid * fluid_mixed = new FluidMixed(fluid, volume_fraction, mix_method_, distr_evolution_);

  // Deep copy taken by constructor of FluidMixed, hence delete fluid here:
  for(size_t i = 0; i < n_fluids; ++i)
    delete fluid[i];

  return fluid_mixed;
}

bool
DistributionsFluidMix::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsFluidMix::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}
