
#include "rplib/fluid.h"

#include "rplib/fluidmix.h"
#include "rplib/distributionsfluidmix.h"

#include "rplib/distributionwithtrend.h"

#include "src/definitions.h"

DistributionsFluidMix::DistributionsFluidMix(std::vector< DistributionsFluid * >          & distr_fluid,
                                             std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                             DEMTools::MixMethod                            mix_method)
: DistributionsFluid()
{
  assert(distr_fluid.size() == distr_vol_frac.size());
  distr_fluid_      = distr_fluid;
  distr_vol_frac_   = distr_vol_frac;
  mix_method_       = mix_method;
}

DistributionsFluidMix::~DistributionsFluidMix()
{
}

Fluid *
DistributionsFluidMix::GenerateSample(const std::vector<double> & trend_params) const
{

  size_t n_fluids = distr_fluid_.size();

  std::vector<double> u(n_fluids, RMISSING);
  for(size_t i=0; i<n_fluids; i++) {
    if(distr_vol_frac_[i] != NULL)
      u[i] = NRLib::Random::Unif01();
  }

  std::vector<Fluid*> fluid_samples(n_fluids);

  for(size_t i = 0; i < n_fluids; ++i)
    fluid_samples[i] = distr_fluid_[i]->GenerateSample(trend_params);

  Fluid * fluid_mixed = GetSample(u, trend_params, fluid_samples);

  // Deep copy taken by constructor of FluidMixed, hence delete fluid here:
  for(size_t i = 0; i < n_fluids; ++i)
    delete fluid_samples[i];

  return fluid_mixed;
}

Fluid *
DistributionsFluidMix::GetSample(const std::vector<double>  & u,
                                 const std::vector<double>  & trend_params,
                                 const std::vector<Fluid *> & fluid_samples) const
{

  size_t n_fluids = fluid_samples.size();

  std::vector<double> volume_fraction(n_fluids, 0.0);

  size_t missing_index = n_fluids;

  for(size_t i = 0; i < n_fluids; ++i) {
    if(u[i] != RMISSING)
      volume_fraction[i] = distr_vol_frac_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    else
      missing_index = i;
  }

  if (missing_index != n_fluids) {
    double sum = 0.0;

    for (size_t i = 0; i < volume_fraction.size(); ++i)
      sum += volume_fraction[i];

    volume_fraction[missing_index] = 1.0 - sum;
  }

  Fluid * fluid_mixed = new FluidMix(fluid_samples, volume_fraction, u, mix_method_);

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

Fluid *
DistributionsFluidMix::UpdateSample(const std::vector< double > & /*corr*/,
                                    const Fluid                 & /*fluid*/) const {

  return NULL;
}
