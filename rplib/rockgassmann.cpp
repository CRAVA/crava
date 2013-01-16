#include "rplib/rockgassmann.h"

#include "rplib/dryrock.h"
#include "rplib/fluid.h"

#include <cassert>

RockGassmann::RockGassmann(const Fluid                * fluid,
                           const DryRock              * dryrock)
: Rock()
{
  // Deep copy of dryrock and fluid:
  fluid_   = fluid->Clone();
  dryrock_ = dryrock->Clone();

  ComputeSeismicAndElasticParams();

}

RockGassmann::RockGassmann() : Rock()
{
   k_ = mu_ = 0;
}

RockGassmann::~RockGassmann()
{
  delete dryrock_;
  delete fluid_;
}

RockGassmann& RockGassmann::operator=(const RockGassmann& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    k_                       = rhs.k_;
    mu_                      = rhs.mu_;

    delete dryrock_;
    delete fluid_;
    dryrock_ = rhs.dryrock_->Clone();
    fluid_   = rhs.fluid_->Clone();
  }
  return *this;
}

Rock *
RockGassmann::Clone() const {
  // Provide base class variables.
  RockGassmann * r = new RockGassmann(*this);

  // Provide variables specific to RockGassmann.
  r->dryrock_                 = this->dryrock_->Clone();          // Deep copy.
  r->fluid_                   = this->fluid_->Clone();          // Deep copy.
  r->k_                       = this->k_;
  r->mu_                      = this->mu_;

  return r;
}

void
RockGassmann::GetElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

Rock *
RockGassmann::Evolve(const std::vector<int>         & /*delta_time*/,
                     const std::vector< Rock * >    & /*rock*/) const {

  return new RockGassmann(*this); //function will be removed later
}

void
RockGassmann::SetPorosity(double porosity) {
  dryrock_->SetTotalPorosity(porosity);
  ComputeSeismicAndElasticParams();
}

void
RockGassmann::ComputeSeismicAndElasticParams() {
  double k_dry, mu_dry, rho_dry;
  dryrock_->GetElasticParams(k_dry, mu_dry, rho_dry);

  double k_fluid, rho_fluid;
  fluid_->GetElasticParams(k_fluid, rho_fluid);

  mu_ = mu_dry; // mu and porosity are assumed unchanged in Gassmann model.

  double total_porosity   = dryrock_->GetTotalPorosity();
  double mineral_moduli_k = dryrock_->GetMineralModuliK();

  double temp1 = (1 - k_dry/mineral_moduli_k)*(1 - k_dry/mineral_moduli_k);

  k_ = k_dry + temp1/(total_porosity/k_fluid + (1 - total_porosity)/mineral_moduli_k + k_dry/(mineral_moduli_k*mineral_moduli_k));

  rho_ = (1.0 - total_porosity)*rho_dry + total_porosity*rho_fluid;
  DEMTools::CalcSeismicParamsFromElasticParams(k_, mu_, rho_, vp_, vs_);
}

