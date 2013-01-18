#include "rplib/rockdem.h"

#include "rplib/solid.h"
#include "rplib/fluid.h"

#include <cassert>

RockDEM::RockDEM(const Solid                         * solid,
                 const std::vector<Fluid*>           & fluid,
                 const std::vector<double>           & inclusion_spectrum,
                 const std::vector<double>           & inclusion_concentration,
                 const std::vector<double>           & u)
: Rock()
{
  u_ = u; // u contains independent samples used in quantiles of (inclusion_spectrum,inclusion_concentration), u.back() is for porosity.

  // Deep copy of solid and fluid:
  solid_ = solid->Clone();
  Clone(fluid);

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;

  ComputeSeismicAndElasticParams();

}

RockDEM::RockDEM() : Rock()
{
   k_ = mu_ = 0;
}

RockDEM::~RockDEM()
{
  delete solid_;
  DeleteInclusion();
}

RockDEM& RockDEM::operator=(const RockDEM& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    k_                       = rhs.k_;
    mu_                      = rhs.mu_;
    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;

    delete solid_;
    DeleteInclusion();
    solid_ = rhs.solid_->Clone();
    Clone(rhs.fluid_);

  }
  return *this;
}

Rock *
RockDEM::Clone() const {
  // Provide base class variables.
  RockDEM * r = new RockDEM(*this);

  // Provide variables specific to RockDEM.
  r->solid_                   = this->solid_->Clone();          // Deep copy.
  r->Clone(this->fluid_);                                       // Deep copy.
  r->inclusion_spectrum_      = this->inclusion_spectrum_;
  r->inclusion_concentration_ = this->inclusion_concentration_;
  r->k_                       = this->k_;
  r->mu_                      = this->mu_;

  return r;
}

void
RockDEM::GetElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

void
RockDEM::SetPorosity(double porosity) {
  // the following assumes inclusions are fluids
  double poro_old = 1.0 - inclusion_concentration_[0];
  double new_volume_fraction_fluid = porosity;
  double new_volume_fraction_solid = 1.0 - new_volume_fraction_fluid;

  inclusion_concentration_[0] = new_volume_fraction_solid;

  for (size_t i = 1; i < inclusion_concentration_.size(); ++i)
    inclusion_concentration_[i] = (porosity * inclusion_concentration_[i]) / poro_old;

  ComputeSeismicAndElasticParams();
}

void
RockDEM::ComputeSeismicAndElasticParams() {
  std::vector<double> fluid_rho, fluid_k, fluid_mu;

  fluid_rho.resize(fluid_.size());
  fluid_k.resize(fluid_.size());
  fluid_mu.resize(fluid_.size(), 0.0);

  for (size_t i = 0; i < fluid_.size(); ++i)
    fluid_[i]->GetElasticParams(fluid_k[i], fluid_rho[i]);

  double solid_rho, solid_k, solid_mu;
  solid_->GetElasticParams(solid_k, solid_mu, solid_rho);

  { //Calculation of effective density if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> fluid_host_inc_rho(fluid_.size() + 1);
    fluid_host_inc_rho[0] = solid_rho;
    std::copy(fluid_rho.begin(), fluid_rho.end(), fluid_host_inc_rho.begin() + 1);
    rho_  = DEMTools::CalcEffectiveDensity(fluid_host_inc_rho, inclusion_concentration_);
  }

  //remove host from inclusion vector
  std::vector<double> inc_conc(fluid_.size());
  std::copy(inclusion_concentration_.begin() + 1, inclusion_concentration_.end(), inc_conc.begin());

  DEMTools::CalcEffectiveBulkAndShearModulus(fluid_k,
                                             fluid_mu,
                                             inclusion_spectrum_,
                                             inc_conc,
                                             solid_k,
                                             solid_mu,
                                             k_,
                                             mu_);

  DEMTools::CalcSeismicParamsFromElasticParams(k_, mu_, rho_, vp_, vs_);
}

void
RockDEM::Clone(const std::vector< Fluid* > & fluid_in)
{
  fluid_.resize(fluid_in.size());
  for (size_t i = 0; i < fluid_in.size(); ++i)
    fluid_[i] = fluid_in[i]->Clone();

}

void
RockDEM::DeleteInclusion()
{
  for (size_t i = 0; i < fluid_.size(); ++i)
    delete fluid_[i];

}

