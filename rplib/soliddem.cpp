#include "rplib/soliddem.h"

#include <cassert>

SolidDEM::SolidDEM(const Solid                         * solid,
                   const Solid                         * solid_inc,
                   const std::vector<double>           & inclusion_spectrum,
                   const std::vector<double>           & inclusion_concentration,
                   const std::vector<double>           & u)
: Solid()
{
  u_ = u; // u contains independent samples used in quantiles of (inclusion_spectrum,inclusion_concentration)

  // Deep copy of solid and solid_inc:
  solid_ = solid->Clone();
  solid_inc_ = solid_inc->Clone();

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;

  ComputeElasticParams();

}

SolidDEM::SolidDEM() : Solid()
{
  rho_ = k_ = mu_ = 0;
}

SolidDEM::~SolidDEM()
{
  delete solid_;
  delete solid_inc_;
}

SolidDEM& SolidDEM::operator=(const SolidDEM& rhs)
{
  if (this != &rhs) {
    Solid::operator=(rhs);

    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;

    delete solid_;
    delete solid_inc_;
    solid_ = rhs.solid_->Clone();
    solid_inc_ = rhs.solid_inc_->Clone();

  }
  return *this;
}

Solid *
SolidDEM::Clone() const {
  // Provide base class variables.
  SolidDEM * r = new SolidDEM(*this);

  // Provide variables specific to SolidDEM.
  r->solid_                   = this->solid_->Clone();          // Deep copy.
  r->solid_inc_                   = this->solid_inc_->Clone();          // Deep copy.
  r->inclusion_spectrum_      = this->inclusion_spectrum_;
  r->inclusion_concentration_ = this->inclusion_concentration_;

  return r;
}

Solid *
SolidDEM::Evolve(const std::vector<int>                & delta_time,
                 const std::vector< const Solid * >    & solids) const {

  size_t n_solids = solids.size();
  std::vector< const SolidDEM * > solid_incl(n_solids);
  std::vector< const Solid * > solid(n_solids);
  std::vector< const Solid * > solid_inc(n_solids);
  for (size_t i = 0; i < n_solids; ++i) {
    solid_incl[i] = dynamic_cast<const SolidDEM*>(solids[i]);
    assert(solid_incl[i] != NULL);
    solid[i]      = solid_incl[i]->GetSolidHost();
    solid_inc[i]  = solid_incl[i]->GetSolidInclusion();
  }
  Solid * solid_new     = solid_->Evolve(delta_time, solid);
  Solid * solid_inc_new = solid_inc_->Evolve(delta_time, solid_inc);

  // Change the assignment of the following three variables when a time develop model has been defined.
  std::vector<double> inclusion_spectrum      = inclusion_spectrum_;
  std::vector<double> inclusion_concentration = inclusion_concentration_;

  Solid * solid_new2 = new SolidDEM(solid_new,
                                    solid_inc_new,
                                    inclusion_spectrum,
                                    inclusion_concentration,
                                    u_);

  // Deep copy taken by constructor of SolidDEM, hence delete
  // solid_new and solid_inc_new here:
  delete solid_new;
  delete solid_inc_new;

  return solid_new2;
}

void
SolidDEM::ComputeElasticParams() {

  double solid_inc_rho, solid_inc_k, solid_inc_mu;
  solid_inc_->GetElasticParams(solid_inc_k, solid_inc_mu, solid_inc_rho);

  double solid_rho, solid_k, solid_mu;
  solid_->GetElasticParams(solid_k, solid_mu, solid_rho);

  { // calculation of effective density
    std::vector<double> rho;
    rho.push_back(solid_inc_rho);
    rho.push_back(solid_rho);

    std::vector<double> volume_fraction = inclusion_concentration_;

    double host_vol_frac = 0;
    for (size_t i = 0; i < inclusion_concentration_.size(); ++i)
      host_vol_frac += inclusion_concentration_[i];

    host_vol_frac = 1.0 - host_vol_frac;

    volume_fraction.push_back(host_vol_frac);
    rho_  = DEMTools::CalcEffectiveDensity(rho, volume_fraction);
  }

  std::vector<double> inclusion_k   =  std::vector<double>(inclusion_spectrum_.size(), solid_inc_k);
  std::vector<double> inclusion_mu  = std::vector<double>(inclusion_spectrum_.size(), solid_inc_mu);

  DEMTools::CalcEffectiveBulkAndShearModulus(inclusion_k,
                                             inclusion_mu,
                                             inclusion_spectrum_,
                                             inclusion_concentration_,
                                             solid_k,
                                             solid_mu,
                                             k_,
                                             mu_);
}

