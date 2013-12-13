#include "rplib/soliddem.h"

#include <cassert>
#include <algorithm>

SolidDEM::SolidDEM(const Solid                         * solid,
                   const std::vector< Solid* >         & solid_inc,
                   const std::vector<double>           & inclusion_spectrum,
                   const std::vector<double>           & inclusion_concentration,
                   const std::vector<double>           & u)
: Solid()
{
  u_ = u; // u contains independent samples used in quantiles of (inclusion_spectrum,inclusion_concentration)

  // Deep copy of solid and solid_inc:
  solid_ = solid->Clone();
  Clone(solid_inc);

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;

  ComputeElasticParams();

}

SolidDEM::SolidDEM()
: Solid(),
  solid_inc_(0)
{
  solid_ = NULL;
}

SolidDEM::~SolidDEM()
{
  delete solid_;

  DeleteInclusion();
}

SolidDEM& SolidDEM::operator=(const SolidDEM& rhs)
{
  if (this != &rhs) {
    Solid::operator=(rhs);

    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;

    delete solid_;
    DeleteInclusion();

    solid_ = rhs.solid_->Clone();
    Clone(rhs.solid_inc_);

  }
  return *this;
}

Solid *
SolidDEM::Clone() const {
  // Provide base class variables.
  SolidDEM * r = new SolidDEM(*this);

  // Provide variables specific to SolidDEM.
  r->solid_                   = this->solid_->Clone();          // Deep copy.
  r->Clone(this->solid_inc_);                                   // Deep copy

  r->inclusion_spectrum_      = this->inclusion_spectrum_;
  r->inclusion_concentration_ = this->inclusion_concentration_;

  return r;
}

void
SolidDEM::ComputeElasticParams() {

  std::vector<double> solid_inc_rho, solid_inc_k, solid_inc_mu;

  solid_inc_rho.resize(solid_inc_.size());
  solid_inc_k.resize(solid_inc_.size());
  solid_inc_mu.resize(solid_inc_.size());

  for (size_t i = 0; i < solid_inc_.size(); ++i)
    solid_inc_[i]->GetElasticParams(solid_inc_k[i], solid_inc_mu[i], solid_inc_rho[i]);

  double solid_rho, solid_k, solid_mu;
  solid_->GetElasticParams(solid_k, solid_mu, solid_rho);

  { //Calculation of effective density if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> solid_host_inc_rho(solid_inc_.size() + 1);
    solid_host_inc_rho[0] = solid_rho;
    std::copy(solid_inc_rho.begin(), solid_inc_rho.end(), solid_host_inc_rho.begin() + 1);
    rho_  = DEMTools::CalcEffectiveDensity(solid_host_inc_rho, inclusion_concentration_);
  }

  //remove host from inclusion vector
  std::vector<double> inc_conc(solid_inc_.size());
  std::copy(inclusion_concentration_.begin()+1, inclusion_concentration_.end(), inc_conc.begin());

  DEMTools::CalcEffectiveBulkAndShearModulus(solid_inc_k,
                                             solid_inc_mu,
                                             inclusion_spectrum_,
                                             inc_conc,
                                             solid_k,
                                             solid_mu,
                                             k_,
                                             mu_);
}


void
SolidDEM::Clone(const std::vector< Solid* > & solid_in)
{
  solid_inc_.resize(solid_in.size());
  for (size_t i = 0; i < solid_in.size(); ++i)
    solid_inc_[i] = solid_in[i]->Clone();

}

void
SolidDEM::DeleteInclusion()
{
  for (size_t i = 0; i < solid_inc_.size(); ++i)
    delete solid_inc_[i];

}

