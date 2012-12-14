#include "rplib/dryrockdem.h"

#include <cassert>
#include <algorithm>

DryRockDEM::DryRockDEM(const DryRock                       * dryrock,
                       const std::vector< DryRock* >       & dryrock_inc,
                       const std::vector<double>           & inclusion_spectrum,
                       const std::vector<double>           & inclusion_concentration,
                       const std::vector<double>           & u)
: DryRock()
{
  u_ = u; // u contains independent samples used in quantiles of (inclusion_spectrum,inclusion_concentration)

  // Deep copy of dryrock and dryrock_inc:
  dryrock_ = dryrock->Clone();
  Clone(dryrock_inc);

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;

  ComputeElasticParams();

}

DryRockDEM::DryRockDEM() : DryRock()
{
  rho_ = k_ = mu_ = 0;
}

DryRockDEM::~DryRockDEM()
{
  delete dryrock_;

  DeleteInclusion();
}

DryRockDEM& DryRockDEM::operator=(const DryRockDEM& rhs)
{
  if (this != &rhs) {
    DryRock::operator=(rhs);

    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;

    delete dryrock_;
    DeleteInclusion();

    dryrock_ = rhs.dryrock_->Clone();
    Clone(rhs.dryrock_inc_);

  }
  return *this;
}

DryRock *
DryRockDEM::Clone() const {
  // Provide base class variables.
  DryRockDEM * r = new DryRockDEM(*this);

  // Provide variables specific to DryRockDEM.
  r->dryrock_                   = this->dryrock_->Clone();          // Deep copy.
  r->Clone(this->dryrock_inc_);                                   // Deep copy

  r->inclusion_spectrum_      = this->inclusion_spectrum_;
  r->inclusion_concentration_ = this->inclusion_concentration_;

  return r;
}

DryRock *
DryRockDEM::Evolve(const std::vector<int>                & /*delta_time*/,
                   const std::vector< const DryRock * >  & /*dryrocks*/) const
{
  return new DryRockDEM(*this);

  //old code, worked when dryrock_inc_ was a single DryRock, not working now. Evolve when model is defined.
  //size_t n_dryrocks = dryrocks.size();
  //std::vector< const DryRockDEM * > dryrock_incl(n_dryrocks);
  //std::vector< const DryRock * > dryrock(n_dryrocks);
  //std::vector< const DryRock * > dryrock_inc(n_dryrocks);
  //for (size_t i = 0; i < n_dryrocks; ++i) {
  //  dryrock_incl[i] = dynamic_cast<const DryRockDEM*>(dryrocks[i]);
  //  assert(dryrock_incl[i] != NULL);
  //  dryrock[i]      = dryrock_incl[i]->GetDryRockHost();
  //  dryrock_inc[i]  = dryrock_incl[i]->GetDryRockInclusion();
  //}
  //DryRock * dryrock_new     = dryrock_->Evolve(delta_time, dryrock);
  //DryRock * dryrock_inc_new = dryrock_inc_->Evolve(delta_time, dryrock_inc);

  //// Change the assignment of the following three variables when a time develop model has been defined.
  //std::vector<double> inclusion_spectrum      = inclusion_spectrum_;
  //std::vector<double> inclusion_concentration = inclusion_concentration_;

  //DryRock * dryrock_new2 = new DryRockDEM(dryrock_new,
  //                                  dryrock_inc_new,
  //                                  inclusion_spectrum,
  //                                  inclusion_concentration,
  //                                  u_);

  //// Deep copy taken by constructor of DryRockDEM, hence delete
  //// dryrock_new and dryrock_inc_new here:
  //delete dryrock_new;
  //delete dryrock_inc_new;

  //return dryrock_new2;
}

void
DryRockDEM::ComputeElasticParams() {

  std::vector<double> dryrock_inc_rho, dryrock_inc_k, dryrock_inc_mu;

  dryrock_inc_rho.resize(dryrock_inc_.size());
  dryrock_inc_k.resize(dryrock_inc_.size());
  dryrock_inc_mu.resize(dryrock_inc_.size());

  for (size_t i = 0; i < dryrock_inc_.size(); ++i)
    dryrock_inc_[i]->GetElasticParams(dryrock_inc_k[i], dryrock_inc_mu[i], dryrock_inc_rho[i]); //NBNB fjellvoll marit what to to here?

  double dryrock_rho, dryrock_k, dryrock_mu;
  dryrock_->GetElasticParams(dryrock_k, dryrock_mu, dryrock_rho);

  { //Calculation of effective density if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> dryrock_host_inc_rho(dryrock_inc_.size() + 1);
    dryrock_host_inc_rho[0] = dryrock_rho;
    std::copy(dryrock_inc_rho.begin(), dryrock_inc_rho.end(), dryrock_host_inc_rho.begin() + 1);
    rho_  = DEMTools::CalcEffectiveDensity(dryrock_host_inc_rho, inclusion_concentration_);
  }

  //remove host from inclusion vector
  std::vector<double> inc_conc(dryrock_inc_.size());
  std::copy(inclusion_concentration_.begin()+1, inclusion_concentration_.end(), inc_conc.begin());

  DEMTools::CalcEffectiveBulkAndShearModulus(dryrock_inc_k,
                                             dryrock_inc_mu,
                                             inclusion_spectrum_,
                                             inc_conc,
                                             dryrock_k,
                                             dryrock_mu,
                                             k_,
                                             mu_);
}


void
DryRockDEM::Clone(const std::vector< DryRock* > & dryrock_in)
{
  dryrock_inc_.resize(dryrock_in.size());
  for (size_t i = 0; i < dryrock_in.size(); ++i)
    dryrock_inc_[i] = dryrock_in[i]->Clone();

}

void
DryRockDEM::DeleteInclusion()
{
  for (size_t i = 0; i < dryrock_inc_.size(); ++i)
    delete dryrock_inc_[i];

}

