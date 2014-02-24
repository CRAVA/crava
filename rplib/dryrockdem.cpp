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

DryRockDEM::DryRockDEM()
: DryRock(),
  dryrock_inc_(0)
{
  dryrock_ = NULL;
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

void
DryRockDEM::ComputeElasticParams() {

  std::vector<double> dryrock_inc_rho, dryrock_inc_k, dryrock_inc_mu;
  std::vector<double> total_porosity(dryrock_inc_.size()), mineral_moduli_k(dryrock_inc_.size());

  dryrock_inc_rho.resize(dryrock_inc_.size());
  dryrock_inc_k.resize(dryrock_inc_.size());
  dryrock_inc_mu.resize(dryrock_inc_.size());


  for (size_t i = 0; i < dryrock_inc_.size(); ++i) {
    dryrock_inc_[i]->GetElasticParams(dryrock_inc_k[i], dryrock_inc_mu[i], dryrock_inc_rho[i]);
    total_porosity[i]   = dryrock_inc_[i]->GetTotalPorosity();
    mineral_moduli_k[i] = dryrock_inc_[i]->GetMineralModuliK();
  }

  double dryrock_rho, dryrock_k, dryrock_mu, dryrock_poro, dryrock_mineral_moduli_k;
  dryrock_->GetElasticParams(dryrock_k, dryrock_mu, dryrock_rho);
  dryrock_poro             = dryrock_->GetTotalPorosity();
  dryrock_mineral_moduli_k = dryrock_->GetMineralModuliK();

  { //Calculation of effective density if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> dryrock_host_inc_rho(dryrock_inc_.size() + 1);
    dryrock_host_inc_rho[0] = dryrock_rho;
    std::copy(dryrock_inc_rho.begin(), dryrock_inc_rho.end(), dryrock_host_inc_rho.begin() + 1);
    rho_  = DEMTools::CalcEffectiveDensity(dryrock_host_inc_rho, inclusion_concentration_);
  }

  { //Calculation of effective total porosity if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> dryrock_host_inc_poro(dryrock_inc_.size() + 1);
    dryrock_host_inc_poro[0] = dryrock_poro;
    std::copy(total_porosity.begin(), total_porosity.end(), dryrock_host_inc_poro.begin() + 1);
    total_porosity_  = DEMTools::CalcEffectivePorosity(dryrock_host_inc_poro, inclusion_concentration_);
  }

  { //Calculation of effective mineral moduli k if inclusion_concentration_ is filled with inclusions AND host.
    std::vector<double> dryrock_host_inc_moduli_k(dryrock_inc_.size() + 1);
    dryrock_host_inc_moduli_k[0] = dryrock_mineral_moduli_k;
    std::copy(mineral_moduli_k.begin(), mineral_moduli_k.end(), dryrock_host_inc_moduli_k.begin() + 1);
    mineral_moduli_k_  = DEMTools::CalcEffectiveElasticModuliUsingHill(dryrock_host_inc_moduli_k, inclusion_concentration_);
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

