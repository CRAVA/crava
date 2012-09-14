#include "rplib/soliddem.h"

#include "rplib/solidmixed.h"

#include <cassert>

SolidDEM::SolidDEM(const Solid                         * solid,
                   const Solid                         * solid_inc,
                   const std::vector<double>           & inclusion_spectrum,
                   const std::vector<double>           & inclusion_concentration,
                   double                                porosity,
                   DistributionsSolidDEMEvolution      * distr_evolution)
: Solid()
{
  // Deep copy of solid and solid_inc:
  solid_ = solid->Clone();
  solid_inc_ = solid_inc->Clone();

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;
  porosity_                = porosity;
  distr_evolution_         = distr_evolution;

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

    k_                       = rhs.k_;
    mu_                      = rhs.mu_;
    rho_                     = rhs.rho_;
    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;
    porosity_                = rhs.porosity_;
    distr_evolution_         = rhs.distr_evolution_;

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
  r->porosity_                = this->porosity_;
  r->distr_evolution_         = this->distr_evolution_;         // Not deep copy.
  r->rho_                     = this->rho_;
  r->k_                       = this->k_;
  r->mu_                      = this->mu_;

  return r;
}

void
SolidDEM::ComputeElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
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
  double  porosity                            = porosity_;

  Solid * solid_new2 = new SolidDEM(solid_new,
                                    solid_inc_new,
                                    inclusion_spectrum,
                                    inclusion_concentration,
                                    porosity,
                                    distr_evolution_);

  // Deep copy taken by constructor of SolidDEM, hence delete
  // solid_new and solid_inc_new here:
  delete solid_new;
  delete solid_inc_new;

  return solid_new2;
}

void
SolidDEM::SetPorosity(double porosity) {
  porosity_ = porosity;
  ComputeElasticParams();
}

void
SolidDEM::ComputeElasticParams() {

  double solid_inc_rho, solid_inc_k, solid_inc_mu;
  solid_inc_->ComputeElasticParams(solid_inc_k, solid_inc_mu, solid_inc_rho);

  double solid_rho, solid_k, solid_mu;
  solid_->ComputeElasticParams(solid_k, solid_mu, solid_rho);

  std::vector<double> rho;
  rho.push_back(solid_inc_rho);
  rho.push_back(solid_rho);
  std::vector<double> poro;
  poro.push_back(porosity_);
  poro.push_back(1.0 - porosity_);
  rho_  = DEMTools::CalcEffectiveDensity(rho, poro);

  std::vector<double> inclusion_k   =  std::vector<double>(inclusion_spectrum_.size(), solid_inc_k);
  std::vector<double> inclusion_mu  = std::vector<double>(inclusion_spectrum_.size(), solid_inc_mu);
  std::vector<double> conc = inclusion_concentration_; // inclusion concentration scaled by porosity
  for (size_t i = 0; i < conc.size(); i++)
   conc[i] *= porosity_;

  DEMTools::CalcEffectiveBulkAndShearModulus(inclusion_k,
                                             inclusion_mu,
                                             inclusion_spectrum_,
                                             conc,
                                             solid_k,
                                             solid_mu,
                                             k_,
                                             mu_);
}

