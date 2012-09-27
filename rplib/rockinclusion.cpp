#include "rplib/rockinclusion.h"

RockInclusion::RockInclusion(const Solid                         * solid,
                             const Fluid                         * fluid,
                             const std::vector<double>           & inclusion_spectrum,
                             const std::vector<double>           & inclusion_concentration,
                             double                                porosity,
                             DistributionsRockInclusionEvolution * distr_evolution)
: Rock()
{
  // Deep copy of solid and fluid:
  solid_ = solid->Clone();
  fluid_ = fluid->Clone();

  inclusion_spectrum_      = inclusion_spectrum;
  inclusion_concentration_ = inclusion_concentration;
  porosity_                = porosity;
  distr_evolution_         = distr_evolution;

  ComputeSeismicAndElasticParams();

}

RockInclusion::RockInclusion() : Rock()
{
   k_ = mu_ = 0;
}

RockInclusion::~RockInclusion()
{
  delete solid_;
  delete fluid_;
}

RockInclusion& RockInclusion::operator=(const RockInclusion& rhs)
{
  if (this != &rhs) {
    Rock::operator=(rhs);

    k_                       = rhs.k_;
    mu_                      = rhs.mu_;
    inclusion_spectrum_      = rhs.inclusion_spectrum_;
    inclusion_concentration_ = rhs.inclusion_concentration_;
    porosity_                = rhs.porosity_;
    distr_evolution_         = rhs.distr_evolution_;

    delete solid_;
    delete fluid_;
    solid_ = rhs.solid_->Clone();
    fluid_ = rhs.fluid_->Clone();

  }
  return *this;
}

Rock *
RockInclusion::Clone() const {
  // Provide base class variables.
  RockInclusion * r = new RockInclusion(*this);

  // Provide variables specific to RockInclusion.
  r->solid_                   = this->solid_->Clone();          // Deep copy.
  r->fluid_                   = this->fluid_->Clone();          // Deep copy.
  r->inclusion_spectrum_      = this->inclusion_spectrum_;
  r->inclusion_concentration_ = this->inclusion_concentration_;
  r->porosity_                = this->porosity_;
  r->distr_evolution_         = this->distr_evolution_;         // Not deep copy
  r->k_                       = this->k_;
  r->mu_                      = this->mu_;

  return r;
}

void
RockInclusion::GetElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

Rock *
RockInclusion::Evolve(const std::vector<int>         & delta_time,
                      const std::vector< Rock * >    & rock) const {

  size_t n_rocks = rock.size();
  std::vector< RockInclusion * > rock_incl(n_rocks);
  std::vector< const Solid * > solid(n_rocks);
  std::vector< const Fluid * > fluid(n_rocks);
  for (size_t i = 0; i < n_rocks; ++i) {
    rock_incl[i] = dynamic_cast<RockInclusion*>(rock[i]);
    assert(rock_incl[i] != NULL);
    solid[i] = rock_incl[i]->GetSolid();
    fluid[i] = rock_incl[i]->GetFluid();
  }
  Solid * solid_new = solid_->Evolve(delta_time, solid);
  Fluid * fluid_new = fluid_->Evolve(delta_time, fluid);

  // Change the assignment of the following three variables when a time develop model has been defined.
  std::vector<double> inclusion_spectrum      = inclusion_spectrum_;
  std::vector<double> inclusion_concentration = inclusion_concentration_;
  double  porosity                            = porosity_;

  Rock * rock_new = new RockInclusion(solid_new,
                                      fluid_new,
                                      inclusion_spectrum,
                                      inclusion_concentration,
                                      porosity,
                                      distr_evolution_);

  // Deep copy taken by constructor of RockInclusion, hence delete
  // solid_new and fluid_new here:
  delete solid_new;
  delete fluid_new;

  return rock_new;
}

void
RockInclusion::SetPorosity(double porosity) {
  porosity_ = porosity;
  ComputeSeismicAndElasticParams();
}

void
RockInclusion::ComputeSeismicAndElasticParams() {

  double fluid_rho, fluid_k;
  fluid_->GetElasticParams(fluid_k, fluid_rho);

  double solid_rho, solid_k, solid_mu;
  solid_->GetElasticParams(solid_k, solid_mu, solid_rho);

  std::vector<double> rho;
  rho.push_back(fluid_rho);
  rho.push_back(solid_rho);
  std::vector<double> poro;
  poro.push_back(porosity_);
  poro.push_back(1.0 - porosity_);
  rho_  = DEMTools::CalcEffectiveDensity(rho, poro);

  std::vector<double> inclusion_k   =  std::vector<double>(inclusion_spectrum_.size(), fluid_k);
  std::vector<double> inclusion_mu  = std::vector<double>(inclusion_spectrum_.size(), 0.0);
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
  DEMTools::CalcSeismicParamsFromElasticParams(k_, mu_, rho_, vp_, vs_);
}

