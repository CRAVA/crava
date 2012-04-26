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

  double fluid_rho, fluid_k;
  fluid_->GetElasticParams(fluid_k, fluid_rho);

  double solid_rho, solid_k, solid_mu;
  solid_->ComputeElasticParams(solid_k, solid_mu, solid_rho);

  rho_  = DEMTools::CalcEffectiveDensity(fluid_rho, porosity_, solid_rho);

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

RockInclusion::RockInclusion() : Rock()
{
  vp_ = vs_ = rho_ = k_ = mu_ = 0;
}

RockInclusion::~RockInclusion()
{
  delete solid_;
  delete fluid_;
}

void
RockInclusion::ComputeSeismicParams(double & vp, double & vs, double & rho) const
{
  // Computation done in constructor.
  vp  = vp_;
  vs  = vs_;
  rho = rho_;
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
  std::vector< Solid * > solid(n_rocks);
  std::vector< Fluid * > fluid(n_rocks);
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

