#ifndef RPLIB_ROCKINCLUSION_H
#define RPLIB_ROCKINCLUSION_H

#include "rplib/rock.h"
#include "rplib/distributionsgeochemicaldem.h"
#include "rplib/solidmixed.h"
#include "rplib/fluidmixed.h"
#include "rplib/demmodelling.h"

#include <vector>

class RockInclusion : public Rock {
public:

  // Parallel classes are DistributionsGeochemicalRockInclusion and DistributionsRockInclusion.
  RockInclusion(const SolidMixed& solid_mix,
                const FluidMixed& fluid_mix,
                const std::vector<double>& bulk_modulus,
                const std::vector<double>& shear_modulus,
                const std::vector<double>& aspect_ratio,
                const std::vector<double>& concentration,
                double porosity)
  : Rock()
  {
    double fluid_rho, fluid_k;
    fluid_mix.GetElasticParams(fluid_k, fluid_rho);

    double solid_rho, solid_k, solid_mu;
    solid_mix.ComputeElasticParams(solid_k, solid_mu, solid_rho);

    rho_  = DEMTools::CalcEffectiveDensity(fluid_rho, porosity, solid_rho);

    std::vector<double> conc = concentration;

    DEMTools::CalcEffectiveBulkAndShearModulus(bulk_modulus,
                                               shear_modulus,
                                               aspect_ratio,
                                               conc,
                                               solid_k,
                                               solid_mu,
                                               k_,
                                               mu_);

  }

  RockInclusion() : Rock() {
    vp_ = vs_ = rho_ = k_ = mu_ = 0;
  }

  virtual ~RockInclusion(){}

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const {
    vp  = vp_;
    vs  = vs_;
    rho = rho_;
  }

  void GetElasticParams(double & k, double & mu, double & rho) const {
    k   = k_;
    mu  = mu_;
    rho = rho_;
  }

  virtual Rock * Evolve(const std::vector<int>         & /*delta_time*/,
                        const std::vector< Rock * >    & /*rock*/,
                        const DistributionsSaturation  * /*dist_sat*/,
                        const DistributionsGeochemical * /*dist_geochem*/) const {

    Rock * new_rock = new RockInclusion;

    return new_rock;
  }

private:
  double vp_, vs_, rho_;
  double k_, mu_;
};

#endif
