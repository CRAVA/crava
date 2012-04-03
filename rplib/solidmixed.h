#ifndef RPLIB_SOLIDMIXED_H
#define RPLIB_SOLIDMIXED_H

#include <vector>
#include <numeric>

#include "rplib/solid.h"
#include "rplib/mineral.h"
#include "rplib/demmodelling.h"

#include "nrlib/exception/exception.hpp"

class SolidMixed : public Solid {
public:

  SolidMixed(const std::vector<Solid*>& solid,
             const std::vector<double>&   volume_fraction)
  : Solid() {

    k_   = mu_ = rho_ = 0;

    if (solid.size() != volume_fraction.size())
      throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
    else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
      throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
    else if (solid.size() != 2) //NBNB fjellvoll this is only tmp.
      throw NRLib::Exception("Invalid arguments:Number of volume fractions is different than 2.");
    else {
      std::vector<double> k(solid.size()), mu(solid.size()), rho(solid.size());
      for (size_t i = 0; i < solid.size(); i++)
        solid[i]->ComputeElasticParams(k[i], mu[i], rho[i]);

      //NBNB fjellvoll only pair-mix allowed at the moment
      k_      = DEMTools::CalcEffectiveElasticModuliUsingHill(k[0], volume_fraction[0], k[1]);
      mu_     = DEMTools::CalcEffectiveElasticModuliUsingHill(mu[0], volume_fraction[0], mu[1]);
      rho_    = DEMTools::CalcEffectiveDensity(rho[0], volume_fraction[0], rho[1]);
    }
  }

  virtual ~SolidMixed(){}

  virtual void ComputeElasticParams(double & k, double & mu, double & rho) const {
    k   = k_;
    mu  = mu_;
    rho = rho_;
  }

  virtual Solid * Evolve(const std::vector<int>             & /*delta_time*/,
                         const std::vector< Solid * >       & /*fluid*/,
                         const DistributionsSolidEvolution  * /*dist_fluid_evolve*/) const {
    return new SolidMixed(*this);
  }

private:
  double k_;
  double mu_;
  double rho_;
};

#endif
