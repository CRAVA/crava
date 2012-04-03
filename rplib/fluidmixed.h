#ifndef RPLIB_FLUIDMIXED_H
#define RPLIB_FLUIDMIXED_H

#include "rplib/fluid.h"
#include "rplib/demmodelling.h"

#include <cassert>
#include <vector>
#include <numeric>

class FluidMixed : public Fluid {
public:

  // Parallel classes are DistributionsFluidMixedEvolution and DistributionsFluidMixed.
  FluidMixed(const    std::vector<double>& k,
             const    std::vector<double>& rho,
             const    std::vector<double>& volume_fraction)
  : Fluid() {
    k_ = rho_ = 0;
    if (k.size() != volume_fraction.size() || rho.size() != volume_fraction.size())
      throw NRLib::Exception("Invalid arguments:Number of properties are different from number of volume fractions.");
    else if (std::accumulate(volume_fraction.begin(), volume_fraction.end(), 0.0) > 1.0) //NBNB fjellvoll possible to give warning to user and then rescale
      throw NRLib::Exception("Invalid arguments:Sum of volume fractions > 1.0");
    else if (k.size() != 2) //NBNB fjellvoll this is only tmp.
      throw NRLib::Exception("Invalid arguments:Number of volume fractions is different than 2.");
    else {
      //NBNB fjellvoll only pair-mix allowed at the moment
      k_    = DEMTools::CalcEffectiveElasticModuliUsingReuss(k[0], volume_fraction[0], k[1]);     // homogeneous
      rho_  = DEMTools::CalcEffectiveDensity(rho[0], volume_fraction[0], rho[1]);
    }
  }

  virtual ~FluidMixed(){}

  virtual void ComputeElasticParams(const double   /*temp*/,
                                    const double   /*pore_pressure*/,
                                    double       & k,
                                    double       & rho) const {
    k = k_;
    rho = rho_;
  }

  virtual Fluid * Evolve(const std::vector<int>             & /*delta_time*/,
                         const std::vector< Fluid * >       & /*fluid*/,
                         const DistributionsFluidEvolution  * /*dist_fluid_evolve*/) const {
    //const DistributionsFluidMixedEvolution * dist_fluid_mix_evolve = dynamic_cast<const DistributionsFluidMixedEvolution*>(dist_fluid_evolve);
    //assert(dist_c02_evolve != NULL);
    //assert(delta_time.size() == fluid.size() + 1);

    //// Temporary implementation that simply makes a copy, but illustrates the possible use of dist_fluid_mix_evolve:
    //double example_param = dist_fluid_mix_evolve->Sample();
    //example_param += 1.0; //FAKE to make avoid compiler warnings
    std::vector<double> k(2,0);
    std::vector<double> rho(2,0);
    std::vector<double> volume_fraction(2,0);
    Fluid * new_fluid = new FluidMixed(k, rho, volume_fraction);
    return new_fluid;
  }

private:
  double k_;
  double rho_;

};

#endif
