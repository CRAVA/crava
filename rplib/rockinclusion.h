#ifndef RPLIB_ROCKINCLUSION_H
#define RPLIB_ROCKINCLUSION_H

#include "rplib/rock.h"
#include "rplib/solidmixed.h"
#include "rplib/fluidmixed.h"
#include "rplib/distributionsrockinclusionevolution.h"
#include "rplib/demmodelling.h"

#include <vector>

class RockInclusion : public Rock {
public:

  RockInclusion(const Solid                         * solid,
                const Fluid                         * fluid,
                const std::vector<double>           & inclusion_spectrum,
                const std::vector<double>           & inclusion_concentration,
                double                                porosity,
                DistributionsRockInclusionEvolution * distr_evolution = NULL);

  RockInclusion();

  virtual ~RockInclusion();

  virtual Rock * Clone() const;

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const;

  void GetElasticParams(double & k, double & mu, double & rho) const;

  const Solid *  GetSolid() const {return solid_;}
  const Fluid *  GetFluid() const {return fluid_;}

  virtual Rock * Evolve(const std::vector<int>         & delta_time,
                        const std::vector< Rock * >    & rock) const;

private:
  //Copy constructor for getting base class variables , used by Clone:
  RockInclusion(const RockInclusion & rhs) : Rock(rhs) {}

  Solid                               * solid_; // Owned and deleted by this class.
  Fluid                               * fluid_; // Owned and deleted by this class.
  std::vector<double>                   inclusion_spectrum_;
  std::vector<double>                   inclusion_concentration_;
  double                                porosity_;
  DistributionsRockInclusionEvolution * distr_evolution_; // Pointer to external object.

  double vp_, vs_, rho_;
  double k_, mu_;
};

#endif
