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

  virtual                               ~RockInclusion();

  // Assignment operator.
  RockInclusion                       & operator=(const RockInclusion& rhs);

  virtual Rock                        * Clone()                                                 const;

  void                                  GetElasticParams(double & k, double & mu, double & rho) const;

  const Solid                         * GetSolid()                                              const {return solid_;}
  const Fluid                         * GetFluid()                                              const {return fluid_;}

  virtual Rock                        * Evolve(const std::vector<int>         & delta_time,
                                               const std::vector< Rock * >    & rock)           const;

  virtual void                          SetPorosity(double porosity);

private:
                                        //Copy constructor for getting base class variables , used by Clone:
                                        RockInclusion(const RockInclusion & rhs) : Rock(rhs) {}

                                        // Calculate elastic and seismic parameters, to be
                                        // used whenever new information is sent to class.
  void                                  ComputeSeismicAndElasticParams();

  Solid                               * solid_; // Owned and deleted by this class.
  Fluid                               * fluid_; // Owned and deleted by this class.
  std::vector<double>                   inclusion_spectrum_;
  std::vector<double>                   inclusion_concentration_;
  double                                porosity_;
  DistributionsRockInclusionEvolution * distr_evolution_; // Pointer to external object.

  double k_, mu_;
};

#endif
