#ifndef RPLIB_ROCK_GASSMANN_H
#define RPLIB_ROCK_GASSMANN_H

#include "rplib/rock.h"
#include "rplib/demmodelling.h"

#include <vector>

class DryRock;
class Fluid;

class RockGassmann : public Rock {
public:

RockGassmann(const Fluid                * fluid,
             const DryRock              * dryrock);

RockGassmann();

  virtual                               ~RockGassmann();

  // Assignment operator.
  RockGassmann                        & operator=(const RockGassmann& rhs);

  virtual Rock                        * Clone()                                                 const;

  void                                  GetElasticParams(double & k, double & mu, double & rho) const;

  const Fluid                         * GetFluid()                                              const { return fluid_; }
  const DryRock                       * GetDryRock()                                            const { return dryrock_; }

  virtual Rock                        * Evolve(const std::vector<int>         & delta_time,
                                               const std::vector< Rock * >    & rock)           const;

  virtual void                          SetPorosity(double porosity);

private:
                                        //Copy constructor for getting base class variables , used by Clone:
                                        RockGassmann(const RockGassmann & rhs) : Rock(rhs) {}

                                        // Calculate elastic and seismic parameters, to be
                                        // used whenever new information is sent to class.
  void                                  ComputeSeismicAndElasticParams();

  Fluid                               * fluid_; // Owned and deleted by this class.
  DryRock                             * dryrock_;

  double                                k_, mu_;
};

#endif
