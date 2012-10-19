#ifndef RPLIB_ROCK_DEM_H
#define RPLIB_ROCK_DEM_H

#include "rplib/rock.h"
#include "rplib/demmodelling.h"

#include <vector>

class Solid;
class Fluid;

class RockDEM : public Rock {
public:

RockDEM(const Solid                         * solid,
        const std::vector<Fluid*>           & fluid,
        const std::vector<double>           & inclusion_spectrum,
        const std::vector<double>           & inclusion_concentration,
        const std::vector<double>           & u);

RockDEM();

  virtual                               ~RockDEM();

  // Assignment operator.
  RockDEM                             & operator=(const RockDEM& rhs);

  virtual Rock                        * Clone()                                                 const;

  void                                  GetElasticParams(double & k, double & mu, double & rho) const;

  const Solid                         * GetSolid()                                              const {return solid_;}
  const std::vector<Fluid*>           & GetFluid()                                              const {return fluid_;}
  const Fluid                         * GetFluid(size_t i)                                      const {return fluid_[i];} // no error checking on valid index range

  virtual Rock                        * Evolve(const std::vector<int>         & delta_time,
                                               const std::vector< Rock * >    & rock)           const;

  virtual void                          SetPorosity(double porosity);

private:
                                        //Copy constructor for getting base class variables , used by Clone:
                                        RockDEM(const RockDEM & rhs) : Rock(rhs) {}

                                        // Calculate elastic and seismic parameters, to be
                                        // used whenever new information is sent to class.
  void                                  ComputeSeismicAndElasticParams();
  void                                  Clone(const std::vector< Fluid* > & fluid_in);
  void                                  DeleteInclusion();

  Solid                               * solid_; // Owned and deleted by this class.
  std::vector<Fluid*>                   fluid_; // Owned and deleted by this class.
  std::vector<double>                   inclusion_spectrum_;
  std::vector<double>                   inclusion_concentration_;

  double k_, mu_;
};

#endif
