#ifndef RPLIB_ROCK_TABULATED_VELOCITY_H
#define RPLIB_ROCK_TABULATED_VELOCITY_H

#include "rplib/rock.h"

#include <vector>

class RockTabulatedVelocity : public Rock {
public:

  RockTabulatedVelocity(double               vp,
                        double               vs,
                        double               density,
                        std::vector<double>  u);

  RockTabulatedVelocity();

  RockTabulatedVelocity(const RockTabulatedVelocity & old_rock);

  virtual ~RockTabulatedVelocity();

  virtual Rock                        * Clone()                                                 const;

  virtual Rock                        * Evolve(const std::vector<int>         & delta_time,
                                               const std::vector< Rock * >    & rock)           const;

  virtual double                        GetPorosity()                                           const { return 0; }

  virtual void                          SetPorosity(double porosity);

};

#endif
