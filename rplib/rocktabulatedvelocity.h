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

  virtual void                          SetPorosity(double porosity);

};

#endif
