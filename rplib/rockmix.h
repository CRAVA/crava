#ifndef RPLIB_ROCKMIX_H
#define RPLIB_ROCKMIX_H

#include <vector>

#include "rplib/rock.h"
#include "rplib/demmodelling.h"

//This file contains two classes RockMixOfRock and RockMixOfSolidAndFluid.

//-------------------------------------- RockMixOfRock ---------------------------------------------------------

class RockMixOfRock : public Rock {
public:

  RockMixOfRock(const std::vector<Rock*>  &   rock,
                const std::vector<double> &   volume_fraction,
                const std::vector<double> &   u,
                DEMTools::MixMethod           mix_method);

  virtual ~RockMixOfRock();

  // Assignment operator.
  RockMixOfRock             & operator=(const RockMixOfRock& rhs);

  virtual Rock              * Clone() const;

  virtual Rock *              Evolve(const std::vector<int>              & delta_time,
                                     const std::vector< Rock * >         & rock)                         const;

  virtual void                SetPorosity(double porosity);


private:
  //Copy constructor for getting base class variables , used by Clone:
  RockMixOfRock(const RockMixOfRock & rhs) : Rock(rhs) {}

  void                        ComputeSeismicVariables();

  std::vector<Rock*>          rock_;           // Owned and deleted by this class.
  std::vector<double>         volume_fraction_;
  DEMTools::MixMethod         mix_method_;
};

//----------------------------------- RockMixOfSolidAndFluid -------------------------------------------------

class Solid;
class Fluid;

class RockMixOfSolidAndFluid : public Rock {
public:

  RockMixOfSolidAndFluid(const std::vector<Solid*> &   solid,
                         const std::vector<Fluid*> &   fluid,
                         const std::vector<double> &   volume_fraction_solid,
                         const std::vector<double> &   volume_fraction_fluid,
                         const std::vector<double> &   u,
                         DEMTools::MixMethod           mix_method);

  virtual ~RockMixOfSolidAndFluid();

  // Assignment operator.
  RockMixOfSolidAndFluid              & operator=(const RockMixOfSolidAndFluid& rhs);

  virtual Rock                        * Clone()                                                                    const;

  virtual Rock *                        Evolve(const std::vector<int>              & delta_time,
                                               const std::vector< Rock * >         & /*rock*/)                     const;

  virtual void                          SetPorosity(double porosity);

private:
  //Copy constructor for getting base class variables , used by Clone:
  RockMixOfSolidAndFluid(const RockMixOfSolidAndFluid & rhs) : Rock(rhs) {}

  void                                  ComputeSeismicVariables();

  std::vector< Solid* >                 solid_;             // Owned and deleted by this class.
  std::vector< Fluid* >                 fluid_;             // Owned and deleted by this class.
  std::vector<double>                   volume_fraction_solid_;
  std::vector<double>                   volume_fraction_fluid_;
  DEMTools::MixMethod                   mix_method_;
};

#endif
