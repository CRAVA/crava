#ifndef RPLIB_ROCKMIX_H
#define RPLIB_ROCKMIX_H

#include <vector>

#include "rplib/rock.h"
#include "rplib/demmodelling.h"



class RockMixOfRock : public Rock {
public:

  RockMixOfRock(const std::vector<Rock*>  &   rock,
                const std::vector<double> &   volume_fraction,
                DEMTools::MixMethod           mix_method);

  virtual ~RockMixOfRock();

  // Assignment operator.
  RockMixOfRock             & operator=(const RockMixOfRock& rhs);

  virtual Rock              * Clone() const;

  virtual Rock *              Evolve(const std::vector<int>              & delta_time,
                                     const std::vector< Rock * >         & rock)                         const;

  virtual double              GetPorosity()                                                              const;

  virtual void                SetPorosity(double porosity);


private:
  //Copy constructor for getting base class variables , used by Clone:
  RockMixOfRock(const RockMixOfRock & rhs) : Rock(rhs) {}

  std::vector<Rock*>          rock_;           // Owned and deleted by this class.
  std::vector<double>         volume_fraction_;
  DEMTools::MixMethod         mix_method_;
  double                      porosity_;
};

#endif
