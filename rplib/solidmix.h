#ifndef RPLIB_SOLIDMIX_H
#define RPLIB_SOLIDMIX_H

#include <vector>

#include "rplib/solid.h"
#include "rplib/demmodelling.h"

class SolidMix : public Solid {
public:

  SolidMix(const std::vector<Solid*>      & solid,
           const std::vector<double>      & volume_fraction,
           const std::vector<double>      & u,
           DEMTools::MixMethod              mix_method);


  virtual ~SolidMix();

  // Assignment operator.
  SolidMix                        & operator=(const SolidMix& rhs);

  virtual Solid                   * Clone() const;

  virtual Solid                   * Evolve(const std::vector<int>             & delta_time,
                                           const std::vector< const Solid * > & solid) const;

private:
  //Copy constructor for getting base class variables , used by Clone:
  SolidMix(const SolidMix & rhs) : Solid(rhs) {}

  std::vector<Solid*>               solid_;           // Owned and deleted by this class.
  std::vector<double>               volume_fraction_;
  DEMTools::MixMethod               mix_method_;
};

#endif
