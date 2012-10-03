#ifndef RPLIB_SOLIDMIXED_H
#define RPLIB_SOLIDMIXED_H

#include <vector>
#include <numeric>

#include "rplib/solid.h"
#include "rplib/distributionssolidmixevolution.h"
#include "rplib/demmodelling.h"

class SolidMixed : public Solid {
public:

  SolidMixed(const std::vector<Solid*>      & solid,
             const std::vector<double>      & volume_fraction,
             DEMTools::MixMethod              mix_method,
             DistributionsSolidMixEvolution * distr_evolution = NULL);


  virtual ~SolidMixed();

  // Assignment operator.
  SolidMixed                      & operator=(const SolidMixed& rhs);

  virtual Solid                   * Clone() const;

  virtual Solid                   * Evolve(const std::vector<int>             & delta_time,
                                           const std::vector< const Solid * > & solid) const;

private:
  //Copy constructor for getting base class variables , used by Clone:
  SolidMixed(const SolidMixed & rhs) : Solid(rhs) {}

  std::vector<Solid*>               solid_;           // Owned and deleted by this class.
  std::vector<double>               volume_fraction_;
  DEMTools::MixMethod               mix_method_;
  DistributionsSolidMixEvolution  * distr_evolution_; // Pointer to external object.
};

#endif
