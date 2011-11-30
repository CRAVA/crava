#ifndef ROCK_H
#define ROCK_H

#include <assert.h>

#include "rplib/fluid.h"
#include "rplib/distributionssaturation.h"
#include "rplib/distributionsgeochemical.h"


// Abstract rock class.
// Each derived class has parallel classes derived from DistributionsRockT0 and DistributionsGeochemical.
class Rock {
public:
  Rock(const std::vector<double> & param, const std::vector<double> & saturation);
  virtual ~Rock();

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const = 0;

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock * Evolve(const std::vector<int>         & delta_time, 
                        const std::vector< Rock * >    & rock, 
                        const DistributionsSaturation  * dist_sat,
                        const DistributionsGeochemical * dist_geochem) const = 0;
  
  void SetSaturation(std::vector<double> saturation) {saturation_ = saturation;}

  static void SetFluid(std::vector<Fluid *> fluid) {fluid_ = fluid;}  
  static std::vector<Fluid *> fluid_;  // The same objects for all Rocks, fixed fluid properties..

protected:
  const std::vector<double> param_;    // Sampled rock physics parameters, apart from saturation and fluid properties.
  std::vector<double> saturation_;     // Sampled saturations.
};

#endif
