#ifndef ROCK_H
#define ROCK_H

#include <assert.h>
#include "rplib/distributionssaturation.h"
#include "rplib/distributionsgeochemical.h"


// Abstract rock class.
// Each derived class has parallel classes derived from DistributionsRockT0 and DistributionsGeochemical.
class Rock {
public:
  Rock();
  virtual ~Rock();

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const = 0;

  // Rock is an abstract class, hence pointer must be used in Evolve.  
  // Allocated memory (using new) MUST be deleted by caller.
  // Derived class Evolve implementation should always start with casting and assert.
  // Example from derived class RockDEM:
          // DistributionsGeochemicalDEM const * dist_geochem_dem 
          //      = dynamic_cast<const DistributionsGeochemicalDEM*>(dist_geochem);
          // assert(dist_geochem_dem != NULL);
          // assert(delta_time.size() == rock.size() + 1);
  virtual Rock * Evolve(const std::vector<int>         & delta_time,
                        const std::vector< Rock * >    & rock,
                        const DistributionsSaturation  * dist_sat,
                        const DistributionsGeochemical * dist_geochem) const = 0;
  

protected:

};

#endif
