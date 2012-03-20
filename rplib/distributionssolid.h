#ifndef DISTRIBUTIONSSOLID_H
#define DISTRIBUTIONSSOLID_H

#include "rplib/solid.h"

// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolid {
public:

  DistributionsSolid(){}

  virtual ~DistributionsSolid(){}

  // Solid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Solid * GenerateSample() const = 0;
};

#endif
