#ifndef DISTRIBUTIONSBRINEEVOLUTION_H
#define DISTRIBUTIONSBRINEEVOLUTION_H

#include "rplib/distributionsfluidevolution.h"

// Parallel classes are Brine and DistributionsBrineT0.
class DistributionsBrineEvolution : public DistributionsFluidEvolution {
public:

  DistributionsBrineEvolution() : DistributionsFluidEvolution(){}

  virtual ~DistributionsBrineEvolution(){}

  double Sample(){return 1.0:}  // FAKE

};

#endif
