#ifndef DISTRIBUTIONSFLUIDEVOLUTION_H
#define DISTRIBUTIONSFLUIDEVOLUTION_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.

class DistributionsFluidEvolution {
public:

  DistributionsFluidEvolution(){}

  virtual ~DistributionsFluidEvolution(){}

  // Functions for sampling the distribution functions and returning the sampled values
  // must be declared and defined for each derived class.
  // These functions are to be used by the Evolve function of the parallel, derived Fluid class.

};

#endif
