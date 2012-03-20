#ifndef DISTRIBUTIONSSOLIDEVOLUTION_H
#define DISTRIBUTIONSSOLIDEVOLUTION_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.

class DistributionsSolidEvolution {
public:

  DistributionsSolidEvolution(){}

  virtual ~DistributionsSolidEvolution(){}

  // Functions for sampling the distribution functions and returning the sampled values
  // must be declared and defined for each derived class.
  // These functions are to be used by the Evolve function of the parallel, derived Solid class.

};

#endif
