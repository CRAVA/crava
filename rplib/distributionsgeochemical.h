#ifndef DISTRIBUTIONSGEOCHEMICAL_H
#define DISTRIBUTIONSGEOCHEMICAL_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of geochemical effects.
// One derived class for each rock physics model, the latter specified in a parallel, derived Rock class.

class DistributionsGeochemical {
public:

  DistributionsGeochemical(){}

  virtual ~DistributionsGeochemical(){}

  // Functions for sampling the distribution functions and returning the sampled values
  // must be declared and defined for each derived class.
  // These functions are to be used by the Evolve function of the parallel, derived Fluid class.

};

#endif
