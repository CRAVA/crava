#ifndef DISTRIBUTIONSFLUIDEVOLUTION_H
#define DISTRIBUTIONSFLUIDEVOLUTION_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.

class DistributionsFluidEvolution {
public:

  DistributionsFluidEvolution(){}

  virtual ~DistributionsFluidEvolution(){}

  virtual void GetParameters(std::vector<double> & param_fluid_evolve) const = 0; //Resample the distributions and return the new values.

};

#endif
