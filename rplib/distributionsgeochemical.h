#ifndef DISTRIBUTIONSGEOCHEMICAL_H
#define DISTRIBUTIONSGEOCHEMICAL_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of geochemical effects.
// One derived class for each rock physics model, the latter specified in a parallel, derived Rock class.

class DistributionsGeochemical {
public:

  DistributionsGeochemical(){}

  virtual ~DistributionsGeochemical(){}

  virtual void GetParameters(std::vector<double> & param_geochem) const = 0; //Resample the distributions and return the new sample.

};

#endif
