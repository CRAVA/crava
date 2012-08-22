#ifndef DISTRIBUTIONS_DRY_ROCK_H
#define DISTRIBUTIONS_DRY_ROCK_H

#include "rplib/dryrock.h"

// Abstract class for holding all t = 0 distribution functions for DryRock parameters.
// One derived class for each DryRock model, the latter specified in a parallel, derived DryRock class.
// The class must be able to produce an object of the specific DryRock class.
class DistributionsDryRock {
public:

  //NBNB: Vet ikke hvoran dry rock skal håndteres. Har derfor laget DistributionsDryRock som ikke nødvendigvis skal brukes

  DistributionsDryRock(){}

  virtual ~DistributionsDryRock(){}

  // DryRock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual DryRock * GenerateSample(const std::vector<double> & /*trend_params*/) const = 0;

  virtual bool                  HasDistribution()                                        const = 0;

  virtual std::vector<bool>     HasTrend()                                               const = 0;
};

#endif
