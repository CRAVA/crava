#ifndef DISTRIBUTIONSGEOCHEMICALDEM_H
#define DISTRIBUTIONSGEOCHEMICALDEM_H

#include "rplib/distributionsgeochemical.h"

// Parallel classes are RockDEM and DistributionsRockT0DEM.
class DistributionsGeochemicalDEM : public DistributionsGeochemical{
public:

  DistributionsGeochemicalDEM() : DistributionsGeochemical() {}

  virtual ~DistributionsGeochemicalDEM(){}

  void GetParameters(std::vector<double> & param_geochem) const {
    param_geochem.resize(1);  //FAKE!
    param_geochem[0] = 1.0;   //FAKE!
  }

};

#endif
