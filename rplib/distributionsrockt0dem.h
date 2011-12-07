#ifndef DISTRIBUTIONSROCKT0DEM_H
#define DISTRIBUTIONSROCKT0DEM_H

#include "rplib/rockdem.h"
#include "rplib/distributionsrockt0.h"

// Parallel classes are RockDEM and DistributionsGeochemicalDEM.
class DistributionsRockT0DEM : public DistributionsRockT0 {
public:

  DistributionsRockT0DEM() : DistributionsRockT0(){}

  virtual ~DistributionsRockT0DEM(){}

  virtual Rock  * GenerateSample(const std::vector<double> & trend_params) const {
    std::vector<double> param(6);      //FAKE
    if (trend_params.size() > 0)       //FAKE
      param[0] -= trend_params[0];     //FAKE
    std::vector<double> saturation(2); //FAKE
    saturation[0] = 1;                 //FAKE
    Rock * rock = new RockDEM(param, saturation);
    return rock;
  }

};

#endif
