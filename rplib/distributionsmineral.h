#ifndef DISTRIBUTIONSMINERAL_H
#define DISTRIBUTIONSMINERAL_H

#include "rplib/distributionssolid.h"
#include "rplib/mineral.h"

class DistributionsMineral : public DistributionsSolid {
public:

  DistributionsMineral() : DistributionsSolid(){}

  virtual ~DistributionsMineral(){}

  virtual Solid * GenerateSample() const {
    double k   = 30.0;     //FAKE
    double mu  = 10.0;     //FAKE
    double rho = 0.05;     //FAKE
    Solid * solid = new Mineral(k, mu, rho);
    return solid;
  }
};

#endif
