#ifndef DISTRIBUTIONSSATURATIONEXMPL_H
#define DISTRIBUTIONSSATURATIONEXMPL_H

#include "rplib/distributionssaturation.h"


class DistributionsSaturationExmpl: public DistributionsSaturation {
public:

  DistributionsSaturationExmpl(const std::vector<Fluid *> & fluid, const std::vector<double> & p) 
    : DistributionsSaturation(fluid), p_(p){}

  virtual ~DistributionsSaturationExmpl(){}

  virtual void GetParameters(std::vector<double> & param_sat) const {
    param_sat.resize(1);  //FAKE!
    param_sat[0] = 0.9;   //FAKE!
  }

private:
  const std::vector<double> & p_;
};

#endif
