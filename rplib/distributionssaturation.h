#ifndef DISTRIBUTIONSSATURATION_H
#define DISTRIBUTIONSSATURATION_H


// Abstract class for holding distribution functions for
// parameters needed for time evolution of saturation.
class DistributionsSaturation {
public:

  // Must use the same fluid vector as set by Rock::SetFluid(fluid).
  DistributionsSaturation(const std::vector<Fluid *> & fluid) : fluid_(fluid){}

  virtual ~DistributionsSaturation(){}

  virtual void GetParameters(std::vector<double> & param_sat) const = 0; //Resample the distributions and return the new sample.

protected:
  const std::vector<Fluid *> fluid_;
};

#endif
