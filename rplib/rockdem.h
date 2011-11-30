#ifndef ROCKDEM_H
#define ROCKDEM_H

#include "rplib/rock.h"

class RockDEM : public Rock {
public:

  // Parallel classes are DistributionsGeochemicalDEM and DistributionsRockT0DEM.
  RockDEM(const std::vector<double> & param, const std::vector<double> & saturation) 
  : Rock(param, saturation)
  {   
  }

  virtual ~RockDEM(){}

  virtual void ComputeSeismicParams(double & vp, double & vs, double & rho) const {
    vp = 5.0;  //FAKE
    vs = 4.0;  //FAKE
    rho = 3.0; //FAKE
  }

  virtual Rock * Evolve(const std::vector<int>         & delta_time, 
                        const std::vector< Rock * >    & rock, 
                        const DistributionsSaturation  * dist_sat,
                        const DistributionsGeochemical * dist_geochem) const {
    
    assert(delta_time.size() == rock.size() + 1);
    
    std::vector<double> param_geochem;
    dist_geochem->GetParameters(param_geochem);
    std::vector<double> param_sat;
    dist_sat->GetParameters(param_sat);
    
    std::vector<double> evolved_param(param_.size(), param_[0] + param_geochem[0]);             //FAKE
    std::vector<double> evolved_saturation(saturation_.size(), saturation_[0] * param_sat[0]);  //FAKE
    Rock * new_rock = new RockDEM(evolved_param, evolved_saturation);
    return new_rock;
  }
};

#endif
