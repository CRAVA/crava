#ifndef ROCKDEM_H
#define ROCKDEM_H

#include "rplib/rock.h"

class RockDEM : public Rock {
public:

  // Parallel classes are DistributionsGeochemicalDEM and DistributionsRockT0DEM.
  RockDEM(const double par_dem, const std::vector<double> & saturation)
  : Rock(saturation)
  {
    par_dem_ = par_dem;
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

    double evolved_param_dem = par_dem_ + param_geochem[0];             //FAKE
    std::vector<double> evolved_saturation(saturation_.size(), saturation_[0] * param_sat[0]);  //FAKE
    Rock * new_rock = new RockDEM(evolved_param_dem, evolved_saturation);
    return new_rock;
  }

private:
  double par_dem_; // Example, to be substituted with real DEM parameters.
};

#endif
