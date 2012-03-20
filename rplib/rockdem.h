#ifndef ROCKDEM_H
#define ROCKDEM_H

#include "rplib/rock.h"
#include "rplib/distributionsgeochemicaldem.h"

class RockDEM : public Rock {
public:

  // Parallel classes are DistributionsGeochemicalDEM and DistributionsRockT0DEM.
  RockDEM(const double par_dem)
  : Rock()
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
    DistributionsGeochemicalDEM const * dist_geochem_dem = dynamic_cast<const DistributionsGeochemicalDEM*>(dist_geochem);
    assert(dist_geochem_dem != NULL);
    assert(delta_time.size() == rock.size() + 1);

    //FAKE, TEMPORARY:
    std::vector<double> param_geochem_dem;
    dist_geochem_dem->GetParameters(param_geochem_dem);
    double evolved_param_dem = par_dem_ + param_geochem_dem[0];
    Rock * new_rock = new RockDEM(evolved_param_dem);

    return new_rock;
  }

private:
  double par_dem_; // Example, to be substituted with real DEM parameters.
};

#endif
