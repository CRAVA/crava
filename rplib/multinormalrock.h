#ifndef RPLIB_MULTINORMALROCK_H
#define RPLIB_MULTINORMALROCK_H

#include "rplib/rock.h"

class MultiNormalRock : public Rock {
public:

  // Parallel classes are
  MultiNormalRock(const std::vector<double>&  param,
                  const std::vector<double>&  saturation);

  virtual ~MultiNormalRock();

  virtual void ComputeSeismicParams(double& vp, double& vs, double& rho) const;

  virtual Rock * Evolve(const std::vector<int>         & /*delta_time*/,
                        const std::vector< Rock * >    & /*rock*/,
                        const DistributionsSaturation  * /*dist_sat*/,
                        const DistributionsGeochemical * /*dist_geochem*/) const { return new MultiNormalRock(*this); }


private:
  std::vector<double>  param_; //vp, vs, rho

};

#endif
