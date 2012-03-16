#include "rplib/multinormalrock.h"

MultiNormalRock::MultiNormalRock(const std::vector<double>&  param,
                                 const std::vector<double>&  saturation)
: Rock(saturation)
{
  param_ = param;
}


MultiNormalRock::~MultiNormalRock()
{
}

void
MultiNormalRock::ComputeSeismicParams(double& vp, double& vs, double& rho) const {
  vp  = param_[0];
  vs  = param_[1];
  rho = param_[2];
}
