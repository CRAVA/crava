#ifndef RPLIB_DISTRIBUTIONSROCK_GASSMANN_H
#define RPLIB_DISTRIBUTIONSROCK_GASSMANN_H

#include "rplib/distributionsrock.h"

//#include "nrlib/grid/grid2d.hpp"

class DistributionsDryRock;
class DistributionsFluid;
class DryRock;
class Fluid;

class DistributionsRockGassmann : public DistributionsRock {
public:

 DistributionsRockGassmann(const DistributionsDryRock            * distr_dryrock,
                           const DistributionsFluid              * distr_fluid,
                           const std::vector<double>             & s_min,
                           const std::vector<double>             & s_max);

 DistributionsRockGassmann(const DistributionsRockGassmann & dist);

  virtual                                        ~DistributionsRockGassmann();

  virtual DistributionsRock                    * Clone() const;

  virtual Rock                                 * UpdateSample(double                      corr_param,
                                                              bool                        param_is_time,
                                                              const std::vector<double> & trend,
                                                              const Rock                * sample);

  virtual bool                                   HasDistribution() const;

  virtual std::vector<bool>                      HasTrend() const;

  virtual bool                                   GetIsOkForBounding()                                                  const { return false; }

private:
  virtual Rock                                 * GenerateSamplePrivate(const std::vector<double> & trend_params);

  Rock                                         * GetSample(const DryRock              * dryrock,
                                                           const Fluid                * fluid);

  DistributionsDryRock                         * distr_dryrock_;
  DistributionsFluid                           * distr_fluid_;

};

#endif
