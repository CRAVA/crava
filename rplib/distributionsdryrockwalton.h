#ifndef RPLIB_DISTRIBUTIONS_DRYROCK_WALTON_H
#define RPLIB_DISTRIBUTIONS_DRYROCK_WALTON_H

#include "rplib/distributionsdryrock.h"

class DistributionWithTrend;
class DistributionsSolid;
class Solid;

class DistributionsDryRockWalton : public DistributionsDryRock {
public:

  DistributionsDryRockWalton(DistributionsSolid                         * distr_solid,
                             DistributionWithTrend                      * distr_friction_weight,
                             DistributionWithTrend                      * distr_pressure,
                             DistributionWithTrend                      * distr_porosity,
                             DistributionWithTrend                      * distr_coord_number,
                             std::vector<double>                        & alpha);

  DistributionsDryRockWalton(const DistributionsDryRockWalton & dist);

  virtual                       ~DistributionsDryRockWalton();

  virtual DistributionsDryRock  * Clone() const;

  virtual DryRock               * GenerateSample(const std::vector<double> & trend_params);

  virtual bool                    HasDistribution() const;

  virtual std::vector<bool>       HasTrend() const;

  virtual DryRock               * UpdateSample(double                      corr_param,
                                               bool                        param_is_time,
                                               const std::vector<double> & trend,
                                               const DryRock             * sample);
protected:

private:
   DryRock                       * GetSample(const std::vector<double>    & u,
                                             const std::vector<double>    & trend_params,
                                             const Solid                  * solid);

   DistributionsSolid            * distr_solid_;
   DistributionWithTrend         * distr_friction_weight_;
   DistributionWithTrend         * distr_pressure_;
   DistributionWithTrend         * distr_porosity_;
   DistributionWithTrend         * distr_coord_number_;
};

#endif
