#ifndef DISTRIBUTIONS_DRY_ROCK_H
#define DISTRIBUTIONS_DRY_ROCK_H

#include "rplib/dryrock.h"

#include <vector>

// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived DryRock class.
// The class must be able to produce an object of the specific DryRock class.
class DistributionsDryRock {
public:

  DistributionsDryRock(){}

  virtual ~DistributionsDryRock(){}

  virtual DistributionsDryRock  * Clone()                                                      const = 0;

  virtual DryRock *               GenerateSample(const std::vector<double> & /*trend_params*/)       = 0;

  std::vector< DryRock* >         GenerateWellSample(const  std::vector<double> & trend_params,
                                                     double                       corr);

  virtual bool                    HasDistribution()                                            const = 0;

  virtual std::vector<bool>       HasTrend()                                                   const = 0;

  DryRock *                       EvolveSample(double           time,
                                               const DryRock &  dryrock)
                                  {
                                    const std::vector<double> trend(2);
                                    return UpdateSample(time, true, trend, &dryrock);
                                  }

  virtual DryRock *               UpdateSample(double                      corr_param,
                                               bool                        param_is_time,
                                               const std::vector<double> & trend,
                                               const DryRock             * sample)                   = 0;

protected:
  std::vector< double >           alpha_;
};

#endif
