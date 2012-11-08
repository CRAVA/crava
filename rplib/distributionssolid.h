#ifndef RPLIB_DISTRIBUTIONS_SOLID_H
#define RPLIB_DISTRIBUTIONS_SOLID_H

#include "rplib/solid.h"

// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolid {
public:

                                DistributionsSolid(){}

  virtual                       ~DistributionsSolid(){}

  virtual Solid *               GenerateSample(const std::vector<double> & /*trend_params*/) const = 0;

  std::vector< Solid* >         GenerateWellSample(const  std::vector<double> & trend_params,
                                                   double                       corr)        const;

  virtual bool                  HasDistribution()                                            const = 0;

  virtual std::vector<bool>     HasTrend()                                                   const = 0;

  Solid *                       EvolveSample(double         time,
                                             const Solid &  solid)                           const;

  virtual Solid *               UpdateSample(double                      corr_param,
                                             bool                        param_is_time,
                                             const std::vector<double> & trend,
                                             const Solid               * sample)            const = 0;

  const std::vector<double>     GetAlpha()                                                  const { return alpha_; }

protected:
  std::vector< double >         alpha_;
};

#endif
