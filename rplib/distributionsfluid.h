#ifndef RPLIB_DISTRIBUTIONS_FLUID_H
#define RPLIB_DISTRIBUTIONS_FLUID_H

#include "rplib/fluid.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluid {
public:

                                DistributionsFluid(){}

  virtual                       ~DistributionsFluid(){}

  virtual DistributionsFluid  * Clone()                                                       const = 0;

  virtual Fluid *               GenerateSample(const std::vector<double> & /*trend_params*/)  const = 0;

  std::vector< Fluid* >         GenerateWellSample(const  std::vector<double> & trend_params,
                                                   double                       corr)         const;

  virtual bool                  HasDistribution()                                             const = 0;

  virtual std::vector<bool>     HasTrend()                                                    const = 0;

  virtual Fluid *               UpdateSample(double                      corr_param,
                                             bool                        param_is_time,
                                             const std::vector<double> & trend,
                                             const Fluid               * sample)             const = 0;

  Fluid *                       EvolveSample(double         time,
                                             const Fluid &  fluid)                            const;

protected:
  std::vector< double >         alpha_;

};

#endif
