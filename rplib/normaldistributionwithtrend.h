#ifndef RPLIB_NORMALDISTRIBUTIONWITHTREND_H
#define RPLIB_NORMALDISTRIBUTIONWITHTREND_H

#include "nrlib/random/distribution.hpp"

#include "rplib/distributionwithtrend.h"

namespace NRLib {
  class Trend;
}

class NormalDistributionWithTrend : public DistributionWithTrend {
 public:

   NormalDistributionWithTrend();

   NormalDistributionWithTrend(const NRLib::Trend * mean,
                               const NRLib::Trend * var,
                               int                  shared);

   NormalDistributionWithTrend(const NormalDistributionWithTrend & dist);

   virtual ~NormalDistributionWithTrend();

   virtual DistributionWithTrend    * Clone() const                           { return new NormalDistributionWithTrend(*this) ;}

   virtual bool                       GetIsDistribution() const               { return(true)                                  ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                       ;}



   virtual double                     ReSample(double s1, double s2);
   virtual double                     GetQuantileValue(double u, double s1, double s2);
   virtual double                     GetMeanValue(double s1, double s2);
   virtual double                     GetVarianceValue(double s1, double s2);

 private:

  const NRLib::Distribution<double> * gaussian_;
  const NRLib::Trend                * mean_;
  const NRLib::Trend                * var_;
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used
};
#endif
