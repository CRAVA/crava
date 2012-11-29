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
                               bool                 shared);

   NormalDistributionWithTrend(const NormalDistributionWithTrend & dist);

   virtual ~NormalDistributionWithTrend();

   virtual DistributionWithTrend    * Clone() const                           { return new NormalDistributionWithTrend(*this) ;}

   virtual bool                       GetIsShared() const                     { return(is_shared_)                            ;}
   virtual bool                       GetIsDistribution() const               { return(true)                                  ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                       ;}

   virtual double                     ReSample(double s1, double s2) const;
   virtual double                     GetQuantileValue(double u, double s1, double s2) const;

 private:

  const NRLib::Distribution<double> * gaussian_;
  const NRLib::Trend                * mean_;
  const NRLib::Trend                * var_;
  const bool                          is_shared_;       // Use is_shared_ like in DistributionWithTrendStorage to know if we have a reservoir variable.
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used


};
#endif
