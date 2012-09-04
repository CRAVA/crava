#ifndef RPLIB_DELTADISTRIBUTIONWITHTREND_H
#define RPLIB_DELTADISTRIBUTIONWITHTREND_H

#include "nrlib/trend/trend.hpp"
#include "nrlib/random/distribution.hpp"

#include "rplib/distributionwithtrend.h"

class DeltaDistributionWithTrend : public DistributionWithTrend {
 public:

   DeltaDistributionWithTrend();

   DeltaDistributionWithTrend(const NRLib::Trend * mean,
                              bool                 shared);

   virtual ~DeltaDistributionWithTrend();

   virtual bool                       GetIsShared() const                     { return(is_shared_)                          ;}
   virtual bool                       GetIsDistribution() const               { return(false)                               ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                     ;}

   virtual double                     ReSample(double s1, double s2) const;

 private:

  const NRLib::Distribution<double> * dirac_;
  const NRLib::Trend                * mean_;
  const bool                          is_shared_;       // Use is_shared_ like in DistributionWithTrendStorage to know if we have a reservoir variable.
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used


};
#endif
