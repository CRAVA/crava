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

   virtual bool                       GetIsShared() const                     { return(share_level_ > None)                   ;}
   virtual bool                       GetIsDistribution() const               { return(true)                                  ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                       ;}

   //Triggers resampling for share_level_ <= level_. Not necessary for share_level_ = 0/None
   virtual void                       TriggerNewSample(int level)             {if(share_level_<=level)
                                                                                 resample_ = true; }

   virtual double                     ReSample(double s1, double s2);
   virtual double                     GetQuantileValue(double u, double s1, double s2);

 private:

  const NRLib::Distribution<double> * gaussian_;
  const NRLib::Trend                * mean_;
  const NRLib::Trend                * var_;
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used
  const int                           share_level_;      // Use like in DistributionWithTrendStorage to know if we have a reservoir variable.
  double                              current_u_;        // Quantile of current sample.
  bool                                resample_;         // If false, and share_level_ > 0, reuse current_u_

};
#endif
