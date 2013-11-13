#ifndef RPLIB_DELTADISTRIBUTIONWITHTREND_H
#define RPLIB_DELTADISTRIBUTIONWITHTREND_H

#include "rplib/distributionwithtrend.h"

namespace NRLib {
  class Trend;
  template <typename T>
  class Distribution;
}

class DeltaDistributionWithTrend : public DistributionWithTrend {
 public:

   DeltaDistributionWithTrend();

   DeltaDistributionWithTrend(const NRLib::Trend * mean,
                              int                  shared);

   DeltaDistributionWithTrend(const DeltaDistributionWithTrend & dist);

   virtual ~DeltaDistributionWithTrend();

   virtual DistributionWithTrend    * Clone() const                           { return new DeltaDistributionWithTrend(*this) ;}

   virtual bool                       GetIsShared() const                     { return(share_level_ > None)                  ;}
   virtual bool                       GetIsDistribution() const               { return(false)                                ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                      ;}

   //Triggers resampling for share_level_ <= level_. Not necessary for share_level_ = 0/None
   virtual void                       TriggerNewSample(int level)             {if(share_level_<=level)
                                                                                 resample_ = true; }

   virtual double                     ReSample(double s1, double s2);
   virtual double                     GetQuantileValue(double u, double s1, double s2);
   virtual double                     GetMeanValue(double s1, double s2);
   virtual double                     GetVarianceValue(double s1, double s2);

 private:

  const NRLib::Distribution<double> * dirac_;
  const NRLib::Trend                * mean_;
  std::vector<bool>                   use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used
 // const int                           share_level_;      // Use like in DistributionWithTrendStorage to know if we have a reservoir variable.
 // double                              current_u_;        // Quantile of current sample.
 // bool                                resample_;         // If false, and share_level_ > 0, reuse current_u_

};
#endif
