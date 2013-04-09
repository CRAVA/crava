#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include <vector>

class DistributionWithTrend {
 public:
   DistributionWithTrend();
   DistributionWithTrend(const int shareLevel,bool reSample);
   DistributionWithTrend(const int shareLevel,double currentU,bool reSample);

   virtual ~DistributionWithTrend();
   double                            GetCurrentSample(const std::vector<double> & trend_params);
   

   virtual DistributionWithTrend    * Clone()                                             const = 0;
   bool                               GetIsShared() const                                 { return(share_level_ > None)                   ;}
   virtual bool                       GetIsDistribution()                                 const = 0;
   virtual std::vector<bool>          GetUseTrendCube()                                   const = 0;

   //Triggers resampling for share_level_ <= level_. Not necessary for share_level_ = 0/None
   void                       TriggerNewSample(int level)             {if(share_level_<=level)
                                                                                 resample_ = true; }

   virtual double                     ReSample(double s1, double s2)                            = 0;
   virtual double                     GetQuantileValue(double u, double s1, double s2)          = 0;
  
   

   void                               FindUseTrendCube(std::vector<bool> & use_trend_cube,
                                                       int                 dim,
                                                       int                 reference);

   enum                               ShareLevel {None, SingleSample, Full}; //Note: New levels should be inserted between SingleSample and Full.
protected:
  const int                           share_level_;      // Use like in DistributionWithTrendStorage to know if we have a reservoir variable.
  double                              current_u_;        // Quantile of current sample.
  bool                                resample_;         // If false, and share_level_ > 0, reuse current_u_

};
#endif
