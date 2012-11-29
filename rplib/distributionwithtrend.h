#ifndef RPLIB_DISTRIBUTIONWITHTREND_H
#define RPLIB_DISTRIBUTIONWITHTREND_H

#include <vector>

class DistributionWithTrend {
 public:
   DistributionWithTrend();

   virtual ~DistributionWithTrend();

   virtual DistributionWithTrend    * Clone()                                             const = 0;

   virtual bool                       GetIsShared()                                       const = 0;
   virtual bool                       GetIsDistribution()                                 const = 0;
   virtual std::vector<bool>          GetUseTrendCube()                                   const = 0;

   virtual double                     ReSample(double s1, double s2)                      const = 0;
   virtual double                     GetQuantileValue(double u, double s1, double s2)    const = 0;

   void                               FindUseTrendCube(std::vector<bool> & use_trend_cube,
                                                       int                 dim,
                                                       int                 reference);

};
#endif
