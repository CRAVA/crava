#ifndef RPLIB_BETAENDMASSDISTRIBUTIONWITHTREND_H
#define RPLIB_BETAENDMASSDISTRIBUTIONWITHTREND_H

#include "rplib/distributionwithtrend.h"

namespace NRLib {
  class Trend;
  template <typename T>
  class Grid2D;
}

class BetaEndMassDistributionWithTrend : public DistributionWithTrend {
 public:

   BetaEndMassDistributionWithTrend();

   BetaEndMassDistributionWithTrend(const NRLib::Trend * mean,
                                    const NRLib::Trend * var,
                                    const double       & lower_limit,
                                    const double       & upper_limit,
                                    const double       & lower_probability,
                                    const double       & upper_probability,
                                    int                  shared);

   BetaEndMassDistributionWithTrend(const BetaEndMassDistributionWithTrend & dist);

   virtual ~BetaEndMassDistributionWithTrend();

   virtual DistributionWithTrend    * Clone() const                           { return new BetaEndMassDistributionWithTrend(*this) ;}

   virtual bool                       GetIsShared() const                     { return(share_level_ > None)                 ;}
   virtual bool                       GetIsDistribution() const               { return(true)                                ;}
   virtual std::vector<bool>          GetUseTrendCube() const                 { return(use_trend_cube_)                     ;}

   //Triggers resampling for share_level_ <= level_. Not necessary for share_level_ = 0/None
   virtual void                       TriggerNewSample(int level)             {if(share_level_<=level)
                                                                                 resample_ = true; }

   virtual double                     ReSample(double s1, double s2);
   virtual double                     GetQuantileValue(double u, double s1, double s2);

private:
  NRLib::Grid2D<NRLib::Distribution<double> *> * beta_endmass_distribution_;
  const NRLib::Trend                           * mean_;
  const NRLib::Trend                           * var_;
  std::vector<bool>                              use_trend_cube_;   // First element true if first trend cube is used, second true if second is used, and both true if both are used
  int                                            ni_;
  int                                            nj_;
  std::vector<double>                            mean_sampling_;
  std::vector<double>                            var_sampling_;
  int                                            n_samples_mean_;
  int                                            n_samples_var_;
};
#endif
