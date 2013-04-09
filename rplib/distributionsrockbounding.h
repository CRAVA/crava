#ifndef RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H
#define RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H

#include "rplib/distributionsrock.h"

class Rock;
class DistributionWithTrend;
class Tabulated;

class DistributionsRockBounding : public DistributionsRock {
public:

  DistributionsRockBounding(DistributionsRock           * upper_rock,
                            DistributionsRock           * lower_rock,
                            DistributionWithTrend       * porosity,
                            DistributionWithTrend       * bulk_weight,
                            DistributionWithTrend       * shear_weight,
                            double                        correlation_weights,
                            const std::vector<double>   & alpha,
                            const std::vector<double>   & s_min,
                            const std::vector<double>   & s_max);

  DistributionsRockBounding(const DistributionsRockBounding & dist);

  virtual ~DistributionsRockBounding();

  virtual DistributionsRock        * Clone() const;

  virtual Rock                     * UpdateSample(double                      corr_param,
                                                  bool                        param_is_time,
                                                  const std::vector<double> & trend,
                                                  const Rock                * sample);

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const;

  virtual bool                       GetIsOkForBounding() const;

private:
  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSamplePrivate(const std::vector<double> & trend_params);

  Rock *                             GetSample(const std::vector<double> & u,
                                               const std::vector<double> & trend_params,
                                               const Rock                * sample_upper_rock,
                                               const Rock                * sample_lower_rock);

  DistributionsRock     * upper_rock_;
  DistributionsRock     * lower_rock_;
  DistributionWithTrend * porosity_;
  DistributionWithTrend * K_weight_;
  DistributionWithTrend * G_weight_;
  double                  correlation_weights_;
  Tabulated             * tabulated_;
};

#endif
