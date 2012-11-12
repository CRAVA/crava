#ifndef RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H
#define RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H

#include "rplib/distributionsrock.h"

class Rock;
class DistributionWithTrend;
class Tabulated;

class DistributionsRockBounding : public DistributionsRock {
public:

  DistributionsRockBounding(const DistributionsRock     * upper_rock,
                            const DistributionsRock     * lower_rock,
                            const DistributionWithTrend * porosity,
                            const DistributionWithTrend * bulk_weight,
                            const DistributionWithTrend * p_wave_weight,
                            double                        correlation_weights,
                            std::vector<double>         & alpha);

  DistributionsRockBounding(const DistributionsRockBounding & dist);

  virtual ~DistributionsRockBounding();

  virtual DistributionsRock        * Clone() const;

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                     * UpdateSample(double                      corr_param,
                                                  bool                        param_is_time,
                                                  const std::vector<double> & trend,
                                                  const Rock                * sample)       const;

  virtual std::vector<double>        GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>      GetCovariance(const std::vector<double> & trend_params)  const;

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const;

  virtual bool                       GetIsOkForBounding() const;

private:

  Rock *                             GetSample(const std::vector<double> & u,
                                               const std::vector<double> & trend_params,
                                               const Rock                * sample_upper_rock,
                                               const Rock                * sample_lower_rock) const;

  const DistributionsRock     * upper_rock_;
  const DistributionsRock     * lower_rock_;
  const DistributionWithTrend * porosity_;
  const DistributionWithTrend * K_weight_;
  const DistributionWithTrend * M_weight_;
  double                        correlation_weights_;
  Tabulated                   * tabulated_;
};

#endif
