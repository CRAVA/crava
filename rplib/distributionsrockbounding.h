#ifndef RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H
#define RPLIB_DISTRIBUTIONS_ROCK_BOUNDING_H

#include "rplib/rock.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"

class DistributionsRockBounding : public DistributionsRock {
public:

  DistributionsRockBounding(const DistributionsRock     * upper_rock,
                            const DistributionsRock     * lower_rock,
                            const DistributionWithTrend * porosity,
                            const DistributionWithTrend * bulk_weight,
                            const DistributionWithTrend * p_wave_weight,
                            double                        correlation_weights);

  virtual ~DistributionsRockBounding();

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                     * UpdateSample(const std::vector<double> & /*corr*/,
                                                  const Rock                & /*rock*/)       const { return NULL; }

  virtual std::vector<double>        GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>      GetCovariance(const std::vector<double> & trend_params)  const;

  virtual Pdf3D                    * GeneratePdf() const;

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const;

private:
  const DistributionsRock     * upper_rock_;
  const DistributionsRock     * lower_rock_;
  const DistributionWithTrend * porosity_;
  const DistributionWithTrend * K_weight_;
  const DistributionWithTrend * M_weight_;
  double                        correlation_weights_;
  Tabulated                   * tabulated_;
  bool                          has_distribution_;
  std::vector<bool>             has_trend_;
};

#endif
