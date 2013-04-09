#ifndef RPLIB_DISTRIBUTIONS_DRYROCK_DEM_H
#define RPLIB_DISTRIBUTIONS_DRYROCK_DEM_H

#include "rplib/distributionsdryrock.h"

class DistributionWithTrend;

class DistributionsDryRockDEM : public DistributionsDryRock {
public:

  DistributionsDryRockDEM(DistributionsDryRock                         * distr_dryrock,
                          std::vector<DistributionsDryRock*>           & distr_dryrock_inc,
                          std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                          std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                          std::vector<double>                          & alpha);

  DistributionsDryRockDEM(const DistributionsDryRockDEM & dist);

  virtual                       ~DistributionsDryRockDEM();

  virtual DistributionsDryRock  * Clone() const;

  virtual DryRock               * GenerateSample(const std::vector<double> & trend_params);

  virtual bool                    HasDistribution() const;

  virtual std::vector<bool>       HasTrend() const;

  virtual DryRock               * UpdateSample(double                      corr_param,
                                               bool                        param_is_time,
                                               const std::vector<double> & trend,
                                               const DryRock             * sample);
protected:

private:
  DryRock                       * GetSample(const std::vector<double>    & u,
                                            const std::vector<double>    & trend_params,
                                            const DryRock                * dryrock,
                                            const std::vector< DryRock* >& dryrock_inc);

  DistributionsDryRock                           * distr_dryrock_;              // Pointer to external object.
  std::vector< DistributionsDryRock*>              distr_dryrock_inc_;          // Pointer to external object.
  std::vector< DistributionWithTrend * >           distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< DistributionWithTrend * >           distr_incl_concentration_; // Pointers to external objects.
};

#endif
