#ifndef RPLIB_DISTRIBUTIONS_DRYROCK_MIX_H
#define RPLIB_DISTRIBUTIONS_DRYROCK_MIX_H

#include "rplib/distributionsdryrock.h"

#include "rplib/demmodelling.h"

class DistributionWithTrend;

class DistributionsDryRockMix : public DistributionsDryRock {
public:

  DistributionsDryRockMix(std::vector< DistributionsDryRock * >        & distr_dryrock,
                          std::vector< DistributionWithTrend * >       & distr_vol_frac,
                          DEMTools::MixMethod                            mix_method,
                          std::vector<double>                          & alpha);

  DistributionsDryRockMix(const DistributionsDryRockMix & dist);

  virtual ~DistributionsDryRockMix();

  virtual DistributionsDryRock                    * Clone() const;

  virtual DryRock                                 * GenerateSample(const std::vector<double> & trend_params);

  virtual bool                                      HasDistribution() const;

  virtual std::vector<bool>                         HasTrend() const;

  virtual DryRock *                                 UpdateSample(double                      corr_param,
                                                                 bool                        param_is_time,
                                                                 const std::vector<double> & trend,
                                                                 const DryRock             * sample);
protected:

private:

  DryRock                                         * GetSample(const std::vector<double>    & u,
                                                              const std::vector<double>    & trend_params,
                                                              const std::vector<DryRock *> & dryrock_samples);

  std::vector< DistributionsDryRock * >           distr_dryrock_;     // Pointers to external objects.
  std::vector< DistributionWithTrend * >          distr_vol_frac_;  // Pointers to external objects.
  DEMTools::MixMethod                             mix_method_;
};

#endif
