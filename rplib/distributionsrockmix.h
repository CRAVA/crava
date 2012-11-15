#ifndef RPLIB_DISTRIBUTIONS_ROCK_MIX_H
#define RPLIB_DISTRIBUTIONS_ROCK_MIX_H

#include "rplib/distributionsrock.h"
#include "rplib/demmodelling.h"

class Rock;
class DistributionWithTrend;

//This file contains two classes DistributionsRockMixOfRock and DistributionsRockMixOfSolidAndFluid.

//-------------------------------------- DistributionsRockMixOfRock ---------------------------------------------------------

class DistributionsRockMixOfRock : public DistributionsRock {
public:

  DistributionsRockMixOfRock(const std::vector< DistributionsRock * >           & distr_rock,
                             const std::vector< DistributionWithTrend * >       & distr_vol_frac,
                             DEMTools::MixMethod                                  mix_method,
                             const std::vector<double>                          & alpha,
                             const std::vector<double>                          & s_min,
                             const std::vector<double>                          & s_max);

  DistributionsRockMixOfRock(const DistributionsRockMixOfRock & dist);

  virtual ~DistributionsRockMixOfRock();

  virtual DistributionsRock             * Clone() const;

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                          * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                          * UpdateSample(double                      corr_param,
                                                       bool                        param_is_time,
                                                       const std::vector<double> & trend,
                                                       const Rock                * sample)       const;

  virtual std::vector<double>             GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>           GetCovariance(const std::vector<double> & trend_params)  const;

  virtual bool                            HasDistribution()                                        const;

  virtual std::vector<bool>               HasTrend()                                               const;

  virtual bool                            GetIsOkForBounding()                                     const { return false; }

private:

  Rock                                  * GetSample(const std::vector<double> & u,
                                                    const std::vector<double> & trend_params,
                                                    const std::vector<Rock *> & sample_rock) const;

  std::vector< DistributionsRock * >      distr_rock_;
  std::vector< DistributionWithTrend * >  distr_vol_frac_;
  DEMTools::MixMethod                     mix_method_;
};

//-------------------------------------- DistributionsRockMixOfSolidAndFluid ---------------------------------------------------------

class DistributionsSolid;
class DistributionsFluid;
class Solid;
class Fluid;


class DistributionsRockMixOfSolidAndFluid : public DistributionsRock {
public:

  DistributionsRockMixOfSolidAndFluid(const std::vector< DistributionsSolid * >           & distr_solid,
                                      const std::vector< DistributionsFluid * >           & distr_fluid,
                                      const std::vector< DistributionWithTrend * >        & distr_vol_frac_solid,
                                      const std::vector< DistributionWithTrend * >        & distr_vol_frac_fluid,
                                      DEMTools::MixMethod                                   mix_method,
                                      const std::vector<double>                           & alpha,
                                      const std::vector<double>                           & s_min,
                                      const std::vector<double>                           & s_max);

  DistributionsRockMixOfSolidAndFluid(const DistributionsRockMixOfSolidAndFluid & dist);

  virtual ~DistributionsRockMixOfSolidAndFluid();

  virtual DistributionsRock             * Clone() const;

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                          * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                          * UpdateSample(double                      corr_param,
                                                       bool                        param_is_time,
                                                       const std::vector<double> & trend,
                                                       const Rock                * sample)       const;

  virtual std::vector<double>             GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>           GetCovariance(const std::vector<double> & trend_params)  const;

  virtual bool                            HasDistribution()                                        const;

  virtual std::vector<bool>               HasTrend()                                               const;

  virtual bool                            GetIsOkForBounding()                                     const;

private:
  Rock                                  * GetSample(const std::vector<double>  & u,
                                                    const std::vector<double>  & trend_params,
                                                    const std::vector<Solid *> & solid_sample,
                                                    const std::vector<Fluid *> & fluid_sample) const;

  std::vector< DistributionsSolid * >     distr_solid_;
  std::vector< DistributionsFluid * >     distr_fluid_;
  std::vector< DistributionWithTrend * >  distr_vol_frac_solid_;
  std::vector< DistributionWithTrend * >  distr_vol_frac_fluid_;
  DEMTools::MixMethod                     mix_method_;
};

#endif
