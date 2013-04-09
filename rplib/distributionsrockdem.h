#ifndef RPLIB_DISTRIBUTIONSROCKDEM_H
#define RPLIB_DISTRIBUTIONSROCKDEM_H

#include "rplib/distributionsrock.h"

#include "nrlib/grid/grid2d.hpp"


class DistributionWithTrend;
class DistributionsSolid;
class DistributionsFluid;
class Solid;
class Fluid;

//NBNB fjellvoll not finished yet
//TODO:covariance and expectations functions and variables.

class DistributionsRockDEM : public DistributionsRock {
public:

 DistributionsRockDEM(DistributionsSolid                           * distr_solid,
                      std::vector< DistributionsFluid *>           & distr_fluid,
                      std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                      std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                      const std::vector<double>                    & alpha,
                      const std::vector<double>                    & s_min,
                      const std::vector<double>                    & s_max);

 DistributionsRockDEM(const DistributionsRockDEM & dist);

  virtual                                        ~DistributionsRockDEM();

  virtual DistributionsRock                    * Clone() const;

  virtual Rock                                 * UpdateSample(double                      corr_param,
                                                              bool                        param_is_time,
                                                              const std::vector<double> & trend,
                                                              const Rock                * sample);

  virtual bool                                   HasDistribution() const;

  virtual std::vector<bool>                      HasTrend() const;

  virtual bool                                   GetIsOkForBounding()                                                  const { return false; }

private:
  virtual Rock                                 * GenerateSamplePrivate(const std::vector<double> & trend_params);

  void                                           SampleVpVsRhoExpectationAndCovariance(std::vector<double>   & expectation,
                                                                                       NRLib::Grid2D<double> & covariance);

  Rock                                         * GetSample(const std::vector<double>  & u,
                                                           const std::vector<double>  & trend_params,
                                                           const Solid                * solid,
                                                           const std::vector< Fluid *>& fluid);

  DistributionsSolid                           * distr_solid_;              // Pointer to external object.
  std::vector< DistributionsFluid *>             distr_fluid_;              // Pointer to external object.
  std::vector< DistributionWithTrend * >         distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< DistributionWithTrend * >         distr_incl_concentration_; // Pointers to external objects.
};

#endif
