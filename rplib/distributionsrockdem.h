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
                      std::vector<double>                          & alpha);

  virtual                                        ~DistributionsRockDEM();

  virtual Rock                                 * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                                 * UpdateSample(double                      corr_param,
                                                              bool                        param_is_time,
                                                              const std::vector<double> & trend,
                                                              const Rock                * sample) const;

  virtual Pdf3D                                * GeneratePdf(void) const; // Returns NULL.

  virtual bool                                   HasDistribution() const;

  virtual std::vector<bool>                      HasTrend() const;

  virtual bool                                   GetIsOkForBounding()                                                  const { return false; }

  virtual std::vector<double>                    GetExpectation(const std::vector<double> & /*trend_params*/) const { return expectation_old_; }

  virtual NRLib::Grid2D<double>                  GetCovariance(const std::vector<double> & /*trend_params*/)  const { return covariance_old_; }

private:
  void                                           SampleVpVsRhoExpectationAndCovariance(std::vector<double>   & expectation,
                                                                                       NRLib::Grid2D<double> & covariance);

  Rock                                         * GetSample(const std::vector<double>  & u,
                                                           const std::vector<double>  & trend_params,
                                                           const Solid                * solid,
                                                           const std::vector< Fluid *>& fluid) const;

  DistributionsSolid                           * distr_solid_;              // Pointer to external object.
  std::vector< DistributionsFluid *>             distr_fluid_;              // Pointer to external object.
  std::vector< DistributionWithTrend * >         distr_incl_spectrum_;      // Pointers to external objects.
  std::vector< DistributionWithTrend * >         distr_incl_concentration_; // Pointers to external objects.

  //NBNB fjellvoll dummy to be removed soon
  std::vector<double>                            expectation_old_;
  NRLib::Grid2D<double>                          covariance_old_;

};

#endif
