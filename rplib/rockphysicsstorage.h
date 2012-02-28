#ifndef RPLIB_ROCK_PHYSICS_STORAGE_HPP
#define RPLIB_ROCK_PHYSICS_STORAGE_HPP

#include "rplib/distributionsrockt0.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"


class RockPhysicsStorage {
public:
  RockPhysicsStorage();

  virtual ~RockPhysicsStorage();

  virtual DistributionsRockT0 * GenerateRockPhysics(const std::string                      & /*path*/,
                                                    const std::vector<std::string>         & /*trend_cube_parameters*/,
                                                    const std::vector<std::vector<float> > & /*trend_cube_sampling*/,
                                                    std::string                            & /*errTxt*/) const = 0;
};

class GaussianRockPhysicsStorage : public RockPhysicsStorage {
public:
  GaussianRockPhysicsStorage(NRLib::TrendStorage * mean_vp,
                             NRLib::TrendStorage * mean_vs,
                             NRLib::TrendStorage * mean_density,
                             NRLib::TrendStorage * variance_vp,
                             NRLib::TrendStorage * variance_vs,
                             NRLib::TrendStorage * variance_density,
                             NRLib::TrendStorage * correlation_vp_vs,
                             NRLib::TrendStorage * correlation_vp_density,
                             NRLib::TrendStorage * correlation_vs_density);

  virtual ~GaussianRockPhysicsStorage();

  virtual DistributionsRockT0 * GenerateRockPhysics(const std::string                      & path,
                                                    const std::vector<std::string>         & trend_cube_parameters,
                                                    const std::vector<std::vector<float> > & trend_cube_sampling,
                                                    std::string                            & errTxt) const;

private:
  int                                 FindNewGridDimension(const std::vector<NRLib::Trend *> trender) const;

  std::vector<int>                    FindNewGridSize(const std::vector<NRLib::Trend *> trender,
                                                      const int                         new_dim) const;

  std::vector<std::vector<double> >   ExpandGrids1D(const std::vector<NRLib::Trend *> trender,
                                                    const std::vector<int>         &  size) const;

  std::vector<NRLib::Grid2D<double> > ExpandGrids2D(const std::vector<NRLib::Trend *> trender,
                                                    const std::vector<int>         &  size) const;

  void LogTransformExpectationAndCovariance(NRLib::Trend *  mean1,
                                            NRLib::Trend *  mean2,
                                            NRLib::Trend *  cov,
                                            NRLib::Trend *& log_mean,
                                            NRLib::Trend *& log_cov,
                                            bool          & diagonal_element) const;

  void CalculateCovarianceFromCorrelation(NRLib::Trend *  corr,
                                          NRLib::Trend *  var1,
                                          NRLib::Trend *  var2,
                                          NRLib::Trend *& cov) const;

  void LogTransformExpectation(const double & expectation,
                               const double & variance,
                               double       & mu) const;

  void LogTransformCovariance(const double & expectation1,
                              const double & expectation2,
                              const double & covariance,
                              double       & s2) const;

  void CalculateCovariance(const double & corr,
                           const double & var1,
                           const double & var2,
                           double       & cov) const;

  NRLib::TrendStorage * mean_vp_;
  NRLib::TrendStorage * mean_vs_;
  NRLib::TrendStorage * mean_density_;
  NRLib::TrendStorage * variance_vp_;
  NRLib::TrendStorage * variance_vs_;
  NRLib::TrendStorage * variance_density_;
  NRLib::TrendStorage * correlation_vp_vs_;
  NRLib::TrendStorage * correlation_vp_density_;
  NRLib::TrendStorage * correlation_vs_density_;
};

#endif
