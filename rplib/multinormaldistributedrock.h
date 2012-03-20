#ifndef RPLIB_MULTINORMALDISTRIBUTEDROCK_H
#define RPLIB_MULTINORMALDISTRIBUTEDROCK_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <rplib/distributionsrock.h>
#include <rplib/trinormalwith2dtrend.h>

class Rock;
class Pdf3D;

class MultiNormalDistributedRock : public DistributionsRock {
public:

  MultiNormalDistributedRock(NRLib::Trend * mean_vp,
                             NRLib::Trend * mean_vs,
                             NRLib::Trend * mean_density,
                             NRLib::Trend * variance_vp,
                             NRLib::Trend * variance_vs,
                             NRLib::Trend * variance_density,
                             NRLib::Trend * correlation_vp_vs,
                             NRLib::Trend * correlation_vp_density,
                             NRLib::Trend * correlation_vs_density);

  virtual ~MultiNormalDistributedRock();

  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual std::vector<double>        GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double>      GetCovariance(const std::vector<double> & trend_params) const;

  virtual Pdf3D                    * GeneratePdf() const;

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


  TriNormalWith2DTrend * mult_normal_distr_;


};

#endif
