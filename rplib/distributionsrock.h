#ifndef DISTRIBUTIONSROCK_H
#define DISTRIBUTIONSROCK_H

#include "rplib/rock.h"
#include <nrlib/flens/nrlib_flens.hpp>
#include "nrlib/grid/grid2d.hpp"


// Abstract class for holding all t = 0 distribution functions for rock physics parameters and saturation.
// One derived class for each rock physics model, the latter specified in a parallel, derived Rock class.
// The class must be able to produce an object of the specific Rock class.
class DistributionsRock {
public:

                                        DistributionsRock(){}

  virtual                               ~DistributionsRock(){}

  virtual DistributionsRock           * Clone()                                                           const = 0;

  virtual Rock                        * GenerateSample(const std::vector<double> & trend_params)          const = 0;

  void                                  GenerateWellSample(double                 corr,
                                                           std::vector<double> &  vp,
                                                           std::vector<double> &  vs,
                                                           std::vector<double> &  rho,
                                                           std::vector<double> &  trend_params)           const;

  Rock                                * EvolveSample(double       time,
                                                     const Rock & rock)                                   const;

  std::vector<double>                   GetExpectation(const std::vector<double> & trend_params)         const;

  NRLib::Grid2D<double>                 GetCovariance(const std::vector<double> & trend_params)          const;

  virtual bool                          HasDistribution()                                                 const = 0;

  virtual std::vector<bool>             HasTrend()                                                        const = 0;

  virtual bool                          GetIsOkForBounding()                                              const = 0;

  virtual Rock                        * UpdateSample(double                      corr_param,
                                                     bool                        param_is_time,
                                                     const std::vector<double> & trend,
                                                     const Rock                * sample)    const = 0;

protected:
                                        //This function should be called last step in constructor
                                        //for all children classes.
  void                                  SetupExpectationAndCovariances(NRLib::Grid2D<std::vector<double> >   & expectation,
                                                                       NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                                                       std::vector<double>                   & tabulated_s0,
                                                                       std::vector<double>                   & tabulated_s1,
                                                                       const std::vector<double>             & s_min,
                                                                       const std::vector<double>             & s_max);

  void                                  FindTabulatedTrendParams(std::vector<double>       & tabulated_s0,
                                                                 std::vector<double>       & tabulated_s1,
                                                                 const std::vector<bool>   & has_trend,
                                                                 const std::vector<double> & s_min,
                                                                 const std::vector<double> & s_max,
                                                                 const size_t                n);

  void                                  SetupTrendMesh(NRLib::Grid2D<std::vector<double> > & trend_params,
                                                       const std::vector<double>           & t1,
                                                       const std::vector<double>           & t2);

  void                                  SampleTrendValues(std::vector<double> & s,
                                                          const double        & s_min,
                                                          const double        & s_max);

  double                                FindExpectation(const std::vector<double> & p);

  double                                FindCovariance(const std::vector<double> & p,
                                                       const double                mup,
                                                       const std::vector<double> & q,
                                                       const double                muq);

  double                                FindInterpolationStartIndex(const std::vector<double> & tabulated_s,
                                                                    const double                s) const;

  void                                  FindInterpolationWeights(double       & w00,
                                                                 double       & w10,
                                                                 double       & w01,
                                                                 double       & w11,
                                                                 const double   di,
                                                                 const double   dj) const;

  void                                  InterpolateExpectation(std::vector<double>                       & mean,
                                                               const NRLib::Grid2D<std::vector<double> > & expectation,
                                                               const double                                w00,
                                                               const double                                w10,
                                                               const double                                w01,
                                                               const double                                w11,
                                                               const double                                di,
                                                               const double                                dj,
                                                               const size_t                                m,
                                                               const size_t                                n,
                                                               const size_t                                p) const;

  void                                  InterpolateCovariance(NRLib::Grid2D<double>                       & cov,
                                                              const NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                                              const double                                  w00,
                                                              const double                                  w10,
                                                              const double                                  w01,
                                                              const double                                  w11,
                                                              const double                                  di,
                                                              const double                                  dj,
                                                              const size_t                                  m,
                                                              const size_t                                  n,
                                                              const size_t                                  p,
                                                              const size_t                                  q) const;

  NRLib::Grid2D<std::vector<double> >   expectation_;
  NRLib::Grid2D<NRLib::Grid2D<double> > covariance_;

  std::vector<double>                   tabulated_s0_;
  std::vector<double>                   tabulated_s1_;

  std::vector<double>                   alpha_;

  std::vector<double>                   s_min_;
  std::vector<double>                   s_max_;
};

#endif
