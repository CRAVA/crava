#ifndef DISTRIBUTIONSROCK_H
#define DISTRIBUTIONSROCK_H

#include "nrlib/grid/grid2d.hpp"
#include "rplib/distributionwithtrend.h"

class Rock;

// Abstract class for holding all t = 0 distribution functions for rock physics parameters and saturation.
// One derived class for each rock physics model, the latter specified in a parallel, derived Rock class.
// The class must be able to produce an object of the specific Rock class.
class DistributionsRock {
public:

                                        DistributionsRock(){resampling_level_ = DistributionWithTrend::SingleSample;}

  virtual                               ~DistributionsRock(){}

  virtual DistributionsRock           * Clone()                                                           const = 0;

  Rock                                * GenerateSample(const std::vector<double> & trend_params);
  Rock                                * GenerateSampleAndReservoirVariables(const std::vector<double> & trend_params, std::vector<double> &resVar );

  void                                  GenerateWellSample(double                 corr,
                                                           std::vector<double> &  vp,
                                                           std::vector<double> &  vs,
                                                           std::vector<double> &  rho,
                                                           std::vector<double> &  trend_params);

  Rock                                * EvolveSample(double       time,
                                                     const Rock & rock);
  
  Rock                                * EvolveSampleAndReservoirVaribles(double             time,
                                                                         const              Rock & rock,
                                                                         std::vector<double> &resVar );

  std::vector<double>                   GetLogExpectation(const std::vector<double> & trend_params)       const;

  NRLib::Grid2D<double>                 GetLogCovariance(const std::vector<double> & trend_params)        const;

  const std::vector<double>           & GetMeanLogExpectation()                                           const { return mean_log_expectation_ ;}

  const NRLib::Grid2D<double>         & GetMeanLogCovariance()                                            const { return mean_log_covariance_  ;}

  virtual bool                          HasDistribution()                                                 const = 0;

  virtual std::vector<bool>             HasTrend()                                                        const = 0;

  virtual bool                          GetIsOkForBounding()                                              const = 0;

  virtual Rock                        * UpdateSample(double                      corr_param,
                                                     bool                        param_is_time,
                                                     const std::vector<double> & trend,
                                                     const Rock                * sample) = 0;

  //Top level objects (those accessed from outside the rock physics model) need more parameters set, so call this.
  void                                  CompleteTopLevelObject(std::vector<DistributionWithTrend *> res_var);

  void                                  SetResamplingLevel(int level) {resampling_level_ = level;}

protected:
  //Since there is a common start of generate sample here, that is the public function.
  //The public function then calls this overloaded function to get the specific object.
  virtual Rock                        * GenerateSamplePrivate(const std::vector<double> & trend_params) = 0;
                                        //This function should be called last step in constructor
                                        //for all children classes.
  void                                  SetupExpectationAndCovariances();
                                                                     /*NRLib::Grid2D<std::vector<double> >   & expectation,
                                                                       NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                                                       std::vector<double>                   & tabulated_s0,
                                                                       std::vector<double>                   & tabulated_s1,
                                                                       const std::vector<double>             & s_min,
                                                                       const std::vector<double>             & s_max);*/

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

  void                                  CheckOrResetS(double                    & s,
                                                      const std::vector<double> & tabulated_s) const;

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
                                                               const size_t                                i0,
                                                               const size_t                                j0,
                                                               const bool                                  do10,
                                                               const bool                                  do01,
                                                               const bool                                  do11,
                                                               const size_t                                p) const;

  void                                  InterpolateCovariance(NRLib::Grid2D<double>                       & cov,
                                                              const NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                                              const double                                  w00,
                                                              const double                                  w10,
                                                              const double                                  w01,
                                                              const double                                  w11,
                                                              const size_t                                  i0,
                                                              const size_t                                  j0,
                                                              const bool                                    do10,
                                                              const bool                                    do01,
                                                              const bool                                    do11,
                                                              const size_t                                  p,
                                                              const size_t                                  q) const;

  NRLib::Grid2D<std::vector<double> >   expectation_;       // Expectation of log(vp,vs,density)
  NRLib::Grid2D<NRLib::Grid2D<double> > covariance_;        // Covariance of log(vp,vs,density)

  std::vector<double>                   mean_log_expectation_;
  NRLib::Grid2D<double>                 mean_log_covariance_;

  std::vector<double>                   tabulated_s0_;      // Tabulated trend values {s0_min, ... , s0_max}
  std::vector<double>                   tabulated_s1_;      // Tabulated trend values {s1_min, ... , s1_max}

  std::vector<double>                   alpha_;

  std::vector<double>                   s_min_;
  std::vector<double>                   s_max_;

  std::vector<DistributionWithTrend *> reservoir_variables_;    //Note: These are stored elsewhere, no delete responsibility
  int                                  resampling_level_;
};

#endif
