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

  virtual Rock                        * GenerateSample(const std::vector<double> & trend_params)          const = 0;

  void                                  GenerateWellSample(double                 corr,
                                                           std::vector<double> &  vp,
                                                           std::vector<double> &  vs,
                                                           std::vector<double> &  rho,
                                                           std::vector<double> &  trend_params)           const;

  Rock                                * EvolveSample(double       time,
                                                     const Rock & rock)                                   const;

  virtual bool                          HasDistribution()                                                 const = 0;

  virtual std::vector<bool>             HasTrend()                                                        const = 0;

  virtual bool                          GetIsOkForBounding()                                              const = 0;

  //NBNB fjellvoll the first two functions should probably survive the new format.
  std::vector<double>                   GetExpectations(const std::vector<double> & trend_params)         const    ;

  NRLib::Grid2D<double>                 GetCovariances(const std::vector<double> & trend_params)          const    ;

  virtual std::vector<double>           GetExpectation(const std::vector<double> & /*trend_params*/)      const = 0;

  virtual NRLib::Grid2D<double>         GetCovariance(const std::vector<double> & /*trend_params*/)       const = 0;

  virtual Rock                        * UpdateSample(double                      corr_param,
                                                     bool                        param_is_time,
                                                     const std::vector<double> & trend,
                                                     const Rock                * sample)    const = 0;

protected:
                                        //This function should be called last step in constructor
                                        //for all children classes.
  void                                  SetupExpectationAndCovariances(const std::vector<double> & s_min,
                                                                       const std::vector<double> & s_max);

  NRLib::Grid2D< std::vector<double> >  expectation_;
  NRLib::Grid2D< std::vector<double> >  covariance_;

  std::vector<double>                   s_min_;
  std::vector<double>                   s_max_;

  std::vector<double>                   alpha_;

};

#endif
