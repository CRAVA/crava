#ifndef DISTRIBUTIONSROCKT0_H
#define DISTRIBUTIONSROCKT0_H

#include "rplib/rock.h"
#include "rplib/pdf3d.h"
#include <nrlib/flens/nrlib_flens.hpp>
#include "nrlib/grid/grid2d.hpp"


// Abstract class for holding all t = 0 distribution functions for rock physics parameters and saturation.
// One derived class for each rock physics model, the latter specified in a parallel, derived Rock class.
// The class must be able to produce an object of the specific Rock class.
class DistributionsRockT0 {
public:

  DistributionsRockT0(){}

  virtual ~DistributionsRockT0(){}

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                * GenerateSample(const std::vector<double> & trend_params) const = 0;

  virtual std::vector<double>   GetExpectation(const std::vector<double> & trend_params) const = 0;

  virtual NRLib::Grid2D<double> GetCovariance(const std::vector<double> & trend_params)  const = 0;

  virtual Pdf3D               * GeneratePdf(void)                                        const = 0;


};

#endif
