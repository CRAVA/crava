#ifndef RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H
#define RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>
#include <nrlib/grid/grid2d.hpp>
#include <rplib/distributionsrockt0.h>
#include <rplib/trinormalwith2dtrend.h>

class Rock;
class Pdf3D;

class MultiNormalDistributedRockT0 : public DistributionsRockT0 {
public:

  MultiNormalDistributedRockT0(const TriNormalWith2DTrend & mult_normal_distr);

  virtual ~MultiNormalDistributedRockT0();

  virtual Rock  * GenerateSample(const std::vector<double> & trend_params) const;

  virtual std::vector<double> GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Grid2D<double> GetCovariance(const std::vector<double> & trend_params) const;

  virtual Pdf3D * GeneratePdf(void) const;

private:
  TriNormalWith2DTrend mult_normal_distr_;


};

#endif
