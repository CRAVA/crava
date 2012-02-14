#ifndef RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H
#define RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>
#include <rplib/distributionsrockt0.h>

class Rock;
class MultiNormalWithTrend;

class MultiNormalDistributedRockT0 : public DistributionsRockT0 {
public:

  MultiNormalDistributedRockT0(const MultiNormalWithTrend& mult_normal_distr);

  virtual ~MultiNormalDistributedRockT0();

  virtual Rock  * GenerateSample(const std::vector<double> & trend_params) const;

  virtual std::vector<double> GetExpectation(const std::vector<double> & trend_params) const;

  virtual NRLib::Matrix       GetCovariance(const std::vector<double> & trend_params) const;


private:
  const MultiNormalWithTrend& mult_normal_distr_;


};

#endif
