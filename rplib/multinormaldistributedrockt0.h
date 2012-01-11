#ifndef RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H
#define RPLIB_MULTINORMALDISTRIBUTEDROCKT0_H

#include <vector>

class Rock;
class MultiNormalWithTrend;

class MultiNormalDistributedRockT0 {
public:

  MultiNormalDistributedRockT0(const MultiNormalWithTrend& mult_normal_distr);

  virtual ~MultiNormalDistributedRockT0();

  virtual Rock  * GenerateSample(const std::vector<double> & trend_params) const;

private:
  const MultiNormalWithTrend& mult_normal_distr_;


};

#endif
