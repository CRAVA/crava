#ifndef RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP

#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsRockStorage {
public:
  DistributionsRockStorage();

  virtual ~DistributionsRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & /*path*/,
                                                        const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                        const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                        std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//

class TabulatedRockStorage : public DistributionsRockStorage {
public:
  TabulatedRockStorage(DistributionWithTrendStorage * vp,
                       DistributionWithTrendStorage * vs,
                       DistributionWithTrendStorage * density,
                       DistributionWithTrendStorage * correlation_vp_vs,
                       DistributionWithTrendStorage * correlation_vp_density,
                       DistributionWithTrendStorage * correlation_vs_density);

  virtual ~TabulatedRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * vs_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_vp_vs_;
  DistributionWithTrendStorage * correlation_vp_density_;
  DistributionWithTrendStorage * correlation_vs_density_;
};

#endif
