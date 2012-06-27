#ifndef RPLIB_DISTRIBUTIONS_SOLID_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_SOLID_STORAGE_HPP

#include "rplib/distributionssolid.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsSolidStorage {
public:
  DistributionsSolidStorage();

  virtual ~DistributionsSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                       & /*path*/,
                                                          const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                          const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                          std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//

class TabulatedSolidStorage : public DistributionsSolidStorage {
public:
  TabulatedSolidStorage(DistributionWithTrendStorage * vp,
                        DistributionWithTrendStorage * vs,
                        DistributionWithTrendStorage * density,
                        DistributionWithTrendStorage * correlation_vp_vs,
                        DistributionWithTrendStorage * correlation_vp_density,
                        DistributionWithTrendStorage * correlation_vs_density);

  virtual ~TabulatedSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                       & path,
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
