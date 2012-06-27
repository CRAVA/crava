#ifndef RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP

#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsFluidStorage {
public:
  DistributionsFluidStorage();

  virtual ~DistributionsFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & /*path*/,
                                                          const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                          const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                          std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//

class TabulatedFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedFluidStorage(DistributionWithTrendStorage * vp,
                        DistributionWithTrendStorage * density,
                        DistributionWithTrendStorage * correlation_vp_density);

  virtual ~TabulatedFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_vp_density_;
};
/*class TabulatedFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedFluidStorage(DistributionWithTrendStorage * vp,
                        DistributionWithTrendStorage * vs,
                        DistributionWithTrendStorage * density,
                        DistributionWithTrendStorage * correlation_vp_vs,
                        DistributionWithTrendStorage * correlation_vp_density,
                        DistributionWithTrendStorage * correlation_vs_density);

  virtual ~TabulatedFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
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
};*/
#endif
