#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsfluidtabulated.h"
#include "rplib/distributionwithtrendstorage.h"


DistributionsFluidStorage::DistributionsFluidStorage()
{
}

DistributionsFluidStorage::~DistributionsFluidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedFluidStorage::TabulatedFluidStorage(DistributionWithTrendStorage * vp,
                                             DistributionWithTrendStorage * density,
                                             DistributionWithTrendStorage * correlation_vp_density)
: vp_(vp),
  density_(density),
  correlation_vp_density_(correlation_vp_density)
{
}

TabulatedFluidStorage::~TabulatedFluidStorage()
{
  delete vp_;
  delete density_;
  delete correlation_vp_density_;
}

DistributionsFluid *
TabulatedFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                  std::string                             & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_vp_density_dist_with_trend = correlation_vp_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = new DistributionsFluidTabulated(vp_dist_with_trend,
                                                               density_dist_with_trend,
                                                               corr_vp_density_dist_with_trend);

  return(fluid);
}

/*//----------------------------------------------------------------------------------//
TabulatedFluidStorage::TabulatedFluidStorage(DistributionWithTrendStorage * vp,
                                             DistributionWithTrendStorage * vs,
                                             DistributionWithTrendStorage * density,
                                             DistributionWithTrendStorage * correlation_vp_vs,
                                             DistributionWithTrendStorage * correlation_vp_density,
                                             DistributionWithTrendStorage * correlation_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density)
{
}

TabulatedFluidStorage::~TabulatedFluidStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsFluid *
TabulatedFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                  std::string                             & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * vs_dist_with_trend              = vs_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_vp_vs_dist_with_trend      = correlation_vp_vs_     ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_vp_density_dist_with_trend = correlation_vp_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_vs_density_dist_with_trend = correlation_vs_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = new DistributionsFluidTabulated(vp_dist_with_trend,
                                                               vs_dist_with_trend,
                                                               density_dist_with_trend,
                                                               corr_vp_vs_dist_with_trend,
                                                               corr_vp_density_dist_with_trend,
                                                               corr_vs_density_dist_with_trend);

  return(fluid);
}*/

