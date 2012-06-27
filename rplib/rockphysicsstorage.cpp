#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/rockphysicsstorage.h"
#include "rplib/distributionsrock.h"
#include "rplib/trinormalwith2dtrend.h"
#include "rplib/multinormaldistributedrock.h"
#include "rplib/distributionwithtrendstorage.h"


RockPhysicsStorage::RockPhysicsStorage()
{
}

RockPhysicsStorage::~RockPhysicsStorage()
{
}

//----------------------------------------------------------------------------------//
GaussianRockPhysicsStorage::GaussianRockPhysicsStorage(NRLib::TrendStorage * mean_vp,
                                                       NRLib::TrendStorage * mean_vs,
                                                       NRLib::TrendStorage * mean_density,
                                                       NRLib::TrendStorage * variance_vp,
                                                       NRLib::TrendStorage * variance_vs,
                                                       NRLib::TrendStorage * variance_density,
                                                       NRLib::TrendStorage * correlation_vp_vs,
                                                       NRLib::TrendStorage * correlation_vp_density,
                                                       NRLib::TrendStorage * correlation_vs_density)
: mean_vp_(mean_vp),
  mean_vs_(mean_vs),
  mean_density_(mean_density),
  variance_vp_(variance_vp),
  variance_vs_(variance_vs),
  variance_density_(variance_density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density)
{
}

GaussianRockPhysicsStorage::~GaussianRockPhysicsStorage()
{
  delete mean_vp_;
  delete mean_vs_;
  delete mean_density_;
  delete variance_vp_;
  delete variance_vs_;
  delete variance_density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsRock *
GaussianRockPhysicsStorage::GenerateRockPhysics(const std::string                       & path,
                                                const std::vector<std::string>          & trend_cube_parameters,
                                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                                std::string                             & errTxt) const
{
  NRLib::Trend * mean_vp_trend                = mean_vp_               ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * mean_vs_trend                = mean_vs_               ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * mean_density_trend           = mean_density_          ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * variance_vp_trend            = variance_vp_           ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * variance_vs_trend            = variance_vs_           ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * variance_density_trend       = variance_density_      ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * correlation_vp_vs_trend      = correlation_vp_vs_     ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * correlation_vp_density_trend = correlation_vp_density_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * correlation_vs_density_trend = correlation_vs_density_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

  DistributionsRock * rock = new MultiNormalDistributedRock(mean_vp_trend,
                                                            mean_vs_trend,
                                                            mean_density_trend,
                                                            variance_vp_trend,
                                                            variance_vs_trend,
                                                            variance_density_trend,
                                                            correlation_vp_vs_trend,
                                                            correlation_vp_density_trend,
                                                            correlation_vs_density_trend);


  delete mean_vp_trend;
  delete mean_vs_trend;
  delete mean_density_trend;
  delete variance_vp_trend;
  delete variance_vs_trend;
  delete variance_density_trend;
  delete correlation_vp_vs_trend;
  delete correlation_vp_density_trend;
  delete correlation_vs_density_trend;

  return(rock);
}

/*DistributionsFluidTabulated::DistributionsFluidTabulated(const DistributionWithTrend * vp,
                                                         const DistributionWithTrend * vs,
                                                         const DistributionWithTrend * density,
                                                         const DistributionWithTrend * corr_vp_vs,
                                                         const DistributionWithTrend * corr_vp_density,
                                                         const DistributionWithTrend * corr_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  corr_vp_vs_(corr_vp_vs),
  corr_vp_density_(corr_vp_density),
  corr_vs_density_(corr_vs_density)
{
}
*/
