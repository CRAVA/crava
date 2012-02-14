#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/rockphysicsstorage.h"
#include "rplib/distributionsrockt0.h"
#include "rplib/multinormalwithtrend.h"
#include "rplib/multinormaldistributedrockt0.h"

RockPhysicsStorage::RockPhysicsStorage()
{
}

RockPhysicsStorage::~RockPhysicsStorage()
{
}

GaussianRockPhysicsStorage::GaussianRockPhysicsStorage(TrendStorage *mean_vp,
                                                       TrendStorage *mean_vs,
                                                       TrendStorage *mean_density,
                                                       TrendStorage *variance_vp,
                                                       TrendStorage *variance_vs,
                                                       TrendStorage *variance_density,
                                                       TrendStorage *correlation_vp_vs,
                                                       TrendStorage *correlation_vp_density,
                                                       TrendStorage *correlation_vs_density)
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
}

DistributionsRockT0 *
GaussianRockPhysicsStorage::GenerateRockPhysics(const std::string & path,
                                                std::string       & errTxt) const
{
  Trend * mean_vp_trend                = mean_vp_               ->GenerateTrend(path,errTxt);
  Trend * mean_vs_trend                = mean_vs_               ->GenerateTrend(path,errTxt);
  Trend * mean_density_trend           = mean_density_          ->GenerateTrend(path,errTxt);
  Trend * variance_vp_trend            = variance_vp_           ->GenerateTrend(path,errTxt);
  Trend * variance_vs_trend            = variance_vs_           ->GenerateTrend(path,errTxt);
  Trend * variance_density_trend       = variance_density_      ->GenerateTrend(path,errTxt);
  Trend * correlation_vp_vs_trend      = correlation_vp_vs_     ->GenerateTrend(path,errTxt);
  Trend * correlation_vp_density_trend = correlation_vp_density_->GenerateTrend(path,errTxt);
  Trend * correlation_vs_density_trend = correlation_vs_density_->GenerateTrend(path,errTxt);

  NRLib::Normal vp0, vs0, density0;
  NRLib::Grid2D<Trend*> cov(3,3,NULL);

  cov(0,0) = variance_vp_trend;
  cov(0,1) = correlation_vp_vs_trend;
  cov(0,2) = correlation_vp_density_trend;
  cov(1,0) = correlation_vp_vs_trend;
  cov(1,1) = variance_vs_trend;
  cov(1,2) = correlation_vs_density_trend;
  cov(2,0) = correlation_vp_density_trend;
  cov(2,1) = correlation_vs_density_trend;
  cov(2,2) = variance_density_trend;

  MultiNormalWithTrend multi(vp0,
                             vs0,
                             density0,
                             mean_vp_trend,
                             mean_vs_trend,
                             mean_density_trend,
                             cov);

  DistributionsRockT0 * rock = new MultiNormalDistributedRockT0(multi);

  return(rock);
}
