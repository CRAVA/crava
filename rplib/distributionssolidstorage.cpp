#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionssolidtabulated.h"
#include "rplib/distributionssolidinclusion.h"
#include "rplib/distributionwithtrendstorage.h"


DistributionsSolidStorage::DistributionsSolidStorage()
{
}

DistributionsSolidStorage::~DistributionsSolidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedSolidStorage::TabulatedSolidStorage(DistributionWithTrendStorage * vp,
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

TabulatedSolidStorage::~TabulatedSolidStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsSolid *
TabulatedSolidStorage::GenerateDistributionsSolid(const std::string                       & path,
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

  DistributionsSolid * solid = new DistributionsSolidTabulated(vp_dist_with_trend,
                                                               vs_dist_with_trend,
                                                               density_dist_with_trend,
                                                               corr_vp_vs_dist_with_trend,
                                                               corr_vp_density_dist_with_trend,
                                                               corr_vs_density_dist_with_trend);

  return(solid);
}

//----------------------------------------------------------------------------------//

DEMSolidStorage::DEMSolidStorage(std::string                                 host_label,
                                 DistributionWithTrendStorage *              host_volume_fraction,
                                 std::vector<std::string>                    inclusion_label,
                                 std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                                 std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMSolidStorage::~DEMSolidStorage()
{
}

DistributionsSolid *
DEMSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                            const std::vector<std::string>          & /*trend_cube_parameters*/,
                                            const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                            std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsSolid * solid = NULL; //new DistributionsSolidInclusion();
  return(solid);
}
