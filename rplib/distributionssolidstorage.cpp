#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionssolidtabulatedvelocity.h"
#include "rplib/distributionssolidinclusion.h"
#include "rplib/distributionwithtrendstorage.h"


DistributionsSolidStorage::DistributionsSolidStorage()
{
}

DistributionsSolidStorage::~DistributionsSolidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedVelocitySolidStorage::TabulatedVelocitySolidStorage(DistributionWithTrendStorage * vp,
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

TabulatedVelocitySolidStorage::~TabulatedVelocitySolidStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsSolid *
TabulatedVelocitySolidStorage::GenerateDistributionsSolid(const std::string                       & path,
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

  DistributionsSolid * solid = new DistributionsSolidTabulatedVelocity(vp_dist_with_trend,
                                                                       vs_dist_with_trend,
                                                                       density_dist_with_trend,
                                                                       corr_vp_vs_dist_with_trend,
                                                                       corr_vp_density_dist_with_trend,
                                                                       corr_vs_density_dist_with_trend);

  return(solid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusSolidStorage::TabulatedModulusSolidStorage(DistributionWithTrendStorage * bulk_modulus,
                                                           DistributionWithTrendStorage * shear_modulus,
                                                           DistributionWithTrendStorage * density,
                                                           DistributionWithTrendStorage * correlation_bulk_shear,
                                                           DistributionWithTrendStorage * correlation_bulk_density,
                                                           DistributionWithTrendStorage * correlation_shear_density)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  correlation_bulk_shear_(correlation_bulk_shear),
  correlation_bulk_density_(correlation_bulk_density),
  correlation_shear_density_(correlation_shear_density)
{
}

TabulatedModulusSolidStorage::~TabulatedModulusSolidStorage()
{
  delete bulk_modulus_;
  delete shear_modulus_;
  delete density_;
  delete correlation_bulk_shear_;
  delete correlation_bulk_density_;
  delete correlation_shear_density_;
}

DistributionsSolid *
TabulatedModulusSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                                         const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                         const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                         std::string                             & /*errTxt*/) const
{
  /*const DistributionWithTrend * bulk_dist_with_trend               = bulk_modulus_              ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * shear_dist_with_trend              = shear_modulus_             ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend            = density_                   ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_bulk_shear_dist_with_trend    = correlation_bulk_shear_    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_bulk_density_dist_with_trend  = correlation_bulk_density_  ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_shear_density_dist_with_trend = correlation_shear_density_ ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsSolid * solid = new DistributionsModulusSolidTabulated(bulk_dist_with_trend,
                                                                      shear_dist_with_trend,
                                                                      density_dist_with_trend,
                                                                      corr_bulk_shear_dist_with_trend,
                                                                      corr_bulk_density_dist_with_trend,
                                                                      corr_shear_density_dist_with_trend);*/
  DistributionsSolid * solid = NULL;

  return(solid);
}

//----------------------------------------------------------------------------------//
ReussSolidStorage::ReussSolidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussSolidStorage::~ReussSolidStorage()
{
}

DistributionsSolid *
ReussSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsSolid * solid = NULL; //new DistributionsSolidReuss();
  return(solid);
}

//----------------------------------------------------------------------------------//
VoigtSolidStorage::VoigtSolidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtSolidStorage::~VoigtSolidStorage()
{
}

DistributionsSolid *
VoigtSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsSolid * solid = NULL; //new DistributionsSolidVoigt();
  return(solid);
}

//----------------------------------------------------------------------------------//
HillSolidStorage::HillSolidStorage(std::vector<std::string>                    constituent_label,
                                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillSolidStorage::~HillSolidStorage()
{
}

DistributionsSolid *
HillSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                             std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsSolid * solid = NULL; //new DistributionsSolidHill();
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
