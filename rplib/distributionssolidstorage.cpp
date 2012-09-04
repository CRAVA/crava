#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionssolidtabulatedvelocity.h"
#include "rplib/distributionssolidinclusion.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"


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

  if(solid == NULL)
    errTxt += "The tabulated model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
                                                         std::string                             & errTxt) const
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

  if(solid == NULL)
    errTxt += "The tabulated model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
ReussSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsSolid * solid = NULL; //new DistributionsSolidReuss();

  if(solid == NULL)
    errTxt += "The Reuss model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
VoigtSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsSolid * solid = NULL; //new DistributionsSolidVoigt();

  if(solid == NULL)
    errTxt += "The Voigt model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
HillSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                             std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsSolid * solid = NULL; //new DistributionsSolidHill();

  if(solid == NULL)
    errTxt += "The Hill model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(solid);
}

//----------------------------------------------------------------------------------//

DEMSolidStorage::DEMSolidStorage(std::string                                 host_label,
                                 DistributionWithTrendStorage *              host_volume_fraction,
                                 DistributionWithTrendStorage *              host_aspect_ratio,
                                 std::vector<std::string>                    inclusion_label,
                                 std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                                 std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  host_aspect_ratio_(host_aspect_ratio),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMSolidStorage::~DEMSolidStorage()
{
  if(host_volume_fraction_->GetIsShared() == false)
    delete host_volume_fraction_;

  if(host_aspect_ratio_->GetIsShared() == false)
    delete host_aspect_ratio_;

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_.size()); i++) {
    if(inclusion_volume_fraction_[i]->GetIsShared() == false)
      delete inclusion_volume_fraction_[i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_.size()); i++) {
    if(inclusion_aspect_ratio_[i]->GetIsShared() == false)
      delete inclusion_aspect_ratio_[i];
  }
}

DistributionsSolid *
DEMSolidStorage::GenerateDistributionsSolid(const std::string                       & /*path*/,
                                            const std::vector<std::string>          & /*trend_cube_parameters*/,
                                            const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                            std::string                             & errTxt) const
{
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());

  std::vector<DistributionWithTrendStorage *> volume_fractions(n_inclusions + 1);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=1; i<n_inclusions+1; i++)
    volume_fractions[i] = inclusion_volume_fraction_[i-1];

  CheckVolumeConsistency(volume_fractions, errTxt);

  DistributionsSolid * solid = NULL; //new DistributionsSolidInclusion();

  if(solid == NULL)
    errTxt += "The DEM model has not been implemented yet for solids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(solid);
}
