#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"


DistributionsDryRockStorage::DistributionsDryRockStorage()
{
}

DistributionsDryRockStorage::~DistributionsDryRockStorage()
{
}
//----------------------------------------------------------------------------------//

TabulatedVelocityDryRockStorage::TabulatedVelocityDryRockStorage(DistributionWithTrendStorage * vp,
                                                                 DistributionWithTrendStorage * vs,
                                                                 DistributionWithTrendStorage * density,
                                                                 DistributionWithTrendStorage * correlation_vp_vs,
                                                                 DistributionWithTrendStorage * correlation_vp_density,
                                                                 DistributionWithTrendStorage * correlation_vs_density,
                                                                 DistributionWithTrendStorage * total_porosity,
                                                                 std::string                    moduli)
: vp_(vp),
  vs_(vs),
  density_(density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

TabulatedVelocityDryRockStorage::~TabulatedVelocityDryRockStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
TabulatedVelocityDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                              std::string                             & errTxt) const
{
  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();

  if(dry_rock == NULL)
    errTxt += "The tabulated model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
TabulatedModulusDryRockStorage::TabulatedModulusDryRockStorage(DistributionWithTrendStorage * bulk_modulus,
                                                               DistributionWithTrendStorage * shear_modulus,
                                                               DistributionWithTrendStorage * density,
                                                               DistributionWithTrendStorage * correlation_bulk_shear,
                                                               DistributionWithTrendStorage * correlation_bulk_density,
                                                               DistributionWithTrendStorage * correlation_shear_density,
                                                               DistributionWithTrendStorage * total_porosity,
                                                               std::string                    moduli)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  correlation_bulk_shear_(correlation_bulk_shear),
  correlation_bulk_density_(correlation_bulk_density),
  correlation_shear_density_(correlation_shear_density),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

TabulatedModulusDryRockStorage::~TabulatedModulusDryRockStorage()
{
  delete bulk_modulus_;
  delete shear_modulus_;
  delete density_;
  delete correlation_bulk_shear_;
  delete correlation_bulk_density_;
  delete correlation_shear_density_;

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
TabulatedModulusDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                             std::string                             & errTxt) const
{
  DistributionsDryRock * dry_rock = NULL;

  if(dry_rock == NULL)
    errTxt += "The tabulated model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
ReussDryRockStorage::ReussDryRockStorage(std::vector<std::string>                    constituent_label,
                                         std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                                         DistributionWithTrendStorage *              total_porosity,
                                         std::string                                 moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

ReussDryRockStorage::~ReussDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
ReussDryRockStorage::GenerateDistributionsDryRock(const std::string                       & path,
                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                  std::string                             & errTxt) const
{
  std::vector<double> volume = getVolume(constituent_volume_fraction_, path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockReuss();

  if(dry_rock == NULL)
    errTxt += "The Reuss model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
VoigtDryRockStorage::VoigtDryRockStorage(std::vector<std::string>                    constituent_label,
                                         std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                                         DistributionWithTrendStorage *              total_porosity,
                                         std::string                                 moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

VoigtDryRockStorage::~VoigtDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
VoigtDryRockStorage::GenerateDistributionsDryRock(const std::string                       & path,
                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                  std::string                             & errTxt) const
{
  std::vector<double> volume = getVolume(constituent_volume_fraction_, path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockVoigt();

  if(dry_rock == NULL)
    errTxt += "The Voigt model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
HillDryRockStorage::HillDryRockStorage(std::vector<std::string>                    constituent_label,
                                       std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                                       DistributionWithTrendStorage *              total_porosity,
                                       std::string                                 moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

HillDryRockStorage::~HillDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
HillDryRockStorage::GenerateDistributionsDryRock(const std::string                       & path,
                                                 const std::vector<std::string>          & trend_cube_parameters,
                                                 const std::vector<std::vector<double> > & trend_cube_sampling,
                                                 std::string                             & errTxt) const
{
  std::vector<double> volume = getVolume(constituent_volume_fraction_, path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockHill();

  if(dry_rock == NULL)
    errTxt += "The Hill model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}

//----------------------------------------------------------------------------------//

DEMDryRockStorage::DEMDryRockStorage(std::string                                 host_label,
                                     DistributionWithTrendStorage *              host_volume_fraction,
                                     std::vector<std::string>                    inclusion_label,
                                     std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                                     std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio,
                                     DistributionWithTrendStorage *              total_porosity,
                                     std::string                                 moduli)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

DEMDryRockStorage::~DEMDryRockStorage()
{
  if(host_volume_fraction_->GetIsSheared() == false)
    delete host_volume_fraction_;

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_.size()); i++) {
    if(inclusion_volume_fraction_[i]->GetIsSheared() == false)
      delete inclusion_volume_fraction_[i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_.size()); i++) {
    if(inclusion_aspect_ratio_[i]->GetIsSheared() == false)
      delete inclusion_aspect_ratio_[i];
  }

  if(total_porosity_->GetIsSheared() == false)
    delete total_porosity_;
}

DistributionsDryRock *
DEMDryRockStorage::GenerateDistributionsDryRock(const std::string                       & path,
                                                const std::vector<std::string>          & trend_cube_parameters,
                                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                                std::string                             & errTxt) const
{
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());

  std::vector<DistributionWithTrendStorage *> volume_fractions(n_inclusions + 1);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=1; i<n_inclusions+1; i++)
    volume_fractions[i] = inclusion_volume_fraction_[i-1];

  std::vector<double> volume = getVolume(volume_fractions, path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();

  if(dry_rock == NULL)
    errTxt += "The DEM model has not been implemented yet for dry-rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(dry_rock);
}
