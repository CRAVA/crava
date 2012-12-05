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

TabulatedVelocityDryRockStorage::TabulatedVelocityDryRockStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                                 std::vector<DistributionWithTrendStorage *> vs,
                                                                 std::vector<DistributionWithTrendStorage *> density,
                                                                 std::vector<double>                         correlation_vp_vs,
                                                                 std::vector<double>                         correlation_vp_density,
                                                                 std::vector<double>                         correlation_vs_density,
                                                                 std::vector<DistributionWithTrendStorage *> total_porosity,
                                                                 std::string                                 moduli)
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
  if(vp_[0]->GetIsShared() == false)
    delete vp_[0];

  if(vs_[0]->GetIsShared() == false)
    delete vs_[0];

  if(density_[0]->GetIsShared() == false)
    delete density_[0];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
TabulatedVelocityDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                              std::string                             & errTxt) const
{
  std::vector<DistributionsDryRock *> dry_rock(1, NULL);

  if(dry_rock[0] == NULL)
    errTxt += "The tabulated model has not been implemented yet for dry-rocks\n";

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
TabulatedModulusDryRockStorage::TabulatedModulusDryRockStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                               std::vector<DistributionWithTrendStorage *> shear_modulus,
                                                               std::vector<DistributionWithTrendStorage *> density,
                                                               std::vector<double>                         correlation_bulk_shear,
                                                               std::vector<double>                         correlation_bulk_density,
                                                               std::vector<double>                         correlation_shear_density,
                                                               std::vector<DistributionWithTrendStorage *> total_porosity,
                                                               std::string                                 moduli)
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

 if(bulk_modulus_[0]->GetIsShared() == false)
    delete bulk_modulus_[0];

  if(shear_modulus_[0]->GetIsShared() == false)
    delete shear_modulus_[0];

  if(density_[0]->GetIsShared() == false)
    delete density_[0];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
TabulatedModulusDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                             std::string                             & errTxt) const
{
  std::vector<DistributionsDryRock *> dry_rock(1, NULL);

  if(dry_rock[0] == NULL)
    errTxt += "The tabulated model has not been implemented yet for dry-rocks\n";

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
ReussDryRockStorage::ReussDryRockStorage(std::vector<std::string >                                 constituent_label,
                                         std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                                         std::vector<DistributionWithTrendStorage *>               total_porosity,
                                         std::string                                               moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

ReussDryRockStorage::~ReussDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++)
    delete constituent_volume_fraction_[0][i];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
ReussDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                  const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                  const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                  std::string                             & errTxt) const
{
  std::vector<DistributionsDryRock *> dry_rock(1, NULL);

  if(dry_rock[0] == NULL)
    errTxt += "The Reuss model has not been implemented yet for dry-rocks\n";

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
VoigtDryRockStorage::VoigtDryRockStorage(std::vector<std::string>                                  constituent_label,
                                         std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                                         std::vector<DistributionWithTrendStorage *>               total_porosity,
                                         std::string                                               moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

VoigtDryRockStorage::~VoigtDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++)
    delete constituent_volume_fraction_[0][i];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
VoigtDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                  const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                  const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                  std::string                             & errTxt) const
{
  std::vector<DistributionsDryRock *> dry_rock(1, NULL);

  if(dry_rock[0] == NULL)
    errTxt += "The Voigt model has not been implemented yet for dry-rocks\n";

  return(dry_rock);
}

//----------------------------------------------------------------------------------//
HillDryRockStorage::HillDryRockStorage(std::vector<std::string>                                  constituent_label,
                                       std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                                       std::vector<DistributionWithTrendStorage *>               total_porosity,
                                       std::string                                               moduli)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction),
  total_porosity_(total_porosity),
  mineral_moduli_(moduli)
{
}

HillDryRockStorage::~HillDryRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++)
    delete constituent_volume_fraction_[0][i];

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
HillDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                 const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                 const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                 std::string                             & errTxt) const
{
  std::vector<DistributionsDryRock *> dry_rock(1, NULL);
  if(dry_rock[0] == NULL)
    errTxt += "The Hill model has not been implemented yet for dry-rocks\n";

  return(dry_rock);
}

//----------------------------------------------------------------------------------//

DEMDryRockStorage::DEMDryRockStorage(std::string                                               host_label,
                                     std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                                     std::vector<std::string>                                  inclusion_label,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio,
                                     std::vector<DistributionWithTrendStorage *>               total_porosity,
                                     std::string                                               moduli)
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
  if(host_volume_fraction_[0]->GetIsShared() == false)
    delete host_volume_fraction_[0];

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_[0].size()); i++) {
    if(inclusion_volume_fraction_[0][i]->GetIsShared() == false)
      delete inclusion_volume_fraction_[0][i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_[0].size()); i++) {
    if(inclusion_aspect_ratio_[0][i]->GetIsShared() == false)
      delete inclusion_aspect_ratio_[0][i];
  }

  if(total_porosity_[0]->GetIsShared() == false)
    delete total_porosity_[0];
}

std::vector<DistributionsDryRock *>
DEMDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                std::string                             & errTxt) const
{

  std::vector<DistributionsDryRock *> dry_rock(1, NULL);

  if(dry_rock[0] == NULL)
    errTxt += "The DEM model has not been implemented yet for dry-rocks\n";
  return(dry_rock);
}
