#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"


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
}

DistributionsDryRock *
TabulatedVelocityDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                              std::string                             & /*errTxt*/) const
{
  //Gj�r sjekk p� om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass p� at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();
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
}

DistributionsDryRock *
TabulatedModulusDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                             std::string                             & /*errTxt*/) const
{
  DistributionsDryRock * dry_rock = NULL;

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
}

DistributionsDryRock *
ReussDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                  const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                  const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                  std::string                             & /*errTxt*/) const
{
  //Gj�r sjekk p� om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass p� at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockReuss();
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
}

DistributionsDryRock *
VoigtDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                  const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                  const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                  std::string                             & /*errTxt*/) const
{
  //Gj�r sjekk p� om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass p� at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockVoigt();
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
}

DistributionsDryRock *
HillDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                 const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                 const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                 std::string                             & /*errTxt*/) const
{
  //Gj�r sjekk p� om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass p� at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockHill();
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
}

DistributionsDryRock *
DEMDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                std::string                             & /*errTxt*/) const
{
  //Gj�r sjekk p� om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass p� at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();
  return(dry_rock);
}
