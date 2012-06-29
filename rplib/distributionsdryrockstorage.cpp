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

TabulatedDryRockStorage::TabulatedDryRockStorage(DistributionWithTrendStorage * vp,
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

TabulatedDryRockStorage::~TabulatedDryRockStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsDryRock *
TabulatedDryRockStorage::GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                      const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                      const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                      std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();
  return(dry_rock);
}

//----------------------------------------------------------------------------------//

DEMDryRockStorage::DEMDryRockStorage(std::string                                 host_label,
                                     DistributionWithTrendStorage *              host_volume_fraction,
                                     std::vector<std::string>                    inclusion_label,
                                     std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                                     std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio ,
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
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsDryRock * dry_rock = NULL; //new DistributionsDryRockInclusion();
  return(dry_rock);
}
