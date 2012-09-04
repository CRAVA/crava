#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsfluidtabulatedvelocity.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"


DistributionsFluidStorage::DistributionsFluidStorage()
{
}

DistributionsFluidStorage::~DistributionsFluidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedVelocityFluidStorage::TabulatedVelocityFluidStorage(DistributionWithTrendStorage * vp,
                                                             DistributionWithTrendStorage * density,
                                                             DistributionWithTrendStorage * correlation_vp_density)
: vp_(vp),
  density_(density),
  correlation_vp_density_(correlation_vp_density)
{
}

TabulatedVelocityFluidStorage::~TabulatedVelocityFluidStorage()
{
  delete vp_;
  delete density_;
  delete correlation_vp_density_;
}

DistributionsFluid *
TabulatedVelocityFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_vp_density_dist_with_trend = correlation_vp_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = new DistributionsFluidTabulatedVelocity(vp_dist_with_trend,
                                                                       density_dist_with_trend,
                                                                       corr_vp_density_dist_with_trend);

  if(fluid == NULL)
    errTxt += "The tabulated model has not been implemented yet for fluids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusFluidStorage::TabulatedModulusFluidStorage(DistributionWithTrendStorage * bulk_modulus,
                                                           DistributionWithTrendStorage * density,
                                                           DistributionWithTrendStorage * correlation_bulk_density)
: bulk_modulus_(bulk_modulus),
  density_(density),
  correlation_bulk_density_(correlation_bulk_density)
{
}

TabulatedModulusFluidStorage::~TabulatedModulusFluidStorage()
{
  delete bulk_modulus_;
  delete density_;
  delete correlation_bulk_density_;
}

DistributionsFluid *
TabulatedModulusFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                         std::string                             & errTxt) const
{
  const DistributionWithTrend * bulk_dist_with_trend           = bulk_modulus_            ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend        = density_                 ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * corr_bulk_dens_dist_with_trend = correlation_bulk_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = NULL;
  //Make new DistributionsFluidTabulatedModulus(bulk_dist_with_trend, density_dist_with_trend, bulk_dist_with_trend)

  //Delete these pointers in new class:
  delete bulk_dist_with_trend;
  delete density_dist_with_trend;
  delete corr_bulk_dens_dist_with_trend;

  if(fluid == NULL)
    errTxt += "The tabulated model has not been implemented yet for fluids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}

//----------------------------------------------------------------------------------//
ReussFluidStorage::ReussFluidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussFluidStorage::~ReussFluidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];
}

DistributionsFluid *
ReussFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  //Make new DistributionsFluidReuss(constituent_label, volume);
  DistributionsFluid * fluid = NULL; //new DistributionsFluidReuss(constituent_label, volume);

  if(fluid == NULL)
    errTxt += "The Reuss model has not been implemented yet for fluids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}

//----------------------------------------------------------------------------------//
VoigtFluidStorage::VoigtFluidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtFluidStorage::~VoigtFluidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];
}

DistributionsFluid *
VoigtFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsFluid * fluid = NULL; //new DistributionsFluidVoigt();

  if(fluid == NULL)
    errTxt += "The Voigt model has not been implemented yet for fluids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}

//----------------------------------------------------------------------------------//
HillFluidStorage::HillFluidStorage(std::vector<std::string>                    constituent_label,
                                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillFluidStorage::~HillFluidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++)
    delete constituent_volume_fraction_[i];
}

DistributionsFluid *
HillFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                             std::string                             & errTxt) const
{
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsFluid * fluid = NULL; //new DistributionsFluidHill();

  if(fluid == NULL)
    errTxt += "The Hill model has not been implemented yet for fluids\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}

//----------------------------------------------------------------------------------//
BatzleWangFluidStorage::BatzleWangFluidStorage(DistributionWithTrendStorage * pore_pressure,
                                               DistributionWithTrendStorage * temperature,
                                               DistributionWithTrendStorage * salinity)
: pore_pressure_(pore_pressure),
  temperature_(temperature),
  salinity_(salinity)
{
}

BatzleWangFluidStorage::~BatzleWangFluidStorage()
{
  if(pore_pressure_->GetIsShared() == false)
    delete pore_pressure_;

  if(temperature_->GetIsShared() == false)
    delete temperature_;

  if(salinity_->GetIsShared() == false)
    delete salinity_;
}

DistributionsFluid *
BatzleWangFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                                   const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                   const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                   std::string                             & errTxt) const
{
  DistributionsFluid * fluid = NULL; //new DistributionsFluidBatzleWang();

  if(fluid == NULL)
    errTxt += "The Batzle-Wang model has not been implemented yet\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(fluid);
}
