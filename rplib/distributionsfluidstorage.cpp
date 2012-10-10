#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsfluidtabulatedvelocity.h"
#include "rplib/distributionsfluidtabulatedmodulus.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionsbrine.h"
#include "rplib/distributionsco2.h"

class DistributionsBrineEvolution;
class DistributionsCO2Evolution;


DistributionsFluidStorage::DistributionsFluidStorage()
{
}

DistributionsFluidStorage::~DistributionsFluidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedVelocityFluidStorage::TabulatedVelocityFluidStorage(DistributionWithTrendStorage * vp,
                                                             DistributionWithTrendStorage * density,
                                                             double                         correlation_vp_density)
: vp_(vp),
  density_(density),
  correlation_vp_density_(correlation_vp_density)
{
}

TabulatedVelocityFluidStorage::~TabulatedVelocityFluidStorage()
{
  if(vp_->GetIsShared() == false)
    delete vp_;
  if(density_->GetIsShared() == false)
    delete density_;
}

DistributionsFluid *
TabulatedVelocityFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = new DistributionsFluidTabulatedVelocity(vp_dist_with_trend,
                                                                       density_dist_with_trend,
                                                                       correlation_vp_density_);

  return(fluid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusFluidStorage::TabulatedModulusFluidStorage(DistributionWithTrendStorage * bulk_modulus,
                                                           DistributionWithTrendStorage * density,
                                                           double                         correlation_bulk_density)
: bulk_modulus_(bulk_modulus),
  density_(density),
  correlation_bulk_density_(correlation_bulk_density)
{
}

TabulatedModulusFluidStorage::~TabulatedModulusFluidStorage()
{
  if(bulk_modulus_->GetIsShared() == false)
    delete bulk_modulus_;
  if(density_->GetIsShared() == false)
    delete density_;
}

DistributionsFluid *
TabulatedModulusFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                         std::string                             & errTxt) const
{
  const DistributionWithTrend * bulk_dist_with_trend    = bulk_modulus_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend = density_     ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsFluid * fluid = new DistributionsFluidTabulatedModulus(bulk_dist_with_trend,
                                                                      density_dist_with_trend,
                                                                      correlation_bulk_density_);

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
  //CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

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
  //CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

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
  //CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

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
BatzleWangFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                   const std::vector<std::string>          & trend_cube_parameters,
                                                   const std::vector<std::vector<double> > & trend_cube_sampling,
                                                   std::string                             & errTxt) const
{
  DistributionsFluid          *                fluid = NULL;

  const DistributionWithTrend * pore_pressure_distr  = pore_pressure_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * temperature_distr    = temperature_  ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * salinity_distr       = salinity_     ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);


  //NBNB fjellvoll //NBNB marit CO2 is not finished yet in XML interface
  bool is_co2 = false;

  if (is_co2)
    fluid       = new DistributionsCO2(temperature_distr, pore_pressure_distr);
  else
    fluid       = new DistributionsBrine(temperature_distr, pore_pressure_distr, salinity_distr);

  return(fluid);
}
