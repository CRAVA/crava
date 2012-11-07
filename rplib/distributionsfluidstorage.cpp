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
#include "rplib/distributionsfluidbatzlewang.h"
#include "rplib/distributionsco2.h"


DistributionsFluidStorage::DistributionsFluidStorage()
{
}

DistributionsFluidStorage::~DistributionsFluidStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedVelocityFluidStorage::TabulatedVelocityFluidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                             std::vector<DistributionWithTrendStorage *> density,
                                                             double                                      correlation_vp_density)
: vp_(vp),
  density_(density),
  correlation_vp_density_(correlation_vp_density)
{
}

TabulatedVelocityFluidStorage::~TabulatedVelocityFluidStorage()
{
  if(vp_[0]->GetIsShared() == false)
    delete vp_[0];
  if(density_[0]->GetIsShared() == false)
    delete density_[0];
}

DistributionsFluid *
TabulatedVelocityFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const
{
  std::vector<double> alpha(2);
  alpha[0] = vp_[0]     ->GetOneYearCorrelation();
  alpha[1] = density_[0]->GetOneYearCorrelation();

  int n_vintages_vp      = static_cast<int>(vp_.size());
  int n_vintages_density = static_cast<int>(density_.size());

  int max_vintage = std::max(n_vintages_vp, n_vintages_density);

  std::vector<DistributionsFluid *>    dist_fluid(max_vintage, NULL);
  std::vector<DistributionWithTrend *> vp_dist_with_trend(max_vintage, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(max_vintage, NULL);

  for(int i=0; i<max_vintage; i++) {
    if(i < n_vintages_vp)
      vp_dist_with_trend[i] = vp_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      vp_dist_with_trend[i] = vp_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    DistributionsFluid * fluid = new DistributionsFluidTabulatedVelocity(vp_dist_with_trend[i],
                                                                         density_dist_with_trend[i],
                                                                         correlation_vp_density_,
                                                                         alpha);

    dist_fluid[i] = fluid;


  }

  return(dist_fluid[0]);
}

//----------------------------------------------------------------------------------//
TabulatedModulusFluidStorage::TabulatedModulusFluidStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                           std::vector<DistributionWithTrendStorage *> density,
                                                           double                         correlation_bulk_density)
: bulk_modulus_(bulk_modulus),
  density_(density),
  correlation_bulk_density_(correlation_bulk_density)
{
}

TabulatedModulusFluidStorage::~TabulatedModulusFluidStorage()
{
  if(bulk_modulus_[0]->GetIsShared() == false)
    delete bulk_modulus_[0];
  if(density_[0]->GetIsShared() == false)
    delete density_[0];
}

DistributionsFluid *
TabulatedModulusFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                         const std::vector<std::string>          & trend_cube_parameters,
                                                         const std::vector<std::vector<double> > & trend_cube_sampling,
                                                         std::string                             & errTxt) const
{
  std::vector<double> alpha(2);
  alpha[0] = bulk_modulus_[0]->GetOneYearCorrelation();
  alpha[1] = density_[0]     ->GetOneYearCorrelation();

  int n_vintages_bulk    = static_cast<int>(bulk_modulus_.size());
  int n_vintages_density = static_cast<int>(density_.size());

  int max_vintage = std::max(n_vintages_bulk, n_vintages_density);

  std::vector<DistributionsFluid *>    dist_fluid(max_vintage, NULL);
  std::vector<DistributionWithTrend *> bulk_dist_with_trend(max_vintage, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(max_vintage, NULL);

  for(int i=0; i<max_vintage; i++) {
    if(i < n_vintages_bulk)
      bulk_dist_with_trend[i] = bulk_modulus_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      bulk_dist_with_trend[i] = bulk_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    DistributionsFluid * fluid = new DistributionsFluidTabulatedModulus(bulk_dist_with_trend[i],
                                                                        density_dist_with_trend[i],
                                                                        correlation_bulk_density_,
                                                                        alpha);

    dist_fluid[i] = fluid;


  }

  return(dist_fluid[0]);
}

//----------------------------------------------------------------------------------//
ReussFluidStorage::ReussFluidStorage(std::vector<std::string>                                  constituent_label,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussFluidStorage::~ReussFluidStorage()
{
  for(size_t i=0; i<constituent_volume_fraction_[0].size(); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == true)
      delete constituent_volume_fraction_[0][i];
  }
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
VoigtFluidStorage::VoigtFluidStorage(std::vector<std::string>                                  constituent_label,
                                     std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtFluidStorage::~VoigtFluidStorage()
{
  for(size_t i=0; i<constituent_volume_fraction_[0].size(); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == true)
      delete constituent_volume_fraction_[0][i];
  }
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
HillFluidStorage::HillFluidStorage(std::vector<std::string>                                  constituent_label,
                                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillFluidStorage::~HillFluidStorage()
{
  for(size_t i=0; i<constituent_volume_fraction_[0].size(); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == true)
      delete constituent_volume_fraction_[0][i];
  }
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
BatzleWangFluidStorage::BatzleWangFluidStorage(std::vector<DistributionWithTrendStorage *> pore_pressure,
                                               std::vector<DistributionWithTrendStorage *> temperature,
                                               std::vector<DistributionWithTrendStorage *> salinity)
: pore_pressure_(pore_pressure),
  temperature_(temperature),
  salinity_(salinity)
{
}

BatzleWangFluidStorage::~BatzleWangFluidStorage()
{
  if(pore_pressure_[0]->GetIsShared() == false)
    delete pore_pressure_[0];

  if(temperature_[0]->GetIsShared() == false)
    delete temperature_[0];

  if(salinity_[0]->GetIsShared() == false)
    delete salinity_[0];
}

DistributionsFluid *
BatzleWangFluidStorage::GenerateDistributionsFluid(const std::string                       & path,
                                                   const std::vector<std::string>          & trend_cube_parameters,
                                                   const std::vector<std::vector<double> > & trend_cube_sampling,
                                                   std::string                             & errTxt) const
{
  std::vector<double> alpha(3);
  alpha[0] = pore_pressure_[0]->GetOneYearCorrelation();
  alpha[1] = temperature_[0]  ->GetOneYearCorrelation();
  alpha[2] = salinity_[0]     ->GetOneYearCorrelation();

  int n_vintages_pressure    = static_cast<int>(pore_pressure_.size());
  int n_vintages_temperature = static_cast<int>(temperature_.size());
  int n_vintages_salinity    = static_cast<int>(salinity_.size());

  int max_vintage;
  max_vintage = std::max(n_vintages_pressure, n_vintages_temperature);
  max_vintage = std::max(max_vintage, n_vintages_salinity);

  std::vector<DistributionsFluid *>    dist_fluid(max_vintage, NULL);
  std::vector<DistributionWithTrend *> pressure_dist_with_trend(max_vintage, NULL);
  std::vector<DistributionWithTrend *> temperature_dist_with_trend(max_vintage, NULL);
  std::vector<DistributionWithTrend *> salinity_dist_with_trend(max_vintage, NULL);

  for(int i=0; i<max_vintage; i++) {
    if(i < n_vintages_pressure)
      pressure_dist_with_trend[i] = pore_pressure_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      pressure_dist_with_trend[i] = pressure_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_temperature)
      temperature_dist_with_trend[i] = temperature_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      temperature_dist_with_trend[i] = temperature_dist_with_trend[i-1]->Clone();

    //NBNB fjellvoll //NBNB marit CO2 is not finished yet in XML interface
    bool is_co2 = false;
    DistributionsFluid * fluid = NULL;

    if (is_co2) {
      std::vector<double> dummy_alpha(2,1);
      fluid       = new DistributionsCO2(temperature_dist_with_trend[i], pressure_dist_with_trend[i], dummy_alpha);
    }
    else
      fluid       = new DistributionsFluidBatzleWang(temperature_dist_with_trend[i], pressure_dist_with_trend[i], salinity_dist_with_trend[i], alpha);


    dist_fluid[i] = fluid;


  }


  DistributionsFluid * dummy = NULL;
  return(dummy);

}
