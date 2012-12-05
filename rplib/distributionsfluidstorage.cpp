#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/iotools/stringtools.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsfluidmix.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsfluidtabulated.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionsfluidbatzlewang.h"
#include "rplib/distributionsfluidco2.h"
#include "rplib/demmodelling.h"


DistributionsFluidStorage::DistributionsFluidStorage()
{
}

DistributionsFluidStorage::~DistributionsFluidStorage()
{
}

std::vector<DistributionsFluid *>
DistributionsFluidStorage::CreateDistributionsFluidMix(const int                                                       & n_vintages,
                                                       const std::string                                               & path,
                                                       const std::vector<std::string>                                  & trend_cube_parameters,
                                                       const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                       const std::map<std::string, DistributionsFluidStorage *>        & model_fluid_storage,
                                                       const std::vector<std::string>                                  & constituent_label,
                                                       const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                       DEMTools::MixMethod                                               mix_method,
                                                       std::string                                                     & errTxt) const
{
  int n_constituents = static_cast<int>(constituent_label.size());

  std::vector<int> n_vintages_constit(n_constituents);
  for(int i=0; i<n_constituents; i++)
    n_vintages_constit[i] = static_cast<int>(constituent_volume_fraction[i].size());

  std::vector<double> alpha(n_constituents);
  for(int i=0; i<n_constituents; i++) {
    if(constituent_volume_fraction[i][0] != NULL)
      alpha[i] = constituent_volume_fraction[i][0]->GetOneYearCorrelation();
    else
      alpha[i] = 1;
  }

  std::vector<std::vector<DistributionsFluid *> > distr_fluid(n_vintages);
  for(int i=0; i<n_vintages; i++)
    distr_fluid[i].resize(n_constituents, NULL);

  for (int s = 0; s < n_constituents; s++) {
    std::vector<DistributionsFluid *> distr_fluid_all_vintages = ReadFluid(n_vintages,
                                                                           constituent_label[s],
                                                                           path,
                                                                           trend_cube_parameters,
                                                                           trend_cube_sampling,
                                                                           model_fluid_storage,
                                                                           errTxt);

    for(int i=0; i<n_vintages; i++) {
      if(i < static_cast<int>(distr_fluid_all_vintages.size()))
        distr_fluid[i][s] = distr_fluid_all_vintages[i];
      else
        distr_fluid[i][s] = distr_fluid[i-1][s]->Clone();
    }
  }

  std::vector<DistributionsFluid *>                  final_dist_fluid(n_vintages, NULL);
  std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);

  for(int i=0; i<n_vintages; i++)
    all_volume_fractions[i].resize(n_constituents, NULL);

  for(int i=0; i<n_vintages; i++) {

    for (int s=0; s<n_constituents; s++) {

      if(i < n_vintages_constit[s]) {
        if(constituent_volume_fraction[s][i] != NULL)
          all_volume_fractions[i][s] = constituent_volume_fraction[s][i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
      }
      else {
        if(all_volume_fractions[i-1][s] != NULL)
          all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
      }
    }

    CheckVolumeConsistency(all_volume_fractions[i], errTxt);

    if (errTxt == "")
      final_dist_fluid[i] = new DistributionsFluidMix(alpha, distr_fluid[i], all_volume_fractions[i], mix_method);
  }

  return(final_dist_fluid);
}

//----------------------------------------------------------------------------------//
TabulatedVelocityFluidStorage::TabulatedVelocityFluidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                             std::vector<DistributionWithTrendStorage *> density,
                                                             std::vector<double>                         correlation_vp_density)
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

std::vector<DistributionsFluid *>
TabulatedVelocityFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                          const std::string                                              & path,
                                                          const std::vector<std::string>                                 & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsFluidStorage *>       & /*model_fluid_storage*/,
                                                          std::string                                                    & errTxt) const
{
  std::vector<double> alpha(2);
  alpha[0] = vp_[0]     ->GetOneYearCorrelation();
  alpha[1] = density_[0]->GetOneYearCorrelation();

  int n_vintages_vp         = static_cast<int>(vp_.size());
  int n_vintages_density    = static_cast<int>(density_.size());
  int n_vintages_vp_density = static_cast<int>(correlation_vp_density_.size());

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vp_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);

  std::vector<double> corr_vp_density;

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_vp)
      vp_dist_with_trend[i] = vp_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      vp_dist_with_trend[i] = vp_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_vp_density)
      corr_vp_density.push_back(correlation_vp_density_[i]);
    else
      corr_vp_density.push_back(correlation_vp_density_[i-1]);
  }

   for(int i=0; i<n_vintages; i++) {
    DistributionsFluid * fluid = new DistributionsFluidTabulated(vp_dist_with_trend[i],
                                                                 density_dist_with_trend[i],
                                                                 corr_vp_density[i],
                                                                 DEMTools::Velocity,
                                                                 alpha);

    dist_fluid[i] = fluid;
  }

  return(dist_fluid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusFluidStorage::TabulatedModulusFluidStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                           std::vector<DistributionWithTrendStorage *> density,
                                                           std::vector<double>                         correlation_bulk_density)
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

std::vector<DistributionsFluid *>
TabulatedModulusFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                         const std::string                                              & path,
                                                         const std::vector<std::string>                                 & trend_cube_parameters,
                                                         const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                         const std::map<std::string, DistributionsFluidStorage *>       & /*model_fluid_storage*/,
                                                         std::string                                                    & errTxt) const
{
  std::vector<double> alpha(2);
  alpha[0] = bulk_modulus_[0]->GetOneYearCorrelation();
  alpha[1] = density_[0]     ->GetOneYearCorrelation();

  int n_vintages_bulk         = static_cast<int>(bulk_modulus_.size());
  int n_vintages_density      = static_cast<int>(density_.size());
  int n_vintages_bulk_density = static_cast<int>(correlation_bulk_density_.size());

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> bulk_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);

  std::vector<double> corr_bulk_density;

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_bulk)
      bulk_dist_with_trend[i] = bulk_modulus_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      bulk_dist_with_trend[i] = bulk_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    double lower_mega = 1.0e+5;
    double upper_mega = 1.0e+8;
    double test_bulk  = bulk_dist_with_trend[0]->ReSample(0,0);
    if(test_bulk < lower_mega || test_bulk > upper_mega)
      errTxt += "Bulk modulus need to be given in megaPa\n";

    if(i < n_vintages_bulk_density)
      corr_bulk_density.push_back(correlation_bulk_density_[i]);
    else
      corr_bulk_density.push_back(correlation_bulk_density_[i-1]);
  }

  for(int i=0; i<n_vintages; i++) {
    DistributionsFluid * fluid = new DistributionsFluidTabulated(bulk_dist_with_trend[i],
                                                                 density_dist_with_trend[i],
                                                                 corr_bulk_density[i],
                                                                 DEMTools::Modulus,
                                                                 alpha);

    dist_fluid[i] = fluid;
  }

  return(dist_fluid);
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

std::vector<DistributionsFluid *>
ReussFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                              const std::string                                              & path,
                                              const std::vector<std::string>                                 & trend_cube_parameters,
                                              const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                              const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                              std::string                                                    & errTxt) const
{
  std::vector<DistributionsFluid *> fluid = CreateDistributionsFluidMix(n_vintages,
                                                                        path,
                                                                        trend_cube_parameters,
                                                                        trend_cube_sampling,
                                                                        model_fluid_storage,
                                                                        constituent_label_,
                                                                        constituent_volume_fraction_,
                                                                        DEMTools::Reuss,
                                                                        errTxt);
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

std::vector<DistributionsFluid *>
VoigtFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                              const std::string                                              & path,
                                              const std::vector<std::string>                                 & trend_cube_parameters,
                                              const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                              const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                              std::string                                                    & errTxt) const
{
  std::vector<DistributionsFluid *> fluid = CreateDistributionsFluidMix(n_vintages,
                                                                        path,
                                                                        trend_cube_parameters,
                                                                        trend_cube_sampling,
                                                                        model_fluid_storage,
                                                                        constituent_label_,
                                                                        constituent_volume_fraction_,
                                                                        DEMTools::Voigt,
                                                                        errTxt);
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

std::vector<DistributionsFluid *>
HillFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                             const std::string                                              & path,
                                             const std::vector<std::string>                                 & trend_cube_parameters,
                                             const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                             const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                             std::string                                                    & errTxt) const
{
  std::vector<DistributionsFluid *> fluid = CreateDistributionsFluidMix(n_vintages,
                                                                        path,
                                                                        trend_cube_parameters,
                                                                        trend_cube_sampling,
                                                                        model_fluid_storage,
                                                                        constituent_label_,
                                                                        constituent_volume_fraction_,
                                                                        DEMTools::Hill,
                                                                        errTxt);
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

std::vector<DistributionsFluid *>
BatzleWangFluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                   const std::string                                              & path,
                                                   const std::vector<std::string>                                 & trend_cube_parameters,
                                                   const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                   const std::map<std::string, DistributionsFluidStorage *>       & /*model_fluid_storage*/,
                                                   std::string                                                    & errTxt) const
{
  std::vector<double> alpha(3);
  alpha[0] = salinity_[0]     ->GetOneYearCorrelation();
  alpha[1] = temperature_[0]  ->GetOneYearCorrelation();
  alpha[2] = pore_pressure_[0]->GetOneYearCorrelation();

  int n_vintages_pressure    = static_cast<int>(pore_pressure_.size());
  int n_vintages_temperature = static_cast<int>(temperature_.size());
  int n_vintages_salinity    = static_cast<int>(salinity_.size());

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> pressure_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> temperature_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> salinity_dist_with_trend(n_vintages, NULL);

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_pressure)
      pressure_dist_with_trend[i] = pore_pressure_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      pressure_dist_with_trend[i] = pressure_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_temperature)
      temperature_dist_with_trend[i] = temperature_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      temperature_dist_with_trend[i] = temperature_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_salinity)
      salinity_dist_with_trend[i] = salinity_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      salinity_dist_with_trend[i] = salinity_dist_with_trend[i-1]->Clone();

    //NBNB fjellvoll //NBNB marit CO2 is not finished yet in XML interface
    bool is_co2 = false;
    DistributionsFluid * fluid = NULL;

    if (is_co2) {
      std::vector<double> alpha2(2);
      std::copy(alpha.begin() + 1, alpha.end(), alpha2.begin());
      fluid       = new DistributionsFluidCO2(temperature_dist_with_trend[i], pressure_dist_with_trend[i], alpha2);
    }
    else
      fluid       = new DistributionsFluidBatzleWang(temperature_dist_with_trend[i], pressure_dist_with_trend[i], salinity_dist_with_trend[i], alpha);


    dist_fluid[i] = fluid;


  }

  return(dist_fluid);
}
