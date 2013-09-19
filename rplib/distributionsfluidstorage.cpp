#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

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

  const std::vector<std::vector<float> >  dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;

  for(int i=0; i<n_constituents; i++)
    CheckValuesInZeroOne(constituent_volume_fraction[i], "volume-fraction", path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, errTxt);

  std::vector<DistributionsFluid *>                  final_dist_fluid(n_vintages, NULL);
  std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);

  for(int i=0; i<n_vintages; i++)
    all_volume_fractions[i].resize(n_constituents, NULL);

  for(int i=0; i<n_vintages; i++) {

    for (int s=0; s<n_constituents; s++) {

      if(i < n_vintages_constit[s]) {
        if(constituent_volume_fraction[s][i] != NULL)
          all_volume_fractions[i][s] = constituent_volume_fraction[s][i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, errTxt);
      }
      else {
        if(all_volume_fractions[i-1][s] != NULL)
          all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
      }
    }

    CheckVolumeConsistency(all_volume_fractions[i], errTxt);
  }

  if (errTxt == "") {
    for(int i=0; i<n_vintages; i++)
      final_dist_fluid[i] = new DistributionsFluidMix(alpha, distr_fluid[i], all_volume_fractions[i], mix_method);

    for(int i=0; i<n_vintages; i++) {
      for(size_t s=0; s<distr_fluid[i].size(); s++)
        delete distr_fluid[i][s];

      for(size_t s=0; s<all_volume_fractions[i].size(); s++) {
        if(all_volume_fractions[i][s] != NULL) {
          if(all_volume_fractions[i][s]->GetIsShared() == false)
            delete all_volume_fractions[i][s];
        }
      }
    }
  }

  return(final_dist_fluid);
}

//----------------------------------------------------------------------------------//
TabulatedVelocityFluidStorage::TabulatedVelocityFluidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                             std::vector<DistributionWithTrendStorage *> density,
                                                             std::vector<DistributionWithTrendStorage *> correlation_vp_density)
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

  if(correlation_vp_density_[0]->GetIsShared() == false)
    delete correlation_vp_density_[0];
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

  std::string tmpErrTxt = "";

  for(int i=0; i<n_vintages_vp; i++) {
    if(vp_[i]->GetEstimate() == true)
      tmpErrTxt += "Vp can not be estimated from wells\n";
  }
  for(int i=0; i<n_vintages_density; i++) {
    if(density_[i]->GetEstimate() == true)
      tmpErrTxt += "Density can not be estimated from wells\n";
  }

  std::vector<double> corr_vp_density(n_vintages, 0);
  for(int i=0; i<n_vintages_vp_density; i++) {
    if(correlation_vp_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-vp-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_vp_density_[i], "correlation", corr_vp_density[i], errTxt);

    if(corr_vp_density[i] > 1 || corr_vp_density[i] < -1)
        errTxt += "<correlation-vp-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  const std::vector<std::vector<float> >  dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vp_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);


  if(tmpErrTxt == "") {

    for(int i=0; i<n_vintages; i++) {

      if(i < n_vintages_vp)
        vp_dist_with_trend[i] = vp_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        vp_dist_with_trend[i] = vp_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_density)
        density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

     if(i >= n_vintages_vp_density)
        corr_vp_density[i] = corr_vp_density[i-1];
    }
  }

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      DistributionsFluid * fluid = new DistributionsFluidTabulated(vp_dist_with_trend[i],
                                                                   density_dist_with_trend[i],
                                                                   corr_vp_density[i],
                                                                   DEMTools::Velocity,
                                                                   alpha);

      dist_fluid[i] = fluid;
    }

    for(int i=0; i<n_vintages; i++) {
      if(vp_dist_with_trend[i]->GetIsShared() == false)
        delete vp_dist_with_trend[i];
      if(density_dist_with_trend[i]->GetIsShared() == false)
        delete density_dist_with_trend[i];
    }
  }

  else{
    errTxt += "Problems with the Tabulated rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_fluid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusFluidStorage::TabulatedModulusFluidStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                           std::vector<DistributionWithTrendStorage *> density,
                                                           std::vector<DistributionWithTrendStorage *> correlation_bulk_density)
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

  if(correlation_bulk_density_[0]->GetIsShared() == false)
    delete correlation_bulk_density_[0];

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

  std::string tmpErrTxt = "";

  for(int i=0; i<n_vintages_bulk; i++) {
    if(bulk_modulus_[i]->GetEstimate() == true)
      tmpErrTxt += "Bulk modulus can not be estimated from wells\n";
  }
  for(int i=0; i<n_vintages_density; i++) {
    if(density_[i]->GetEstimate() == true)
      tmpErrTxt += "Density can not be estimated from wells\n";
  }

  std::vector<double> corr_bulk_density(n_vintages, 0);
  for(int i=0; i<n_vintages_bulk_density; i++) {
    if(correlation_bulk_density_[i]->GetEstimate() == true)
      tmpErrTxt += "<correlation-bulk-density> can not be estimated from wells\n";
    else
      FindDoubleValueFromDistributionWithTrend(correlation_bulk_density_[i], "correlation", corr_bulk_density[i], errTxt);

    if(corr_bulk_density[i] > 1 || corr_bulk_density[i] < -1)
        errTxt += "<correlation-bulk-density> should be in the interval [-1,1] in the tabulated model\n";
  }

  const std::vector<std::vector<float> >  dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> bulk_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_bulk)
        bulk_dist_with_trend[i] = bulk_modulus_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        bulk_dist_with_trend[i] = bulk_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_density)
        density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

      double lower_mega_fluid = 1.0e+4;
      double upper_mega_fluid = 1.0e+7;
      double test_bulk  = bulk_dist_with_trend[0]->ReSample(0,0);
      if(test_bulk < lower_mega_fluid || test_bulk > upper_mega_fluid)
        tmpErrTxt += "Bulk modulus need to be given in kPa\n";

      if(i >= n_vintages_bulk_density)
        corr_bulk_density.push_back(corr_bulk_density[i-1]);
    }
  }

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      DistributionsFluid * fluid = new DistributionsFluidTabulated(bulk_dist_with_trend[i],
                                                                   density_dist_with_trend[i],
                                                                   corr_bulk_density[i],
                                                                   DEMTools::Modulus,
                                                                   alpha);

      dist_fluid[i] = fluid;
    }

    for(int i=0; i<n_vintages; i++) {
      if(bulk_dist_with_trend[i]->GetIsShared() == false)
        delete bulk_dist_with_trend[i];
      if(density_dist_with_trend[i]->GetIsShared() == false)
        delete density_dist_with_trend[i];
    }
  }

  else{
    errTxt += "Problems with the Tabulated rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
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
    if(constituent_volume_fraction_[0][i] != NULL) {
      if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
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
  std::vector<DistributionsFluid *> fluid;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    fluid = CreateDistributionsFluidMix(n_vintages,
                                        path,
                                        trend_cube_parameters,
                                        trend_cube_sampling,
                                        model_fluid_storage,
                                        constituent_label_,
                                        constituent_volume_fraction_,
                                        DEMTools::Reuss,
                                        tmpErrTxt);

  else {
    errTxt += "\nProblems with the Reuss rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

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
    if(constituent_volume_fraction_[0][i] != NULL) {
      if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
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
  std::vector<DistributionsFluid *> fluid;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    fluid = CreateDistributionsFluidMix(n_vintages,
                                        path,
                                        trend_cube_parameters,
                                        trend_cube_sampling,
                                        model_fluid_storage,
                                        constituent_label_,
                                        constituent_volume_fraction_,
                                        DEMTools::Voigt,
                                        tmpErrTxt);

  else {
    errTxt += "\nProblems with the Voigt rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

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
    if(constituent_volume_fraction_[0][i] != NULL) {
      if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
        delete constituent_volume_fraction_[0][i];
    }
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
  std::vector<DistributionsFluid *> fluid;

  std::string tmpErrTxt = "";

  for(size_t i=0; i<constituent_volume_fraction_.size(); i++) {
    for(size_t j=0; j<constituent_volume_fraction_[i].size(); j++) {
      if(constituent_volume_fraction_[i][j] != NULL && constituent_volume_fraction_[i][j]->GetEstimate() == true) {
        tmpErrTxt += "The volume fractions can not be estimated from wells\n";
        break;
      }
    }
    if(tmpErrTxt != "")
      break;
  }

  if(tmpErrTxt == "")
    fluid = CreateDistributionsFluidMix(n_vintages,
                                        path,
                                        trend_cube_parameters,
                                        trend_cube_sampling,
                                        model_fluid_storage,
                                        constituent_label_,
                                        constituent_volume_fraction_,
                                        DEMTools::Hill,
                                        tmpErrTxt);

  else {
    errTxt += "\nProblems with the Hill rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

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
  std::string tmpErrTxt = "";

  std::vector<double> alpha(3);
  alpha[0] = salinity_[0]     ->GetOneYearCorrelation();
  alpha[1] = temperature_[0]  ->GetOneYearCorrelation();
  alpha[2] = pore_pressure_[0]->GetOneYearCorrelation();

  for(size_t i=0; i<pore_pressure_.size(); i++) {
    if(pore_pressure_[i]->GetEstimate() == true)
      tmpErrTxt += "Pore pressure can not be estimated from wells for a <solid>\n";
  }
  for(size_t i=0; i<temperature_.size(); i++) {
    if(temperature_[i]->GetEstimate() == true)
      tmpErrTxt += "Temperature can not be estimated from wells for a <solid>\n";
  }
  for(size_t i=0; i<salinity_.size(); i++) {
    if(salinity_[i]->GetEstimate() == true)
      tmpErrTxt += "Salinity can not be estimatedfrom wells for a <solid>\n";
  }

  int n_vintages_pressure    = static_cast<int>(pore_pressure_.size());
  int n_vintages_temperature = static_cast<int>(temperature_.size());
  int n_vintages_salinity    = static_cast<int>(salinity_.size());

  const std::vector<std::vector<float> >  dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> pressure_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> temperature_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> salinity_dist_with_trend(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_pressure)
        pressure_dist_with_trend[i] = pore_pressure_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        pressure_dist_with_trend[i] = pressure_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_temperature)
        temperature_dist_with_trend[i] = temperature_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        temperature_dist_with_trend[i] = temperature_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_salinity)
        salinity_dist_with_trend[i] = salinity_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        salinity_dist_with_trend[i] = salinity_dist_with_trend[i-1]->Clone();

    }
  }
  DistributionsFluid * fluid = NULL;

  if(tmpErrTxt == "") {

    for(int i=0; i<n_vintages; i++) {
      fluid       = new DistributionsFluidBatzleWang(temperature_dist_with_trend[i], pressure_dist_with_trend[i], salinity_dist_with_trend[i], alpha);
      dist_fluid[i] = fluid;
    }

    for(int i=0; i<n_vintages; i++) {
      if(temperature_dist_with_trend[i]->GetIsShared() == false)
        delete temperature_dist_with_trend[i];
      if(pressure_dist_with_trend[i]->GetIsShared() == false)
        delete pressure_dist_with_trend[i];
      if(salinity_dist_with_trend[i]->GetIsShared() == false)
        delete salinity_dist_with_trend[i];
    }
  }

  else {
    errTxt += "\nProblems with the Batzle-Wang rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_fluid);
}


//----------------------------------------------------------------------------------//
CO2FluidStorage::CO2FluidStorage(std::vector<DistributionWithTrendStorage *> pressure,
                                 std::vector<DistributionWithTrendStorage *> temperature)
: pressure_(pressure),
  temperature_(temperature)
{
}

CO2FluidStorage::~CO2FluidStorage()
{
  if(pressure_[0]->GetIsShared() == false)
    delete pressure_[0];

  if(temperature_[0]->GetIsShared() == false)
    delete temperature_[0];
}

std::vector<DistributionsFluid *>
CO2FluidStorage::GenerateDistributionsFluid(const int                                                      & n_vintages,
                                            const std::string                                              & path,
                                            const std::vector<std::string>                                 & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                            const std::map<std::string, DistributionsFluidStorage *>       & /*model_fluid_storage*/,
                                            std::string                                                    & errTxt) const
{
  std::string tmpErrTxt = "";

  std::vector<double> alpha(2);
  alpha[0] = temperature_[0]  ->GetOneYearCorrelation();
  alpha[1] = pressure_[0]->GetOneYearCorrelation();

  for(size_t i=0; i<pressure_.size(); i++) {
    if(pressure_[i]->GetEstimate() == true)
      tmpErrTxt += "Pressure can not be estimated from wells for a <solid>\n";
  }
  for(size_t i=0; i<temperature_.size(); i++) {
    if(temperature_[i]->GetEstimate() == true)
      tmpErrTxt += "Temperature can not be estimated from wells for a <solid>\n";
  }

  int n_vintages_pressure    = static_cast<int>(pressure_.size());
  int n_vintages_temperature = static_cast<int>(temperature_.size());

  const std::vector<std::vector<float> >  dummy_blocked_logs;
  const std::vector<std::vector<double> > dummy_s1;
  const std::vector<std::vector<double> > dummy_s2;

  std::vector<DistributionsFluid *>    dist_fluid(n_vintages, NULL);
  std::vector<DistributionWithTrend *> pressure_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> temperature_dist_with_trend(n_vintages, NULL);

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_pressure)
        pressure_dist_with_trend[i] = pressure_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        pressure_dist_with_trend[i] = pressure_dist_with_trend[i-1]->Clone();

      if(i < n_vintages_temperature)
        temperature_dist_with_trend[i] = temperature_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, dummy_blocked_logs, dummy_s1, dummy_s2, tmpErrTxt);
      else
        temperature_dist_with_trend[i] = temperature_dist_with_trend[i-1]->Clone();

    }
  }

  DistributionsFluid * fluid = NULL;

  if(tmpErrTxt == "") {
    for(int i=0; i<n_vintages; i++) {
      fluid         = new DistributionsFluidCO2(temperature_dist_with_trend[i], pressure_dist_with_trend[i], alpha);
      dist_fluid[i] = fluid;
    }

    for(int i=0; i<n_vintages; i++) {
      if(temperature_dist_with_trend[i]->GetIsShared() == false)
        delete temperature_dist_with_trend[i];
      if(pressure_dist_with_trend[i]->GetIsShared() == false)
        delete pressure_dist_with_trend[i];
    }
  }

  else {
    errTxt += "\nProblems with the CO2 rock physics model for <fluid>:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_fluid);
}
