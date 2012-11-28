#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"

#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsfluidstorage.h"

#include "rplib/distributionsrocktabulated.h"
#include "rplib/distributionsrockbounding.h"
#include "rplib/distributionsrockmix.h"
#include "rplib/distributionsrockdem.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"

#include "src/modelsettings.h"

#include <typeinfo>

DistributionsRockStorage::DistributionsRockStorage()
{
}

DistributionsRockStorage::~DistributionsRockStorage()
{
}

std::vector<DistributionsRock *>
DistributionsRockStorage::CreateDistributionsRockMix(const int                                                       & n_vintages,
                                                     const std::string                                               & path,
                                                     const std::vector<std::string>                                  & trend_cube_parameters,
                                                     const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                     const std::vector<std::string>                                  & constituent_label,
                                                     const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                     const std::map<std::string, DistributionsRockStorage *>         & model_rock_storage,
                                                     const std::map<std::string, DistributionsSolidStorage *>        & model_solid_storage,
                                                     const std::map<std::string, DistributionsDryRockStorage *>      & model_dry_rock_storage,
                                                     const std::map<std::string, DistributionsFluidStorage *>        & model_fluid_storage,
                                                     DEMTools::MixMethod                                               mix_method,
                                                     std::string                                                     & errTxt) const
{
  std::string tmpErrTxt = "";

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

  std::vector<double> s_min;
  std::vector<double> s_max;

  FindSMinMax(trend_cube_sampling, s_min, s_max);

  std::vector<DistributionsRock *>                   final_dist_rock(n_vintages, NULL);
  std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);

  for(int i=0; i<n_vintages; i++)
    all_volume_fractions[i].resize(n_constituents, NULL);

  for(int i=0; i<n_vintages; i++) {

    for (int s=0; s<n_constituents; s++) {

      if(i < n_vintages_constit[s]) {
        if(constituent_volume_fraction[s][i] != NULL)
          all_volume_fractions[i][s] = constituent_volume_fraction[s][i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, tmpErrTxt);
      }
      else {
        if(all_volume_fractions[i-1][s] != NULL)
          all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
      }
    }

    CheckVolumeConsistency(all_volume_fractions[i], errTxt);
  }

  bool mix_rock  = false;
  bool mix_fluid = false;
  bool mix_solid = false;

  std::vector<int> constituent_type(n_constituents);

  FindMixTypesForRock(constituent_label,
                      n_constituents,
                      model_rock_storage,
                      model_solid_storage,
                      model_dry_rock_storage,
                      model_fluid_storage,
                      mix_rock,
                      mix_solid,
                      mix_fluid,
                      constituent_type,
                      tmpErrTxt);


  if(mix_rock == true) {

    std::vector<std::vector<DistributionsRock *> > distr_rock(n_vintages);
    for(int i=0; i<n_vintages; i++)
      distr_rock[i].resize(n_constituents, NULL);

    for(int s=0; s<n_constituents; s++) {
      std::vector<DistributionsRock *> distr_rock_all_vintages = ReadRock(n_vintages,
                                                                          constituent_label[s],
                                                                          path,
                                                                          trend_cube_parameters,
                                                                          trend_cube_sampling,
                                                                          model_rock_storage,
                                                                          model_solid_storage,
                                                                          model_dry_rock_storage,
                                                                          model_fluid_storage,
                                                                          tmpErrTxt);

      int n_vintages_rock = static_cast<int>(distr_rock_all_vintages.size());

      for(int i=0; i<n_vintages; i++) {
        if(i < n_vintages_rock)
          distr_rock[i][s] = distr_rock_all_vintages[i];
        else
          distr_rock[i][s] = distr_rock[i-1][s]->Clone();
      }
    }

    for(int i=0; i<n_vintages; i++)
      final_dist_rock[i] = new DistributionsRockMixOfRock(distr_rock[i],
                                                          all_volume_fractions[i],
                                                          mix_method,
                                                          alpha,
                                                          s_min,
                                                          s_max);
  }

  else if(mix_fluid == true && mix_solid == true) {
    std::vector<std::vector<DistributionsFluid *> >    distr_fluid(n_vintages);
    std::vector<std::vector<DistributionsSolid *> >    distr_solid(n_vintages);
    std::vector<std::vector<DistributionWithTrend *> > fluid_volume_fractions(n_vintages);
    std::vector<std::vector<DistributionWithTrend *> > solid_volume_fractions(n_vintages);
    std::vector<double>                                solid_alpha;
    std::vector<double>                                fluid_alpha;

    int n_fluids = 0;
    int n_solids = 0;

    for(int s=0; s<n_constituents; s++) {

      if(constituent_type[s] == ModelSettings::FLUID) {

        std::vector<DistributionsFluid *> constit_fluid_all_vintages = ReadFluid(n_vintages,
                                                                                 constituent_label[s],
                                                                                 path,
                                                                                 trend_cube_parameters,
                                                                                 trend_cube_sampling,
                                                                                 model_fluid_storage,
                                                                                 errTxt);

        int n_vintages_fluid = static_cast<int>(constit_fluid_all_vintages.size());

        for(int i=0; i<n_vintages; i++)
          distr_fluid[i].resize(n_fluids+1);

        for(int i=0; i<n_vintages; i++) {
          if(i < n_vintages_fluid)
            distr_fluid[i][n_fluids] = constit_fluid_all_vintages[i];
          else
            distr_fluid[i][n_fluids] = distr_fluid[i-1][n_fluids]->Clone();
        }

        for(int i=0; i<n_vintages; i++)
          fluid_volume_fractions[i].push_back(all_volume_fractions[i][s]);

        fluid_alpha.push_back(alpha[s]);

        n_fluids++;

      }


      else {

        std::vector<DistributionsSolid *> constit_solid_all_vintages = ReadSolid(n_vintages,
                                                                                 constituent_label[s],
                                                                                 path,
                                                                                 trend_cube_parameters,
                                                                                 trend_cube_sampling,
                                                                                 model_solid_storage,
                                                                                 errTxt);

        int n_vintages_solid = static_cast<int>(constit_solid_all_vintages.size());

        for(int i=0; i<n_vintages_solid; i++)
          distr_solid[i].resize(n_solids+1);

        for(int i=0; i<n_vintages; i++) {
          if(i < n_vintages_solid)
            distr_solid[i][n_solids] = constit_solid_all_vintages[i];
          else
            distr_solid[i][n_solids] = distr_solid[i-1][n_solids]->Clone();
        }

        for(int i=0; i<n_vintages; i++)
          solid_volume_fractions[i].push_back(all_volume_fractions[i][s]);

        solid_alpha.push_back(alpha[s]);

        n_solids++;

      }

    }

    // Reuse alpha, but now it is sorted
    for(int i=0; i<n_solids; i++)
      alpha[i] = solid_alpha[i];
    for(int i=0; i<n_fluids; i++)
      alpha[i+n_solids] = fluid_alpha[i];

    for(int i=0; i<n_vintages; i++)
      final_dist_rock[i] = new DistributionsRockMixOfSolidAndFluid(distr_solid[i],
                                                                   distr_fluid[i],
                                                                   solid_volume_fractions[i],
                                                                   fluid_volume_fractions[i],
                                                                   mix_method,
                                                                   alpha,
                                                                   s_min,
                                                                   s_max);
  }

  return(final_dist_rock);
}

//----------------------------------------------------------------------------------//
TabulatedVelocityRockStorage::TabulatedVelocityRockStorage(std::vector<DistributionWithTrendStorage *> vp,
                                                           std::vector<DistributionWithTrendStorage *> vs,
                                                           std::vector<DistributionWithTrendStorage *> density,
                                                           double                                      correlation_vp_vs,
                                                           double                                      correlation_vp_density,
                                                           double                                      correlation_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density)
{
}

TabulatedVelocityRockStorage::~TabulatedVelocityRockStorage()
{
  if(vp_[0]->GetIsShared() == false)
    delete vp_[0];
  if(vs_[0]->GetIsShared() == false)
    delete vs_[0];
  if(density_[0]->GetIsShared() == false)
    delete density_[0];
}

std::vector<DistributionsRock *>
TabulatedVelocityRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                                        const std::string                                           & path,
                                                        const std::vector<std::string>                              & trend_cube_parameters,
                                                        const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                        const std::map<std::string, DistributionsRockStorage *>     & /*model_rock_storage*/,
                                                        const std::map<std::string, DistributionsSolidStorage *>    & /*model_solid_storage*/,
                                                        const std::map<std::string, DistributionsDryRockStorage *>  & /*model_dry_rock_storage*/,
                                                        const std::map<std::string, DistributionsFluidStorage *>    & /*model_fluid_storage*/,
                                                        std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Tabulated' rock physics model\n");

  std::vector<double> alpha(3);
  alpha[0] = vp_[0]     ->GetOneYearCorrelation();
  alpha[1] = vs_[0]     ->GetOneYearCorrelation();
  alpha[2] = density_[0]->GetOneYearCorrelation();

  std::vector<double> s_min;
  std::vector<double> s_max;

  FindSMinMax(trend_cube_sampling, s_min, s_max);

  int n_vintages_vp      = static_cast<int>(vp_.size());
  int n_vintages_vs      = static_cast<int>(vs_.size());
  int n_vintages_density = static_cast<int>(density_.size());

  std::vector<DistributionsRock *>     dist_rock(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vp_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> vs_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_vp)
      vp_dist_with_trend[i] = vp_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      vp_dist_with_trend[i] = vp_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_vs)
      vs_dist_with_trend[i] = vs_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      vs_dist_with_trend[i] = vs_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    DistributionsRock * rock = new DistributionsRockTabulated(vp_dist_with_trend[i],
                                                              vs_dist_with_trend[i],
                                                              density_dist_with_trend[i],
                                                              correlation_vp_vs_,
                                                              correlation_vp_density_,
                                                              correlation_vs_density_,
                                                              DEMTools::Velocity,
                                                              alpha,
                                                              s_min,
                                                              s_max);

    dist_rock[i] = rock;


  }

  return(dist_rock);

}

//----------------------------------------------------------------------------------//
TabulatedModulusRockStorage::TabulatedModulusRockStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                                                         std::vector<DistributionWithTrendStorage *> shear_modulus,
                                                         std::vector<DistributionWithTrendStorage *> density,
                                                         double                                      correlation_bulk_shear,
                                                         double                                      correlation_bulk_density,
                                                         double                                      correlation_shear_density)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  correlation_bulk_shear_(correlation_bulk_shear),
  correlation_bulk_density_(correlation_bulk_density),
  correlation_shear_density_(correlation_shear_density)
{
}

TabulatedModulusRockStorage::~TabulatedModulusRockStorage()
{
 if(bulk_modulus_[0]->GetIsShared() == false)
    delete bulk_modulus_[0];
  if(shear_modulus_[0]->GetIsShared() == false)
    delete shear_modulus_[0];
  if(density_[0]->GetIsShared() == false)
    delete density_[0];
}

std::vector<DistributionsRock *>
TabulatedModulusRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                                       const std::string                                           & path,
                                                       const std::vector<std::string>                              & trend_cube_parameters,
                                                       const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                       const std::map<std::string, DistributionsRockStorage *>     & /*model_rock_storage*/,
                                                       const std::map<std::string, DistributionsSolidStorage *>    & /*model_solid_storage*/,
                                                       const std::map<std::string, DistributionsDryRockStorage *>  & /*model_dry_rock_storage*/,
                                                       const std::map<std::string, DistributionsFluidStorage *>    & /*model_fluid_storage*/,
                                                       std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Tabulated' rock physics model\n");

  std::vector<double> alpha(3);
  alpha[0] = bulk_modulus_[0] ->GetOneYearCorrelation();
  alpha[1] = shear_modulus_[0]->GetOneYearCorrelation();
  alpha[2] = density_[0]      ->GetOneYearCorrelation();

  std::vector<double> s_min;
  std::vector<double> s_max;

  FindSMinMax(trend_cube_sampling, s_min, s_max);

  int n_vintages_bulk    = static_cast<int>(bulk_modulus_.size());
  int n_vintages_shear   = static_cast<int>(shear_modulus_.size());
  int n_vintages_density = static_cast<int>(density_.size());

  std::vector<DistributionsRock *>     dist_rock(n_vintages, NULL);
  std::vector<DistributionWithTrend *> bulk_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> shear_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> density_dist_with_trend(n_vintages, NULL);

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_bulk)
      bulk_dist_with_trend[i] = bulk_modulus_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      bulk_dist_with_trend[i] = bulk_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_shear)
      shear_dist_with_trend[i] = shear_modulus_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      shear_dist_with_trend[i] = shear_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_density)
      density_dist_with_trend[i] = density_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    else
      density_dist_with_trend[i] = density_dist_with_trend[i-1]->Clone();

    DistributionsRock * rock = new DistributionsRockTabulated(bulk_dist_with_trend[i],
                                                              shear_dist_with_trend[i],
                                                              density_dist_with_trend[i],
                                                              correlation_bulk_shear_,
                                                              correlation_bulk_density_,
                                                              correlation_shear_density_,
                                                              DEMTools::Modulus,
                                                              alpha,
                                                              s_min,
                                                              s_max);

    dist_rock[i] = rock;
  }

  return(dist_rock);
}
//----------------------------------------------------------------------------------//
ReussRockStorage::ReussRockStorage(std::vector<std::string>                                  constituent_label,
                                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussRockStorage::~ReussRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
      delete constituent_volume_fraction_[0][i];
  }
}

std::vector<DistributionsRock *>
ReussRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                            const std::string                                           & path,
                                            const std::vector<std::string>                              & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                            const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                            const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                            const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                            const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                            std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Reuss' rock physics model\n");

  std::vector<DistributionsRock *> rock;

  std::string tmpErrTxt = "";

  rock =   CreateDistributionsRockMix(n_vintages,
                                      path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      DEMTools::Reuss,
                                      tmpErrTxt);

  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Reuss rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
VoigtRockStorage::VoigtRockStorage(std::vector<std::string>                                  constituent_label,
                                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtRockStorage::~VoigtRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
      delete constituent_volume_fraction_[0][i];
  }
}

std::vector<DistributionsRock *>
VoigtRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                            const std::string                                           & path,
                                            const std::vector<std::string>                              & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                            const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                            const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                            const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                            const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                            std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Voigt' rock physics model\n");

  std::vector<DistributionsRock *> rock;

  std::string tmpErrTxt = "";

  rock =   CreateDistributionsRockMix(n_vintages,
                                      path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      DEMTools::Voigt,
                                      tmpErrTxt);


  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Voigt rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
HillRockStorage::HillRockStorage(std::vector<std::string>                                  constituent_label,
                                 std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillRockStorage::~HillRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_[0].size()); i++) {
    if(constituent_volume_fraction_[0][i]->GetIsShared() == false)
      delete constituent_volume_fraction_[0][i];
  }
}

std::vector<DistributionsRock *>
HillRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                           const std::string                                           & path,
                                           const std::vector<std::string>                              & trend_cube_parameters,
                                           const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                           const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                           const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                           const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                           const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                           std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Hill' rock physics model\n");

  std::vector<DistributionsRock *> rock;

  std::string tmpErrTxt = "";

  rock =   CreateDistributionsRockMix(n_vintages,
                                      path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      DEMTools::Hill,
                                      tmpErrTxt);

  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Hill rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
DEMRockStorage::DEMRockStorage(std::string                                               host_label,
                               std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                               std::vector<std::string>                                  inclusion_label,
                               std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                               std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMRockStorage::~DEMRockStorage()
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
}

std::vector<DistributionsRock *>
DEMRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                          const std::string                                           & path,
                                          const std::vector<std::string>                              & trend_cube_parameters,
                                          const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                          const std::map<std::string, DistributionsRockStorage *>     & /*model_rock_storage*/,
                                          const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                          const std::map<std::string, DistributionsDryRockStorage *>  & /*model_dry_rock_storage*/,
                                          const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                          std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'DEM' rock physics model\n");

    // Remember: Host info is included first in constituent vectors
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());
  int n_constituents = n_inclusions + 1;

  std::vector<std::vector<DistributionWithTrendStorage *> > volume_fractions(n_constituents);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=0; i<n_inclusions; i++)
    volume_fractions[i+1] = inclusion_volume_fraction_[i];

  // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions
  std::vector<double> alpha(n_inclusions + n_constituents);

  for(int i=0; i<n_inclusions; i++)
    alpha[i] = inclusion_aspect_ratio_[i][0]->GetOneYearCorrelation();

  for(int i=0; i<n_constituents; i++) {
    if(volume_fractions[i][0] != NULL)
      alpha[i+n_inclusions] = volume_fractions[i][0]->GetOneYearCorrelation();
    else
      alpha[i+n_inclusions] = 1;
  }

  std::vector<double> s_min;
  std::vector<double> s_max;

  FindSMinMax(trend_cube_sampling, s_min, s_max);

  std::vector<int> n_vintages_aspect(n_constituents);
  for(int i=0; i<n_inclusions; i++)
    n_vintages_aspect[i] = static_cast<int>(inclusion_aspect_ratio_[i].size());

  std::vector<int> n_vintages_volume(n_constituents);
  for(int i=0; i<n_constituents; i++)
    n_vintages_volume[i] = static_cast<int>(volume_fractions[i].size());

 //Read host label
  std::vector<DistributionsSolid *> final_distr_solid(n_vintages);

  std::vector<DistributionsSolid *> distr_solid = ReadSolid(n_vintages,
                                                            host_label_,
                                                            path,
                                                            trend_cube_parameters,
                                                            trend_cube_sampling,
                                                            model_solid_storage,
                                                            errTxt);

  int n_vintages_solid = static_cast<int>(distr_solid.size());

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_solid)
      final_distr_solid[i] = distr_solid[i];
    else
      final_distr_solid[i] = final_distr_solid[i-1]->Clone();
  }

  //Read inclusion label
  std::vector<std::vector<DistributionsFluid *> > final_distr_fluid_inc(n_vintages);
  for(int i=0; i<n_vintages; i++)
    final_distr_fluid_inc[i].resize(n_inclusions, NULL);

  for (int s = 0; s < n_inclusions; s++) {
    std::vector<DistributionsFluid *> distr_fluid_all_vintages = ReadFluid(n_vintages,
                                                                           inclusion_label_[s],
                                                                           path,
                                                                           trend_cube_parameters,
                                                                           trend_cube_sampling,
                                                                           model_fluid_storage,
                                                                           errTxt);

    int n_vintages_inc = static_cast<int>(distr_fluid_all_vintages.size());

    for(int i=0; i<n_vintages; i++) {
      if(i < n_vintages_inc)
        final_distr_fluid_inc[i][s] = distr_fluid_all_vintages[i];
      else
        final_distr_fluid_inc[i][s] = final_distr_fluid_inc[i-1][s]->Clone();
    }
  }

  std::vector<DistributionsRock *>                   final_dist_rock(n_vintages, NULL);
  std::vector<std::vector<DistributionWithTrend *> > all_volume_fractions(n_vintages);
  std::vector<std::vector<DistributionWithTrend *> > all_aspect_ratios(n_vintages);

  for(int i=0; i<n_vintages; i++) {
    all_volume_fractions[i].resize(n_constituents, NULL);
    all_aspect_ratios[i].resize(n_inclusions, NULL);
  }

  for(int i=0; i<n_vintages; i++) {

    for (int s = 0; s < n_inclusions; s++) {

      if(i < n_vintages_aspect[s]) {
        if(inclusion_aspect_ratio_[s][i] != NULL)
          all_aspect_ratios[i][s] = inclusion_aspect_ratio_[s][i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
      }
      else
        all_aspect_ratios[i][s] = all_aspect_ratios[i-1][s]->Clone();
    }

    for (int s = 0; s < n_constituents; s++) {

      if(i < n_vintages_volume[s]) {
        if(volume_fractions[s][i] != NULL)
          all_volume_fractions[i][s] = volume_fractions[s][i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
      }
      else
        all_volume_fractions[i][s] = all_volume_fractions[i-1][s]->Clone();
    }

    CheckVolumeConsistency(all_volume_fractions[i], errTxt);

    if (errTxt == "")
      final_dist_rock[i] = new DistributionsRockDEM(final_distr_solid[i],
                                                    final_distr_fluid_inc[i],
                                                    all_aspect_ratios[i],
                                                    all_volume_fractions[i],
                                                    alpha,
                                                    s_min,
                                                    s_max);
  }

  return(final_dist_rock);
}

//----------------------------------------------------------------------------------//
GassmannRockStorage::GassmannRockStorage(std::string dry_rock,
                                         std::string fluid)
: dry_rock_(dry_rock),
  fluid_(fluid)
{
}

GassmannRockStorage::~GassmannRockStorage()
{
}

std::vector<DistributionsRock *>
GassmannRockStorage::GenerateDistributionsRock(const int                                                   & /*n_vintages*/,
                                               const std::string                                           & /*path*/,
                                               const std::vector<std::string>                              & /*trend_cube_parameters*/,
                                               const std::vector<std::vector<double> >                     & /*trend_cube_sampling*/,
                                               const std::map<std::string, DistributionsRockStorage *>     & /*model_rock_storage*/,
                                               const std::map<std::string, DistributionsSolidStorage *>    & /*model_solid_storage*/,
                                               const std::map<std::string, DistributionsDryRockStorage *>  & /*model_dry_rock_storage*/,
                                               const std::map<std::string, DistributionsFluidStorage *>    & /*model_fluid_storage*/,
                                               std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Gassmann' rock physics model\n");

  std::vector<DistributionsRock *> rock(1, NULL);

  if(rock[0] == NULL)
    errTxt += "The Gassmann model has not been implemented yet for rocks\n";

  return(rock);
}

//----------------------------------------------------------------------------------//
BoundingRockStorage::BoundingRockStorage(std::string                                 upper_rock,
                                         std::string                                 lower_rock,
                                         std::vector<DistributionWithTrendStorage *> porosity,
                                         std::vector<DistributionWithTrendStorage *> bulk_weight,
                                         std::vector<DistributionWithTrendStorage *> p_wave_weight,
                                         double                                      correlation_weights)
: upper_rock_(upper_rock),
  lower_rock_(lower_rock),
  porosity_(porosity),
  bulk_weight_(bulk_weight),
  p_wave_weight_(p_wave_weight),
  correlation_weights_(correlation_weights)
{
}

BoundingRockStorage::~BoundingRockStorage()
{
  if(bulk_weight_[0]->GetIsShared() == false)
    delete bulk_weight_[0];
  if(p_wave_weight_[0]->GetIsShared() == false)
    delete p_wave_weight_[0];
}

std::vector<DistributionsRock *>
BoundingRockStorage::GenerateDistributionsRock(const int                                                   & n_vintages,
                                               const std::string                                           & path,
                                               const std::vector<std::string>                              & trend_cube_parameters,
                                               const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                               const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                               const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                               const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                               const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                               std::string                                                 & errTxt) const
{
  LogKit::LogFormatted(LogKit::Low," Generating 'Bounding' rock physics model\n");

  std::string tmpErrTxt = "";

  int n_vintages_porosity      = static_cast<int>(porosity_.size());
  int n_vintages_bulk_weight   = static_cast<int>(bulk_weight_.size());
  int n_vintages_p_wave_weight = static_cast<int>(p_wave_weight_.size());

  std::vector<double> alpha(3);
  alpha[0] = porosity_[0]     ->GetOneYearCorrelation();
  alpha[1] = bulk_weight_[0]  ->GetOneYearCorrelation();
  alpha[2] = p_wave_weight_[0]->GetOneYearCorrelation();

  std::vector<double> s_min;
  std::vector<double> s_max;

  FindSMinMax(trend_cube_sampling, s_min, s_max);

  std::vector<DistributionsRock *> final_distr_upper_rock(n_vintages);
  std::vector<DistributionsRock *> distr_upper_rock;

  LogKit::LogFormatted(LogKit::Low,"Upper bound:\n");
  distr_upper_rock = ReadRock(n_vintages,
                              upper_rock_,
                              path,
                              trend_cube_parameters,
                              trend_cube_sampling,
                              model_rock_storage,
                              model_solid_storage,
                              model_dry_rock_storage,
                              model_fluid_storage,
                              tmpErrTxt);

  for(int i=0; i<n_vintages; i++) {
    if(i < static_cast<int>(distr_upper_rock.size()))
      final_distr_upper_rock[i] = distr_upper_rock[i];
    else
      final_distr_upper_rock[i] = final_distr_upper_rock[i-1]->Clone();
  }

  if(distr_upper_rock[0]->GetIsOkForBounding() == false) {
    tmpErrTxt += "The upper bound does not follow the requirements for the bounding model.\n";
    tmpErrTxt += " The solid and fluid being mix need to be tabulated where the variables don't have distributions nor trends\n";
  }

  std::map<std::string, DistributionsRockStorage *>::const_iterator m_all = model_rock_storage.find(upper_rock_);
  if (m_all != model_rock_storage.end()) {
    if(typeid(*(m_all->second)) != typeid(VoigtRockStorage))
      tmpErrTxt += "The upper bound in the Bounding rock physics model needs to follow a Voigt model\n";
  }

  std::vector<DistributionsRock *> final_distr_lower_rock(n_vintages);
  std::vector<DistributionsRock *> distr_lower_rock;

  LogKit::LogFormatted(LogKit::Low,"Lower bound:\n");
  distr_lower_rock = ReadRock(n_vintages,
                              lower_rock_,
                              path,
                              trend_cube_parameters,
                              trend_cube_sampling,
                              model_rock_storage,
                              model_solid_storage,
                              model_dry_rock_storage,
                              model_fluid_storage,
                              tmpErrTxt);

  for(int i=0; i<n_vintages; i++) {
    if(i < static_cast<int>(distr_lower_rock.size()))
      final_distr_lower_rock[i] = distr_lower_rock[i];
    else
      final_distr_lower_rock[i] = final_distr_lower_rock[i-1]->Clone();
  }

  if(distr_lower_rock[0]->GetIsOkForBounding() == false) {
    tmpErrTxt += "The lower bound does not follow the requirements for the bounding model.\n";
    tmpErrTxt += " The solid and fluid being mix need to be tabulated where the variables don't have distributions nor trends\n";
  }

   std::map<std::string, DistributionsRockStorage *>::const_iterator n_all = model_rock_storage.find(lower_rock_);
  if (n_all != model_rock_storage.end()) {
    if(typeid(*(n_all->second)) != typeid(ReussRockStorage))
      tmpErrTxt += "The lower bound in the Bounding rock physics model needs to follow a Reuss model\n";
  }

  std::vector<DistributionsRock *>     dist_rock(n_vintages, NULL);
  std::vector<DistributionWithTrend *> porosity_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> bulk_weight_dist_with_trend(n_vintages, NULL);
  std::vector<DistributionWithTrend *> p_wave_weight_dist_with_trend(n_vintages, NULL);

  for(int i=0; i<n_vintages; i++) {
    if(i < n_vintages_porosity)
      porosity_dist_with_trend[i] = porosity_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, tmpErrTxt);
    else
      porosity_dist_with_trend[i] = porosity_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_bulk_weight)
      bulk_weight_dist_with_trend[i] = bulk_weight_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, tmpErrTxt);
    else
      bulk_weight_dist_with_trend[i] = bulk_weight_dist_with_trend[i-1]->Clone();

    if(i < n_vintages_p_wave_weight)
      p_wave_weight_dist_with_trend[i] = p_wave_weight_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, tmpErrTxt);
    else
      p_wave_weight_dist_with_trend[i] = p_wave_weight_dist_with_trend[i-1]->Clone();

    DistributionsRock * rock = new DistributionsRockBounding(final_distr_upper_rock[i],
                                                             final_distr_lower_rock[i],
                                                             porosity_dist_with_trend[i],
                                                             bulk_weight_dist_with_trend[i],
                                                             p_wave_weight_dist_with_trend[i],
                                                             correlation_weights_,
                                                             alpha,
                                                             s_min,
                                                             s_max);

    dist_rock[i] = rock;


  }

  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Bounding rock physics model:\n";
    errTxt += tmpErrTxt;
  }

  return(dist_rock);

}

