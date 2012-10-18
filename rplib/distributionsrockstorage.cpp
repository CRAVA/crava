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

#include "rplib/distributionsrocktabulatedvelocity.h"
#include "rplib/distributionsrocktabulatedmodulus.h"
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

DistributionsRock *
DistributionsRockStorage::CreateDistributionsRockMix(const std::string                                          & path,
                                                     const std::vector<std::string>                             & trend_cube_parameters,
                                                     const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                                     const std::vector<std::string>                             & constituent_label,
                                                     const std::vector<DistributionWithTrendStorage *>          & constituent_volume_fraction,
                                                     const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                                                     const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                                     const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                                                     const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                                     std::map<std::string, DistributionsRock *>                 & rock_distribution,
                                                     std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                                     std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
                                                     std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                                     DEMTools::MixMethod                                          mix_method,
                                                     std::string                                                & errTxt) const
{
  DistributionsRock * rock = NULL;

  std::string tmpErrTxt = "";

  int n_constituents = static_cast<int>(constituent_label.size());

  std::vector<DistributionWithTrend *> distr_volume_fraction(n_constituents);
  for(int i=0; i<n_constituents; i++) {
    if(constituent_volume_fraction[i] == NULL)
      distr_volume_fraction[i] = NULL;
    else
      distr_volume_fraction[i] = constituent_volume_fraction[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  }

  CheckVolumeConsistency(distr_volume_fraction, tmpErrTxt);

  std::vector<DistributionsRock *> distributions_rock(n_constituents);
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
    std::vector<DistributionsRock *> distr_rock(n_constituents);
    for(int i=0; i<n_constituents; i++) {
      distr_rock[i] = ReadRock(constituent_label[i],
                               path,
                               trend_cube_parameters,
                               trend_cube_sampling,
                               model_rock_storage,
                               model_solid_storage,
                               model_dry_rock_storage,
                               model_fluid_storage,
                               rock_distribution,
                               solid_distribution,
                               dry_rock_distribution,
                               fluid_distribution,
                               tmpErrTxt);
    }

    rock = new DistributionsRockMixOfRock(distr_rock, distr_volume_fraction, mix_method);
  }
  else if(mix_fluid == true && mix_solid == true) {
    std::vector<DistributionsFluid *>    distr_fluid;
    std::vector<DistributionsSolid *>    distr_solid;
    std::vector<DistributionWithTrend *> fluid_volume_fractions;
    std::vector<DistributionWithTrend *> solid_volume_fractions;

    for(int i=0; i<n_constituents; i++) {

      if(constituent_type[i] == ModelSettings::FLUID) {
        DistributionsFluid * constit_fluid = NULL;

        constit_fluid = ReadFluid(constituent_label[i],
                                  path,
                                  trend_cube_parameters,
                                  trend_cube_sampling,
                                  model_fluid_storage,
                                  fluid_distribution,
                                  errTxt);

        distr_fluid.push_back(constit_fluid);

        fluid_volume_fractions.push_back(distr_volume_fraction[i]);
      }
      else {
        DistributionsSolid * constit_solid = NULL;

        constit_solid = ReadSolid(constituent_label[i],
                                  path,
                                  trend_cube_parameters,
                                  trend_cube_sampling,
                                  model_solid_storage,
                                  solid_distribution,
                                  errTxt);

        distr_solid.push_back(constit_solid);

        solid_volume_fractions.push_back(distr_volume_fraction[i]);
      }
    }

    rock = new DistributionsRockMixOfSolidAndFluid(distr_solid, distr_fluid, solid_volume_fractions, fluid_volume_fractions, mix_method);
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
TabulatedVelocityRockStorage::TabulatedVelocityRockStorage(DistributionWithTrendStorage * vp,
                                                           DistributionWithTrendStorage * vs,
                                                           DistributionWithTrendStorage * density,
                                                           double                         correlation_vp_vs,
                                                           double                         correlation_vp_density,
                                                           double                         correlation_vs_density)
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
  if(vp_->GetIsShared() == false)
    delete vp_;
  if(vs_->GetIsShared() == false)
    delete vs_;
  if(density_->GetIsShared() == false)
    delete density_;
}

DistributionsRock *
TabulatedVelocityRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                                        const std::vector<std::string>                             & trend_cube_parameters,
                                                        const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                                        const std::map<std::string, DistributionsRockStorage *>    & /*model_rock_storage*/,
                                                        const std::map<std::string, DistributionsSolidStorage *>   & /*model_solid_storage*/,
                                                        const std::map<std::string, DistributionsDryRockStorage *> & /*model_dry_rock_storage*/,
                                                        const std::map<std::string, DistributionsFluidStorage *>   & /*model_fluid_storage*/,
                                                        std::map<std::string, DistributionsRock *>                 & /*rock_distribution*/,
                                                        std::map<std::string, DistributionsSolid *>                & /*solid_distribution*/,
                                                        std::map<std::string, DistributionsDryRock *>              & /*dry_rock_distribution*/,
                                                        std::map<std::string, DistributionsFluid *>                & /*fluid_distribution*/,
                                                        std::string                                                & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * vs_dist_with_trend              = vs_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsRock * rock = new DistributionsRockTabulatedVelocity(vp_dist_with_trend,
                                                                    vs_dist_with_trend,
                                                                    density_dist_with_trend,
                                                                    correlation_vp_vs_,
                                                                    correlation_vp_density_,
                                                                    correlation_vs_density_);

  return(rock);
}

//----------------------------------------------------------------------------------//
TabulatedModulusRockStorage::TabulatedModulusRockStorage(DistributionWithTrendStorage * bulk_modulus,
                                                         DistributionWithTrendStorage * shear_modulus,
                                                         DistributionWithTrendStorage * density,
                                                         double                         correlation_bulk_shear,
                                                         double                         correlation_bulk_density,
                                                         double                         correlation_shear_density)
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
  delete bulk_modulus_;
  delete shear_modulus_;
  delete density_;
}

DistributionsRock *
TabulatedModulusRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                                       const std::vector<std::string>                             & trend_cube_parameters,
                                                       const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                                       const std::map<std::string, DistributionsRockStorage *>    & /*model_rock_storage*/,
                                                       const std::map<std::string, DistributionsSolidStorage *>   & /*model_solid_storage*/,
                                                       const std::map<std::string, DistributionsDryRockStorage *> & /*model_dry_rock_storage*/,
                                                       const std::map<std::string, DistributionsFluidStorage *>   & /*model_fluid_storage*/,
                                                       std::map<std::string, DistributionsRock *>                 & /*rock_distribution*/,
                                                       std::map<std::string, DistributionsSolid *>                & /*solid_distribution*/,
                                                       std::map<std::string, DistributionsDryRock *>              & /*dry_rock_distribution*/,
                                                       std::map<std::string, DistributionsFluid *>                & /*fluid_distribution*/,
                                                       std::string                                                & errTxt) const
{

  const DistributionWithTrend * bulk_with_trend    = bulk_modulus_ ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * shear_with_trend   = shear_modulus_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_with_trend = density_      ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsRock * rock = new DistributionsRockTabulatedModulus(bulk_with_trend,
                                                                   shear_with_trend,
                                                                   density_with_trend,
                                                                   correlation_bulk_shear_,
                                                                   correlation_bulk_density_,
                                                                   correlation_shear_density_);

  return(rock);
}
//----------------------------------------------------------------------------------//
ReussRockStorage::ReussRockStorage(std::vector<std::string>                    constituent_label,
                                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussRockStorage::~ReussRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
ReussRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                            const std::vector<std::string>                             & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                            const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                                            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                            const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                                            const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                            std::map<std::string, DistributionsRock *>                 & rock_distribution,
                                            std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                            std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
                                            std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                            std::string                                                & errTxt) const
{

  DistributionsRock * rock = NULL;

  std::string tmpErrTxt    = "";

  rock =   CreateDistributionsRockMix(path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      rock_distribution,
                                      solid_distribution,
                                      dry_rock_distribution,
                                      fluid_distribution,
                                      DEMTools::Reuss,
                                      tmpErrTxt);

  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Reuss rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
VoigtRockStorage::VoigtRockStorage(std::vector<std::string>                    constituent_label,
                                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtRockStorage::~VoigtRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
VoigtRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                            const std::vector<std::string>                             & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                            const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                                            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                            const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                                            const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                            std::map<std::string, DistributionsRock *>                 & rock_distribution,
                                            std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                            std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
                                            std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                            std::string                                                & errTxt) const
{

  DistributionsRock * rock = NULL;

  std::string tmpErrTxt    = "";

  rock =   CreateDistributionsRockMix(path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      rock_distribution,
                                      solid_distribution,
                                      dry_rock_distribution,
                                      fluid_distribution,
                                      DEMTools::Voigt,
                                      tmpErrTxt);


  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Voigt rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
HillRockStorage::HillRockStorage(std::vector<std::string>                    constituent_label,
                                 std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillRockStorage::~HillRockStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
HillRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                            const std::vector<std::string>                             & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                            const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                                            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                            const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                                            const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                            std::map<std::string, DistributionsRock *>                 & rock_distribution,
                                            std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                            std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
                                            std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                            std::string                                                & errTxt) const
{

  DistributionsRock * rock = NULL;

  std::string tmpErrTxt    = "";

  rock =   CreateDistributionsRockMix(path,
                                      trend_cube_parameters,
                                      trend_cube_sampling,
                                      constituent_label_,
                                      constituent_volume_fraction_,
                                      model_rock_storage,
                                      model_solid_storage,
                                      model_dry_rock_storage,
                                      model_fluid_storage,
                                      rock_distribution,
                                      solid_distribution,
                                      dry_rock_distribution,
                                      fluid_distribution,
                                      DEMTools::Hill,
                                      tmpErrTxt);

  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Hill rock physics model for <rock>:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
DEMRockStorage::DEMRockStorage(std::string                                 host_label,
                               DistributionWithTrendStorage *              host_volume_fraction,
                               DistributionWithTrendStorage *              host_aspect_ratio,
                               std::vector<std::string>                    inclusion_label,
                               std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                               std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  host_aspect_ratio_(host_aspect_ratio),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMRockStorage::~DEMRockStorage()
{
  if(host_volume_fraction_->GetIsShared() == false)
    delete host_volume_fraction_;

  if(host_aspect_ratio_->GetIsShared() == false)
    delete host_aspect_ratio_;

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_.size()); i++) {
    if(inclusion_volume_fraction_[i]->GetIsShared() == false)
      delete inclusion_volume_fraction_[i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_.size()); i++) {
    if(inclusion_aspect_ratio_[i]->GetIsShared() == false)
      delete inclusion_aspect_ratio_[i];
  }
}

DistributionsRock *
DEMRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                          const std::vector<std::string>                             & trend_cube_parameters,
                                          const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                          const std::map<std::string, DistributionsRockStorage *>    & /*model_rock_storage*/,
                                          const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                          const std::map<std::string, DistributionsDryRockStorage *> & /*model_dry_rock_storage*/,
                                          const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                          std::map<std::string, DistributionsRock *>                 & /*rock_distribution*/,
                                          std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                          std::map<std::string, DistributionsDryRock *>              & /*dry_rock_distribution*/,
                                          std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                          std::string                                                & errTxt) const
{
  // Remember: Host info is included first in inclusion vectors
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());

  std::vector<DistributionWithTrendStorage *> volume_fractions(n_inclusions + 1);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=1; i<n_inclusions+1; i++)
    volume_fractions[i] = inclusion_volume_fraction_[i-1];

  DistributionsRock * rock = NULL;

  std::vector< DistributionWithTrend *> inclusion_volume_fraction_distr(inclusion_volume_fraction_.size()+1, NULL);
  std::vector< DistributionWithTrend *> inclusion_aspect_ratio_distr(inclusion_aspect_ratio_.size()+1, NULL);

  inclusion_volume_fraction_distr[0]  = host_volume_fraction_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  inclusion_aspect_ratio_distr[0]     = host_aspect_ratio_   ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt); //NBNB fjellvoll host har altså aspect ratio.

  for (size_t i = 1; i < inclusion_volume_fraction_.size(); ++i)
    inclusion_volume_fraction_distr[i] = inclusion_volume_fraction_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  for (size_t i = 1; i < inclusion_aspect_ratio_.size(); ++i)
    inclusion_aspect_ratio_distr[i] = inclusion_aspect_ratio_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  //Read host label
  DistributionsSolid  * final_distr_solid = NULL;
  final_distr_solid = ReadSolid(host_label_,path, trend_cube_parameters, trend_cube_sampling, model_solid_storage, solid_distribution, errTxt);

  //Read inclusion label
  std::vector< DistributionsFluid* > final_distr_fluid_inc;
  size_t s;
  for (s = 0; s != inclusion_label_.size(); ++s) {
    DistributionsFluid * incl_fluid;
    incl_fluid = ReadFluid(inclusion_label_[s],path, trend_cube_parameters, trend_cube_sampling, model_fluid_storage, fluid_distribution, errTxt);
    final_distr_fluid_inc.push_back(incl_fluid);
  }

  //Questions //NBNB fjellvoll //NBNB marit
  //1.Do we support more than one inclusion?

  //CheckVolumeConsistency(distr_porosity, errTxt); //Fix when questions are solved

  if (errTxt == "") {
    rock = new DistributionsRockDEM(final_distr_solid,
                                    final_distr_fluid_inc[0],         //tmp solution only single inclusion supported
                                    inclusion_aspect_ratio_distr,
                                    inclusion_volume_fraction_distr);
  }

  return(rock);
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

DistributionsRock *
GassmannRockStorage::GenerateDistributionsRock(const std::string                                          & /*path*/,
                                               const std::vector<std::string>                             & /*trend_cube_parameters*/,
                                               const std::vector<std::vector<double> >                    & /*trend_cube_sampling*/,
                                               const std::map<std::string, DistributionsRockStorage *>    & /*model_rock_storage*/,
                                               const std::map<std::string, DistributionsSolidStorage *>   & /*model_solid_storage*/,
                                               const std::map<std::string, DistributionsDryRockStorage *> & /*model_dry_rock_storage*/,
                                               const std::map<std::string, DistributionsFluidStorage *>   & /*model_fluid_storage*/,
                                               std::map<std::string, DistributionsRock *>                 & /*rock_distribution*/,
                                               std::map<std::string, DistributionsSolid *>                & /*solid_distribution*/,
                                               std::map<std::string, DistributionsDryRock *>              & /*dry_rock_distribution*/,
                                               std::map<std::string, DistributionsFluid *>                & /*fluid_distribution*/,
                                               std::string                                                & errTxt) const
{
  DistributionsRock * rock = NULL; //new DistributionsRockGassmann();

  if(rock == NULL)
    errTxt += "The Gassmann model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(rock);
}

//----------------------------------------------------------------------------------//
BoundingRockStorage::BoundingRockStorage(std::string                    upper_rock,
                                         std::string                    lower_rock,
                                         DistributionWithTrendStorage * porosity,
                                         DistributionWithTrendStorage * bulk_weight,
                                         DistributionWithTrendStorage * p_wave_weight,
                                         double                         correlation_weights)
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
  if(bulk_weight_->GetIsShared() == false)
    delete bulk_weight_;
  if(p_wave_weight_->GetIsShared() == false)
    delete p_wave_weight_;
}

DistributionsRock *
BoundingRockStorage::GenerateDistributionsRock(const std::string                                          & path,
                                               const std::vector<std::string>                             & trend_cube_parameters,
                                               const std::vector<std::vector<double> >                    & trend_cube_sampling,
                                               const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                                               const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                                               const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                                               const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                                               std::map<std::string, DistributionsRock *>                 & rock_distribution,
                                               std::map<std::string, DistributionsSolid *>                & solid_distribution,
                                               std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
                                               std::map<std::string, DistributionsFluid *>                & fluid_distribution,
                                               std::string                                                & errTxt) const
{
  std::string tmpErrTxt = "";

  DistributionsRock   * distr_upper_rock = NULL;

  distr_upper_rock = ReadRock(upper_rock_,
                              path,
                              trend_cube_parameters,
                              trend_cube_sampling,
                              model_rock_storage,
                              model_solid_storage,
                              model_dry_rock_storage,
                              model_fluid_storage,
                              rock_distribution,
                              solid_distribution,
                              dry_rock_distribution,
                              fluid_distribution,
                              tmpErrTxt);

  if(distr_upper_rock->GetIsOkForBounding() == false) {
    tmpErrTxt += "The upper bound does not follow the requirements for the bounding model.\n";
    tmpErrTxt += " The solid and fluid being mix need to be tabulated where the variables don't have distributions nor trends\n";
  }

  DistributionsRock   * distr_lower_rock = NULL;

  distr_lower_rock = ReadRock(lower_rock_,
                              path,
                              trend_cube_parameters,
                              trend_cube_sampling,
                              model_rock_storage,
                              model_solid_storage,
                              model_dry_rock_storage,
                              model_fluid_storage,
                              rock_distribution,
                              solid_distribution,
                              dry_rock_distribution,
                              fluid_distribution,
                              tmpErrTxt);

  if(distr_lower_rock->GetIsOkForBounding() == false) {
    tmpErrTxt += "The lower bound does not follow the requirements for the bounding model.\n";
    tmpErrTxt += " The solid and fluid being mix need to be tabulated where the variables don't have distributions nor trends\n";
  }

  std::string upper_type = typeid(distr_upper_rock).name(); //Marit: Denne testen fungerer ikke
  if(upper_type != "class DistributionsRock *")
    tmpErrTxt += "The upper bound in the Bounding rock physics model needs to follow a Voigt model\n";

  std::string lower_type = typeid(distr_lower_rock).name(); //Marit: Denne testen fungerer ikke.
  if(lower_type != "class DistributionsRock *")
    tmpErrTxt += "The lower bound in the Bounding rock physics model needs to follow a Reuss model\n";

  DistributionWithTrend * distr_porosity      = porosity_     ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  DistributionWithTrend * distr_bulk_weight   = bulk_weight_  ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  DistributionWithTrend * distr_p_wave_weight = p_wave_weight_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsRock * rock = new DistributionsRockBounding(distr_upper_rock, distr_lower_rock, distr_porosity, distr_bulk_weight, distr_p_wave_weight, correlation_weights_);


  if(tmpErrTxt != "") {
    errTxt += "\nProblems with the Bounding rock physics model:\n";
    errTxt += tmpErrTxt;
  }

  return(rock);
}

