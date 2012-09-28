#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
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
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"

#include <typeinfo>

DistributionsRockStorage::DistributionsRockStorage()
{
}

DistributionsRockStorage::~DistributionsRockStorage()
{
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
ReussRockStorage::GenerateDistributionsRock(const std::string                                          & /*path*/,
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
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsRock * rock = NULL; //new DistributionsRockReuss();

  if(rock == NULL)
    errTxt += "The Reuss model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
VoigtRockStorage::GenerateDistributionsRock(const std::string                                          & /*path*/,
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
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsRock * rock = NULL; //new DistributionsRockVoigt();

  if(rock == NULL)
    errTxt += "The Voigt model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
HillRockStorage::GenerateDistributionsRock(const std::string                                          & /*path*/,
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
  CheckVolumeConsistency(constituent_volume_fraction_, errTxt);

  DistributionsRock * rock = NULL; //new DistributionsRockHill();

  if(rock == NULL)
    errTxt += "The Hill model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
DEMRockStorage::GenerateDistributionsRock(const std::string                                          & /*path*/,
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
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());

  std::vector<DistributionWithTrendStorage *> volume_fractions(n_inclusions + 1);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=1; i<n_inclusions+1; i++)
    volume_fractions[i] = inclusion_volume_fraction_[i-1];

  CheckVolumeConsistency(volume_fractions, errTxt);

  DistributionsRock * rock = NULL; //new DistributionsRockInclusion();

  if(rock == NULL)
    errTxt += "The DEM model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
                                         DistributionWithTrendStorage * bulk_weight,
                                         DistributionWithTrendStorage * p_wave_weight,
                                         double                         correlation_weights)
: upper_rock_(upper_rock),
  lower_rock_(lower_rock),
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

  std::string upper_type = typeid(distr_upper_rock).name();
  if(upper_type != "voigt")
    tmpErrTxt += "The upper bound in the Bounding rock physics model needs follow a Voigt model";

  std::string lower_type = typeid(distr_lower_rock).name();
  if(lower_type != "reuss")
    tmpErrTxt += "The lower bound in the Bounding rock physics model needs follow a Reuss model";

  //Sjekk at porosity følger fordeling
  //Sjekk at porosity er delt reservoirvariabel

  DistributionWithTrend * distr_bulk_weight    = bulk_weight_  ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  DistributionWithTrend * distr_p_wave_weight = p_wave_weight_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  if(distr_bulk_weight->GetIsDistribution() == false)
    tmpErrTxt += "The bulk-modulus-weights need to follow a distribution in the Bounding rock physics model\n";

  if(distr_p_wave_weight->GetIsDistribution() == false)
    tmpErrTxt += "The p-wave-weights need to follow a distribution in the Bounding rock physics model\n";


  //DistributionsRock * rock = NULL; //new DistributionsRockBounding();


  if(tmpErrTxt != "") {
    errTxt += "Problems with the Bounding rock physics model:\n";
    errTxt += tmpErrTxt;
  }

  return(distr_upper_rock);
}
