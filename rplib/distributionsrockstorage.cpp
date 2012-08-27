#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionsrocktabulatedvelocity.h"
#include "rplib/distributionsrocktrinormal.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsstoragekit.h"


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
                                                           DistributionWithTrendStorage * correlation_vp_vs,
                                                           DistributionWithTrendStorage * correlation_vp_density,
                                                           DistributionWithTrendStorage * correlation_vs_density)
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
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsRock *
TabulatedVelocityRockStorage::GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const
{
  bool is_gaussian = false;


  if(vp_                    ->GetIsGaussian() == true &&
     vs_                    ->GetIsGaussian() == true &&
     density_               ->GetIsGaussian() == true &&
     correlation_vp_vs_     ->GetIsGaussian() == false &&
     correlation_vp_density_->GetIsGaussian() == false &&
     correlation_vs_density_->GetIsGaussian() == false) {

       is_gaussian = true;
  }

  DistributionsRock * rock = NULL;

  if(is_gaussian) {
    const NRLib::TrendStorage * mean_vp                = vp_                    ->CloneMean();
    const NRLib::TrendStorage * mean_vs                = vs_                    ->CloneMean();
    const NRLib::TrendStorage * mean_density           = density_               ->CloneMean();
    const NRLib::TrendStorage * variance_vp            = vp_                    ->CloneVariance();
    const NRLib::TrendStorage * variance_vs            = vs_                    ->CloneVariance();
    const NRLib::TrendStorage * variance_density       = density_               ->CloneVariance();
    const NRLib::TrendStorage * correlation_vp_vs      = correlation_vp_vs_     ->CloneMean();
    const NRLib::TrendStorage * correlation_vp_density = correlation_vp_density_->CloneMean();
    const NRLib::TrendStorage * correlation_vs_density = correlation_vs_density_->CloneMean();



    NRLib::Trend * mean_vp_trend                = mean_vp               ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * mean_vs_trend                = mean_vs               ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * mean_density_trend           = mean_density          ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_vp_trend            = variance_vp           ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_vs_trend            = variance_vs           ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_density_trend       = variance_density      ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * correlation_vp_vs_trend      = correlation_vp_vs     ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * correlation_vp_density_trend = correlation_vp_density->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * correlation_vs_density_trend = correlation_vs_density->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

    rock = new DistributionsRockTriNormal(mean_vp_trend,
                                          mean_vs_trend,
                                          mean_density_trend,
                                          variance_vp_trend,
                                          variance_vs_trend,
                                          variance_density_trend,
                                          correlation_vp_vs_trend,
                                          correlation_vp_density_trend,
                                          correlation_vs_density_trend);

    delete mean_vp;
    delete mean_vs;
    delete mean_density;
    delete variance_vp;
    delete variance_vs;
    delete variance_density;
    delete correlation_vp_vs;
    delete correlation_vp_density;
    delete correlation_vs_density;

    delete mean_vp_trend;
    delete mean_vs_trend;
    delete mean_density_trend;
    delete variance_vp_trend;
    delete variance_vs_trend;
    delete variance_density_trend;
    delete correlation_vp_vs_trend;
    delete correlation_vp_density_trend;
    delete correlation_vs_density_trend;
  }

  else {
    const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    const DistributionWithTrend * vs_dist_with_trend              = vs_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    const DistributionWithTrend * corr_vp_vs_dist_with_trend      = correlation_vp_vs_     ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    const DistributionWithTrend * corr_vp_density_dist_with_trend = correlation_vp_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
    const DistributionWithTrend * corr_vs_density_dist_with_trend = correlation_vs_density_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

    rock = new DistributionsRockTabulatedVelocity(vp_dist_with_trend,
                                                  vs_dist_with_trend,
                                                  density_dist_with_trend,
                                                  corr_vp_vs_dist_with_trend,
                                                  corr_vp_density_dist_with_trend,
                                                  corr_vs_density_dist_with_trend);
  }

  return(rock);
}

//----------------------------------------------------------------------------------//
TabulatedModulusRockStorage::TabulatedModulusRockStorage(DistributionWithTrendStorage * bulk_modulus,
                                                         DistributionWithTrendStorage * shear_modulus,
                                                         DistributionWithTrendStorage * density,
                                                         DistributionWithTrendStorage * correlation_bulk_shear,
                                                         DistributionWithTrendStorage * correlation_bulk_density,
                                                         DistributionWithTrendStorage * correlation_shear_density)
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
  delete correlation_bulk_shear_;
  delete correlation_bulk_density_;
  delete correlation_shear_density_;
}

DistributionsRock *
TabulatedModulusRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                                       const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                       const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                       std::string                             & errTxt) const
{
  DistributionsRock * rock = NULL;

  if(rock == NULL)
    errTxt += "The tabulated moduli model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

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
    if(constituent_volume_fraction_[i]->GetIsSheared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
ReussRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                            const std::vector<std::string>          & /*trend_cube_parameters*/,
                                            const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                            std::string                             & errTxt) const
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
    if(constituent_volume_fraction_[i]->GetIsSheared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
VoigtRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                            const std::vector<std::string>          & /*trend_cube_parameters*/,
                                            const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                            std::string                             & errTxt) const
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
    if(constituent_volume_fraction_[i]->GetIsSheared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsRock *
HillRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                           const std::vector<std::string>          & /*trend_cube_parameters*/,
                                           const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                           std::string                             & errTxt) const
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
  if(host_volume_fraction_->GetIsSheared() == false)
    delete host_volume_fraction_;

  if(host_aspect_ratio_->GetIsSheared() == false)
    delete host_aspect_ratio_;

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_.size()); i++) {
    if(inclusion_volume_fraction_[i]->GetIsSheared() == false)
      delete inclusion_volume_fraction_[i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_.size()); i++) {
    if(inclusion_aspect_ratio_[i]->GetIsSheared() == false)
      delete inclusion_aspect_ratio_[i];
  }
}

DistributionsRock *
DEMRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                          const std::vector<std::string>          & /*trend_cube_parameters*/,
                                          const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                          std::string                             & errTxt) const
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
GassmannRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                               const std::vector<std::string>          & /*trend_cube_parameters*/,
                                               const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                               std::string                             & errTxt) const
{
  DistributionsRock * rock = NULL; //new DistributionsRockGassmann();

  if(rock == NULL)
    errTxt += "The Gassmann model has not been implemented yet for rocks\n"; //Marit: Denne feilmeldingen fjernes når modellen er implementert

  return(rock);
}
