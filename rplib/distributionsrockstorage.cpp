#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionsrocktabulated.h"
#include "rplib/distributionsrocktrinormal.h"
#include "rplib/distributionwithtrendstorage.h"


DistributionsRockStorage::DistributionsRockStorage()
{
}

DistributionsRockStorage::~DistributionsRockStorage()
{
}

//----------------------------------------------------------------------------------//
TabulatedRockStorage::TabulatedRockStorage(DistributionWithTrendStorage * vp,
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

TabulatedRockStorage::~TabulatedRockStorage()
{
  delete vp_;
  delete vs_;
  delete density_;
  delete correlation_vp_vs_;
  delete correlation_vp_density_;
  delete correlation_vs_density_;
}

DistributionsRock *
TabulatedRockStorage::GenerateDistributionsRock(const std::string                       & path,
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

    rock = new DistributionsRockTabulated(vp_dist_with_trend,
                                                              vs_dist_with_trend,
                                                              density_dist_with_trend,
                                                              corr_vp_vs_dist_with_trend,
                                                              corr_vp_density_dist_with_trend,
                                                              corr_vs_density_dist_with_trend);
  }

  return(rock);
}

//----------------------------------------------------------------------------------//

DEMRockStorage::DEMRockStorage(std::string                                 host_label,
                               DistributionWithTrendStorage *              host_volume_fraction,
                               std::vector<std::string>                    inclusion_label,
                               std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction)
{
}

DEMRockStorage::~DEMRockStorage()
{
}

DistributionsRock *
DEMRockStorage::GenerateDistributionsRock(const std::string                       & /*path*/,
                                          const std::vector<std::string>          & /*trend_cube_parameters*/,
                                          const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                          std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsRock * rock = NULL; //new DistributionsRockInclusion();
  return(rock);
}
