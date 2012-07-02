#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsfluidtabulatedvelocity.h"
#include "rplib/distributionwithtrendstorage.h"


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
TabulatedModulusFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                                         const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                         const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                         std::string                             & /*errTxt*/) const
{
  DistributionsFluid * fluid = NULL;

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
}

DistributionsFluid *
ReussFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsFluid * fluid = NULL; //new DistributionsFluidReuss();
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
}

DistributionsFluid *
VoigtFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                              std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsFluid * fluid = NULL; //new DistributionsFluidVoigt();
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
}

DistributionsFluid *
HillFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                             const std::vector<std::string>          & /*trend_cube_parameters*/,
                                             const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                             std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsFluid * fluid = NULL; //new DistributionsFluidHill();
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
}

DistributionsFluid *
BatzleWangFluidStorage::GenerateDistributionsFluid(const std::string                       & /*path*/,
                                                   const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                   const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                   std::string                             & /*errTxt*/) const
{
  //Gjør sjekk på om volume-fractions er double ved nytt kall. Feilmelding dersom ikke double.
  //Sjekk siste tallet i inclusion_volume_fraction, og pass på at det summeres til 1

  DistributionsFluid * fluid = NULL; //new DistributionsFluidBatzleWang();
  return(fluid);
}
