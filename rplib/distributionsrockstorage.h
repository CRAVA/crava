#ifndef RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsRockStorage {
public:
  DistributionsRockStorage();

  virtual ~DistributionsRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt)            const = 0;

protected:
  std::vector<DistributionsRock *> CreateDistributionsRockMix(const std::string                                               & path,
                                                              const std::vector<std::string>                                  & trend_cube_parameters,
                                                              const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                              const std::vector<std::string>                                  & constituent_label,
                                                              const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                              const std::map<std::string, DistributionsRockStorage *>         & model_rock_storage,
                                                              const std::map<std::string, DistributionsSolidStorage *>        & model_solid_storage,
                                                              const std::map<std::string, DistributionsDryRockStorage *>      & model_dry_rock_storage,
                                                              const std::map<std::string, DistributionsFluidStorage *>        & model_fluid_storage,
                                                              std::map<std::string, std::vector<DistributionsRock *> >        & rock_distribution,
                                                              std::map<std::string, std::vector<DistributionsSolid *> >       & solid_distribution,
                                                              std::map<std::string, std::vector<DistributionsDryRock *> >     & dry_rock_distribution,
                                                              std::map<std::string, std::vector<DistributionsFluid *> >       & fluid_distribution,
                                                              DEMTools::MixMethod                                               mix_method,
                                                              std::string                                                     & errTxt) const;

  };
//----------------------------------------------------------------------------------//
class TabulatedVelocityRockStorage : public DistributionsRockStorage {
public:
  TabulatedVelocityRockStorage(std::vector<DistributionWithTrendStorage *> vp,
                               std::vector<DistributionWithTrendStorage *> vs,
                               std::vector<DistributionWithTrendStorage *> density,
                               double                                      correlation_vp_vs,
                               double                                      correlation_vp_density,
                               double                                      correlation_vs_density);

  virtual ~TabulatedVelocityRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> vp_;
  std::vector<DistributionWithTrendStorage *> vs_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_vp_vs_;
  double                                      correlation_vp_density_;
  double                                      correlation_vs_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusRockStorage : public DistributionsRockStorage {
public:
  TabulatedModulusRockStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                              std::vector<DistributionWithTrendStorage *> shear_modulus,
                              std::vector<DistributionWithTrendStorage *> density,
                              double                                      correlation_bulk_shear,
                              double                                      correlation_bulk_density,
                              double                                      correlation_shear_density);

  virtual ~TabulatedModulusRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> bulk_modulus_;
  std::vector<DistributionWithTrendStorage *> shear_modulus_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_bulk_shear_;
  double                                      correlation_bulk_density_;
  double                                      correlation_shear_density_;
};

//----------------------------------------------------------------------------------//
class ReussRockStorage : public DistributionsRockStorage {
public:
  ReussRockStorage(std::vector<std::string>                                  constituent_label,
                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~ReussRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtRockStorage : public DistributionsRockStorage {
public:
  VoigtRockStorage(std::vector<std::string>                                  constituent_label,
                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~VoigtRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillRockStorage : public DistributionsRockStorage {
public:
  HillRockStorage(std::vector<std::string>                                  constituent_label,
                  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~HillRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class DEMRockStorage : public DistributionsRockStorage {
public:
  DEMRockStorage(std::string                                               host_label,
                 std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                 std::vector<std::string>                                  inclusion_label,
                 std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                 std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio);

  virtual ~DEMRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string                                               host_label_;
  std::vector<DistributionWithTrendStorage *>               host_volume_fraction_;
  std::vector<std::string>                                  inclusion_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio_;
};

//----------------------------------------------------------------------------------//
class GassmannRockStorage : public DistributionsRockStorage {
public:
  GassmannRockStorage(std::string dry_rock,
                      std::string fluid);

  virtual ~GassmannRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string dry_rock_;
  std::string fluid_;
};

//----------------------------------------------------------------------------------//
class BoundingRockStorage : public DistributionsRockStorage {
public:
  BoundingRockStorage(std::string                                 upper_rock,
                      std::string                                 lower_rock,
                      std::vector<DistributionWithTrendStorage *> porosity,
                      std::vector<DistributionWithTrendStorage *> bulk_weight,
                      std::vector<DistributionWithTrendStorage *> p_wave_weight,
                      double                                      correlation_weights);

  virtual ~BoundingRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     std::map<std::string, std::vector<DistributionsRock *> >    & rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsSolid *> >   & solid_distribution,
                                                                     std::map<std::string, std::vector<DistributionsDryRock *> > & dry_rock_distribution,
                                                                     std::map<std::string, std::vector<DistributionsFluid *> >   & fluid_distribution,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string                                 upper_rock_;
  std::string                                 lower_rock_;
  std::vector<DistributionWithTrendStorage *> porosity_;
  std::vector<DistributionWithTrendStorage *> bulk_weight_;
  std::vector<DistributionWithTrendStorage *> p_wave_weight_;
  double                                      correlation_weights_;

};

#endif
