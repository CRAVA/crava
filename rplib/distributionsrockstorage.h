#ifndef RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP

#include "rplib/demmodelling.h"

#include <map>

class DistributionWithTrendStorage;
class DistributionsSolidStorage;
class DistributionsDryRockStorage;
class DistributionsFluidStorage;

class DistributionsRock;
class DistributionsFluid;

class BlockedLogsForRockPhysics;

class DistributionsRockStorage {
public:
  DistributionsRockStorage();

  virtual ~DistributionsRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt)            const = 0;
protected:
  std::vector<DistributionsRock *> CreateDistributionsRockMix(const int                                                       & n_vintages,
                                                              const std::string                                               & path,
                                                              const std::vector<std::string>                                  & trend_cube_parameters,
                                                              const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                              const std::vector<BlockedLogsForRockPhysics *>                  & blockedLogs,
                                                              const std::vector<std::string>                                  & constituent_label,
                                                              const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                              const std::map<std::string, DistributionsRockStorage *>         & model_rock_storage,
                                                              const std::map<std::string, DistributionsSolidStorage *>        & model_solid_storage,
                                                              const std::map<std::string, DistributionsDryRockStorage *>      & model_dry_rock_storage,
                                                              const std::map<std::string, DistributionsFluidStorage *>        & model_fluid_storage,
                                                              const int                                                         output_other,
                                                              DEMTools::MixMethod                                               mix_method,
                                                              std::string                                                     & errTxt) const;

  };
//----------------------------------------------------------------------------------//
class TabulatedVelocityRockStorage : public DistributionsRockStorage {
public:
  TabulatedVelocityRockStorage(std::vector<DistributionWithTrendStorage *> vp,
                               std::vector<DistributionWithTrendStorage *> vs,
                               std::vector<DistributionWithTrendStorage *> density,
                               std::vector<DistributionWithTrendStorage *> correlation_vp_vs,
                               std::vector<DistributionWithTrendStorage *> correlation_vp_density,
                               std::vector<DistributionWithTrendStorage *> correlation_vs_density,
                               std::string                                 rock_name);

  virtual ~TabulatedVelocityRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> vp_;
  std::vector<DistributionWithTrendStorage *> vs_;
  std::vector<DistributionWithTrendStorage *> density_;
  std::vector<DistributionWithTrendStorage *> correlation_vp_vs_;        /// Converted to double
  std::vector<DistributionWithTrendStorage *> correlation_vp_density_;   /// Converted to double
  std::vector<DistributionWithTrendStorage *> correlation_vs_density_;   /// Converted to double
  std::string                                 rock_name_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusRockStorage : public DistributionsRockStorage {
public:
  TabulatedModulusRockStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                              std::vector<DistributionWithTrendStorage *> shear_modulus,
                              std::vector<DistributionWithTrendStorage *> density,
                              std::vector<DistributionWithTrendStorage *> correlation_bulk_shear,
                              std::vector<DistributionWithTrendStorage *> correlation_bulk_density,
                              std::vector<DistributionWithTrendStorage *> correlation_shear_density,
                              std::string                                 rock_name);

  virtual ~TabulatedModulusRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> bulk_modulus_;
  std::vector<DistributionWithTrendStorage *> shear_modulus_;
  std::vector<DistributionWithTrendStorage *> density_;
  std::vector<DistributionWithTrendStorage *> correlation_bulk_shear_;      /// Converted to double
  std::vector<DistributionWithTrendStorage *> correlation_bulk_density_;    /// Converted to double
  std::vector<DistributionWithTrendStorage *> correlation_shear_density_;   /// Converted to double
  std::string                                 rock_name_;
};

//----------------------------------------------------------------------------------//
class ReussRockStorage : public DistributionsRockStorage {
public:
  ReussRockStorage(std::vector<std::string>                                  constituent_label,
                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                   std::string                                               rock_name);

  virtual ~ReussRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
  std::string                                               rock_name_;
};

//----------------------------------------------------------------------------------//
class VoigtRockStorage : public DistributionsRockStorage {
public:
  VoigtRockStorage(std::vector<std::string>                                  constituent_label,
                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                   std::string                                               rock_name);

  virtual ~VoigtRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
  std::string                                               rock_name_;
};

//----------------------------------------------------------------------------------//
class HillRockStorage : public DistributionsRockStorage {
public:
  HillRockStorage(std::vector<std::string>                                  constituent_label,
                  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction,
                  std::string                                               rock_name);

  virtual ~HillRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
  std::string                                               rock_name_;
};

//----------------------------------------------------------------------------------//
class DEMRockStorage : public DistributionsRockStorage {
public:
  DEMRockStorage(std::string                                               host_label,
                 std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                 std::vector<std::string>                                  inclusion_label,
                 std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                 std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio,
                 std::string                                               rock_name);

  virtual ~DEMRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string                                               host_label_;
  std::vector<DistributionWithTrendStorage *>               host_volume_fraction_;
  std::vector<std::string>                                  inclusion_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio_;
  std::string                                               rock_name_;
};

//----------------------------------------------------------------------------------//
class GassmannRockStorage : public DistributionsRockStorage {
public:
  GassmannRockStorage(std::string dry_rock,
                      std::string fluid,
                      std::string rock_name_);

  virtual ~GassmannRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string dry_rock_;
  std::string fluid_;
  std::string rock_name_;
};

//----------------------------------------------------------------------------------//
class BoundingRockStorage : public DistributionsRockStorage {
public:
  BoundingRockStorage(std::string                                 upper_rock,
                      std::string                                 lower_rock,
                      std::vector<DistributionWithTrendStorage *> porosity,
                      std::vector<DistributionWithTrendStorage *> bulk_weight,
                      std::vector<DistributionWithTrendStorage *> shear_weight,
                      double                                      correlation_weights,
                      std::string                                 rock_name);

  virtual ~BoundingRockStorage();

  virtual std::vector<DistributionsRock *> GenerateDistributionsRock(const int                                                   & n_vintages,
                                                                     const std::string                                           & path,
                                                                     const std::vector<std::string>                              & trend_cube_parameters,
                                                                     const std::vector<std::vector<double> >                     & trend_cube_sampling,
                                                                     const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
                                                                     const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
                                                                     const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
                                                                     const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
                                                                     const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
                                                                     const int                                                     output_other,
                                                                     std::string                                                 & errTxt) const;

private:
  std::string                                 upper_rock_;
  std::string                                 lower_rock_;
  std::vector<DistributionWithTrendStorage *> porosity_;
  std::vector<DistributionWithTrendStorage *> bulk_weight_;
  std::vector<DistributionWithTrendStorage *> shear_weight_;
  double                                      correlation_weights_;
  std::string                                 rock_name_;
};

#endif
