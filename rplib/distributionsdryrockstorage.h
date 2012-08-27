#ifndef RPLIB_DISTRIBUTIONS_DRY_ROCK_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_DRY_ROCK_STORAGE_HPP

#include "rplib/distributionsdryrock.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsDryRockStorage {
public:
  DistributionsDryRockStorage();

  virtual ~DistributionsDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & /*path*/,
                                                              const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                              const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                              std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//

class TabulatedVelocityDryRockStorage : public DistributionsDryRockStorage {
public:
  TabulatedVelocityDryRockStorage(DistributionWithTrendStorage * vp,
                                  DistributionWithTrendStorage * vs,
                                  DistributionWithTrendStorage * density,
                                  DistributionWithTrendStorage * correlation_vp_vs,
                                  DistributionWithTrendStorage * correlation_vp_density,
                                  DistributionWithTrendStorage * correlation_vs_density,
                                  DistributionWithTrendStorage * total_porosity,
                                  std::string                    moduli);

  virtual ~TabulatedVelocityDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * vs_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_vp_vs_;
  DistributionWithTrendStorage * correlation_vp_density_;
  DistributionWithTrendStorage * correlation_vs_density_;
  DistributionWithTrendStorage * total_porosity_;
  std::string                    mineral_moduli_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusDryRockStorage : public DistributionsDryRockStorage {
public:
  TabulatedModulusDryRockStorage(DistributionWithTrendStorage * bulk_modulus,
                                 DistributionWithTrendStorage * shear_modulus,
                                 DistributionWithTrendStorage * density,
                                 DistributionWithTrendStorage * correlation_bulk_shear,
                                 DistributionWithTrendStorage * correlation_bulk_density,
                                 DistributionWithTrendStorage * correlation_shear_density,
                                 DistributionWithTrendStorage * total_porosity,
                                 std::string                    moduli);

  virtual ~TabulatedModulusDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * bulk_modulus_;
  DistributionWithTrendStorage * shear_modulus_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_bulk_shear_;
  DistributionWithTrendStorage * correlation_bulk_density_;
  DistributionWithTrendStorage * correlation_shear_density_;
  DistributionWithTrendStorage * total_porosity_;
  std::string                    mineral_moduli_;
};

//----------------------------------------------------------------------------------//
class ReussDryRockStorage : public DistributionsDryRockStorage {
public:
  ReussDryRockStorage(std::vector<std::string>                    constituent_label,
                      std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                      DistributionWithTrendStorage *              total_porosity,
                      std::string                                 moduli);

  virtual ~ReussDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  std::vector<std::string>                      constituent_label_;
  std::vector<DistributionWithTrendStorage *>   constituent_volume_fraction_;
  DistributionWithTrendStorage                * total_porosity_;
  std::string                                   mineral_moduli_;
};

//----------------------------------------------------------------------------------//
class VoigtDryRockStorage : public DistributionsDryRockStorage {
public:
  VoigtDryRockStorage(std::vector<std::string>                    constituent_label,
                      std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                      DistributionWithTrendStorage *              total_porosity,
                      std::string                                 moduli);

  virtual ~VoigtDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  std::vector<std::string>                      constituent_label_;
  std::vector<DistributionWithTrendStorage *>   constituent_volume_fraction_;
  DistributionWithTrendStorage                * total_porosity_;
  std::string                                   mineral_moduli_;
};

//----------------------------------------------------------------------------------//
class HillDryRockStorage : public DistributionsDryRockStorage {
public:
  HillDryRockStorage(std::vector<std::string>                    constituent_label,
                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction,
                     DistributionWithTrendStorage *              total_porosity,
                     std::string                                 moduli);

  virtual ~HillDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  std::vector<std::string>                      constituent_label_;
  std::vector<DistributionWithTrendStorage *>   constituent_volume_fraction_;
  DistributionWithTrendStorage                * total_porosity_;
  std::string                                   mineral_moduli_;
};

//----------------------------------------------------------------------------------//

class DEMDryRockStorage : public DistributionsDryRockStorage {
public:
  DEMDryRockStorage(std::string                                 host_label,
                    DistributionWithTrendStorage *              host_volume_fraction,
                    DistributionWithTrendStorage *              host_aspect_ratio,
                    std::vector<std::string>                    inclusion_label,
                    std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                    std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio,
                    DistributionWithTrendStorage *              total_porosity,
                    std::string                                 moduli);

  virtual ~DEMDryRockStorage();

  virtual DistributionsDryRock * GenerateDistributionsDryRock(const std::string                       & path,
                                                              const std::vector<std::string>          & trend_cube_parameters,
                                                              const std::vector<std::vector<double> > & trend_cube_sampling,
                                                              std::string                             & errTxt) const;

private:
  std::string                                 host_label_;
  DistributionWithTrendStorage *              host_volume_fraction_;
  DistributionWithTrendStorage *              host_aspect_ratio_;
  std::vector<std::string>                    inclusion_label_;
  std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction_;
  std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio_;
  DistributionWithTrendStorage              * total_porosity_;
  std::string                                 mineral_moduli_;
};

#endif
