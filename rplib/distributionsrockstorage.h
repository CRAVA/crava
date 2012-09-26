#ifndef RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_ROCK_STORAGE_HPP

#include "rplib/distributionsrock.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsRockStorage {
public:
  DistributionsRockStorage();

  virtual ~DistributionsRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & /*path*/,
                                                        const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                        const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                        std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//
class TabulatedVelocityRockStorage : public DistributionsRockStorage {
public:
  TabulatedVelocityRockStorage(DistributionWithTrendStorage * vp,
                               DistributionWithTrendStorage * vs,
                               DistributionWithTrendStorage * density,
                               double                         correlation_vp_vs,
                               double                         correlation_vp_density,
                               double                         correlation_vs_density);

  virtual ~TabulatedVelocityRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * vs_;
  DistributionWithTrendStorage * density_;
  double                         correlation_vp_vs_;
  double                         correlation_vp_density_;
  double                         correlation_vs_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusRockStorage : public DistributionsRockStorage {
public:
  TabulatedModulusRockStorage(DistributionWithTrendStorage * bulk_modulus,
                              DistributionWithTrendStorage * shear_modulus,
                              DistributionWithTrendStorage * density,
                              double                         correlation_bulk_shear,
                              double                         correlation_bulk_density,
                              double                         correlation_shear_density);

  virtual ~TabulatedModulusRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * bulk_modulus_;
  DistributionWithTrendStorage * shear_modulus_;
  DistributionWithTrendStorage * density_;
  double                         correlation_bulk_shear_;
  double                         correlation_bulk_density_;
  double                         correlation_shear_density_;
};

//----------------------------------------------------------------------------------//
class ReussRockStorage : public DistributionsRockStorage {
public:
  ReussRockStorage(std::vector<std::string>                    constituent_label,
                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~ReussRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtRockStorage : public DistributionsRockStorage {
public:
  VoigtRockStorage(std::vector<std::string>                    constituent_label,
                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~VoigtRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillRockStorage : public DistributionsRockStorage {
public:
  HillRockStorage(std::vector<std::string>                    constituent_label,
                  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~HillRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class DEMRockStorage : public DistributionsRockStorage {
public:
  DEMRockStorage(std::string                                 host_label,
                 DistributionWithTrendStorage *              host_volume_fraction,
                 DistributionWithTrendStorage *              host_aspect_ratio,
                 std::vector<std::string>                    inclusion_label,
                 std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                 std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio);

  virtual ~DEMRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
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
};

//----------------------------------------------------------------------------------//
class GassmannRockStorage : public DistributionsRockStorage {
public:
  GassmannRockStorage(std::string dry_rock,
                      std::string fluid);

  virtual ~GassmannRockStorage();

  virtual DistributionsRock * GenerateDistributionsRock(const std::string                       & path,
                                                        const std::vector<std::string>          & trend_cube_parameters,
                                                        const std::vector<std::vector<double> > & trend_cube_sampling,
                                                        std::string                             & errTxt) const;

private:
  std::string dry_rock_;
  std::string fluid_;
};

#endif
