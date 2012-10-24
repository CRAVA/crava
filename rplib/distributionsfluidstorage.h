#ifndef RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP

#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrendstorage.h"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"

class DistributionsFluidStorage {
public:
  DistributionsFluidStorage();

  virtual ~DistributionsFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & /*path*/,
                                                          const std::vector<std::string>          & /*trend_cube_parameters*/,
                                                          const std::vector<std::vector<double> > & /*trend_cube_sampling*/,
                                                          std::string                             & /*errTxt*/)                    const = 0;
};

//----------------------------------------------------------------------------------//

class TabulatedVelocityFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedVelocityFluidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                std::vector<DistributionWithTrendStorage *> density,
                                double                                      correlation_vp_density);

  virtual ~TabulatedVelocityFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> vp_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_vp_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedModulusFluidStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                               std::vector<DistributionWithTrendStorage *> density,
                               double                                      correlation_bulk_density);

  virtual ~TabulatedModulusFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> bulk_modulus_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_bulk_density_;
};

//----------------------------------------------------------------------------------//
class ReussFluidStorage : public DistributionsFluidStorage {
public:
  ReussFluidStorage(std::vector<std::string>                                  constituent_label,
                    std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~ReussFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtFluidStorage : public DistributionsFluidStorage {
public:
  VoigtFluidStorage(std::vector<std::string>                                  constituent_label,
                    std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~VoigtFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillFluidStorage : public DistributionsFluidStorage {
public:
  HillFluidStorage(std::vector<std::string>                                  constituent_label,
                   std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~HillFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class BatzleWangFluidStorage : public DistributionsFluidStorage {
public:
  BatzleWangFluidStorage(std::vector<DistributionWithTrendStorage *> pore_pressure,
                         std::vector<DistributionWithTrendStorage *> temperature,
                         std::vector<DistributionWithTrendStorage *> salinity);

  virtual ~BatzleWangFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> pore_pressure_;
  std::vector<DistributionWithTrendStorage *> temperature_;
  std::vector<DistributionWithTrendStorage *> salinity_;
};

#endif
