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
  TabulatedVelocityFluidStorage(DistributionWithTrendStorage * vp,
                                DistributionWithTrendStorage * density,
                                DistributionWithTrendStorage * correlation_vp_density);

  virtual ~TabulatedVelocityFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_vp_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedModulusFluidStorage(DistributionWithTrendStorage * bulk_modulus,
                               DistributionWithTrendStorage * density,
                               DistributionWithTrendStorage * correlation_bulk_density);

  virtual ~TabulatedModulusFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * bulk_modulus_;
  DistributionWithTrendStorage * density_;
  DistributionWithTrendStorage * correlation_bulk_density_;
};

//----------------------------------------------------------------------------------//
class ReussFluidStorage : public DistributionsFluidStorage {
public:
  ReussFluidStorage(std::vector<std::string>                    constituent_label,
                    std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~ReussFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtFluidStorage : public DistributionsFluidStorage {
public:
  VoigtFluidStorage(std::vector<std::string>                    constituent_label,
                    std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~VoigtFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillFluidStorage : public DistributionsFluidStorage {
public:
  HillFluidStorage(std::vector<std::string>                    constituent_label,
                    std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~HillFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class BatzleWangFluidStorage : public DistributionsFluidStorage {
public:
  BatzleWangFluidStorage(DistributionWithTrendStorage * pore_pressure,
                         DistributionWithTrendStorage * temperature,
                         DistributionWithTrendStorage * salinity);

  virtual ~BatzleWangFluidStorage();

  virtual DistributionsFluid * GenerateDistributionsFluid(const std::string                       & path,
                                                          const std::vector<std::string>          & trend_cube_parameters,
                                                          const std::vector<std::vector<double> > & trend_cube_sampling,
                                                          std::string                             & errTxt) const;

private:
  DistributionWithTrendStorage * pore_pressure_;
  DistributionWithTrendStorage * temperature_;
  DistributionWithTrendStorage * salinity_;
};

#endif
