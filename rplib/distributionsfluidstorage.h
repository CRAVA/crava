#ifndef RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_FLUID_STORAGE_HPP

#include "rplib/demmodelling.h"

#include <map>

class DistributionWithTrendStorage;
class DistributionsFluid;

class DistributionsFluidStorage {
public:
  DistributionsFluidStorage();

  virtual ~DistributionsFluidStorage();

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & /*n_vintages*/,
                                                                       const std::string                                              & /*path*/,
                                                                       const std::vector<std::string>                                 & /*trend_cube_parameters*/,
                                                                       const std::vector<std::vector<double> >                        & /*trend_cube_sampling*/,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & /*model_fluid_storage*/,
                                                                       std::string                                                    & /*errTxt*/) const = 0;

protected:
  std::vector<DistributionsFluid *>         CreateDistributionsFluidMix(const int                                                       & n_vintages,
                                                                        const std::string                                               & path,
                                                                        const std::vector<std::string>                                  & trend_cube_parameters,
                                                                        const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                                        const std::map<std::string, DistributionsFluidStorage *>        & model_fluid_storage,
                                                                        const std::vector<std::string>                                  & constituent_label,
                                                                        const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                                        DEMTools::MixMethod                                               mix_method,
                                                                        std::string                                                     & errTxt) const;
};

//----------------------------------------------------------------------------------//

class TabulatedVelocityFluidStorage : public DistributionsFluidStorage {
public:
  TabulatedVelocityFluidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                std::vector<DistributionWithTrendStorage *> density,
                                double                                      correlation_vp_density);

  virtual ~TabulatedVelocityFluidStorage();

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

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

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

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

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

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

 virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

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

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

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

  virtual std::vector<DistributionsFluid *> GenerateDistributionsFluid(const int                                                      & n_vintages,
                                                                       const std::string                                              & path,
                                                                       const std::vector<std::string>                                 & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                        & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsFluidStorage *>       & model_fluid_storage,
                                                                       std::string                                                    & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> pore_pressure_;
  std::vector<DistributionWithTrendStorage *> temperature_;
  std::vector<DistributionWithTrendStorage *> salinity_;
};

#endif
