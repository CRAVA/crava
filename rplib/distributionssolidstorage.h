#ifndef RPLIB_DISTRIBUTIONS_SOLID_STORAGE_HPP
#define RPLIB_DISTRIBUTIONS_SOLID_STORAGE_HPP

#include "rplib/demmodelling.h"

class DistributionsSolid;
class DistributionWithTrendStorage;

#include <map>

class DistributionsSolidStorage {
public:
  DistributionsSolidStorage();

  virtual ~DistributionsSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const = 0;
protected:
  DistributionsSolid         * CreateDistributionsSolidMix(const std::string                                        & path,
                                                           const std::vector<std::string>                           & trend_cube_parameters,
                                                           const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                           const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                           const std::vector<std::string>                           & constituent_label,
                                                           const std::vector<DistributionWithTrendStorage *>        & constituent_volume_fraction,
                                                           std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                           DEMTools::MixMethod                                        mix_method,
                                                           std::string                                              & errTxt) const;
};

//----------------------------------------------------------------------------------//
class TabulatedVelocitySolidStorage : public DistributionsSolidStorage {
public:
  TabulatedVelocitySolidStorage(DistributionWithTrendStorage * vp,
                                DistributionWithTrendStorage * vs,
                                DistributionWithTrendStorage * density,
                                double                         correlation_vp_vs,
                                double                         correlation_vp_density,
                                double                         correlation_vs_density);

  virtual ~TabulatedVelocitySolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const;

private:
  DistributionWithTrendStorage * vp_;
  DistributionWithTrendStorage * vs_;
  DistributionWithTrendStorage * density_;
  double                         correlation_vp_vs_;
  double                         correlation_vp_density_;
  double                         correlation_vs_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusSolidStorage : public DistributionsSolidStorage {
public:
  TabulatedModulusSolidStorage(DistributionWithTrendStorage * bulk_modulus,
                               DistributionWithTrendStorage * shear_modulus,
                               DistributionWithTrendStorage * density,
                               double                         correlation_bulk_shear,
                               double                         correlation_bulk_density,
                               double                         correlation_shear_density);

  virtual ~TabulatedModulusSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const ;

private:
  DistributionWithTrendStorage * bulk_modulus_;
  DistributionWithTrendStorage * shear_modulus_;
  DistributionWithTrendStorage * density_;
  double                         correlation_bulk_shear_;
  double                         correlation_bulk_density_;
  double                         correlation_shear_density_;
};

//----------------------------------------------------------------------------------//
class ReussSolidStorage : public DistributionsSolidStorage {
public:
  ReussSolidStorage(std::vector<std::string>                    constituent_label,
                    std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~ReussSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtSolidStorage : public DistributionsSolidStorage {
public:
  VoigtSolidStorage(std::vector<std::string>                    constituent_label,
                    std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~VoigtSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillSolidStorage : public DistributionsSolidStorage {
public:
  HillSolidStorage(std::vector<std::string>                    constituent_label,
                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction);

  virtual ~HillSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const;

private:
  std::vector<std::string>                    constituent_label_;
  std::vector<DistributionWithTrendStorage *> constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class DEMSolidStorage : public DistributionsSolidStorage {
public:
  DEMSolidStorage(std::string                                 host_label,
                  DistributionWithTrendStorage *              host_volume_fraction,
                  DistributionWithTrendStorage *              host_aspect_ratio,
                  std::vector<std::string>                    inclusion_label,
                  std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                  std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio);

  virtual ~DEMSolidStorage();

  virtual DistributionsSolid * GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                          std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                          std::string                                              & errTxt) const;

private:
  std::string                                 host_label_;
  DistributionWithTrendStorage *              host_volume_fraction_;
  DistributionWithTrendStorage *              host_aspect_ratio_;
  std::vector<std::string>                    inclusion_label_;
  std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction_;
  std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio_;
};

#endif
