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

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const = 0;
protected:
  std::vector<DistributionsSolid *> CreateDistributionsSolidMix(const std::string                                               & path,
                                                                const std::vector<std::string>                                  & trend_cube_parameters,
                                                                const std::vector<std::vector<double> >                         & trend_cube_sampling,
                                                                const std::map<std::string, DistributionsSolidStorage *>        & model_solid_storage,
                                                                const std::vector<std::string>                                  & constituent_label,
                                                                const std::vector<std::vector<DistributionWithTrendStorage *> > & constituent_volume_fraction,
                                                                std::map<std::string, std::vector<DistributionsSolid *> >       & solid_distribution,
                                                                DEMTools::MixMethod                                               mix_method,
                                                                std::string                                                     & errTxt) const;
};

//----------------------------------------------------------------------------------//
class TabulatedVelocitySolidStorage : public DistributionsSolidStorage {
public:
  TabulatedVelocitySolidStorage(std::vector<DistributionWithTrendStorage *> vp,
                                std::vector<DistributionWithTrendStorage *> vs,
                                std::vector<DistributionWithTrendStorage *> density,
                                double                                      correlation_vp_vs,
                                double                                      correlation_vp_density,
                                double                                      correlation_vs_density);

  virtual ~TabulatedVelocitySolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const;

private:
  std::vector<DistributionWithTrendStorage *> vp_;
  std::vector<DistributionWithTrendStorage *> vs_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_vp_vs_;
  double                                      correlation_vp_density_;
  double                                      correlation_vs_density_;
};

//----------------------------------------------------------------------------------//
class TabulatedModulusSolidStorage : public DistributionsSolidStorage {
public:
  TabulatedModulusSolidStorage(std::vector<DistributionWithTrendStorage *> bulk_modulus,
                               std::vector<DistributionWithTrendStorage *> shear_modulus,
                               std::vector<DistributionWithTrendStorage *> density,
                               double                                      correlation_bulk_shear,
                               double                                      correlation_bulk_density,
                               double                                      correlation_shear_density);

  virtual ~TabulatedModulusSolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const ;

private:
  std::vector<DistributionWithTrendStorage *> bulk_modulus_;
  std::vector<DistributionWithTrendStorage *> shear_modulus_;
  std::vector<DistributionWithTrendStorage *> density_;
  double                                      correlation_bulk_shear_;
  double                                      correlation_bulk_density_;
  double                                      correlation_shear_density_;
};

//----------------------------------------------------------------------------------//
class ReussSolidStorage : public DistributionsSolidStorage {
public:
  ReussSolidStorage(std::vector<std::string>                                  constituent_label,
                    std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~ReussSolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class VoigtSolidStorage : public DistributionsSolidStorage {
public:
  VoigtSolidStorage(std::vector<std::string>                                  constituent_label,
                    std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~VoigtSolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const;

private:
  std::vector<std::string>                                 constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class HillSolidStorage : public DistributionsSolidStorage {
public:
  HillSolidStorage(std::vector<std::string>                   constituent_label,
    std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction);

  virtual ~HillSolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const;

private:
  std::vector<std::string>                                  constituent_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > constituent_volume_fraction_;
};

//----------------------------------------------------------------------------------//
class DEMSolidStorage : public DistributionsSolidStorage {
public:
  DEMSolidStorage(std::string                                               host_label,
                  std::vector<DistributionWithTrendStorage *>               host_volume_fraction,
                  std::vector<std::string>                                  inclusion_label,
                  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction,
                  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio);

  virtual ~DEMSolidStorage();

  virtual std::vector<DistributionsSolid *> GenerateDistributionsSolid(const std::string                                         & path,
                                                                       const std::vector<std::string>                            & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> >                   & trend_cube_sampling,
                                                                       const std::map<std::string, DistributionsSolidStorage *>  & model_solid_storage,
                                                                       std::map<std::string, std::vector<DistributionsSolid *> > & solid_distribution,
                                                                       std::string                                               & errTxt) const;

private:
  std::string                                               host_label_;
  std::vector<DistributionWithTrendStorage *>               host_volume_fraction_;
  std::vector<std::string>                                  inclusion_label_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_volume_fraction_;
  std::vector<std::vector<DistributionWithTrendStorage *> > inclusion_aspect_ratio_;
};

#endif
