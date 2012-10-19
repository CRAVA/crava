#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionssolidtabulatedvelocity.h"
#include "rplib/distributionssoliddem.h"
#include "rplib/distributionssolidmix.h"
#include "rplib/distributionwithtrendstorage.h"

#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionssolidtabulatedmodulus.h"


DistributionsSolidStorage::DistributionsSolidStorage()
{
}

DistributionsSolidStorage::~DistributionsSolidStorage()
{
}

DistributionsSolid         *
DistributionsSolidStorage::CreateDistributionsSolidMix(const std::string                                        & path,
                                                       const std::vector<std::string>                           & trend_cube_parameters,
                                                       const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                       const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                                       const std::vector<std::string>                           & constituent_label,
                                                       const std::vector<DistributionWithTrendStorage *>        & constituent_volume_fraction,
                                                       std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                                       DEMTools::MixMethod                                        mix_method,
                                                       std::string                                              & errTxt) const
{

  DistributionsSolid * solid = NULL;

  std::vector<DistributionsSolid*> final_distr_solid;
  std::vector<DistributionWithTrend*> final_distr_vol_frac;
  size_t s;
  for (s = 0; s != constituent_label.size(); ++s) {
    DistributionsSolid * constit_solid = NULL;
    constit_solid = ReadSolid(constituent_label[s], path, trend_cube_parameters, trend_cube_sampling, model_solid_storage, solid_distribution, errTxt);
    final_distr_solid.push_back(constit_solid);

    if (errTxt == "" && constituent_volume_fraction[s]) {
      const DistributionWithTrend * vol_frac = constituent_volume_fraction[s]->GenerateDistributionWithTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
      final_distr_vol_frac.push_back(const_cast<DistributionWithTrend *>(vol_frac));
    }
    else
      final_distr_vol_frac.push_back(NULL);
  }

  CheckVolumeConsistency(final_distr_vol_frac, errTxt);

  if (errTxt == "") {
    solid = new DistributionsSolidMix(final_distr_solid, final_distr_vol_frac, mix_method);
  }

  return(solid);


}

//----------------------------------------------------------------------------------//
TabulatedVelocitySolidStorage::TabulatedVelocitySolidStorage(DistributionWithTrendStorage * vp,
                                                             DistributionWithTrendStorage * vs,
                                                             DistributionWithTrendStorage * density,
                                                             double                         correlation_vp_vs,
                                                             double                         correlation_vp_density,
                                                             double                         correlation_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density)
{
}

TabulatedVelocitySolidStorage::~TabulatedVelocitySolidStorage()
{
  if(vp_->GetIsShared() == false)
    delete vp_;
  if(vs_->GetIsShared() == false)
    delete vs_;
  if(density_->GetIsShared() == false)
    delete density_;
}

DistributionsSolid *
TabulatedVelocitySolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                                          const std::vector<std::string>                           & trend_cube_parameters,
                                                          const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                          const std::map<std::string, DistributionsSolidStorage *> & /*model_solid_storage*/,
                                                          std::map<std::string, DistributionsSolid *>              & /*solid_distribution*/,
                                                          std::string                                              & errTxt) const
{
  const DistributionWithTrend * vp_dist_with_trend              = vp_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * vs_dist_with_trend              = vs_                    ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend         = density_               ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsSolid * solid = new DistributionsSolidTabulatedVelocity(vp_dist_with_trend,
                                                                       vs_dist_with_trend,
                                                                       density_dist_with_trend,
                                                                       correlation_vp_vs_,
                                                                       correlation_vp_density_,
                                                                       correlation_vs_density_);

  return(solid);
}

//----------------------------------------------------------------------------------//
TabulatedModulusSolidStorage::TabulatedModulusSolidStorage(DistributionWithTrendStorage * bulk_modulus,
                                                           DistributionWithTrendStorage * shear_modulus,
                                                           DistributionWithTrendStorage * density,
                                                           double                         correlation_bulk_shear,
                                                           double                         correlation_bulk_density,
                                                           double                         correlation_shear_density)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  correlation_bulk_shear_(correlation_bulk_shear),
  correlation_bulk_density_(correlation_bulk_density),
  correlation_shear_density_(correlation_shear_density)
{
}

TabulatedModulusSolidStorage::~TabulatedModulusSolidStorage()
{
  if(bulk_modulus_->GetIsShared() == false)
    delete bulk_modulus_;
  if(shear_modulus_->GetIsShared() == false)
    delete shear_modulus_;
  if(density_->GetIsShared()== false)
    delete density_;
}

DistributionsSolid *
TabulatedModulusSolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                                         const std::vector<std::string>                           & trend_cube_parameters,
                                                         const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                                         const std::map<std::string, DistributionsSolidStorage *> & /*model_solid_storage*/,
                                                         std::map<std::string, DistributionsSolid *>              & /*solid_distribution*/,
                                                         std::string                                              & errTxt) const
{
  const DistributionWithTrend * bulk_dist_with_trend               = bulk_modulus_              ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * shear_dist_with_trend              = shear_modulus_             ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);
  const DistributionWithTrend * density_dist_with_trend            = density_                   ->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  DistributionsSolid * solid = new DistributionsSolidTabulatedModulus(bulk_dist_with_trend,
                                                                      shear_dist_with_trend,
                                                                      density_dist_with_trend,
                                                                      correlation_bulk_shear_,
                                                                      correlation_bulk_density_,
                                                                      correlation_shear_density_);

  return(solid);
}

//----------------------------------------------------------------------------------//
ReussSolidStorage::ReussSolidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

ReussSolidStorage::~ReussSolidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i] && constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
ReussSolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                              const std::vector<std::string>                           & trend_cube_parameters,
                                              const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                              const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                              std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                              std::string                                              & errTxt) const
{
  DistributionsSolid * solid = CreateDistributionsSolidMix(path,
                                                           trend_cube_parameters,
                                                           trend_cube_sampling,
                                                           model_solid_storage,
                                                           constituent_label_,
                                                           constituent_volume_fraction_,
                                                           solid_distribution,
                                                           DEMTools::Reuss,
                                                           errTxt);
  return(solid);
}

//----------------------------------------------------------------------------------//
VoigtSolidStorage::VoigtSolidStorage(std::vector<std::string>                    constituent_label,
                                     std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

VoigtSolidStorage::~VoigtSolidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i] && constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
VoigtSolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                              const std::vector<std::string>                           & trend_cube_parameters,
                                              const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                              const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                              std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                              std::string                                              & errTxt) const
{
  DistributionsSolid * solid = CreateDistributionsSolidMix(path,
                                                           trend_cube_parameters,
                                                           trend_cube_sampling,
                                                           model_solid_storage,
                                                           constituent_label_,
                                                           constituent_volume_fraction_,
                                                           solid_distribution,
                                                           DEMTools::Voigt,
                                                           errTxt);
  return(solid);
}

//----------------------------------------------------------------------------------//
HillSolidStorage::HillSolidStorage(std::vector<std::string>                    constituent_label,
                                   std::vector<DistributionWithTrendStorage *> constituent_volume_fraction)
: constituent_label_(constituent_label),
  constituent_volume_fraction_(constituent_volume_fraction)
{
}

HillSolidStorage::~HillSolidStorage()
{
  for(int i=0; i<static_cast<int>(constituent_volume_fraction_.size()); i++) {
    if(constituent_volume_fraction_[i] && constituent_volume_fraction_[i]->GetIsShared() == false)
      delete constituent_volume_fraction_[i];
  }
}

DistributionsSolid *
HillSolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                             const std::vector<std::string>                           & trend_cube_parameters,
                                             const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                             const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                             std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                             std::string                                              & errTxt) const
{
  DistributionsSolid * solid = CreateDistributionsSolidMix(path,
                                                           trend_cube_parameters,
                                                           trend_cube_sampling,
                                                           model_solid_storage,
                                                           constituent_label_,
                                                           constituent_volume_fraction_,
                                                           solid_distribution,
                                                           DEMTools::Hill,
                                                           errTxt);
  return(solid);
}

//----------------------------------------------------------------------------------//

DEMSolidStorage::DEMSolidStorage(std::string                                 host_label,
                                 DistributionWithTrendStorage *              host_volume_fraction,
                                 std::vector<std::string>                    inclusion_label,
                                 std::vector<DistributionWithTrendStorage *> inclusion_volume_fraction,
                                 std::vector<DistributionWithTrendStorage *> inclusion_aspect_ratio)
: host_label_(host_label),
  host_volume_fraction_(host_volume_fraction),
  inclusion_label_(inclusion_label),
  inclusion_volume_fraction_(inclusion_volume_fraction),
  inclusion_aspect_ratio_(inclusion_aspect_ratio)
{
}

DEMSolidStorage::~DEMSolidStorage()
{
  if(host_volume_fraction_->GetIsShared() == false)
    delete host_volume_fraction_;

  for(int i=0; i<static_cast<int>(inclusion_volume_fraction_.size()); i++) {
    if(inclusion_volume_fraction_[i] && inclusion_volume_fraction_[i]->GetIsShared() == false)
      delete inclusion_volume_fraction_[i];
  }

  for(int i=0; i<static_cast<int>(inclusion_aspect_ratio_.size()); i++) {
    if(inclusion_aspect_ratio_[i] && inclusion_aspect_ratio_[i]->GetIsShared() == false)
      delete inclusion_aspect_ratio_[i];
  }
}

DistributionsSolid *
DEMSolidStorage::GenerateDistributionsSolid(const std::string                                        & path,
                                            const std::vector<std::string>                           & trend_cube_parameters,
                                            const std::vector<std::vector<double> >                  & trend_cube_sampling,
                                            const std::map<std::string, DistributionsSolidStorage *> & model_solid_storage,
                                            std::map<std::string, DistributionsSolid *>              & solid_distribution,
                                            std::string                                              & errTxt) const
{
  // Remember: Host info is included first in inclusion vectors
  int n_inclusions = static_cast<int>(inclusion_volume_fraction_.size());

  std::vector<DistributionWithTrendStorage *> volume_fractions(n_inclusions + 1);
  volume_fractions[0] = host_volume_fraction_;

  for(int i=1; i<n_inclusions+1; i++)
    volume_fractions[i] = inclusion_volume_fraction_[i-1];

  DistributionsSolid * solid = NULL;

  std::vector< DistributionWithTrend *> inclusion_volume_fraction_distr(inclusion_volume_fraction_.size()+1, NULL);
  std::vector< DistributionWithTrend *> inclusion_aspect_ratio_distr(inclusion_aspect_ratio_.size(), NULL);

  inclusion_volume_fraction_distr[0]  = host_volume_fraction_->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  for (size_t i = 1; i < inclusion_volume_fraction_.size(); ++i)
    inclusion_volume_fraction_distr[i] = inclusion_volume_fraction_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  for (size_t i = 0; i < inclusion_aspect_ratio_.size(); ++i)
    inclusion_aspect_ratio_distr[i] = inclusion_aspect_ratio_[i]->GenerateDistributionWithTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

  //Read host label
  DistributionsSolid  * final_distr_solid = NULL;
  final_distr_solid = ReadSolid(host_label_,path, trend_cube_parameters, trend_cube_sampling, model_solid_storage, solid_distribution, errTxt);

  //Read inclusion label
  std::vector< DistributionsSolid* > final_distr_solid_inc;
  size_t s;
  for (s = 0; s != inclusion_label_.size(); ++s) {
    DistributionsSolid * incl_solid;
    incl_solid = ReadSolid(inclusion_label_[s],path, trend_cube_parameters, trend_cube_sampling, model_solid_storage, solid_distribution, errTxt);
    final_distr_solid_inc.push_back(incl_solid);
  }

  //CheckVolumeConsistency(distr_porosity, errTxt); //Fix when questions are solved

  if (errTxt == "") {
    solid = new DistributionsSolidDEM(final_distr_solid,
                                      final_distr_solid_inc,
                                      inclusion_aspect_ratio_distr,
                                      inclusion_volume_fraction_distr);
  }

  return(solid);
}
