
#include "src/definitions.h"
#include "src/modelsettings.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionwithtrend.h"

void CheckVolumeConsistency(const std::vector<DistributionWithTrend *> & volume_fraction,
                            std::string                                & errTxt)
{
  int n_constituents = static_cast<int>(volume_fraction.size());

  if(n_constituents > 2) {
    for(int i = 0; i<n_constituents; i++) {
      if(volume_fraction[i]->GetIsDistribution() == true)
        errTxt += "The volume fractions can not be defined by a distribution when more than two constituents are used in a rock physics model\n";
    }
  }

  int n_missing = 0;

  for(int i=0; i<n_constituents; i++) {
    if(volume_fraction[i] == NULL)
      n_missing++;
  }

  if(n_missing == 0)
    errTxt += "One of the volume frations must be unspecified in the rock physics models where elements with corresponding volume frations are given\n";
  else if(n_missing > 1)
    errTxt += "All but one of the volume frations must be defined in the rock physics models where elements with corresponding volume frations are given\n";
}

void FindMixTypesForRock(std::vector<std::string>  constituent_label,
                         int n_constituents,
                         const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                         const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                         const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                         const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                         bool & mix_rock,
                         bool & mix_solid,
                         bool & mix_fluid,
                         std::vector<int> & constituent_type,
                         std::string & tmpErrTxt)
{

  for(int i=0; i<n_constituents; i++) {
    std::map<std::string, DistributionsRockStorage *>::const_iterator m = model_rock_storage.find(constituent_label[i]);
    if(m != model_rock_storage.end()) {
      constituent_type[i] = ModelSettings::ROCK;
      mix_rock = true;
    }
    else {
      std::map<std::string, DistributionsFluidStorage *>::const_iterator m = model_fluid_storage.find(constituent_label[i]);
      if(m != model_fluid_storage.end()) {
        constituent_type[i] = ModelSettings::FLUID;
        mix_fluid = true;
      }
      else {
        std::map<std::string, DistributionsSolidStorage *>::const_iterator m = model_solid_storage.find(constituent_label[i]);
        if(m != model_solid_storage.end()) {
          constituent_type[i] = ModelSettings::SOLID;
          mix_solid = true;
        }
        else {
          std::map<std::string, DistributionsDryRockStorage *>::const_iterator m = model_dry_rock_storage.find(constituent_label[i]);
          if(m != model_dry_rock_storage.end()) {
            constituent_type[i] = ModelSettings::DRY_ROCK;
            tmpErrTxt += "A dry-rock can not be used as constituent for mixing a rock\n";
          }
          else
            tmpErrTxt += "Failed to find label " + constituent_label[i] + "\n";
        }
      }
    }
  }

  if(mix_rock == true && mix_fluid == true && mix_solid == true)
    tmpErrTxt += "Fluids and solids can not be mixed with rocks in the Reuss model\n";
}


DistributionsRock *
ReadRock(const std::string                                          & target_rock,
         const std::string                                          & path,
         const std::vector<std::string>                             & trend_cube_parameters,
         const std::vector<std::vector<double> >                    & trend_cube_sampling,
         const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
         const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
         const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
         const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
         std::map<std::string, DistributionsRock *>                 & rock_distribution,
         std::map<std::string, DistributionsSolid *>                & solid_distribution,
         std::map<std::string, DistributionsDryRock *>              & dry_rock_distribution,
         std::map<std::string, DistributionsFluid *>                & fluid_distribution,
         std::string                                                & errTxt)
{
  DistributionsRock * final_rock = NULL;
  DistributionsRock * rock       = NULL;

  std::map<std::string, DistributionsRock *>::iterator m = rock_distribution.find(target_rock);
  if (m == rock_distribution.end()) { // label not found in rock_distribution map
    std::map<std::string, DistributionsRockStorage *>::const_iterator m_all = model_rock_storage.find(target_rock);
    if (m_all == model_rock_storage.end()) // fatal error
      errTxt += "Failed to find rock label " + target_rock + " requested in the rock physics model\n";
    else { //label found
      DistributionsRockStorage     * storage     = m_all->second;
      rock                                       = storage->GenerateDistributionsRock(path,
                                                                                      trend_cube_parameters,
                                                                                      trend_cube_sampling,
                                                                                      model_rock_storage,
                                                                                      model_solid_storage,
                                                                                      model_dry_rock_storage,
                                                                                      model_fluid_storage,
                                                                                      rock_distribution,
                                                                                      solid_distribution,
                                                                                      dry_rock_distribution,
                                                                                      fluid_distribution,
                                                                                      errTxt);
      rock_distribution[m_all->first]            = rock;
      final_rock                                 = rock;
    }
  }
  else // label found
    final_rock = m->second;

  return(final_rock);

}

DistributionsSolid *
ReadSolid(const std::string                                          & target_solid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
          std::map<std::string, DistributionsSolid *>                & solid_distribution,
          std::string                                                & errTxt)
{

  DistributionsSolid * final_solid = NULL;
  DistributionsSolid * solid       = NULL;

  std::map<std::string, DistributionsSolid *>::iterator m = solid_distribution.find(target_solid);

  if (m == solid_distribution.end()) { // label not found in solid_distribution map
    std::map<std::string, DistributionsSolidStorage *>::const_iterator m_all = model_solid_storage.find(target_solid);

    if (m_all == model_solid_storage.end()) // fatal error
      errTxt += "Failed to find solid label " + target_solid + "\n";

    else { //label found
      DistributionsSolidStorage  * storage = m_all->second;
      solid                                = storage->GenerateDistributionsSolid(path,
                                                                                 trend_cube_parameters,
                                                                                 trend_cube_sampling,
                                                                                 model_solid_storage,
                                                                                 solid_distribution,
                                                                                 errTxt);
      solid_distribution[m_all->first]     = solid;
      final_solid                          = solid;
    }
  }
  else // label found
    final_solid = m->second;

  return(final_solid);
}

DistributionsFluid *
ReadFluid(const std::string                                          & target_fluid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
          std::map<std::string, DistributionsFluid *>                & fluid_distribution,
          std::string                                                & errTxt)
{

  DistributionsFluid * final_fluid = NULL;
  DistributionsFluid * fluid       = NULL;

  std::map<std::string, DistributionsFluid *>::iterator m = fluid_distribution.find(target_fluid);

  if (m == fluid_distribution.end()) { // label not found in fluid_distribution map
    std::map<std::string, DistributionsFluidStorage *>::const_iterator m_all = model_fluid_storage.find(target_fluid);

    if (m_all == model_fluid_storage.end()) // fatal error
      errTxt += "Failed to find fluid label " + target_fluid + "\n";

    else { //label found
      DistributionsFluidStorage  * storage = m_all->second;
      fluid                                = storage->GenerateDistributionsFluid(path,
                                                                                 trend_cube_parameters,
                                                                                 trend_cube_sampling,
                                                                                 errTxt);
      fluid_distribution[m_all->first]     = fluid;
      final_fluid                          = fluid;
    }
  }
  else // label found
    final_fluid = m->second;

  return(final_fluid);
}
