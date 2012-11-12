
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
    errTxt += "One of the volume fracions must be unspecified in the rock physics models where elements with corresponding volume frations are given\n";
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


std::vector<DistributionsRock *>
ReadRock(const std::string                                           & target_rock,
         const std::string                                           & path,
         const std::vector<std::string>                              & trend_cube_parameters,
         const std::vector<std::vector<double> >                     & trend_cube_sampling,
         const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
         const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
         const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
         const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
         std::string                                                 & errTxt)
{
  std::vector<DistributionsRock *> rock;

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
                                                                                    errTxt);
  }

  return(rock);

}

std::vector<DistributionsSolid *>
ReadSolid(const std::string                                          & target_solid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
          std::string                                                & errTxt)
{
  std::vector<DistributionsSolid *> solid;

  std::map<std::string, DistributionsSolidStorage *>::const_iterator m_all = model_solid_storage.find(target_solid);

  if (m_all == model_solid_storage.end()) // fatal error
    errTxt += "Failed to find solid label " + target_solid + "\n";

  else { //label found
    DistributionsSolidStorage  * storage = m_all->second;
    solid                                = storage->GenerateDistributionsSolid(path,
                                                                               trend_cube_parameters,
                                                                               trend_cube_sampling,
                                                                               model_solid_storage,
                                                                               errTxt);
  }

  return(solid);
}

std::vector<DistributionsFluid *>
ReadFluid(const std::string                                          & target_fluid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
          std::string                                                & errTxt)
{

  std::vector<DistributionsFluid *> fluid;

  std::map<std::string, DistributionsFluidStorage *>::const_iterator m_all = model_fluid_storage.find(target_fluid);

  if (m_all == model_fluid_storage.end()) // fatal error
    errTxt += "Failed to find fluid label " + target_fluid + "\n";

  else { //label found
    DistributionsFluidStorage  * storage = m_all->second;
    fluid                                = storage->GenerateDistributionsFluid(path,
                                                                               trend_cube_parameters,
                                                                               trend_cube_sampling,
                                                                               errTxt);
  }

  return(fluid);
}

void
CheckVintageConsistency(const std::vector<int> & vintage_number,
                        std::string            & errTxt)
{
  int n_vintages = static_cast<int>(vintage_number.size());

  if(n_vintages > 1) {
    if(vintage_number[0] < 1)
      errTxt += "The vintage numbers need to be larger than zero in the rock physics model\n";

    int compare = vintage_number[0];
    for(int i=1; i<n_vintages; i++) {
      if(vintage_number[i] <= compare) {
        errTxt += "The vintage numbers need to be given in ascending order\n";
        break;
      }
      else
        compare = vintage_number[i];
    }
  }

  if(errTxt != "")
    errTxt += "Remember that the first vintage given under <reservoir><variable> is given vinage number 1\n";
}
