
#include "src/definitions.h"

#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsstoragekit.h"
#include "rplib/distributionwithtrendstorage.h"

void CheckVolumeConsistency(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
                            std::string                                       & errTxt)
{
  int n_constituents = static_cast<int>(volume_fraction.size());

  /*if(n_constituents > 2) {
    for(int i = 0; i<n_constituents; i++) {
      if(volume_fraction[i]->GetIsDistribution() == true)
        errTxt += "The volume fractions can not be defined by a distribution when more than two constituents are used in a rock physics model\n";
    }
  }*/

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
         std::map<std::string, const DistributionWithTrend *>       & reservoir_variables,
         std::string                                                & errTxt)
{
  DistributionsRock * final_rock = NULL;
  DistributionsRock * rock       = NULL;

  std::map<std::string, DistributionsRock *>::iterator m = rock_distribution.find(target_rock);
  if (m == rock_distribution.end()) { // label not found in rock_distribution map
    std::map<std::string, DistributionsRockStorage *>::const_iterator m_all = model_rock_storage.find(target_rock);
    if (m_all == model_rock_storage.end()) { // fatal error
      errTxt += "Failed to find rock label " + target_rock + " requested in the rock physics model\n";
    }
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
                                                                                      reservoir_variables,
                                                                                      errTxt);
      rock_distribution[m_all->first]            = rock;
      final_rock                                 = rock;
    }
  }
  else // label found
    final_rock = m->second;

  return(final_rock);

}
