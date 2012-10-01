#ifndef RPLIB_DISTRIBUTIONS_STORAGE_KIT_H
#define RPLIB_DISTRIBUTIONS_STORAGE_KIT_H

#include <string>
#include <vector>

#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsrock.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionsfluidstorage.h"

void CheckVolumeConsistency(const std::vector<DistributionWithTrendStorage *> & volume_fraction,
                            std::string                                       & errTxt);

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
         std::string                                                & errTxt);
#endif
