#ifndef RPLIB_DISTRIBUTIONS_STORAGE_KIT_H
#define RPLIB_DISTRIBUTIONS_STORAGE_KIT_H

#include <string>
#include <vector>
#include <map>

class DistributionWithTrend;
class DistributionsRockStorage;
class DistributionsSolidStorage;
class DistributionsDryRockStorage;
class DistributionsFluidStorage;

class DistributionsRock;
class DistributionsSolid;
class DistributionsFluid;
class DistributionsDryRock;


void CheckVolumeConsistency(const std::vector<DistributionWithTrend *> & volume_fraction,
                            std::string                                & errTxt);

void FindMixTypesForRock(std::vector<std::string> constituent_label,
                         int n_constituents,
                         const std::map<std::string, DistributionsRockStorage *>    & model_rock_storage,
                         const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                         const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                         const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
                         bool & mix_rock,
                         bool & mix_solid,
                         bool & mix_fluid,
                         std::vector<int> & constituent_type,
                         std::string & tmpErrTxt);

void FindMixTypesForDryRock(std::vector<std::string> constituent_label,
                            int n_constituents,
                            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                            const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                            bool & mix_dryrock,
                            bool & mix_solid,
                            std::vector<int> & constituent_type,
                            std::string & tmpErrTxt);

std::vector<DistributionsRock *>
ReadRock(const int                                                   & n_vintages,
         const std::string                                           & target_rock,
         const std::string                                           & path,
         const std::vector<std::string>                              & trend_cube_parameters,
         const std::vector<std::vector<double> >                     & trend_cube_sampling,
         const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
         const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
         const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
         const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
         std::string                                                 & errTxt);

std::vector<DistributionsSolid *>
ReadSolid(const int                                                  & n_vintages,
          const std::string                                          & target_solid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
          std::string                                                & errTxt);

std::vector<DistributionsDryRock *>
ReadDryRock(const int                                                  & n_vintages,
            const std::string                                          & target_dryrock,
            const std::string                                          & path,
            const std::vector<std::string>                             & trend_cube_parameters,
            const std::vector<std::vector<double> >                    & trend_cube_sampling,
            const std::map<std::string, DistributionsDryRockStorage *> & model_dryrock_storage,
            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
            std::string                                                & errTxt);

std::vector<DistributionsFluid *>
ReadFluid(const int                                                  & n_vintages,
          const std::string                                          & target_fluid,
          const std::string                                          & path,
          const std::vector<std::string>                             & trend_cube_parameters,
          const std::vector<std::vector<double> >                    & trend_cube_sampling,
          const std::map<std::string, DistributionsFluidStorage *>   & model_fluid_storage,
          std::string                                                & errTxt);

void
FindSMinMax(const std::vector<std::vector<double> > & trend_cube_sampling,
            std::vector<double>                     & s_min,
            std::vector<double>                     & s_max);

void
CheckPositiveDefiniteCorrMatrix(double corr01,
                                double corr02,
                                double corr12,
                                std::string & errTxt);
void
FindSamplingMinMax(const std::vector<std::vector<double> > & trend_cube_sampling,
                   std::vector<double>                     & s_min,
                   std::vector<double>                     & s_max);

void CheckValuesInZeroOne(const std::vector<DistributionWithTrendStorage *> & test_objects,
                          const std::string                                 & type,
                          const std::string                                 & path,
                          const std::vector<std::string>                    & trend_cube_parameters,
                          const std::vector<std::vector<double> >           & trend_cube_sampling,
                          std::string                                       & errTxt);

#endif
