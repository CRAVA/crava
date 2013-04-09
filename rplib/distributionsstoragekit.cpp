
#include "src/definitions.h"
#include "src/modelsettings.h"
#include "nrlib/flens/nrlib_flens.hpp"


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
ReadRock(const int                                                   & n_vintages,
         const std::string                                           & target_rock,
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
    rock                                       = storage->GenerateDistributionsRock(n_vintages,
                                                                                    path,
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
ReadSolid(const int                                                  & n_vintages,
          const std::string                                          & target_solid,
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
    solid                                = storage->GenerateDistributionsSolid(n_vintages,
                                                                               path,
                                                                               trend_cube_parameters,
                                                                               trend_cube_sampling,
                                                                               model_solid_storage,
                                                                               errTxt);
  }

  return(solid);
}

std::vector<DistributionsDryRock *>
ReadDryRock(const int                                                  & n_vintages,
            const std::string                                          & target_dryrock,
            const std::string                                          & path,
            const std::vector<std::string>                             & trend_cube_parameters,
            const std::vector<std::vector<double> >                    & trend_cube_sampling,
            const std::map<std::string, DistributionsDryRockStorage *> & model_dryrock_storage,
            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
            std::string                                                & errTxt)
{
  std::vector<DistributionsDryRock *> dryrock;

  std::map<std::string, DistributionsDryRockStorage *>::const_iterator m_all = model_dryrock_storage.find(target_dryrock);

  if (m_all == model_dryrock_storage.end()) // fatal error
    errTxt += "Failed to find dryrock label " + target_dryrock + "\n";

  else { //label found
    DistributionsDryRockStorage  * storage = m_all->second;
    dryrock                                = storage->GenerateDistributionsDryRock(n_vintages,
                                                                                   path,
                                                                                   trend_cube_parameters,
                                                                                   trend_cube_sampling,
                                                                                   model_dryrock_storage,
                                                                                   model_solid_storage,
                                                                                   errTxt);
  }

  return(dryrock);
}

std::vector<DistributionsFluid *>
ReadFluid(const int                                                  & n_vintages,
          const std::string                                          & target_fluid,
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
    fluid                                = storage->GenerateDistributionsFluid(n_vintages,
                                                                               path,
                                                                               trend_cube_parameters,
                                                                               trend_cube_sampling,
                                                                               model_fluid_storage,
                                                                               errTxt);
  }

  return(fluid);
}

void
FindSMinMax(const std::vector<std::vector<double> > & trend_cube_sampling,
            std::vector<double>                     & s_min,
            std::vector<double>                     & s_max)
{
  s_min.resize(2, 0.0);
  s_max.resize(2, 0.0);

  int n_trend_cubes = static_cast<int>(trend_cube_sampling.size());

  for(int i=0; i<n_trend_cubes; i++) {
    s_min[i] = 0;
    s_max[i] = trend_cube_sampling[i][trend_cube_sampling[0].size()-1] - trend_cube_sampling[i][0];
  }
}

void
CheckPositiveDefiniteCorrMatrix(double corr01, double corr02, double corr12, std::string & errTxt)
{

  NRLib::Matrix corr_matrix(3,3);
  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr01;
  corr_matrix(1,0) = corr01;
  corr_matrix(0,2) = corr02;
  corr_matrix(2,0) = corr02;
  corr_matrix(1,2) = corr12;
  corr_matrix(2,1) = corr12;

  NRLib::Vector eigen_values(3);
  NRLib::Matrix eigen_vectors(3,3);
  NRLib::ComputeEigenVectors(corr_matrix, eigen_values, eigen_vectors);

  bool pos_def = true;
  for( int i=0; i<3; i++) {
    if(eigen_values(i) < 0)
      pos_def = false;
  }

  if(pos_def == false)
    errTxt += "The correlations given in the tabulated rock physics model need to generate a positive definite matrix\n";
}

void FindSamplingMinMax(const std::vector<std::vector<double> > & trend_cube_sampling,
                        std::vector<double>                     & s_min,
                        std::vector<double>                     & s_max)
{
  s_min.resize(2, 0.0);
  s_max.resize(2, 0.0);

  int n_trend_cubes = static_cast<int>(trend_cube_sampling.size());

  for(int i=0; i<n_trend_cubes; i++) {
    s_min[i] = trend_cube_sampling[i][0];
    s_max[i] = trend_cube_sampling[i][trend_cube_sampling[0].size()-1] - trend_cube_sampling[i][0];
  }
}

void CheckValuesInZeroOne(const std::vector<DistributionWithTrendStorage *> & test_objects,
                          const std::string                                 & type,
                          const std::string                                 & path,
                          const std::vector<std::string>                    & trend_cube_parameters,
                          const std::vector<std::vector<double> >           & trend_cube_sampling,
                          std::string                                       & errTxt)
{
  for(size_t i=0; i<test_objects.size(); i++) {

    if(test_objects[i] != NULL) {
      if(typeid(*(test_objects[i])) == typeid(DeltaDistributionWithTrendStorage)) {
        const DeltaDistributionWithTrendStorage * delta = dynamic_cast<const DeltaDistributionWithTrendStorage *>(test_objects[i]);

        NRLib::TrendStorage * mean_storage = delta->CloneMean();
        NRLib::Trend        * mean         = mean_storage->GenerateTrend(path, trend_cube_parameters, trend_cube_sampling, errTxt);

        double min = mean->GetMinValue();
        double max = mean->GetMaxValue();

        if(min < 0 || min > 1 || max > 1 || min > max) {
          errTxt += "The "+type+" must be in [0,1]\n";
          break;
        }

        delete mean_storage;
        delete mean;
      }
      else if(typeid(*(test_objects[i])) == typeid(BetaDistributionWithTrendStorage)) {
        const BetaDistributionWithTrendStorage * beta = dynamic_cast<const BetaDistributionWithTrendStorage *>(test_objects[i]);

        double lower_limit = beta->GetLowerLimit();
        double upper_limit = beta->GetUpperLimit();

        if(lower_limit < 0 || lower_limit > 1 || upper_limit > 1 || lower_limit > upper_limit)
          errTxt += "The limits in the Beta distribution must be in [0,1] for "+type+"\n";
      }
      else if(typeid(*(test_objects[i])) == typeid(BetaEndMassDistributionWithTrendStorage)) {
        const BetaEndMassDistributionWithTrendStorage * beta = dynamic_cast<const BetaEndMassDistributionWithTrendStorage *>(test_objects[i]);

        double lower_limit = beta->GetLowerLimit();
        double upper_limit = beta->GetUpperLimit();

        if(lower_limit < 0 || lower_limit > 1 || upper_limit > 1 || lower_limit > upper_limit)
          errTxt += "The limits in the Beta-end-mass distribution must be in [0,1] for "+type+"\n";
      }
      else
        errTxt += "The "+type+" must be in [0,1]. It should be given by a value or the Beta distribution with limits [0,1]\n";
    }
  }
}
