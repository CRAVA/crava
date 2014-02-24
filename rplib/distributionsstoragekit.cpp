
#include "src/definitions.h"
#include "src/modelsettings.h"
#include "nrlib/flens/nrlib_flens.hpp"

//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
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

//--------------------------------------------------------------//
void FindMixTypesForDryRock(std::vector<std::string> constituent_label,
                            int n_constituents,
                            const std::map<std::string, DistributionsSolidStorage *>   & model_solid_storage,
                            const std::map<std::string, DistributionsDryRockStorage *> & model_dry_rock_storage,
                            bool & mix_dryrock,
                            bool & mix_solid,
                            std::vector<int> & constituent_type,
                            std::string & tmpErrTxt)
{
  for(int i=0; i<n_constituents; i++) {
    std::map<std::string, DistributionsDryRockStorage *>::const_iterator m = model_dry_rock_storage.find(constituent_label[i]);
      if(m != model_dry_rock_storage.end()) {
        constituent_type[i] = ModelSettings::DRY_ROCK;
        mix_dryrock = true;
      }
      else {
        std::map<std::string, DistributionsSolidStorage *>::const_iterator m = model_solid_storage.find(constituent_label[i]);
        if(m != model_solid_storage.end()) {
          constituent_type[i] = ModelSettings::SOLID;
          mix_solid = true;
        }
        else
          tmpErrTxt += "Failed to find label " + constituent_label[i] + "as a dryrock or a solid\n";
      }
  }
}

//--------------------------------------------------------------//
std::vector<DistributionsRock *>
ReadRock(const int                                                   & n_vintages,
         const std::string                                           & target_rock,
         const std::string                                           & path,
         const std::vector<std::string>                              & trend_cube_parameters,
         const std::vector<std::vector<double> >                     & trend_cube_sampling,
         const std::vector<BlockedLogsForRockPhysics *>              & blockedLogs,
         const std::map<std::string, DistributionsRockStorage *>     & model_rock_storage,
         const std::map<std::string, DistributionsSolidStorage *>    & model_solid_storage,
         const std::map<std::string, DistributionsDryRockStorage *>  & model_dry_rock_storage,
         const std::map<std::string, DistributionsFluidStorage *>    & model_fluid_storage,
         const int                                                     output_other,
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
                                                                                    blockedLogs,
                                                                                    model_rock_storage,
                                                                                    model_solid_storage,
                                                                                    model_dry_rock_storage,
                                                                                    model_fluid_storage,
                                                                                    output_other,
                                                                                    errTxt);
  }

  return(rock);

}
//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
void
CheckPositiveDefiniteCorrMatrix(double corr01, double corr02, double corr12, std::string & errTxt)
{

  NRLib::Matrix corr_matrix(3, 3, 0);
  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr01;
  corr_matrix(1,0) = corr01;
  corr_matrix(0,2) = corr02;
  corr_matrix(2,0) = corr02;
  corr_matrix(1,2) = corr12;
  corr_matrix(2,1) = corr12;

  NRLib::Vector eigen_values(3, 0);
  NRLib::Matrix eigen_vectors(3, 3, 0);
  NRLib::ComputeEigenVectors(corr_matrix, eigen_values, eigen_vectors);

  bool pos_def = true;
  for( int i=0; i<3; i++) {
    if(eigen_values(i) < 0)
      pos_def = false;
  }

  if(pos_def == false)
    errTxt += "The correlations given in the tabulated rock physics model need to generate a positive definite matrix\n";
}
//--------------------------------------------------------------//
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
//--------------------------------------------------------------//
void CheckValuesInZeroOne(const std::vector<DistributionWithTrendStorage *> & test_objects,
                          const std::string                                 & type,
                          const std::string                                 & path,
                          const std::vector<std::string>                    & trend_cube_parameters,
                          const std::vector<std::vector<double> >           & trend_cube_sampling,
                          const std::vector<std::vector<float> >            & blocked_logs,
                          std::string                                       & errTxt)
{
  std::vector<std::vector<double> > dummy_s1;
  std::vector<std::vector<double> > dummy_s2;

  NRLib::Trend * mean_trend_dummy = NULL;

  for(size_t i=0; i<test_objects.size(); i++) {

    if(test_objects[i] != NULL) {
      if(typeid(*(test_objects[i])) == typeid(DeltaDistributionWithTrendStorage)) {
        const DeltaDistributionWithTrendStorage * delta = dynamic_cast<const DeltaDistributionWithTrendStorage *>(test_objects[i]);

        NRLib::TrendStorage * mean_storage = delta->CloneMean();
        NRLib::Trend        * mean         = mean_storage->GenerateTrend(path,
                                                                         trend_cube_parameters,
                                                                         trend_cube_sampling,
                                                                         blocked_logs,
                                                                         dummy_s1,
                                                                         dummy_s2,
                                                                         NRLib::TrendStorage::MEAN,
                                                                         mean_trend_dummy,
                                                                         errTxt);

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

void
FindDoubleValueFromDistributionWithTrend(const DistributionWithTrendStorage * dist_with_trend,
                                         std::string                          type,
                                         double                             & value,
                                         std::string                        & errTxt)
{

  if(typeid(*(dist_with_trend)) == typeid(DeltaDistributionWithTrendStorage)) {
    const DeltaDistributionWithTrendStorage * d1 = dynamic_cast<const DeltaDistributionWithTrendStorage *>(dist_with_trend);
    if(typeid((*d1->GetMean())) == typeid(NRLib::TrendConstantStorage)) {
      const NRLib::TrendConstantStorage * t1 = dynamic_cast<const NRLib::TrendConstantStorage *>(d1->GetMean());
      value = t1->GetMean();
    }
  }
  else {
    errTxt += "All "+type+" variables need to be double values, not trends or distributions\n";
  }
}
// ----------------------------------------------------------------------------------------- //
void PreprocessDataForSpearmanCorrelation(const std::vector<std::vector<double> > & s1,
                                          const std::vector<std::vector<double> > & s2,
                                          const std::vector<std::vector<float> >  & blocked_logs_x,
                                          const std::vector<std::vector<float> >  & blocked_logs_y,
                                          DistributionWithTrend                   * distribution_with_trend_x,
                                          DistributionWithTrend                   * distribution_with_trend_y,
                                          std::vector<double>                     & x,
                                          std::vector<double>                     & x_mean,
                                          std::vector<double>                     & x_variance,
                                          std::vector<double>                     & y,
                                          std::vector<double>                     & y_mean,
                                          std::vector<double>                     & y_variance,
                                          std::string                             & errTxt)
{
  x.clear();
  x_mean.clear();
  x_variance.clear();

  y.clear();
  y_mean.clear();
  y_variance.clear();

  size_t nx            = blocked_logs_x.size();
  size_t ny            = blocked_logs_y.size();

  // retriev observations from wells
  if (nx != ny || nx == 0) {
    errTxt += "Error: Incompatible dimensions in EstimateSpearmanCorrelation(...).\n";
  } else {
    for (size_t i = 0; i < nx; i++) {

      size_t nx_i = blocked_logs_x[i].size();
      size_t ny_i = blocked_logs_y[i].size();

      if (nx_i != ny_i) {
        errTxt += "Error: Incompatible dimensions in EstimateSpearmanCorrelation(...).\n";
      } else {
        for (size_t j = 0; j < nx_i; j++) {

          double s1_j         = s1[i][j];
          double s2_j         = s2[i][j];

          double x_log_j      = blocked_logs_x[i][j];
          double x_mean_j     = distribution_with_trend_x->GetMeanValue(s1_j, s2_j);
          double x_variance_j = distribution_with_trend_x->GetVarianceValue(s1_j, s2_j);

          double y_log_j      = blocked_logs_y[i][j];
          double y_mean_j     = distribution_with_trend_y->GetMeanValue(s1_j, s2_j);
          double y_variance_j = distribution_with_trend_y->GetVarianceValue(s1_j, s2_j);

          if (x_log_j != RMISSING && x_mean_j != RMISSING && x_variance_j != RMISSING && y_log_j != RMISSING && y_mean_j != RMISSING && y_variance_j != RMISSING) {
            x.push_back(std::exp(x_log_j));
            x_mean.push_back(x_mean_j);
            x_variance.push_back(x_variance_j);

            y.push_back(std::exp(y_log_j));
            y_mean.push_back(y_mean_j);
            y_variance.push_back(y_variance_j);
          }
        }
      }
    }
  }
}
// ----------------------------------------------------------------------------------------- //
double EstimateSpearmanCorrelation(const std::vector<double> & x,
                                   const std::vector<double> & x_mean,
                                   const std::vector<double> & x_variance,
                                   const std::vector<double> & y,
                                   const std::vector<double> & y_mean,
                                   const std::vector<double> & y_variance,
                                   std::string               & errTxt)

{
  size_t nx = x.size();
  size_t ny = y.size();

  std::vector<std::pair<float, size_t> > x_order_tmp;
  std::vector<std::pair<float, size_t> > y_order_tmp;

  bool removed_missing = false;

  // retriev observations from wells
  size_t l = 0;
  if (nx != ny || nx == 0) {
    errTxt += "Error: Incompatible dimensions in EstimateSpearmanCorrelation(...).\n";
    return(RMISSING);
  } else {
    for (size_t i = 0; i < nx; i++) {
      double x_i          = x[i];
      double x_mean_i     = x_mean[i];
      double x_variance_i = x_variance[i];

      double y_i          = y[i];
      double y_mean_i     = y_mean[i];
      double y_variance_i = y_variance[i];

      if (x_i == RMISSING || x_mean_i == RMISSING || x_variance_i == RMISSING || y_i == RMISSING || y_mean_i == RMISSING || y_variance_i == RMISSING) {
        removed_missing = true;
      } else {
        if (x_variance_i == 0.0) {
          x_variance_i = 1.0;
        }
        if (x_variance_i == 0.0) {
          y_variance_i = 1.0;
        }
        x_order_tmp.push_back(std::make_pair((x_i - x_mean_i)/std::pow(x_variance_i, 0.5), l));
        y_order_tmp.push_back(std::make_pair((y_i - y_mean_i)/std::pow(y_variance_i, 0.5), l));
        l++;
      }
    }
  }


  // rank observations needed to calcualte Spearman correlation
  std::sort(x_order_tmp.begin(), x_order_tmp.end());
  std::sort(y_order_tmp.begin(), y_order_tmp.end());

  std::vector<std::pair<size_t, size_t> > x_rank_tmp(l);
  std::vector<std::pair<size_t, size_t> > y_rank_tmp(l);

  for (size_t i = 0; i < l; i++) {
    x_rank_tmp[i].first  = x_order_tmp[i].second;
    x_rank_tmp[i].second = i;

    y_rank_tmp[i].first  = y_order_tmp[i].second;
    y_rank_tmp[i].second = i;
  }

  std::sort(x_rank_tmp.begin(), x_rank_tmp.end());
  std::sort(y_rank_tmp.begin(), y_rank_tmp.end());

  std::vector<size_t> x_rank(l);
  std::vector<size_t> y_rank(l);

  for (size_t i = 0; i < l; i++) {
    x_rank[i] = x_rank_tmp[i].second + 1;
    y_rank[i] = y_rank_tmp[i].second + 1;
  }

  // calcualte Spearman correlation
  double sum_x  = 0.0;
  double sum_y  = 0.0;
  double sum_xx = 0.0;
  double sum_yy = 0.0;
  double sum_xy = 0.0;

  for (size_t i = 0; i < l; i++) {
    double x_rank_i = static_cast<double>(x_rank[i]);
    double y_rank_i = static_cast<double>(y_rank[i]);

    sum_x           = sum_x  + x_rank_i;
    sum_y           = sum_y  + y_rank_i;
    sum_xx          = sum_xx + x_rank_i*x_rank_i;
    sum_yy          = sum_yy + y_rank_i*y_rank_i;
    sum_xy          = sum_xy + x_rank_i*y_rank_i;
  }

  double rho_numerator           =  l*sum_xy - sum_x*sum_y;
  double rho_denominator_squared = (l*sum_xx - sum_x*sum_x)*(l*sum_yy - sum_y*sum_y);

  return(rho_numerator/std::pow(rho_denominator_squared, 0.5));
}
