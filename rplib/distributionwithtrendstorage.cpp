#include "src/io.h"

#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"
#include "rplib/normaldistributionwithtrend.h"
#include "rplib/betadistributionwithtrend.h"
#include "rplib/betaendmassdistributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trendkit.hpp"
#include "nrlib/iotools/stringtools.hpp"


DistributionWithTrendStorage::DistributionWithTrendStorage()
{
}

DistributionWithTrendStorage::~DistributionWithTrendStorage()
{
}


DistributionWithTrend *
DistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                            const std::vector<std::string>          & trend_cube_parameters,
                                                            const std::vector<std::vector<double> > & trend_cube_sampling,
                                                            std::string                             & errTxt)
{
  const std::vector<std::vector<double> >   dummy_blocked_logs; // Use dummy when not estimatin.
  const std::vector<std::vector<double> >   dummy_s1;
  const std::vector<std::vector<double> >   dummy_s2;
  const int                                 dummy_output_other = -999;
  DistributionWithTrend * result = GenerateDistributionWithTrend(path,
                                                                 trend_cube_parameters,
                                                                 trend_cube_sampling,
                                                                 dummy_blocked_logs,
                                                                 dummy_s1,
                                                                 dummy_s2,
                                                                 dummy_output_other,
                                                                 "dummy",
                                                                 errTxt);
  return(result);
}


//--------------------------------------------------------------//

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage()
: is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                    = NULL;
  distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(double mean,
                                                                     bool   is_shared,
                                                                     bool   estimate)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
  //Use this constructor when mean_ is a double

  mean_                    = new NRLib::TrendConstantStorage(mean, estimate);
  distribution_with_trend_ = NULL;
  estimate_                = estimate;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                     bool                        is_shared)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                    = mean;
  distribution_with_trend_ = NULL;

  if(mean_->GetEstimate() == true)
    estimate_ = true;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(const DeltaDistributionWithTrendStorage & dist)
: is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_),
  estimate_(dist.estimate_)
{
  mean_ = dist.mean_->Clone();

 if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::~DeltaDistributionWithTrendStorage()
{
  delete mean_;

  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
DeltaDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                 const std::vector<std::string>          & trend_cube_parameters,
                                                                 const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                 const std::vector<std::vector<double> > & blocked_logs,
                                                                 const std::vector<std::vector<double> > & s1,
                                                                 const std::vector<std::vector<double> > & s2,
                                                                 const int                                 output_other,
                                                                 std::string                               filename_prefix,
                                                                 std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {  //Make sure shared variables are only generated one time
    NRLib::Trend * mean_trend_dummy = NULL;

    NRLib::Trend * mean_trend = mean_->GenerateTrend(path,
                                                     trend_cube_parameters,
                                                     trend_cube_sampling,
                                                     blocked_logs,
                                                     s1,
                                                     s2,
                                                     NRLib::TrendStorage::MEAN,
                                                     mean_trend_dummy,
                                                     errTxt);

    // write mean and variance surface to file
    if ((output_other & IO::ROCK_PHYSICS_TRENDS) > 0) {
      if (mean_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_mean", trend_cube_sampling, mean_trend);
    }


    int share_level = DistributionWithTrend::None;
    if(is_shared_ == true)
      share_level = DistributionWithTrend::SingleSample;

    if(errTxt == "")
      distribution_with_trend_= new DeltaDistributionWithTrend(mean_trend, share_level);

    delete mean_trend;
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
DeltaDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

//--------------------------------------------------------------//

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage()
: is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                    = NULL;
  variance_                = NULL;
  distribution_with_trend_ = NULL;
}

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                       const NRLib::TrendStorage * variance,
                                                                       bool                        is_shared)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                     = mean->Clone();
  variance_                 = variance->Clone();
  distribution_with_trend_  = NULL;

  if(mean_->GetEstimate() == true || variance_->GetEstimate() == true)
    estimate_ = true;
}

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage(const NormalDistributionWithTrendStorage & dist)
: is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_),
  estimate_(dist.estimate_)
{
  mean_     = dist.mean_    ->Clone();
  variance_ = dist.variance_->Clone();

  if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

NormalDistributionWithTrendStorage::~NormalDistributionWithTrendStorage()
{
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
NormalDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                  const std::vector<std::vector<double> > & blocked_logs,
                                                                  const std::vector<std::vector<double> > & s1,
                                                                  const std::vector<std::vector<double> > & s2,
                                                                  const int                                 output_other,
                                                                  std::string                               filename_prefix,
                                                                  std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {     //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend_dummy = NULL;
    NRLib::Trend * mean_trend       = mean_->GenerateTrend(path,
                                                           trend_cube_parameters,
                                                           trend_cube_sampling,
                                                           blocked_logs,
                                                           s1,
                                                           s2,
                                                           NRLib::TrendStorage::MEAN,
                                                           mean_trend_dummy,
                                                           errTxt);

    NRLib::Trend * variance_trend   = variance_->GenerateTrend(path,
                                                               trend_cube_parameters,
                                                               trend_cube_sampling,
                                                               blocked_logs,
                                                               s1,
                                                               s2,
                                                               NRLib::TrendStorage::VAR,
                                                               mean_trend,
                                                               errTxt);

    // write mean and variance surface to file
    if ((output_other & IO::ROCK_PHYSICS_TRENDS) > 0) {
      if (mean_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_mean", trend_cube_sampling, mean_trend);
      if (variance_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_var",  trend_cube_sampling, variance_trend);
    }

    int share_level = DistributionWithTrend::None;
    if(is_shared_ == true)
      share_level = DistributionWithTrend::SingleSample;

    if(errTxt == "")
      distribution_with_trend_= new NormalDistributionWithTrend(mean_trend, variance_trend, share_level);

    delete mean_trend;
    delete variance_trend;
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
NormalDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}
//--------------------------------------------------------------//

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage()
: lower_limit_(0),
  upper_limit_(0),
  is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                    = NULL;
  variance_                = NULL;
  distribution_with_trend_ = NULL;
}

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                   const NRLib::TrendStorage * variance,
                                                                   const double              & lower_limit,
                                                                   const double              & upper_limit,
                                                                   bool                        is_shared)
: lower_limit_(lower_limit),
  upper_limit_(upper_limit),
  is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                     = mean->Clone();
  variance_                 = variance->Clone();
  distribution_with_trend_  = NULL;

  if(mean_->GetEstimate() == true || variance_->GetEstimate() == true)
    estimate_ = true;
}

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage(const BetaDistributionWithTrendStorage & dist)
: lower_limit_(dist.upper_limit_),
  upper_limit_(dist.lower_limit_),
  is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_),
  estimate_(dist.estimate_)
{
  mean_     = dist.mean_    ->Clone();
  variance_ = dist.variance_->Clone();

 if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

BetaDistributionWithTrendStorage::~BetaDistributionWithTrendStorage()
{
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
BetaDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                const std::vector<std::string>          & trend_cube_parameters,
                                                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                const std::vector<std::vector<double> > & blocked_logs,
                                                                const std::vector<std::vector<double> > & s1,
                                                                const std::vector<std::vector<double> > & s2,
                                                                const int                                 output_other,
                                                                std::string                               filename_prefix,
                                                                std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {     //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend_dummy = NULL;


    NRLib::Trend * mean_trend      = mean_    ->GenerateTrend(path,
                                                              trend_cube_parameters,
                                                              trend_cube_sampling,
                                                              blocked_logs,
                                                              s1,
                                                              s2,
                                                              NRLib::TrendStorage::MEAN,
                                                              mean_trend_dummy,
                                                              errTxt);

    NRLib::Trend * variance_trend  = variance_->GenerateTrend(path,
                                                              trend_cube_parameters,
                                                              trend_cube_sampling,
                                                              blocked_logs,
                                                              s1,
                                                              s2,
                                                              NRLib::TrendStorage::VAR,
                                                              mean_trend,
                                                              errTxt);

    // write mean and variance surface to file
    if ((output_other & IO::ROCK_PHYSICS_TRENDS) > 0) {
      if (mean_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_mean", trend_cube_sampling, mean_trend);
      if (variance_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_var",  trend_cube_sampling, variance_trend);
    }

    CheckBetaConsistency(mean_trend, variance_trend, lower_limit_, upper_limit_, errTxt);

    int share_level = DistributionWithTrend::None;
    if(is_shared_ == true)
      share_level = DistributionWithTrend::SingleSample;

    if(errTxt == "")
      distribution_with_trend_= new BetaDistributionWithTrend(mean_trend, variance_trend, lower_limit_, upper_limit_, share_level);

    delete mean_trend;
    delete variance_trend;
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
BetaDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

void
BetaDistributionWithTrendStorage::CheckBetaConsistency(NRLib::Trend * mean,
                                                       NRLib::Trend * variance,
                                                       double       & lower_limit,
                                                       double       & upper_limit,
                                                       std::string  & errTxt)
{

  bool mean_outside = false;
  if(mean->GetMinValue() < lower_limit || mean->GetMaxValue() > upper_limit) {
    errTxt += "The mean values of the Beta distribution must be in the interval ["+NRLib::ToString(lower_limit)+","+NRLib::ToString(upper_limit)+"]\n";
    mean_outside = true;
  }

  int n_samples_mean;
  int n_samples_var;
  int ni;
  int nj;

  if(mean->GetTrendDimension() == 0) {
    ni = 1;
    n_samples_mean = 1;
  }
  else {
    ni = 2;
    n_samples_mean = 4;
  }

  if(variance->GetTrendDimension() == 0) {
    nj = 1;
    n_samples_var = 1;
  }
  else {
    nj = 2;
    n_samples_var = 4;
  }

  double mean_min = mean->GetMinValue();
  double mean_max = mean->GetMaxValue();

  double var_min = variance->GetMinValue();
  double var_max = variance->GetMaxValue();

  double dx;
  if(ni == 1)
    dx = 0;
  else
    dx = (mean_max - mean_min)/(n_samples_mean - 1);

  double dy;
  if(nj == 1)
    dy = 0;
  else
    dy = (var_max  - var_min) /(n_samples_var - 1);

  std::vector<double> mean_sampling(n_samples_mean);
  for(int i=0; i<n_samples_mean; i++)
    mean_sampling[i] = mean_min + i*dx;

  std::vector<double> var_sampling(n_samples_var);
  for(int j=0; j<n_samples_var; j++)
    var_sampling[j]  = var_min  + j*dy;

  for(int i=0; i<n_samples_mean; i++) {
    for(int j=0; j<n_samples_var; j++) {

      double diff       = upper_limit-lower_limit;
      double moved_mean = (mean_sampling[i] - lower_limit) / diff;
      double moved_var  = var_sampling[j] / std::pow(diff,2);

      double alpha =    moved_mean  * (moved_mean*(1-moved_mean)/moved_var - 1);
      double beta  = (1-moved_mean) * (moved_mean*(1-moved_mean)/moved_var - 1);

      if(alpha <= 0 || beta <= 0) {
        errTxt += "The combination mean="+NRLib::ToString(mean_sampling[i])+", variance="+NRLib::ToString(var_sampling[j])+" in [a,b]=["+NRLib::ToString(lower_limit)+","+NRLib::ToString(upper_limit)+"]\n";
        errTxt += "  provides alpha<=0 and/or beta<=0 in the Beta distribution\n";
        if(mean_outside == true)
          errTxt += "  Check that the upper and lower limits in the Beta distribution are correct\n";
        else
          errTxt += "  The variance should be smaller than "+NRLib::ToString(moved_mean*(1-moved_mean)*std::pow(diff,2))+" for mean="+NRLib::ToString(mean_sampling[i])+"\n";
      }
      if(moved_var > moved_mean*(1-moved_mean)) {
        errTxt += "  The Beta distribution demands that v < e*(1-e)\n";
        errTxt += "  where e = (mean-a)/(b-a) and v = variance/(b-a)^2.\n";
      }
    }
  }
}
//--------------------------------------------------------------//

BetaEndMassDistributionWithTrendStorage::BetaEndMassDistributionWithTrendStorage()
: lower_limit_(0),
  upper_limit_(0),
  lower_probability_(0),
  upper_probability_(0),
  is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                    = NULL;
  variance_                = NULL;
  distribution_with_trend_ = NULL;
}

BetaEndMassDistributionWithTrendStorage::BetaEndMassDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                                 const NRLib::TrendStorage * variance,
                                                                                 const double              & lower_limit,
                                                                                 const double              & upper_limit,
                                                                                 const double              & lower_probability,
                                                                                 const double              & upper_probability,
                                                                                 bool                        is_shared)
: lower_limit_(lower_limit),
  upper_limit_(upper_limit),
  lower_probability_(lower_probability),
  upper_probability_(upper_probability),
  is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0),
  estimate_(false)
{
  mean_                     = mean->Clone();
  variance_                 = variance->Clone();
  distribution_with_trend_  = NULL;

  if(mean_->GetEstimate() == true || variance_->GetEstimate() == true)
    estimate_ = true;
}

BetaEndMassDistributionWithTrendStorage::BetaEndMassDistributionWithTrendStorage(const BetaEndMassDistributionWithTrendStorage & dist)
: lower_limit_(dist.upper_limit_),
  upper_limit_(dist.lower_limit_),
  lower_probability_(dist.lower_probability_),
  upper_probability_(dist.upper_probability_),
  is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_),
  estimate_(dist.estimate_)
{
  mean_     = dist.mean_    ->Clone();
  variance_ = dist.variance_->Clone();

 if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

BetaEndMassDistributionWithTrendStorage::~BetaEndMassDistributionWithTrendStorage()
{
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
BetaEndMassDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                       const std::vector<std::string>          & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                       const std::vector<std::vector<double> > & blocked_logs,
                                                                       const std::vector<std::vector<double> > & s1,
                                                                       const std::vector<std::vector<double> > & s2,
                                                                       const int                                 output_other,
                                                                       std::string                               filename_prefix,
                                                                       std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {     //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend_dummy = NULL;

    NRLib::Trend * mean_trend      = mean_    ->GenerateTrend(path,
                                                              trend_cube_parameters,
                                                              trend_cube_sampling,
                                                              blocked_logs,
                                                              s1,
                                                              s2,
                                                              NRLib::TrendStorage::MEAN,
                                                              mean_trend_dummy,
                                                              errTxt);

    NRLib::Trend * variance_trend  = variance_->GenerateTrend(path,
                                                              trend_cube_parameters,
                                                              trend_cube_sampling,
                                                              blocked_logs,
                                                              s1,
                                                              s2,
                                                              NRLib::TrendStorage::VAR,
                                                              mean_trend,
                                                              errTxt);

    // write mean and variance surface to file
    if ((output_other & IO::ROCK_PHYSICS_TRENDS) > 0) {
      if (mean_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_mean", trend_cube_sampling, mean_trend);
      if (variance_->GetEstimate() == true)
        WriteTrendToFile(filename_prefix + "_var",  trend_cube_sampling, variance_trend);
    }

    BetaDistributionWithTrendStorage::CheckBetaConsistency(mean_trend, variance_trend, lower_limit_, upper_limit_, errTxt);

    int share_level = DistributionWithTrend::None;
    if(is_shared_ == true)
      share_level = DistributionWithTrend::SingleSample;

    if(errTxt == "")
      distribution_with_trend_= new BetaEndMassDistributionWithTrend(mean_trend, variance_trend, lower_limit_, upper_limit_, lower_probability_, upper_probability_, share_level);

    delete mean_trend;
    delete variance_trend;
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
BetaEndMassDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

//--------------------------------------------------------------//

void DistributionWithTrendStorage::WriteTrendToFile(const std::string                       & filename,
                                                    const std::vector<std::vector<double> > & trend_cube_sampling,
                                                    const NRLib::Trend                      * trend) const
{
  if (typeid(*trend)  == typeid(NRLib::Trend1D)) {

    int                 s_ref = trend->GetReference();
    std::vector<double> s_tmp;

    if (0 < s_ref && s_ref < 3) {
      if (s_ref == 1)
        s_tmp = trend_cube_sampling[0];
      else if (s_ref == 2)
        s_tmp = trend_cube_sampling[1];


      std::vector<double> trend_vector(s_tmp.size(), RMISSING);

      int j_dummy = -999;
      int k_dummy = -999;
      for (int i = 0; i < static_cast<int>(s_tmp.size()); ++i) {
        if (s_ref == 1)
          trend_vector[i] = trend->GetTrendElement(i, j_dummy, k_dummy);
        else
          trend_vector[i] = trend->GetTrendElement(j_dummy, i, k_dummy);
      }
      std::string full_fileName = IO::makeFullFileName(IO::PathToRockPhysics(), filename + IO::Suffix1DTrend());
      NRLib::WriteTrend1D(full_fileName, s_tmp, trend_vector);
    }
  }

  else if (typeid(*trend)  == typeid(NRLib::Trend2D)) {

    std::vector<double> s1_tmp   = trend_cube_sampling[0];
    std::vector<double> s2_tmp   = trend_cube_sampling[1];

    NRLib::Grid2D<double> trend_grid2D(s1_tmp.size(), s2_tmp.size(), RMISSING);

    int k_dummy = -999;
    for (int i = 0; i < static_cast<int>(s1_tmp.size()); ++i) {
      for (int j = 0; j < static_cast<int>(s2_tmp.size()); ++j) {
        trend_grid2D(i, j) = trend->GetTrendElement(i, j, k_dummy);
      }
    }

    NRLib::RegularSurface<double> trend_surface(s1_tmp.front(), s2_tmp.front(), s1_tmp.back() - s1_tmp.front(), s2_tmp.back() - s2_tmp.front(), trend_grid2D);

    std::string full_fileName = IO::makeFullFileName(IO::PathToRockPhysics(), filename + IO::SuffixStormBinary()); //Should be able to specify suffix from model file
    trend_surface.WriteToFile(full_fileName);
  }
}
