
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/betawithendmass.hpp"
#include "nrlib/random/uniform.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/betadistributionwithtrend.h"
#include "rplib/betaendmassdistributionwithtrend.h"

BetaEndMassDistributionWithTrend::BetaEndMassDistributionWithTrend()
: DistributionWithTrend()
{
  beta_endmass_distribution_ = NULL;
  mean_                      = NULL;
  var_                       = NULL;
  ni_                        = 0;
  nj_                        = 0;
  n_samples_mean_            = 0;
  n_samples_var_             = 0;
}

BetaEndMassDistributionWithTrend::BetaEndMassDistributionWithTrend(const NRLib::Trend * mean,
                                                                   const NRLib::Trend * var,
                                                                   const double       & lower_limit,
                                                                   const double       & upper_limit,
                                                                   const double       & lower_probability,
                                                                   const double       & upper_probability,
                                                                   int                  shared)
: DistributionWithTrend(shared,true)
{
  mean_ = mean->Clone();
  var_  = var->Clone();

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(use_trend_cube_, mean_->GetTrendDimension(), mean_->GetReference());
  FindUseTrendCube(use_trend_cube_, var_ ->GetTrendDimension(), var_ ->GetReference());

  if(mean_->GetTrendDimension() == 0) {
    ni_ = 1;
    n_samples_mean_ = 1;
  }
  else {
    ni_ = 2;
    n_samples_mean_ = 100;
  }

  if(var_->GetTrendDimension() == 0) {
    nj_ = 1;
    n_samples_var_ = 1;
  }
  else {
    nj_ = 2;
    n_samples_var_ = 100;
  }

  double mean_min = mean_->GetMinValue();
  double mean_max = mean_->GetMaxValue();

  double var_min = var_->GetMinValue();
  double var_max = var_->GetMaxValue();

  double dx;
  if(ni_ == 1)
    dx = 0;
  else
    dx = (mean_max - mean_min)/(n_samples_mean_ - 1);

  double dy;
  if(nj_ == 1)
    dy = 0;
  else
    dy = (var_max  - var_min) /(n_samples_var_ - 1);

  mean_sampling_.resize(n_samples_mean_);
  for(int i=0; i<n_samples_mean_; i++)
    mean_sampling_[i] = mean_min + i*dx;

  var_sampling_.resize(n_samples_var_);
  for(int j=0; j<n_samples_var_; j++)
    var_sampling_[j]  = var_min  + j*dy;

  double a;
  double b;

  beta_endmass_distribution_ = new NRLib::Grid2D<NRLib::Distribution<double> *>(n_samples_mean_, n_samples_var_, NULL);

  for(int i=0; i<n_samples_mean_; i++) {
    for(int j=0; j<n_samples_var_; j++) {

      BetaDistributionWithTrend::CalculateAlpha(mean_sampling_[i], var_sampling_[j], lower_limit, upper_limit, a);
      BetaDistributionWithTrend::CalculateBeta(mean_sampling_[i],  var_sampling_[j], lower_limit, upper_limit, b);

      NRLib::Distribution<double> * dist = NULL;

      if(a > 0 && b > 0)
        dist = new NRLib::BetaWithEndMass(lower_limit, upper_limit, lower_probability, upper_probability, a, b);

      (*beta_endmass_distribution_)(i,j) = dist;

    }
  }
}

BetaEndMassDistributionWithTrend::BetaEndMassDistributionWithTrend(const BetaEndMassDistributionWithTrend & dist)
: DistributionWithTrend(dist.share_level_,dist.current_u_,dist.resample_),
  use_trend_cube_(dist.use_trend_cube_),
  ni_(dist.ni_),
  nj_(dist.nj_),
  mean_sampling_(dist.mean_sampling_),
  var_sampling_(dist.var_sampling_),
  n_samples_mean_(dist.n_samples_mean_),
  n_samples_var_(dist.n_samples_var_)
{
  beta_endmass_distribution_ = new NRLib::Grid2D<NRLib::Distribution<double> *>(n_samples_mean_, n_samples_var_, NULL);

  NRLib::Grid2D<NRLib::Distribution<double> *> * dist_beta = dist.beta_endmass_distribution_;

  for(int i=0; i<n_samples_mean_; i++) {
    for(int j=0; j<n_samples_var_; j++)
      (*beta_endmass_distribution_)(i,j) = (*dist_beta)(i,j)->Clone();
  }

  mean_ = dist.mean_->Clone();
  var_  = dist.var_ ->Clone();
}

BetaEndMassDistributionWithTrend::~BetaEndMassDistributionWithTrend()
{
  for(int i=0; i<n_samples_mean_; i++) {
    for(int j=0; j<n_samples_var_; j++)
      delete (*beta_endmass_distribution_)(i,j);
  }
  delete beta_endmass_distribution_;

  delete mean_;
  delete var_;
}

double
BetaEndMassDistributionWithTrend::ReSample(double s1, double s2)
{
  double uniform = NRLib::Random::Unif01();
  double value   = GetQuantileValue(uniform, s1, s2);

  return value;
}

double
BetaEndMassDistributionWithTrend::GetMeanValue(double s1, double s2) {
  double dummy = 0.0;
  return(mean_->GetValue(s1, s2, dummy));
}
double
BetaEndMassDistributionWithTrend::GetVarianceValue(double s1, double s2) {
  double dummy = 0.0;
  return(var_->GetValue(s1, s2, dummy));
}


double
BetaEndMassDistributionWithTrend::GetQuantileValue(double u, double s1, double s2)
{
  // Want sample from Y(s1, s2) ~ Beta(0, 1, alpha(s1,s2), beta(s1,s2))
  // Generate sample from Z ~ Uniform(0, 1)
  // Calculate X_{ij}~Beta(alpha_i,beta_j) using quantile mapping in locations (i,j) around Y(s1,s2) where the Beta distribution is initialized
  // Calculate y(s1, s2) by interpolating between the generated X_{ij}

  double y;

  if(share_level_ > None && resample_ == false)
    u = current_u_;
  else {
    current_u_ = u;
    resample_ = false;
  }

  if(ni_ == 1 && nj_ == 1)
    y = (*beta_endmass_distribution_)(0,0)->Quantile(u);

  else {
    double dummy = 0;

    double mean_trend = mean_->GetValue(s1, s2, dummy);
    double var_trend  = var_ ->GetValue(s1, s2, dummy);

    int ix;
    for(ix = 0; ix < n_samples_mean_; ix++) {
      if(mean_trend >= mean_sampling_[ix])
        break;
    }

    int jy;
    for(jy = 0; jy < n_samples_var_; jy++) {
      if(var_trend >= var_sampling_[jy])
        break;
    }

    double li;
    if(ni_ == 1)
      li = 0;
    else
      li = mean_sampling_[ix+1] - mean_sampling_[ix];

    double lj;
    if(nj_ == 1)
      lj = 0;
    else
      lj = var_sampling_[jy+1] - var_sampling_[jy];

    NRLib::Grid2D<double> quantile_grid(ni_, nj_, 0);

    for(int i=0; i<ni_; i++) {
      for(int j=0; j<nj_; j++) {
        NRLib::Distribution<double> * dist = (*beta_endmass_distribution_)(ix+i, jy+j);
        quantile_grid(i,j) = dist->Quantile(u);
      }
    }

    NRLib::RegularSurface<double> surf(mean_sampling_[ix], var_sampling_[jy], li, lj, quantile_grid);

    y = surf.GetZ(mean_trend, var_trend);

  }

  return(y);
}
