
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/beta.hpp"
#include "nrlib/random/uniform.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/betadistributionwithtrend.h"

BetaDistributionWithTrend::BetaDistributionWithTrend()
: is_shared_(false)
{
}

BetaDistributionWithTrend::BetaDistributionWithTrend(const NRLib::Trend * mean,
                                                     const NRLib::Trend * var,
                                                     bool                 shared)
: mean_(mean),
  var_(var),
  is_shared_(shared)
{
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

  beta_distribution_ = new NRLib::Grid2D<NRLib::Distribution<double> *>(n_samples_mean_, n_samples_var_, NULL);

  for(int i=0; i<n_samples_mean_; i++) {
    for(int j=0; j<n_samples_var_; j++) {

      CalculateAlpha(mean_sampling_[i], var_sampling_[j], a);
      CalculateBeta(mean_sampling_[i],  var_sampling_[j], b);

      NRLib::Distribution<double> * dist = new NRLib::Beta(0, 1, a, b);

      (*beta_distribution_)(i,j) = dist;

    }
  }
}

BetaDistributionWithTrend::BetaDistributionWithTrend(const BetaDistributionWithTrend & dist)
: beta_distribution_(dist.beta_distribution_),
  mean_(dist.mean_),
  var_(dist.var_),
  is_shared_(dist.is_shared_),
  use_trend_cube_(dist.use_trend_cube_),
  ni_(dist.ni_),
  nj_(dist.nj_),
  mean_sampling_(dist.mean_sampling_),
  var_sampling_(dist.var_sampling_),
  n_samples_mean_(dist.n_samples_mean_),
  n_samples_var_(dist.n_samples_var_)
{
}

BetaDistributionWithTrend::~BetaDistributionWithTrend()
{
  for(int i=0; i<n_samples_mean_; i++) {
    for(int j=0; j<n_samples_var_; j++)
      delete (*beta_distribution_)(i,j);
  }
  delete beta_distribution_;

  delete mean_;
  delete var_;
}

double
BetaDistributionWithTrend::ReSample(double s1, double s2) const
{
  double uniform     = NRLib::Random::Unif01();
  double probability = GetQuantileValue(uniform, s1, s2);

  return(probability);
}

double
BetaDistributionWithTrend::GetQuantileValue(double u, double s1, double s2) const
{
  // Want sample from Y(s1, s2) ~ Beta(0, 1, alpha(s1,s2), beta(s1,s2))
  // Generate sample from Z ~ Uniform(0, 1)
  // Calculate X_{ij}~Beta(alpha_i,beta_j) using quantile mapping in locations (i,j) around Y(s1,s2) where the Beta distribution is initialized
  // Calculate y(s1, s2) by interpolating between the generated X_{ij}

  double y;

  if(ni_ == 1 && nj_ == 1)
    y = (*beta_distribution_)(0,0)->Quantile(u);

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
        NRLib::Distribution<double> * dist = (*beta_distribution_)(ix+i, jy+j);
        quantile_grid(i,j) = dist->Quantile(u);
      }
    }

    NRLib::RegularSurface<double> surf(mean_sampling_[ix], var_sampling_[jy], li, lj, quantile_grid);

    y = surf.GetZ(mean_trend, var_trend);
  }

  return(y);
}

void
BetaDistributionWithTrend::CalculateAlpha(double mean, double var, double & alpha) const
{
  alpha = mean * (mean*(1-mean)/var - 1);
}

void
BetaDistributionWithTrend::CalculateBeta(double mean, double var, double & beta) const
{
  beta = (1-mean) * (mean*(1-mean)/var - 1);
}
