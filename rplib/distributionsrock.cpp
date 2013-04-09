#include "rplib/distributionsrock.h"
#include "rplib/rock.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/statistics/statistics.hpp"
#include <nrlib/flens/nrlib_flens.hpp>


Rock * DistributionsRock::GenerateSampleAndReservoirVariables(const std::vector<double> & trend_params, std::vector<double> &resVar )
{
  size_t nResVar=reservoir_variables_.size();
  resVar.resize(nResVar);
  //Note: If this is not a top level rock, reservoir_variables_ are not set, so the loop is skipped.
  for(size_t i=0;i<nResVar;i++)
    reservoir_variables_[i]->TriggerNewSample(resampling_level_);
  Rock * result = GenerateSamplePrivate(trend_params);
  
  for(size_t i=0;i<nResVar;i++)
    resVar[i]=reservoir_variables_[i]->GetCurrentSample(trend_params); 

  return(result);
}


Rock * DistributionsRock::GenerateSample(const std::vector<double> & trend_params )
{
  size_t nResVar=reservoir_variables_.size();
  //Note: If this is not a top level rock, reservoir_variables_ are not set, so the loop is skipped.
  for(size_t i=0;i<nResVar;i++)
    reservoir_variables_[i]->TriggerNewSample(resampling_level_);
  Rock * result = GenerateSamplePrivate(trend_params);
  return(result);
}





void DistributionsRock::GenerateWellSample(double                 corr,
                                           std::vector<double>  & vp,
                                           std::vector<double>  & vs,
                                           std::vector<double>  & rho,
                                           std::vector<double>  & trend_params)
{
  Rock * rock = GenerateSample(trend_params);
  rock->GetSeismicParams(vp[0],vs[0],rho[0]);

  for(size_t i=1;i<vp.size();i++) {
    // The time parameter in UpdateSample() is always false in this case
    Rock * rock_update = UpdateSample(corr, false, trend_params, rock);
    rock_update->GetSeismicParams(vp[i],vs[i],rho[i]);
    delete rock;
    rock = rock_update;
  }
  delete rock;
}

Rock * DistributionsRock::EvolveSample(double       time,
                                       const Rock & rock)
{
  //Note: If this is not a top level rock, reservoir_variables_ are not set, so the loop is skipped.
  for(size_t i=0;i<reservoir_variables_.size();i++)
    reservoir_variables_[i]->TriggerNewSample(resampling_level_);

  const std::vector<double> trend(2);
  return UpdateSample(time, true, trend, &rock);
}


Rock * DistributionsRock::EvolveSampleAndReservoirVaribles(double       time,
                                       const Rock & rock,
                                       std::vector<double> &resVar )
{
  size_t nResVar=reservoir_variables_.size();
  resVar.resize(nResVar);
  //Note: If this is not a top level rock, reservoir_variables_ are not set, so the loop is skipped.
  for(size_t i=0;i<reservoir_variables_.size();i++)
    reservoir_variables_[i]->TriggerNewSample(resampling_level_);

  const std::vector<double> trend(2);
  Rock * rock_update = UpdateSample(time, true, trend, &rock);

  for(size_t i=0;i<nResVar;i++)
    resVar[i]=reservoir_variables_[i]->GetCurrentSample(trend); 

  return  rock_update;

}



//-----------------------------------------------------------------------------------------------------------
void  DistributionsRock::SetupExpectationAndCovariances()
//-----------------------------------------------------------------------------------------------------------
{
  int n  = 1024*16; // Number of samples generated for each distribution
  int m  =   10; // Number of samples to use when sampling from s_min to s_max

  FindTabulatedTrendParams(tabulated_s0_,
                           tabulated_s1_,
                           HasTrend(),
                           s_min_,
                           s_max_,
                           m);

  int mi = static_cast<int>(tabulated_s0_.size());
  int mj = static_cast<int>(tabulated_s1_.size());

  expectation_.Resize(mi, mj);
  covariance_.Resize(mi, mj);

  NRLib::Grid2D<std::vector<double> > trend_params;

  SetupTrendMesh(trend_params,
                 tabulated_s0_,
                 tabulated_s1_);

  NRLib::Grid2D<double> cov(3,3);
  std::vector<double>   mean(3);
//  std::vector<double>   a(n);
 // std::vector<double>   b(n);
 // std::vector<double>   c(n);

  bool failed = false;

  for (int i = 0 ; i < mi ; i++) {
    for (int j = 0 ; j < mj ; j++) {

      const std::vector<double> & tp = trend_params(i,j); // trend_params = two-dimensional

      NRLib::Vector log_vp(n);
      NRLib::Vector log_vs(n);
      NRLib::Vector log_rho(n);

      double vp;
      double vs;
      double rho;

      for (int k = 0 ; k < n ; k++) {
        Rock * rock = GenerateSample(tp);

        rock->GetSeismicParams(vp, vs, rho);

        log_vp(k) = std::log(vp);
        log_vs(k) = std::log(vs);
        log_rho(k) = std::log(rho);

        delete rock;

        if(vp <= 0 || vs < 0 || rho <=0) {
          failed = true;
          break;
        }
      }

      std::vector<double>   expectation_small(3, 0.0);
      NRLib::Grid2D<double> covariance_small(3, 3, 0.0);

      if(failed == false) {

        std::vector<NRLib::Vector> m(3);

        m[0] = log_vp;
        m[1] = log_vs;
        m[2] = log_rho;

        for (int k = 0; k < 3; k++) {
          expectation_small[k] = NRLib::Mean(m[k]);
          for (int l = k; l < 3; l++) {
            covariance_small(k,l) = NRLib::Cov(m[k], m[l]);
            covariance_small(l,k) = covariance_small(k,l);
          }
        }
      }

      expectation_(i,j) = expectation_small;
      covariance_(i,j) = covariance_small;
    }
  }

  mean_log_expectation_.resize(3, 0);
  mean_log_covariance_.Resize(3, 3, 0);

  if(failed == false) {
    for(int i = 0 ; i < mi ; i++) {
      for(int j = 0 ; j < mj ; j++) {

        std::vector<double>   this_expectation = expectation_(i,j);
        NRLib::Grid2D<double> this_covariance  = covariance_(i,j);

        for(int k=0; k<3; k++)
          mean_log_expectation_[k] += this_expectation[k];

        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++)
            mean_log_covariance_(k,l) += this_covariance(k,l);
        }

      }
    }

    for(int k=0; k<3; k++)
      mean_log_expectation_[k] /= (mi*mj);

    for(int k=0; k<3; k++) {
      for(int l=0; l<3; l++)
        mean_log_covariance_(k,l) /= (mi*mj);
    }
  }
}

//----------------------------------------------------------------------------------------
void DistributionsRock::FindTabulatedTrendParams(std::vector<double>       & tabulated_s0,
                                                 std::vector<double>       & tabulated_s1,
                                                 const std::vector<bool>   & has_trend,
                                                 const std::vector<double> & s_min,
                                                 const std::vector<double> & s_max,
                                                 const size_t                m)
//----------------------------------------------------------------------------------------
{
  bool t1 = has_trend[0];
  bool t2 = has_trend[1];

  std::vector<double> no_trend(1, 0.0);

  if (t1 && !t2) {
    tabulated_s0.resize(m);
    SampleTrendValues(tabulated_s0, s_min[0], s_max[0]);
    tabulated_s1 = no_trend;
  }
  else if (!t1 && t2) {
    tabulated_s1.resize(m);
    SampleTrendValues(tabulated_s1, s_min[1], s_max[1]);
    tabulated_s0 = no_trend;
  }
  else if (t1 && t2) {
    tabulated_s0.resize(m);
    tabulated_s1.resize(m);
    SampleTrendValues(tabulated_s0, s_min[0], s_max[0]);
    SampleTrendValues(tabulated_s1, s_min[1], s_max[1]);
  }
  else {
    tabulated_s0 = no_trend;
    tabulated_s1 = no_trend;
  }

  //for (size_t i = 0 ; i < tabulated_s0.size() ; i++)
  //  std::cout << "tab_s0 " << i << " " << tabulated_s0[i] <<std::endl;
  //for (size_t i = 0 ; i < tabulated_s1.size() ; i++)
  //  std::cout << "tab_s1 " << i << " " << tabulated_s1[i] <<std::endl;
}

//----------------------------------------------------------------------------------------
void DistributionsRock::SetupTrendMesh(NRLib::Grid2D<std::vector<double> > & trend_params,
                                       const std::vector<double>           & t1,
                                       const std::vector<double>           & t2)
//----------------------------------------------------------------------------------------
{
  trend_params.Resize(t1.size(), t2.size());

  std::vector<double> params(2);

  for (size_t i = 0 ; i < t1.size() ; i++) {
    for (size_t j = 0 ; j < t2.size() ; j++) {
      params[0] = t1[i];
      params[1] = t2[j];
      trend_params(i,j) = params;
    }
  }
}


//--------------------------------------------------------------------
void DistributionsRock::SampleTrendValues(std::vector<double> & s,
                                          const double        & s_min,
                                          const double        & s_max)
//--------------------------------------------------------------------
{
  size_t n    = s.size();
  double step = (s_max - s_min)/(n - 1);

  s[0]   = s_min;
  s[n-1] = s_max;

  for (size_t i = 1 ; i < n - 1 ; i++) {
    s[i] = s_min + step*static_cast<double>(i);
  }
}

//------------------------------------------------------------------------------------------------------
std::vector<double> DistributionsRock::GetLogExpectation(const std::vector<double> & trend_params) const
//------------------------------------------------------------------------------------------------------
{
  double s0 = trend_params[0];
  double s1 = trend_params[1];

  CheckOrResetS(s0, tabulated_s0_);
  CheckOrResetS(s1, tabulated_s1_);

  size_t m  = tabulated_s0_.size();
  size_t n  = tabulated_s1_.size();

  if (m == 1 && n == 1) {
    return expectation_(0,0);
  }

  double di  = FindInterpolationStartIndex(tabulated_s0_, s0);
  double dj  = FindInterpolationStartIndex(tabulated_s1_, s1);
  size_t i0  = static_cast<size_t>(floor(di));
  size_t j0  = static_cast<size_t>(floor(dj));

  //printf("\ns0=%7.2f s1=%7.2f\n",s0,s1);
  //printf("tabulated_s0_min=%7.2f  tabulated_s0_max=%7.2f\n",tabulated_s0_[0],tabulated_s0_[m-2]);
  //printf("tabulated_s1_min=%7.2f  tabulated_s1_max=%7.2f\n",tabulated_s1_[0],tabulated_s1_[n-2]);
  //printf("di=%5.2f  dj=%5.2f   i0=%lu  j0=%lu   m=%lu  n=%lu\n",di,dj,i0,j0,m,n);

  double w00 = 0.0;
  double w10 = 0.0;
  double w01 = 0.0;
  double w11 = 0.0;

  FindInterpolationWeights(w00, w10, w01, w11, di, dj);

  bool do10 = m > 1 && i0 < m - 1;
  bool do01 = n > 1 && j0 < n - 1;
  bool do11 = do10 && do01;

  std::vector<double> mean(3);

  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 0); // Interpolate 1st parameter
  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 1); // Interpolate 2nd parameter
  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 2); // Interpolate 3rd parameter

  return mean;
}

//-----------------------------------------------------------------------------------
void DistributionsRock::CheckOrResetS(double                    & s,
                                      const std::vector<double> & tabulated_s) const
//-----------------------------------------------------------------------------------
{
  size_t n = tabulated_s.size();

  if (n > 1) { // There is a trend
    double smin = tabulated_s[0];
    double smax = tabulated_s[n-1];
    if (s < smin) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error,"\nERROR: First trend parameter (%.2f) is smaller than assumed lowest value (%.2f)\n", s, smin);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error,"       Setting trend parameter to lowest value.\n");
      s = smin;
    }
    if (s > smax) {
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error,"\nERROR: First trend parameter (%.2f) is larger than assumed largest value (%.2f)\n", s, smax);
      NRLib::LogKit::LogFormatted(NRLib::LogKit::Error,"       Setting trend parameter to largest value.\n");

      std::cout << std::setprecision(12) << s << " " << smax << std::endl;

      s = smax;
    }
  }
  else {
    s = 0.0;
  }
}


//--------------------------------------------------------------------------------------------
double DistributionsRock::FindInterpolationStartIndex(const std::vector<double> & tabulated_s,
                                                      const double                s) const
//--------------------------------------------------------------------------------------------
{
  double di = 0.0; // Default if there is no trend

  if (tabulated_s.size() > 1) {
    double dx;
    dx = tabulated_s[1] - tabulated_s[0]; // Assumes equally spaced table elements
    di = (s - tabulated_s[0])/dx;
  }
  return di;
}

//-------------------------------------------------------------------------
void DistributionsRock::FindInterpolationWeights(double       & w00,
                                                 double       & w10,
                                                 double       & w01,
                                                 double       & w11,
                                                 const double   di,
                                                 const double   dj) const
//-------------------------------------------------------------------------
{
  //                     w10 - w11
  // Find weights         |     |
  //                     w00 - w01
  //
  double di0 = floor(di);
  double dj0 = floor(dj);
  double u   = di - di0;
  double v   = dj - dj0;

  w00  = (1 - u)*(1 - v);
  w10  =      u *(1 - v);
  w01  = (1 - u)*     v;
  w11  =      u *     v;
}

//-----------------------------------------------------------------------------------------------------
void DistributionsRock::InterpolateExpectation(std::vector<double>                       & mean,
                                               const NRLib::Grid2D<std::vector<double> > & expectation,
                                               const double                                w00,
                                               const double                                w10,
                                               const double                                w01,
                                               const double                                w11,
                                               const size_t                                i0,
                                               const size_t                                j0,
                                               const bool                                  do10,
                                               const bool                                  do01,
                                               const bool                                  do11,
                                               const size_t                                p) const
//-----------------------------------------------------------------------------------------------------
{
  //
  // Find values in cell corners
  //
  double v00 = 0.0;
  double v10 = 0.0;
  double v01 = 0.0;
  double v11 = 0.0;

  v00 = expectation(i0, j0)[p];
  if (do10)
    v10 = expectation(i0 + 1, j0)[p];
  if (do01)
    v01 = expectation(i0, j0 + 1)[p];
  if (do11)
    v11 = expectation(i0 + 1, j0 + 1)[p];

  //
  // Interpolate ...
  //
  mean[p] = w00*v00 + w10*v10 + w01*v01 + w11*v11;

  //printf("i0=%lu  j0=%lu\n",i0,j0);
  //printf("v00 =%7.5f  w00=%5.3f\n",v00,w00);
  //printf("v10 =%7.5f  w10=%5.3f\n",v10,w10);
  //printf("v01 =%7.5f  w01=%5.3f\n",v01,w01);
  //printf("v11 =%7.5f  w11=%5.3f\n",v11,w11);
  //printf("mean=%7.5f\n\n",mean[p]);

  //  exit(1);

}


//----------------------------------------------------------------------------------------------------
NRLib::Grid2D<double> DistributionsRock::GetLogCovariance(const std::vector<double> & trend_params) const
//----------------------------------------------------------------------------------------------------
{
  double s0 = trend_params[0];
  double s1 = trend_params[1];

  CheckOrResetS(s0, tabulated_s0_);
  CheckOrResetS(s1, tabulated_s1_);

  size_t m  = tabulated_s0_.size();
  size_t n  = tabulated_s1_.size();

  if (m == 1 && n == 1) {
    return covariance_(0,0);
  }

  double di  = FindInterpolationStartIndex(tabulated_s0_, s0);
  double dj  = FindInterpolationStartIndex(tabulated_s1_, s1);
  size_t i0  = static_cast<size_t>(floor(di));
  size_t j0  = static_cast<size_t>(floor(dj));

  double w00 = 0.0;
  double w10 = 0.0;
  double w01 = 0.0;
  double w11 = 0.0;

  FindInterpolationWeights(w00, w10, w01, w11, di, dj);

  bool do10 = m > 1 && i0 < m - 1;
  bool do01 = n > 1 && j0 < n - 1;
  bool do11 = do10 && do01;

  NRLib::Grid2D<double> cov(3,3);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 0, 0);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 0, 1);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 0, 2);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 1, 1);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 1, 2);
  InterpolateCovariance(cov, covariance_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 2, 2);

  cov(1, 0) = cov(0, 1);
  cov(2, 0) = cov(0, 2);
  cov(2, 1) = cov(1, 2);

  return cov;
}

//-----------------------------------------------------------------------------------------------------
void DistributionsRock::InterpolateCovariance(NRLib::Grid2D<double>                       & cov,
                                              const NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                              const double                                  w00,
                                              const double                                  w10,
                                              const double                                  w01,
                                              const double                                  w11,
                                              const size_t                                  i0,
                                              const size_t                                  j0,
                                              const bool                                    do10,
                                              const bool                                    do01,
                                              const bool                                    do11,
                                              const size_t                                  p,
                                              const size_t                                  q) const
//-----------------------------------------------------------------------------------------------------
{
  //
  // Find values in cell corners
  //
  double v00 = 0.0;
  double v10 = 0.0;
  double v01 = 0.0;
  double v11 = 0.0;

  v00 = covariance(i0, j0)(p,q);
  if (do10)
    v10 = covariance(i0 + 1, j0)(p,q);
  if (do01)
    v01 = covariance(i0, j0 + 1)(p,q);
  if (do11)
    v11 = covariance(i0 + 1, j0 + 1)(p,q);

  //
  // Interpolate ...
  //
  cov(p,q) = w00*v00 + w10*v10 + w01*v01 + w11*v11;
}


void DistributionsRock::CompleteTopLevelObject(std::vector<DistributionWithTrend *> res_var)
{
  reservoir_variables_ = res_var;
  SetResamplingLevel(DistributionWithTrend::Full);
  SetupExpectationAndCovariances();
}
