#include "rplib/distributionsrock.h"
#include "rplib/rock.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/grid/grid2d.hpp"

void DistributionsRock::GenerateWellSample(double                 corr,
                                           std::vector<double>  & vp,
                                           std::vector<double>  & vs,
                                           std::vector<double>  & rho,
                                           std::vector<double>  & trend_params) const
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
                                       const Rock & rock) const
{
    const std::vector<double> trend(2);
    return UpdateSample(time, true, trend, &rock);
}


//-----------------------------------------------------------------------------------------------------------
void  DistributionsRock::SetupExpectationAndCovariances(NRLib::Grid2D<std::vector<double> >   & expectation,
                                                        NRLib::Grid2D<NRLib::Grid2D<double> > & covariance,
                                                        std::vector<double>                   & tabulated_s0,
                                                        std::vector<double>                   & tabulated_s1,
                                                        const std::vector<double>             & s_min,
                                                        const std::vector<double>             & s_max)
//-----------------------------------------------------------------------------------------------------------
{
  size_t n  = 1024; // Number of samples generated for each distribution
  size_t m  =   10; // Number of samples to use when sampling from s_min to s_max

  FindTabulatedTrendParams(tabulated_s0,
                           tabulated_s1,
                           HasTrend(),
                           s_min,
                           s_max,
                           m);

  size_t mi = tabulated_s0.size();
  size_t mj = tabulated_s1.size();

  expectation.Resize(mi, mj);
  covariance.Resize(mi, mj);

  NRLib::Grid2D<std::vector<double> > trend_params;

  SetupTrendMesh(trend_params,
                 tabulated_s0,
                 tabulated_s1);

  NRLib::Grid2D<double> cov(3,3);
  std::vector<double>   mean(3);
  std::vector<double>   a(n);
  std::vector<double>   b(n);
  std::vector<double>   c(n);

  for (size_t i = 0 ; i < mi ; i++) {
    for (size_t j = 0 ; j < mj ; j++) {

      const std::vector<double> & tp = trend_params(i,j); // trend_params = two-dimensional

      for (size_t k = 0 ; k < n ; k++) {
        Rock * rock = GenerateSample(tp);
        rock->GetSeismicParams(a[k], b[k], c[k]);
        delete rock;
      }

      //
      // Expectation
      //
      mean[0] = FindLogExpectation(a);
      mean[1] = FindLogExpectation(b);
      mean[2] = FindLogExpectation(c);

      expectation(i, j) = mean;

      //
      // Covariance
      //
      cov(0,0) = FindLogCovariance(a,mean[0],a,mean[0]);
      cov(1,1) = FindLogCovariance(b,mean[1],b,mean[1]);
      cov(2,2) = FindLogCovariance(c,mean[2],c,mean[2]);

      cov(0,1) = FindLogCovariance(a,mean[0],b,mean[1]);
      cov(1,0) = cov(0,1);
      cov(0,2) = FindLogCovariance(a,mean[0],c,mean[2]);
      cov(2,0) = cov(0,2);
      cov(1,2) = FindLogCovariance(b,mean[1],c,mean[2]);
      cov(2,1) = cov(1,2);

      covariance(i, j) = cov;

      // Temporary logging
      std::vector<double> s(3);
      s[0] = std::sqrt(cov(0,0));
      s[1] = std::sqrt(cov(1,1));
      s[2] = std::sqrt(cov(2,2));

      /*
      printf("Expectations :  %.6f %.6f %.6f\n", mean[0], mean[1], mean[2]);
      printf("Var          :  %.6f %.6f %.6f\n", cov(0,0), cov(1,1), cov(2,2));
      printf("Std          :  %.6f %.6f %.6f\n\n", s[0], s[1], s[2]);
      printf("cor ab: %.6f\n"  ,cov(0,1)/(s[0]*s[1]));
      printf("cor ac: %.6f\n"  ,cov(0,2)/(s[0]*s[2]));
      printf("cor bc: %.6f\n"  ,cov(1,2)/(s[1]*s[2]));
      */
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

  std::vector<double> no_trend(1, 0.0); // use when first and/or second trend is missing

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

//----------------------------------------------------------------------
double DistributionsRock::FindLogExpectation(const std::vector<double> & p)
//----------------------------------------------------------------------
{
  int    n    = static_cast<int>(p.size());
  double mean = 0.0;
  for (int i = 0 ; i < n ; i++) {
    mean += std::log(p[i]);
  }
  mean /= n;
  return mean;
}

//----------------------------------------------------------------------
double DistributionsRock::FindLogCovariance(const std::vector<double> & p,
                                            const double                mup,
                                            const std::vector<double> & q,
                                            const double                muq)
//----------------------------------------------------------------------
{
  int    n   = static_cast<int>(p.size());
  double cov = 0.0;
  for (int i = 0 ; i < n ; i++) {
    cov += (std::log(p[i]) - mup)*(std::log(q[i]) - muq);
  }
  if (n > 1)
    cov /= n - 1;
  return cov;
}

//---------------------------------------------------------------------------------------------------
std::vector<double> DistributionsRock::GetLogExpectation(const std::vector<double> & trend_params) const
//---------------------------------------------------------------------------------------------------
{
  /* would not compile ....
  if (trend_params[0] < s_min[0]) {
    LogKit::LogFormatted(LogKit::Error,"ERROR: First trend parameter (%.2f) is smaller than assumed lowest value (%.2f)\n", trend_params[0], s_min[0]);
    LogKit::LogFormatted(LogKit::Error,"       Setting trend parameter to lowest value.\n");
  }
  if (trend_params[0] > s_max[0]) {
    LogKit::LogFormatted(LogKit::Error,"ERROR: First trend parameter (%.2f) is larger than assumed largest value (%.2f)\n", trend_params[0], s_max[0]);
    LogKit::LogFormatted(LogKit::Error,"       Setting trend parameter to largest value.\n");
  }
  if (trend_params[1] < s_min[1]) {
    LogKit::LogFormatted(LogKit::Error,"ERROR: First trend parameter (%.2f) is smaller than assumed lowest value (%.2f)\n", trend_params[1], s_min[1]);
    LogKit::LogFormatted(LogKit::Error,"       Setting trend parameter to lowest value.\n");
  }
  if (trend_params[1] > s_max[1]) {
    LogKit::LogFormatted(LogKit::Error,"ERROR: First trend parameter (%.2f) is larger than assumed largest value (%.2f)\n", trend_params[1], s_max[1]);
    LogKit::LogFormatted(LogKit::Error,"       Setting trend parameter to largest value.\n");
  }
  */

  double s0 = trend_params[0];
  double s1 = trend_params[1];

  size_t m  = tabulated_s0_.size();
  size_t n  = tabulated_s1_.size();

  if (m == 1 && n == 1) {
    return expectation_(0,0);
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

  std::vector<double> mean(3);

  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 0); // Interpolate 1st parameter
  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 1); // Interpolate 2nd parameter
  InterpolateExpectation(mean, expectation_, w00, w10, w01, w11, i0, j0, do10, do01, do11, 2); // Interpolate 3rd parameter

  return mean;
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
    di = s/dx;
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

  //printf("di=%5.2f  dj=%5.2f   m=%lu  n=%lu\n",di,dj,m,n);
  //printf("v00 =%7.2f  w00=%5.3f\n",v00,w00);
  //printf("v10 =%7.2f  w10=%5.3f\n",v10,w10);
  //printf("v01 =%7.2f  w01=%5.3f\n",v01,w01);
  //printf("v11 =%7.2f  w11=%5.3f\n",v11,w11);
  //printf("mean=%7.2f\n\n",mean[p]);
}


//----------------------------------------------------------------------------------------------------
NRLib::Grid2D<double> DistributionsRock::GetLogCovariance(const std::vector<double> & trend_params) const
//----------------------------------------------------------------------------------------------------
{
  double s0 = trend_params[0];
  double s1 = trend_params[1];

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
