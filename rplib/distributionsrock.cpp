#include "rplib/distributionsrock.h"
#include "rplib/rock.h"

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

void  DistributionsRock::SetupExpectationAndCovariances(const std::vector<double> & s_min,
                                                        const std::vector<double> & s_max)
{
  size_t n = 1024; // Number of samples generated for each distribution
  size_t m =   10; // Number of samples to use when sampling from s_min to s_max

  std::vector<double>   mean(3);
  NRLib::Grid2D<double> cov(3,3);

  NRLib::Grid2D<std::vector<double> > trend_params;

  FindTrendParams(trend_params,
                  HasTrend(),
                  s_min,
                  s_max,
                  m);

  size_t ni = trend_params.GetNI();
  size_t nj = trend_params.GetNJ();

  expectation_.Resize(ni, nj);
  covariance_.Resize(ni, nj);

  std::vector<double> a(n);
  std::vector<double> b(n);
  std::vector<double> c(n);

  for (size_t i = 0 ; i < ni ; i++) {
    for (size_t j = 0 ; j < nj ; j++) {

      const std::vector<double> & tp = trend_params(i,j); // trend_params = two-dimensional

      for (size_t k = 0 ; k < n ; k++) {
        Rock * rock = GenerateSample(tp);
        rock->GetSeismicParams(a[k], b[k], c[k]);
        delete rock;
      }

      //
      // Expectation
      //
      mean[0] = FindExpectation(a);
      mean[1] = FindExpectation(b);
      mean[2] = FindExpectation(c);

      expectation_(i, j) = mean;

      //
      // Covariance
      //
      cov(0,0) = FindCovariance(a,mean[0],a,mean[0]);
      cov(1,1) = FindCovariance(b,mean[1],b,mean[1]);
      cov(2,2) = FindCovariance(c,mean[2],c,mean[2]);

      cov(0,1) = FindCovariance(a,mean[0],b,mean[1]);
      cov(1,0) = cov(0,1);
      cov(0,2) = FindCovariance(a,mean[0],c,mean[2]);
      cov(2,0) = cov(0,2);
      cov(1,2) = FindCovariance(b,mean[1],c,mean[2]);
      cov(2,1) = cov(1,2);

      covariance_(i, j) = cov;

      // Temporary logging
      std::vector<double> s(3);
      s[0] = std::sqrt(cov(0,0));
      s[1] = std::sqrt(cov(1,1));
      s[2] = std::sqrt(cov(2,2));
      printf("Expectations :  %.6f %.6f %.6f\n", mean[0], mean[1], mean[2]);
      printf("Var          :  %.6f %.6f %.6f\n", cov(0,0), cov(1,1), cov(2,2));
      printf("Std          :  %.6f %.6f %.6f\n\n", s[0], s[1], s[2]);
      printf("cor ab: %.6f\n"  ,cov(0,1)/(s[0]*s[1]));
      printf("cor ac: %.6f\n"  ,cov(0,2)/(s[0]*s[2]));
      printf("cor bc: %.6f\n"  ,cov(1,2)/(s[1]*s[2]));
    }
  }

  exit(1);
}

void DistributionsRock::FindTrendParams(NRLib::Grid2D<std::vector<double> > & trend_params,
                                        const std::vector<bool>             & has_trend,
                                        const std::vector<double>           & s_min,
                                        const std::vector<double>           & s_max,
                                        const size_t                          n)
{
  bool t1 = has_trend[0];
  bool t2 = has_trend[1];

  std::vector<double> s0(1, 0.0); // use when first and/or second trend is missing
  std::vector<double> s1(n);
  std::vector<double> s2(n);

  if (t1 && !t2) {
    SampleTrendValues(s1, s_min[0], s_max[0]);
    SetupTrendMesh(trend_params, s1, s0);
  }
  else if (!t1 && t2) {
    SampleTrendValues(s2, s_min[1], s_max[1]);
    SetupTrendMesh(trend_params, s0, s2);
  }
  else if (t1 && t2) {
    SampleTrendValues(s0, s_min[0], s_max[0]);
    SampleTrendValues(s1, s_min[1], s_max[1]);
    SetupTrendMesh(trend_params, s1, s2);
  }
  else {
    SetupTrendMesh(trend_params, s0, s0);
  }
}

void DistributionsRock::SetupTrendMesh(NRLib::Grid2D<std::vector<double> > & trend_params,
                                       const std::vector<double>           & t1,
                                       const std::vector<double>           & t2)
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


void DistributionsRock::SampleTrendValues(std::vector<double> & s,
                                          const double        & s_min,
                                          const double        & s_max)
{
  size_t n    = s.size();
  double step = (s_max - s_min)/n;

  s[0]   = s_min;
  s[n-1] = s_max;

  for (size_t i = 1 ; i < n - 1 ; i++) {
    s[i] = s_min + step*(i - 1);
  }
}

double DistributionsRock::FindExpectation(const std::vector<double> & p)
{
  int    n    = p.size();
  double mean = 0.0;
  for (int i = 0 ; i < n ; i++) {
    mean += p[i];
  }
  mean /= n;
  return mean;
}

double DistributionsRock::FindCovariance(const std::vector<double> & p,
                                         const double                mup,
                                         const std::vector<double> & q,
                                         const double                muq)
{
  int    n   = p.size();
  double cov = 0.0;
  for (int i = 0 ; i < n ; i++) {
    cov += (p[i] - mup)*(q[i] - muq);
  }
  if (n > 1)
    cov /= n - 1;
  return cov;
}
