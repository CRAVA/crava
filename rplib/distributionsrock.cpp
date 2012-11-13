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
  //  int n = 1024;
  int n = 10000;

  //
  // er det riktig å sette opp en vektor for dette når trend parametrene ikke eksisterer?
  //
  std::vector<double> trend_params(2);
  trend_params[0] = 0.0;
  trend_params[1] = 0.0;

  const std::vector<bool> & has_trend = HasTrend();
  if (has_trend[0]) {
    trend_params[0] = 0.0;    // Sette p1-trend-verdi her.
    if (has_trend[1])
      trend_params[1] = 0.0;  // Sette p2-trend-verdi her.
    else if (has_trend[2])
      trend_params[1] = 0.0;  // Sette p3-trend-verdi her.
  }
  else if (has_trend[1]) {
    trend_params[0] = 0.0;    // Sette p2-trend-verdi her.
    if (has_trend[2])
      trend_params[1] = 0.0;  // Sette p3-trend-verdi her.
  }
  else if (has_trend[2]) {
    trend_params[0] = 0.0;    // Sette p3-trend-verdi her.
  }

  //create_grid()


  std::vector<double> a(n);
  std::vector<double> b(n);
  std::vector<double> c(n);
  for (int i = 0 ; i < n ; i++) {
    Rock * rock = GenerateSample(trend_params);
    rock->GetSeismicParams(a[i], b[i], c[i]);
    delete rock;
  }

  //
  // Expectation
  //
  std::vector<double> mean(3);
  mean[0] = FindExpectation(a);
  mean[1] = FindExpectation(b);
  mean[2] = FindExpectation(c);

  printf("Expectations :  %.6f %.6f %.6f\n\n", mean[0], mean[1], mean[2]);

  //
  // Covariance
  //
  //  NRLib::Grid2D<double> corr_matrix(3,3,0);

  std::vector<double> cov(6);

  cov[0] = FindCovariance(a,mean[0],a,mean[0]);
  cov[1] = FindCovariance(b,mean[1],b,mean[1]);
  cov[2] = FindCovariance(c,mean[2],c,mean[2]);

  cov[3] = FindCovariance(a,mean[0],b,mean[1]);
  cov[4] = FindCovariance(a,mean[0],c,mean[2]);
  cov[5] = FindCovariance(b,mean[1],c,mean[2]);

  std::vector<double> s(3);
  s[0] = std::sqrt(cov[0]);
  s[1] = std::sqrt(cov[1]);
  s[2] = std::sqrt(cov[2]);

  printf("Std :  %.6f %.6f %.6f\n\n", s[0], s[1], s[2]);

  printf("cov :  %.6f %.6f %.6f\n"         ,cov[0],cov[3],cov[4]);
  printf("cov :            %.6f %.6f\n"    ,cov[1],cov[5]);
  printf("cov :                     %.6f\n",cov[2]);

  printf("\n");

  printf("cor aa: %.6f\n"  ,cov[0]/(s[0]*s[0]));
  printf("cor bb: %.6f\n"  ,cov[1]/(s[1]*s[1]));
  printf("cor cc: %.6f\n"  ,cov[2]/(s[2]*s[2]));
  printf("\n");
  printf("cor ab: %.6f\n"  ,cov[3]/(s[0]*s[1]));
  printf("cor ac: %.6f\n"  ,cov[4]/(s[0]*s[2]));
  printf("cor bc: %.6f\n"  ,cov[5]/(s[1]*s[2]));

  exit(1);
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
