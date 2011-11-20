#include "rplib/multinormalwithtrend.h"

#include "rplib/trend.h"

#include "nrlib/exception/exception.hpp"

#include "lib/lib_matr.h"




MultiNormalWithTrend::
MultiNormalWithTrend(const NRLib::Distribution<double>&        vp, 
                     const NRLib::Distribution<double>&        vs, 
                     const NRLib::Distribution<double>&        rho,
                     const Trend&                              mean_trend_vp,
                     const Trend&                              mean_trend_vs,
                     const Trend&                              mean_trend_rho,
                     const std::vector< std::vector<Trend*> >  cov_trend) :
  vp_(vp), vs_(vs), rho_(rho), 
  mean_trend_vp_(mean_trend_vp), mean_trend_vs_(mean_trend_vs), mean_trend_rho_(mean_trend_rho), 
  cov_trend_(cov_trend)
{
  if (cov_trend_.size() != 3)
    throw NRLib::IndexOutOfRange("MultiNormalWithTrend: Number of Cov trend rows is not valid.");

  unsigned int i;
  for (i = 0; i < cov_trend_.size(); ++i) {
    if (cov_trend_[i].size() != 3)
      throw NRLib::IndexOutOfRange("MultiNormalWithTrend: Number of Cov trend columns is not valid.");
  }

}

MultiNormalWithTrend::~MultiNormalWithTrend()
{
  
}

void   
MultiNormalWithTrend::ReSample(double s1, double s2,
                               double& vp, double& vs, double& rho) const {
  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp_.Draw();
  rc[1] = vs_.Draw();
  rc[2] = rho_.Draw();

  // fill cov matrix
  double ** cov_matrix = new double*[3];
  for (int i = 0; i < 3; ++i)
    cov_matrix[i] = new double[3];

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) 
      cov_matrix[j][i] = cov_trend_[j][i]->GetValue(s1, s2);
  }
  
  //cholesky 
  lib_matrCholR(3, cov_matrix); 

  MatrProdTranspCholVecRR(3, cov_matrix, &rc[0] , &rhs[0]); // rhs is given value 

  vp  =  rhs[0] + mean_trend_vp_.GetValue(s1, s2);
  vs  =  rhs[1] + mean_trend_vs_.GetValue(s1, s2);
  rho =  rhs[2] + mean_trend_rho_.GetValue(s1, s2);

  //clean up
  for (int i = 0; i < 3; ++i)
    delete [] cov_matrix[i];
  delete [] cov_matrix;

}

void
MultiNormalWithTrend::
MatrProdTranspCholVecRR(int n, double **mat, double *in_vec,double* out_vec) const 
{
  int i,j;
  for (i = 0; i < n; i++)
  { 
    out_vec[i]=0.0;
    for (j = 0; j <= i; j++)
    {
      out_vec[i] += mat[i][j]*in_vec[j]; 
    }
  }
}

