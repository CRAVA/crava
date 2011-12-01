#include "rplib/multinormalwithtrend.h"

#include "rplib/trend.h"

#include "nrlib/exception/exception.hpp"

#include "lib/lib_matr.h"

#include <cmath>




MultiNormalWithTrend::
MultiNormalWithTrend(const NRLib::Normal&                      vp, 
                     const NRLib::Normal&                      vs, 
                     const NRLib::Normal&                      rho,
                     const Trend&                              mean_trend_vp,
                     const Trend&                              mean_trend_vs,
                     const Trend&                              mean_trend_rho,
                     const std::vector< std::vector<Trend*> >& cov_trend) :
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
  double ** cov_matrix = CreateCovMatrix(s1, s2);
  
  //cholesky 
  lib_matrCholR(3, cov_matrix); 

  MatrProdTranspCholVecRR(3, cov_matrix, &rc[0] , &rhs[0]); // rhs is given value 

  vp  =  rhs[0] + mean_trend_vp_.GetValue(s1, s2);
  vs  =  rhs[1] + mean_trend_vs_.GetValue(s1, s2);
  rho =  rhs[2] + mean_trend_rho_.GetValue(s1, s2);

  DeleteCovMatrix(cov_matrix);

}

void     
MultiNormalWithTrend::ReSample(double s1, double s2,
                               double** cov_matrix_cholesky,
                               double& vp, double& vs, double& rho,
                               bool is_cholesky) const {

  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp_.Draw();
  rc[1] = vs_.Draw();
  rc[2] = rho_.Draw();
  
  //cholesky 
  if (!is_cholesky)
    lib_matrCholR(3, cov_matrix_cholesky); 

  MatrProdTranspCholVecRR(3, cov_matrix_cholesky, &rc[0] , &rhs[0]); // rhs is given value 

  vp  =  rhs[0] + mean_trend_vp_.GetValue(s1, s2);
  vs  =  rhs[1] + mean_trend_vs_.GetValue(s1, s2);
  rho =  rhs[2] + mean_trend_rho_.GetValue(s1, s2);

}

void     
MultiNormalWithTrend::EstimateExpectation(double s1, double s2, int sample_size, 
                                          double& exp_vp, double& exp_vs, double& exp_rho) const {
  
  exp_vp = exp_vs = exp_rho = 0;

  if (sample_size <= 0)
    return;

  // fill cov matrix
  double ** cov_matrix = CreateCovMatrix(s1, s2);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  for (int i = 0; i < sample_size; ++i) {
    double vp, vs, rho;
    ReSample(s1, s2, cov_matrix,
             vp, vs, rho);

    exp_vp  += vp;
    exp_vs  += vs;
    exp_rho += rho;

  }

  exp_vp  /= sample_size;
  exp_vs  /= sample_size;
  exp_rho /= sample_size;

  DeleteCovMatrix(cov_matrix);

}

void     
MultiNormalWithTrend::EstimateExpectationAndVariance(double s1, double s2, int sample_size,
                                                     double& exp_vp, double& exp_vs, double& exp_rho,
                                                     double& var_vp, double& var_vs, double& var_rho) const {
  
  exp_vp = exp_vs = exp_rho = 0;
  var_vp = var_vs = var_rho = 0;

  if (sample_size <= 0)
    return;

  // fill cov matrix
  double ** cov_matrix = CreateCovMatrix(s1, s2);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  for (int i = 0; i < sample_size; ++i) {
    double vp, vs, rho;
    ReSample(s1, s2, cov_matrix,
             vp, vs, rho);

    exp_vp  += vp;
    exp_vs  += vs;
    exp_rho += rho;

    var_vp  += vp*vp;
    var_vs  += vs*vs;
    var_rho += rho*rho;

  }

  exp_vp  /= sample_size;
  exp_vs  /= sample_size;
  exp_rho /= sample_size;

  var_vp   = var_vp/sample_size - exp_vp*exp_vp;
  var_vs   = var_vs/sample_size - exp_vs*exp_vs;
  var_rho  = var_rho/sample_size - exp_rho*exp_rho;

  DeleteCovMatrix(cov_matrix);

}

double** 
MultiNormalWithTrend::CreateEstimateOfCovMatrix(double s1, double s2, int sample_size) const {
  if (sample_size <= 0)
    return NULL;

  double exp_vp, exp_vs, exp_rho;
  double var_vp, var_vs, var_rho;

  exp_vp = exp_vs = exp_rho = 0;
  var_vp = var_vs = var_rho = 0;
  

  // fill cov matrix
  double ** cov_matrix = CreateCovMatrix(s1, s2);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  double vp_vs   = 0;
  double vp_rho  = 0;
  double vs_rho  = 0;

  for (int i = 0; i < sample_size; ++i) {
    double vp = 0, vs = 0, rho = 0;
    ReSample(s1, s2, cov_matrix,
             vp, vs, rho);

    vp_vs    += vp*vs;
    vp_rho   += vp*rho;
    vs_rho   += vs*rho;

    exp_vp  += vp;
    exp_vs  += vs;
    exp_rho += rho;

    var_vp  += vp*vp;
    var_vs  += vs*vs;
    var_rho += rho*rho;
  }

  exp_vp  /= sample_size;
  exp_vs  /= sample_size;
  exp_rho /= sample_size;

  var_vp   = var_vp/sample_size - exp_vp*exp_vp;
  var_vs   = var_vs/sample_size - exp_vs*exp_vs;
  var_rho  = var_rho/sample_size - exp_rho*exp_rho;

  vp_vs     = vp_vs/sample_size - exp_vp*exp_vs;
  vp_rho    = vp_rho/sample_size - exp_vp*exp_rho;
  vs_rho    = vs_rho/sample_size - exp_vs*exp_rho;

  cov_matrix[0][0] = var_vp;  cov_matrix[0][1] = vp_vs;   cov_matrix[0][2] = vp_rho;
  cov_matrix[1][0] = vp_vs;   cov_matrix[1][1] = var_vs;  cov_matrix[1][2] = vs_rho;
  cov_matrix[2][0] = vp_rho;  cov_matrix[2][1] = vs_rho;  cov_matrix[2][2] = var_rho; 

  return cov_matrix;
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

double** 
MultiNormalWithTrend::CreateCovMatrix(double s1, double s2) const {
  // fill cov matrix
  double ** cov_matrix = new double*[3];
  for (int i = 0; i < 3; ++i)
    cov_matrix[i] = new double[3];

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i)
      cov_matrix[j][i] = cov_trend_[j][i]->GetValue(s1, s2);
  }

  return cov_matrix;
}

void     
MultiNormalWithTrend::DeleteCovMatrix(double** cov_matrix) const {
  //clean up
  for (int i = 0; i < 3; ++i)
    delete [] cov_matrix[i];
  delete [] cov_matrix;

}

