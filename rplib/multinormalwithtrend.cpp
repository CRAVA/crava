#include "rplib/multinormalwithtrend.h"

#include "nrlib/trend/trend.hpp"
#include "nrlib/exception/exception.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "lib/lib_matr.h"

#include <cmath>




MultiNormalWithTrend::
MultiNormalWithTrend(const NRLib::Normal&                      vp01, //Marit: Hva brukes de tre første variablene til?
                     const NRLib::Normal&                      vs01,
                     const NRLib::Normal&                      rho01,
                     const Trend*                              mean_trend_vp,
                     const Trend*                              mean_trend_vs,
                     const Trend*                              mean_trend_rho,
                     const NRLib::Grid2D<Trend*>               cov_trend) : //Marit: Vil ta inn verdiene direkte uten grid
  vp01_(vp01), vs01_(vs01), rho01_(rho01),
  mean_trend_vp_(mean_trend_vp), mean_trend_vs_(mean_trend_vs), mean_trend_rho_(mean_trend_rho),
  cov_trend_(cov_trend)
{
  if (cov_trend_.GetNI() != 3 || cov_trend_.GetNJ() != 3)
    throw NRLib::IndexOutOfRange("MultiNormalWithTrend: Size of covariance trend matrix is not valid.");

}

MultiNormalWithTrend::~MultiNormalWithTrend()
{

}

void
MultiNormalWithTrend::ReSample(double s1, double s2,
                               double& vp, double& vs, double& rho) const {
  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp01_.Draw();
  rc[1] = vs01_.Draw();
  rc[2] = rho01_.Draw();

  // fill cov matrix
  double ** cov_matrix = CreateCovMatrix(s1, s2);

  //cholesky
  lib_matrCholR(3, cov_matrix);

  MatrProdTranspCholVecRR(3, cov_matrix, &rc[0] , &rhs[0]); // rhs is given value

  vp  =  rhs[0] + mean_trend_vp_->GetValue(s1, s2);
  vs  =  rhs[1] + mean_trend_vs_->GetValue(s1, s2);
  rho =  rhs[2] + mean_trend_rho_->GetValue(s1, s2);

  DeleteCovMatrix(cov_matrix);

}

void
MultiNormalWithTrend::ReSample(double s1, double s2,
                               double** cov_matrix_cholesky,
                               double& vp, double& vs, double& rho,
                               bool is_cholesky) const {

  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp01_.Draw();
  rc[1] = vs01_.Draw();
  rc[2] = rho01_.Draw();

  //cholesky
  if (!is_cholesky)
    lib_matrCholR(3, cov_matrix_cholesky);

  MatrProdTranspCholVecRR(3, cov_matrix_cholesky, &rc[0] , &rhs[0]); // rhs is given value

  vp  =  rhs[0] + mean_trend_vp_->GetValue(s1, s2);
  vs  =  rhs[1] + mean_trend_vs_->GetValue(s1, s2);
  rho =  rhs[2] + mean_trend_rho_->GetValue(s1, s2);

}

void
MultiNormalWithTrend::GetExpectation(double s1, double s2, std::vector<double>& expectation) const
{
  expectation[0] = mean_trend_vp_->GetValue(s1, s2);
  expectation[1] = mean_trend_vs_->GetValue(s1, s2);
  expectation[2] = mean_trend_rho_->GetValue(s1, s2);
}

void
MultiNormalWithTrend::GetCovariance(double s1, double s2, NRLib::Matrix & covariance) const
{
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i)
      covariance(j, i) = cov_trend_(i, j)->GetValue(s1, s2);
  }
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


void
MultiNormalWithTrend::CalculatePDF(double s1, double s2, double obs_vp, double obs_vs, double obs_rho, float& prob) const
{
  double E_vp  = mean_trend_vp_ ->GetValue(s1,s2);
  double E_vs  = mean_trend_vs_ ->GetValue(s1,s2);
  double E_rho = mean_trend_rho_->GetValue(s1,s2);

  int dim = 3;

  double ** cov_matrix = CreateCovMatrix(s1, s2);

  cov_matrix[0][0] = std::log(1+cov_matrix[0][0]/(E_vp*E_vp));
  cov_matrix[1][1] = std::log(1+cov_matrix[1][1]/(E_vs*E_vs));
  cov_matrix[2][2] = std::log(1+cov_matrix[2][2]/(E_rho*E_rho));

  for(int i=0; i<2; i++) {
    for(int j=i+1; j<3; j++) {
      cov_matrix[i][j] = cov_matrix[i][j]*std::sqrt(cov_matrix[i][i]*cov_matrix[j][j]);
      cov_matrix[j][i] = cov_matrix[i][j];
    }
  }

  E_vp  = std::log(E_vp)  - 0.5*cov_matrix[0][0];
  E_vs  = std::log(E_vs)  - 0.5*cov_matrix[1][1];
  E_rho = std::log(E_rho) - 0.5*cov_matrix[2][2];

  NRLib::Vector diff(dim);
  diff(0) = obs_vp  - E_vp;
  diff(1) = obs_vs  - E_vs;
  diff(2) = obs_rho - E_rho;

  NRLib::Matrix inv_cov_mat(dim,dim);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      inv_cov_mat(i,j) = cov_matrix[i][j];
  NRLib::invert(inv_cov_mat);

  double determinant;
  CalculateDeterminant(cov_matrix, determinant);
  DeleteCovMatrix(cov_matrix);

  NRLib::Vector b(dim);
  b = inv_cov_mat * diff;

  double prod = 0;
  for(int i=0; i<dim; i++)
    prod += diff(i)*b(i);

  double aa = std::exp(-0.5*prod);
  double bb = std::pow(2*NRLib::Pi,1.5)*std::pow(determinant,0.5);

  prob = static_cast<float>(std::exp(std::log(aa)-std::log(bb)));
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

void MultiNormalWithTrend::
CalculateDeterminant(double** cov_matrix, double& determinant) const
{
  determinant = cov_matrix[0][0]*(cov_matrix[2][2]*cov_matrix[1][1]-cov_matrix[2][1]*cov_matrix[1][2])
               -cov_matrix[1][0]*(cov_matrix[2][2]*cov_matrix[0][1]-cov_matrix[2][1]*cov_matrix[0][2])
               +cov_matrix[2][0]*(cov_matrix[1][2]*cov_matrix[0][1]-cov_matrix[1][1]*cov_matrix[0][2]);
}

double**
MultiNormalWithTrend::CreateCovMatrix(double s1, double s2) const {
  // fill cov matrix
  double ** cov_matrix = new double*[3];
  for (int i = 0; i < 3; ++i)
    cov_matrix[i] = new double[3];

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i)
      cov_matrix[j][i] = cov_trend_(i, j)->GetValue(s1, s2);
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

