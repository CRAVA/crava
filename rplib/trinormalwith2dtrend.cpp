#include "rplib/trinormalwith2dtrend.h"

#include "nrlib/trend/trend.hpp"
#include "nrlib/exception/exception.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

#include "lib/lib_matr.h"

#include <cmath>


TriNormalWith2DTrend::TriNormalWith2DTrend(NRLib::Trend               *  mean_trend_vp,
                                           NRLib::Trend               *  mean_trend_vs,
                                           NRLib::Trend               *  mean_trend_rho,
                                           NRLib::Grid2D<NRLib::Trend *> covariance_matrix)
: mean_trend_vp_(mean_trend_vp),
  mean_trend_vs_(mean_trend_vs),
  mean_trend_rho_(mean_trend_rho),
  covariance_matrix_(covariance_matrix)
{
}

TriNormalWith2DTrend::~TriNormalWith2DTrend()
{

}

void
TriNormalWith2DTrend::ReSample(const NRLib::Normal & vp01,
                               const NRLib::Normal & vs01,
                               const NRLib::Normal & rho01,
                               const double        & s1,
                               const double        & s2,
                               double              & vp,
                               double              & vs,
                               double              & rho) const
{
  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp01.Draw();
  rc[1] = vs01.Draw();
  rc[2] = rho01.Draw();

  std::vector<double> expectation(3);
  GetExpectation(s1,s2,expectation);

  // fill cov matrix
  double ** cov_matrix = GetCovMatrix(s1, s2);

  //cholesky
  lib_matrCholR(3, cov_matrix);

  MatrProdTranspCholVecRR(3, cov_matrix, &rc[0] , &rhs[0]); // rhs is given value
  DeleteCovMatrix(cov_matrix);

  vp  =  rhs[0] + expectation[0];
  vs  =  rhs[1] + expectation[1];
  rho =  rhs[2] + expectation[2];

}

void
TriNormalWith2DTrend::ReSample(const NRLib::Normal & vp01,
                               const NRLib::Normal & vs01,
                               const NRLib::Normal & rho01,
                               double             ** cov_matrix_cholesky,
                               double              & E_vp,
                               double              & E_vs,
                               double              & E_rho,
                               double              & vp,
                               double              & vs,
                               double              & rho,
                               const bool          & is_cholesky) const {

  std::vector<double> rhs(3); // right hand side
  std::vector<double> rc(3);  // random component

  rc[0] = vp01.Draw();
  rc[1] = vs01.Draw();
  rc[2] = rho01.Draw();

  //cholesky
  if (!is_cholesky)
    lib_matrCholR(3, cov_matrix_cholesky);

  MatrProdTranspCholVecRR(3, cov_matrix_cholesky, &rc[0] , &rhs[0]); // rhs is given value

  vp  =  rhs[0] + E_vp;
  vs  =  rhs[1] + E_vs;
  rho =  rhs[2] + E_rho;

}

void
TriNormalWith2DTrend::GetExpectation(const double        & s1,
                                     const double        & s2,
                                     std::vector<double> & expectation) const
{
  double dummy = 0;

  expectation[0] = mean_trend_vp_ ->GetValue(s1,s2,dummy);
  expectation[1] = mean_trend_vs_ ->GetValue(s1,s2,dummy);
  expectation[2] = mean_trend_rho_->GetValue(s1,s2,dummy);

}

void
TriNormalWith2DTrend::GetCovariance(const double          & s1,
                                    const double          & s2,
                                    NRLib::Grid2D<double> & covariance) const
{
  double dummy = 0;

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      covariance(i,j) = covariance_matrix_(i,j)->GetValue(s1,s2,dummy);
    }
  }
}

void
TriNormalWith2DTrend::DebugEstimateExpectation(const NRLib::Normal & vp01,
                                               const NRLib::Normal & vs01,
                                               const NRLib::Normal & rho01,
                                               const double        & s1,
                                               const double        & s2,
                                               const int           & sample_size,
                                               double              & exp_vp,
                                               double              & exp_vs,
                                               double              & exp_rho) const
{

  exp_vp = exp_vs = exp_rho = 0;

  if (sample_size <= 0)
    return;

  std::vector<double> expectation(3);
  GetExpectation(s1,s2,expectation);

  // fill cov matrix
  double ** cov_matrix = GetCovMatrix(s1, s2);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  for (int i = 0; i < sample_size; ++i) {
    double vp, vs, rho;
    ReSample(vp01,vs01, rho01,
             cov_matrix,
             expectation[0],
             expectation[1],
             expectation[2],
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
TriNormalWith2DTrend::DebugEstimateExpectationAndVariance(const NRLib::Normal & vp01,
                                                          const NRLib::Normal & vs01,
                                                          const NRLib::Normal & rho01,
                                                          const double        & s1,
                                                          const double        & s2,
                                                          const int           & sample_size,
                                                          double              & exp_vp,
                                                          double              & exp_vs,
                                                          double              & exp_rho,
                                                          double              & var_vp,
                                                          double              & var_vs,
                                                          double              & var_rho) const {

  exp_vp = exp_vs = exp_rho = 0;
  var_vp = var_vs = var_rho = 0;

  if (sample_size <= 0)
    return;

  std::vector<double> expectation(3);
  GetExpectation(s1,s2,expectation);

  // fill cov matrix
  double ** cov_matrix = GetCovMatrix(s1, s2);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  for (int i = 0; i < sample_size; ++i) {
    double vp, vs, rho;
    ReSample(vp01, vs01, rho01,
             cov_matrix,
             expectation[0],
             expectation[1],
             expectation[2],
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
TriNormalWith2DTrend::CalculatePDF(const double & s1,
                                   const double & s2,
                                   const double & obs_vp,
                                   const double & obs_vs,
                                   const double & obs_rho,
                                   double       & prob) const
{
  int dim = 3;

  std::vector<double> expectation(3);
  GetExpectation(s1,s2,expectation);

  NRLib::Vector diff(dim);
  diff(0) = obs_vp  - expectation[0];
  diff(1) = obs_vs  - expectation[1];
  diff(2) = obs_rho - expectation[2];

  double ** cov_matrix = GetCovMatrix(s1, s2);

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

  prob = std::exp(std::log(aa)-std::log(bb));
}

double**
TriNormalWith2DTrend::DebugCreateEstimateOfCovMatrix(const NRLib::Normal & vp01,
                                                     const NRLib::Normal & vs01,
                                                     const NRLib::Normal & rho01,
                                                     const double        & s1,
                                                     const double        & s2,
                                                     const int           & sample_size) const {
  if (sample_size <= 0)
    return NULL;

  double exp_vp, exp_vs, exp_rho;
  double var_vp, var_vs, var_rho;

  exp_vp = exp_vs = exp_rho = 0;
  var_vp = var_vs = var_rho = 0;


  // fill cov matrix
  double ** cov_matrix = GetCovMatrix(s1, s2);

  std::vector<double> expectation(3);
  GetExpectation(s1,s2,expectation);

  // cholesky
  lib_matrCholR(3, cov_matrix);

  double vp_vs   = 0;
  double vp_rho  = 0;
  double vs_rho  = 0;

  for (int i = 0; i < sample_size; ++i) {
    double vp = 0, vs = 0, rho = 0;
    ReSample(vp01, vs01, rho01,
             cov_matrix,
             expectation[0],
             expectation[1],
             expectation[2],
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
TriNormalWith2DTrend::MatrProdTranspCholVecRR(int       n,
                                             double ** mat,
                                             double  * in_vec,
                                             double  * out_vec) const
{
  int i,j;
  for (i = 0; i < n; i++) {
    out_vec[i]=0.0;
    for (j = 0; j <= i; j++)
      out_vec[i] += mat[i][j]*in_vec[j];
  }
}

void TriNormalWith2DTrend::CalculateDeterminant(double ** cov_matrix,
                                                double  & determinant) const
{
  determinant = cov_matrix[0][0]*(cov_matrix[2][2]*cov_matrix[1][1]-cov_matrix[2][1]*cov_matrix[1][2])
               -cov_matrix[1][0]*(cov_matrix[2][2]*cov_matrix[0][1]-cov_matrix[2][1]*cov_matrix[0][2])
               +cov_matrix[2][0]*(cov_matrix[1][2]*cov_matrix[0][1]-cov_matrix[1][1]*cov_matrix[0][2]);
}

double**
TriNormalWith2DTrend::GetCovMatrix(const double & s1,
                                   const double & s2) const
{
  NRLib::Grid2D<double> covariance(3,3,0);
  GetCovariance(s1,s2,covariance);

  double ** cov_matrix = new double*[3];
  for (int i = 0; i < 3; ++i)
    cov_matrix[i] = new double[3];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      cov_matrix[i][j] = covariance(i,j);

  return cov_matrix;
}

void
TriNormalWith2DTrend::DeleteCovMatrix(double ** cov_matrix) const {
  //clean up
  for (int i = 0; i < 3; ++i)
    delete [] cov_matrix[i];
  delete [] cov_matrix;
}
