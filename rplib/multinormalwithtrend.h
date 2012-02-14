#ifndef RPLIB_MULTINORMALWITHTREND_H
#define RPLIB_MULTINORMALWITHTREND_H

#include "nrlib/random/normal.hpp"
#include "nrlib/grid/grid2d.hpp"

#include <vector>
#include <nrlib/flens/nrlib_flens.hpp>

#include <vector>

class Trend;

namespace NRLib {
  template<class A>
  class Grid2D;
}

class MultiNormalWithTrend {
 public:
   MultiNormalWithTrend(const NRLib::Normal&                      vp01,
                        const NRLib::Normal&                      vs01,
                        const NRLib::Normal&                      rho01,
                        const Trend*                              mean_trend_vp,
                        const Trend*                              mean_trend_vs,
                        const Trend*                              mean_trend_rho,
                        const NRLib::Grid2D<Trend*>               cov_trend); //If error, throws NRLib::IndexOutOfRange
                        //const Trend*                              variance_trend_vp,
                        //const Trend*                              variance_trend_vs,
                        //const Trend*                              variance_trend_rho,
                        //const Trend*                              correlation_trend_vp_vs,
                        //const Trend*                              correlation_trend_vp_rho,
                        //const Trend*                              correlation_trend_vs_rho); //If error, throws NRLib::IndexOutOfRange

   virtual ~MultiNormalWithTrend();

   void     ReSample(double s1, double s2,
                     double& vp, double& vs, double& rho) const;

   void     ReSample(double s1, double s2,
                     double** cov_matrix_cholesky,
                     double& vp, double& vs, double& rho,
                     bool is_cholesky = true) const;

   void     GetExpectation(double s1, double s2, std::vector<double>& expectation) const;

   void     GetCovariance(double s1, double s2, NRLib::Matrix & covariance) const;

   void     EstimateExpectation(double s1, double s2, int sample_size,
                                double& exp_vp, double& exp_vs, double& exp_rho) const;

   void     EstimateExpectationAndVariance(double s1, double s2, int sample_size,
                                           double& exp_vp, double& exp_vs, double& exp_rho,
                                           double& var_vp, double& var_vs, double& var_rho) const;

   void     CalculatePDF(double s1, double s2, double obs_vp, double obs_vs, double obs_rho, float& prob) const;

   double** CreateEstimateOfCovMatrix(double s1, double s2, int sample_size) const;
 private:
   void     MatrProdTranspCholVecRR(int n, double **mat, double *in_vec,double* out_vec) const; //NBNB fjellvoll replace with library function?
   double** CreateCovMatrix(double s1, double s2) const;
   void     DeleteCovMatrix(double** cov_matrix) const;
   void     CalculateDeterminant(double** cov_matrix, double& determinant) const;


 private:
  const NRLib::Normal&                      vp01_;            // standardized distribution of vp
  const NRLib::Normal&                      vs01_;            // standardized distribution of vs given vp
  const NRLib::Normal&                      rho01_;           // standardized distribution of rho given vp and vs
  const Trend*                              mean_trend_vp_;   // mean trend of vp
  const Trend*                              mean_trend_vs_;   // mean trend of vs
  const Trend*                              mean_trend_rho_;  // mean trend of rho
  const NRLib::Grid2D<Trend*>               cov_trend_;

};
#endif
