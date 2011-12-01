#ifndef RPLIB_MULTINORMALWITHTREND_H
#define RPLIB_MULTINORMALWITHTREND_H

#include "nrlib/random/normal.hpp"

#include <vector>

class Trend;

class MultiNormalWithTrend {
 public:
   MultiNormalWithTrend(const NRLib::Normal&                      vp, 
                        const NRLib::Normal&                      vs, 
                        const NRLib::Normal&                      rho,
                        const Trend&                              mean_trend_vp,
                        const Trend&                              mean_trend_vs,
                        const Trend&                              mean_trend_rho,
                        const std::vector< std::vector<Trend*> >& cov_trend); //If error, throws NRLib::IndexOutOfRange

   virtual ~MultiNormalWithTrend();
   
   void     ReSample(double s1, double s2, 
                     double& vp, double& vs, double& rho) const; 

   void     ReSample(double s1, double s2,
                     double** cov_matrix_cholesky,
                     double& vp, double& vs, double& rho,
                     bool is_cholesky = true) const;

   void     EstimateExpectation(double s1, double s2, int sample_size,
                                double& exp_vp, double& exp_vs, double& exp_rho) const;

   void     EstimateExpectationAndVariance(double s1, double s2, int sample_size,
                                           double& exp_vp, double& exp_vs, double& exp_rho,
                                           double& var_vp, double& var_vs, double& var_rho) const;

   double** CreateEstimateOfCovMatrix(double s1, double s2, int sample_size) const;
 private:
   void     MatrProdTranspCholVecRR(int n, double **mat, double *in_vec,double* out_vec) const; //NBNB fjellvoll replace with library function?
   double** CreateCovMatrix(double s1, double s2) const;
   void     DeleteCovMatrix(double** cov_matrix) const;
  

 private:
  const NRLib::Normal&                      vp_;              // standardized distribution of vp 
  const NRLib::Normal&                      vs_;              // standardized distribution of vs given vp
  const NRLib::Normal&                      rho_;             // standardized distribution of rho given vp and vs
  const Trend&                              mean_trend_vp_;   // mean trend of vp
  const Trend&                              mean_trend_vs_;   // mean trend of vs
  const Trend&                              mean_trend_rho_;  // mean trend of rho
  const std::vector< std::vector<Trend*> >& cov_trend_;
   
};
#endif
