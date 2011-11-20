#ifndef RPLIB_MULTINORMALWITHTREND_H
#define RPLIB_MULTINORMALWITHTREND_H

#include "nrlib/random/distribution.hpp"

#include <vector>

class Trend;

class MultiNormalWithTrend {
 public:
   MultiNormalWithTrend(const NRLib::Distribution<double>&        vp, 
                        const NRLib::Distribution<double>&        vs, 
                        const NRLib::Distribution<double>&        rho,
                        const Trend&                              mean_trend_vp,
                        const Trend&                              mean_trend_vs,
                        const Trend&                              mean_trend_rho,
                        const std::vector< std::vector<Trend*> >  cov_trend); //If error, throws NRLib::IndexOutOfRange

   virtual ~MultiNormalWithTrend();
   
   void     ReSample(double s1, double s2, 
                     double& vp, double& vs, double& rho) const; 
 private:
   void MatrProdTranspCholVecRR(int n, double **mat, double *in_vec,double* out_vec) const; //NBNB fjellvoll replace with library function?
  

 private:
  const NRLib::Distribution<double>&        vp_;              // standardized distribution of vp 
  const NRLib::Distribution<double>&        vs_;              // standardized distribution of vs given vp
  const NRLib::Distribution<double>&        rho_;             // standardized distribution of rho given vp and vs
  const Trend&                              mean_trend_vp_;   // mean trend of vp
  const Trend&                              mean_trend_vs_;   // mean trend of vs
  const Trend&                              mean_trend_rho_;  // mean trend of rho
  const std::vector< std::vector<Trend*> >& cov_trend_;
   
};
#endif
