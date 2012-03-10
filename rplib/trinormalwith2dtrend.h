#ifndef RPLIB_TRINORMALWITH2DTREND_H
#define RPLIB_TRINORMALWITH2DTREND_H

#include "nrlib/random/normal.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/flens/nrlib_flens.hpp"
#include "nrlib/trend/trend.hpp"

#include <vector>

class Trend;

namespace NRLib {
  template<class A>
  class Grid2D;
}

class TriNormalWith2DTrend {
 public:

   TriNormalWith2DTrend(NRLib::Trend               *  mean_trend_vp,
                        NRLib::Trend               *  mean_trend_vs,
                        NRLib::Trend               *  mean_trend_rho,
                        NRLib::Grid2D<NRLib::Trend *> covariance_matrix);

   TriNormalWith2DTrend();

   virtual ~TriNormalWith2DTrend();

   void     ReSample(const NRLib::Normal & vp01,
                     const NRLib::Normal & vs01,
                     const NRLib::Normal & rho01,
                     const double        & s1,
                     const double        & s2,
                     double              & vp,
                     double              & vs,
                     double              & rho) const;

   void     ReSample(const NRLib::Normal & vp01,
                     const NRLib::Normal & vs01,
                     const NRLib::Normal & rho01,
                     double             ** cov_matrix_cholesky,
                     double              & E_vp,
                     double              & E_vs,
                     double              & E_rho,
                     double              & vp,
                     double              & vs,
                     double              & rho,
                     const bool          & is_cholesky = true) const;

   void     GetExpectation(const double        & s1,
                           const double        & s2,
                           std::vector<double> & expectation) const;

   void     GetCovariance(const double          & s1,
                          const double          & s2,
                          NRLib::Grid2D<double> & covariance) const;

   void     DebugEstimateExpectation(const NRLib::Normal & vp01,
                                     const NRLib::Normal & vs01,
                                     const NRLib::Normal & rho01,
                                     const double        & s1,
                                     const double        & s2,
                                     const int           & sample_size,
                                     double              & exp_vp,
                                     double              & exp_vs,
                                     double              & exp_rho) const;

   void     DebugEstimateExpectationAndVariance(const NRLib::Normal & vp01,
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
                                                double              & var_rho) const;

   void     CalculatePDF(const double & s1,
                         const double & s2,
                         const double & obs_vp,
                         const double & obs_vs,
                         const double & obs_rho,
                         double       & prob) const;

   double** DebugCreateEstimateOfCovMatrix(const NRLib::Normal & vp01,
                                           const NRLib::Normal & vs01,
                                           const NRLib::Normal & rho01,
                                           const double        & s1,
                                           const double        & s2,
                                           const int           & sample_size) const;
 private:
   void      MatrProdTranspCholVecRR(int       n,
                                     double ** mat,
                                     double  * in_vec,
                                     double  * out_vec) const; //NBNB fjellvoll replace with library function?

   double ** GetCovMatrix(const double & s1,
                          const double & s2) const;

   void      DeleteCovMatrix(double** cov_matrix) const;

   void      CalculateDeterminant(double ** cov_matrix,
                                  double  & determinant) const;

 private:
  NRLib::Trend               *  mean_trend_vp_;            // mean trend of vp
  NRLib::Trend               *  mean_trend_vs_;            // mean trend of vs
  NRLib::Trend               *  mean_trend_rho_;           // mean trend of rho
  NRLib::Grid2D<NRLib::Trend *> covariance_matrix_;

};
#endif
