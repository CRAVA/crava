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
   TriNormalWith2DTrend(NRLib::Trend * mean_trend_vp,
                        NRLib::Trend * mean_trend_vs,
                        NRLib::Trend * mean_trend_rho,
                        NRLib::Trend * variance_trend_vp,
                        NRLib::Trend * variance_trend_vs,
                        NRLib::Trend * variance_trend_rho,
                        NRLib::Trend * correlation_trend_vp_vs,
                        NRLib::Trend * correlation_trend_vp_rho,
                        NRLib::Trend * correlation_trend_vs_rho);

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

   void     CalculateExpectation(const double        & s1,
                                 const double        & s2,
                                 std::vector<double> & expectation) const;

   void     CalculateCovariance(const double  & s1,
                                const double  & s2,
                                NRLib::Matrix & covariance) const;

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
                         float        & prob) const;

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

   double ** CalculateCovMatrix(const double & s1,
                                const double & s2) const;

   void      DeleteCovMatrix(double** cov_matrix) const;

   void      CalculateDeterminant(double ** cov_matrix,
                                  double  & determinant) const;

   void      CalculateExpectation(double       & E_vp,
                                  double       & E_vs,
                                  double       & E_rho,
                                  const double & s1,
                                  const double & s2) const;

   void      CalculateExpectation(double       & E_vp,
                                  double       & E_vs,
                                  double       & E_rho,
                                  const double & s1,
                                  const double & s2,
                                  double      ** cov_matrix,
                                  const bool   & is_cholesky = false) const;

 private:
  NRLib::Trend *      mean_trend_vp_;            // mean trend of vp
  NRLib::Trend *      mean_trend_vs_;            // mean trend of vs
  NRLib::Trend *      mean_trend_rho_;           // mean trend of rho
  NRLib::Trend *      variance_trend_vp_;        // variance trend of vp
  NRLib::Trend *      variance_trend_vs_;        // variance trend of vs
  NRLib::Trend *      variance_trend_rho_;       // variance trend of rho
  NRLib::Trend *      correlation_trend_vp_vs_;  // correlation trend of vp,vs
  NRLib::Trend *      correlation_trend_vp_rho_; // corraletion trend of vp,rho
  NRLib::Trend *      correlation_trend_vs_rho_; // correlation trend of vs,rho

};
#endif
