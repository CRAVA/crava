#ifndef RPLIB_DEMMODELLING_H
#define RPLIB_DEMMODELLING_H

#include <vector>

namespace DEMTools {
  enum MixMethod{Hill, Reuss, Voigt};
  enum TabulatedMethod{Modulus, Velocity};

  //list of main functions
  double CalcBulkModulusOfBrineFromTPS(double temperature,
                                       double pressure,
                                       double salinity);

  double CalcDensityOfBrineFromTPS(double temperature,
                                   double pressure,
                                   double salinity);

  void   CalcCo2Prop(double& bulk_modulus,
                     double& density,
                     double ti,
                     double pi);

  double CalcEffectiveElasticModuliUsingHill(const std::vector<double> & prop,
                                             const std::vector<double> & volumefraction);

  double CalcEffectiveElasticModuliUsingReuss(const std::vector<double> & prop,
                                              const std::vector<double> & volumefraction);

  double CalcEffectiveElasticModuliUsingVoigt(const std::vector<double> & prop,
                                              const std::vector<double> & volumefraction);

  double CalcEffectiveDensity(const std::vector<double> & rho,
                              const std::vector<double> & volumefraction);

  double CalcEffectivePorosity(const std::vector<double> & porosity,
                               const std::vector<double> & volumefraction);

  void   CalcEffectiveBulkAndShearModulus(const std::vector<double>&       bulk_modulus,
                                          const std::vector<double>&       shear_modulus,
                                          const std::vector<double>&       aspect_ratio,
                                          std::vector<double>&             concentration,
                                          double                           bulk_modulus_bg,
                                          double                           shear_modulus_bg,
                                          double&                    effective_bulk_modulus,
                                          double&                    effective_shear_modulus);


  //list of helper functions called by the main functions
  double CalcVelocityOfBrineFromTPS(double temperature,
                                    double pressure,
                                    double salinity);

  double CalcVelocityOfWaterFromTP(double temperature,
                                   double pressure);

  double CalcDensityOfWaterFromTP(double temperature,
                                  double pressure);

  void CalcSeismicParamsFromElasticParams(const double & bulk_modulus,
                                          const double & shear_modulus,
                                          const double & density,
                                          double       & vp,
                                          double       & vs);

  void CalcElasticParamsFromSeismicParams(const double & vp,
                                          const double & vs,
                                          const double & density,
                                          double       & bulk_modulus,
                                          double       & shear_modulus);

  //list of debug/testing functions
  void DebugTestCalcEffectiveModulus(double& effective_bulk_modulus,
                                     double& effective_shear_modulus,
                                     double& effective_density);

  void DebugTestCalcEffectiveModulus2(double& effective_bulk_modulus,
                                      double& effective_shear_modulus,
                                      double& effective_density);

  // not working  because Distribution->DistributionWithTrend
  /*void DebugTestCalcEffectiveModulus3(double& effective_bulk_modulus,
                                      double& effective_shear_modulus,
                                      double& effective_density);*/

  void DebugTestCalcEffectiveModulus4(double& effective_bulk_modulus,
                                      double& effective_shear_modulus,
                                      double& effective_density);

  // not working because Distribution->DistributionWithTrend
  //void DebugTestDeletionAndCopying();

  void UpdateU(std::vector<double> & u,
               double                corr_param,
               bool                  param_is_time,
               std::vector<double>   alpha = std::vector<double>(0));
}
#endif

