#ifndef RPLIB_DEMMODELLING_H
#define RPLIB_DEMMODELLING_H

#include <vector>

namespace DEMTools {
  //list of main functions
  double CalcBulkModulusOfBrineFromTPS(double temperature,
                                       double pressure,
                                       double salinity);

  double CalcDensityOfBrineFromTPS(double temperature,
                                   double pressure,
                                   double salinity);

  void   CalcCo2Prop(double& bulk_modulus,
                     double& density,
                     double& velocity,
                     double ti,
                     double pi);

  double CalcEffectiveElasticModuliUsingHill(double property1,
                                             double volumefraction1,
                                             double property2,
                                             double volumefraction2 = -1.0);

  double CalcEffectiveElasticModuliUsingReuss(double property1,
                                              double volumefraction1,
                                              double property2,
                                              double volumefraction2 = -1.0);

  double CalcEffectiveElasticModuliUsingVoigt(double property1,
                                              double volumefraction1,
                                              double property2,
                                              double volumefraction2 = -1.0);

  double CalcEffectiveDensity(double rho1,
                              double volumefraction1,
                              double rho2,
                              double volumefraction2 = -1.0);

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

  //list of debug/testing functions
  void DebugTestCalcEffectiveModulus(double& effective_bulk_modulus,
                                     double& effective_shear_modulus,
                                     double& effective_density);




}
#endif

