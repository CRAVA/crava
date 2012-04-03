#include "rplib/demmodelling.h"

#include "rplib/dem.h"
#include "rplib/mineral.h"
#include "rplib/brine.h"
#include "rplib/co2.h"
#include "rplib/solidmixed.h"
#include "rplib/fluidmixed.h"
#include "rplib/rockinclusion.h"

#include "nrlib/exception/exception.hpp"

#include <cmath>

double
DEMTools::CalcBulkModulusOfBrineFromTPS(double temperature,
                                        double pressure,
                                        double salinity) {


  double vb   = CalcVelocityOfBrineFromTPS(temperature,
                                           pressure,
                                           salinity);

  double rhob = CalcDensityOfBrineFromTPS(temperature,
                                          pressure,
                                          salinity);

  return rhob*(vb*vb)/1E3;
}

double
DEMTools::CalcDensityOfBrineFromTPS(double temperature,
                                    double pressure,
                                    double salinity) {

  return CalcDensityOfWaterFromTP(temperature,pressure) +
        salinity*(0.668 + 0.44*salinity + 1E-6*(300*pressure - 2400*pressure*salinity +
        temperature*(80 + 3*temperature - 3300*salinity - 13*pressure + 47*pressure*salinity)));

}

void
DEMTools::CalcCo2Prop(double& bulk_modulus,
                      double& density,
                      double ti,
                      double pi) {

  double ttmp[] = {1.7000000e+001,  2.7000000e+001,  3.7000000e+001,  4.7000000e+001,  5.7000000e+001,  6.7000000e+001,  7.7000000e+001,  8.7000000e+001,  9.7000000e+001,  1.0700000e+002,  1.1700000e+002,  1.2700000e+002};
  size_t nt     = 12;

  double pmpatmp[] = {1.0005000e-001,  1.0005000e+000,  4.0020000e+000,  7.0035000e+000,  1.0005000e+001,  1.4007000e+001,  2.0010000e+001,  2.5012500e+001,  3.0015000e+001,  4.0020000e+001};
  size_t np        = 10;

  std::vector<double> t(ttmp, ttmp + nt);
  std::vector<double> pmpa(pmpatmp, pmpatmp + np);

  static std::vector< std::vector<double> > co2_bulk(np, std::vector<double>(nt, 0.0));
  /*static std::vector< std::vector<double> > co2_velocity(np, std::vector<double>(nt, 0.0));*/ //Velocity interpolation commented out
  static std::vector< std::vector<double> > co2_density(np, std::vector<double>(nt, 0.0));

  //{ // local scope co2 velocity
  //  double tmp0[] = {2.6400000e+002,  2.6800000e+002,  2.7200000e+002,  2.7600000e+002,  2.8000000e+002,  2.8400000e+002,  2.8800000e+002,  2.9100000e+002,  2.9500000e+002,  2.9900000e+002,  3.0200000e+002,  3.0600000e+002};
  //  co2_velocity[0] = std::vector<double>(tmp0, tmp0 + nt);

  //  double tmp1[] = {2.5600000e+002,  2.6100000e+002,  2.6600000e+002,  2.7000000e+002,  2.7500000e+002,  2.7900000e+002,  2.8400000e+002,  2.8800000e+002,  2.9200000e+002,  2.9600000e+002,  3.0000000e+002,  3.0400000e+002};
  //  co2_velocity[1] = std::vector<double>(tmp1, tmp1 + nt);

  //  double tmp2[] = {2.2000000e+002,  2.3100000e+002,  2.4100000e+002,  2.5100000e+002,  2.5900000e+002,  2.6400000e+002,  2.7100000e+002,  2.7600000e+002,  2.8100000e+002,  2.8700000e+002,  2.9200000e+002,  2.9600000e+002};
  //  co2_velocity[2] = std::vector<double>(tmp2, tmp2 + nt);

  //  double tmp3[] = {3.7900000e+002,  1.8500000e+002,  2.0900000e+002,  2.1900000e+002,  2.3700000e+002,  2.5200000e+002,  2.6300000e+002,  2.7000000e+002,  2.7600000e+002,  2.8200000e+002,  2.8500000e+002,  2.9000000e+002};
  //  co2_velocity[3] = std::vector<double>(tmp3, tmp3 + nt);

  //  double tmp4[] = {4.8300000e+002,  4.1100000e+002,  1.4200000e+002,  1.8300000e+002,  2.1200000e+002,  2.3800000e+002,  2.4900000e+002,  2.5700000e+002,  2.6500000e+002,  2.7200000e+002, 2.7800000e+002,  2.8600000e+002};
  //  co2_velocity[4] = std::vector<double>(tmp4, tmp4 + nt);

  //  double tmp5[] = {5.5400000e+002,  4.9600000e+002,  3.6600000e+002,  2.6900000e+002,  2.4700000e+002,  2.4800000e+002,  2.5500000e+002,  2.6000000e+002,  2.6900000e+002,  2.7800000e+002,  2.8500000e+002,  2.9100000e+002};
  //  co2_velocity[5] = std::vector<double>(tmp5, tmp5 + nt);

  //  double tmp6[] = {6.3700000e+002,  5.8600000e+002, 5.0200000e+002,  4.4100000e+002,  3.8200000e+002,  3.5800000e+002,  3.5300000e+002,  3.4500000e+002,  3.3700000e+002,  3.2800000e+002,  3.2000000e+002,  3.1400000e+002};
  //  co2_velocity[6] = std::vector<double>(tmp6, tmp6 + nt);

  //  double tmp7[] = { 6.8800000e+002,  6.4000000e+002,  5.9100000e+002,  5.3800000e+002,  4.8800000e+002,  4.6000000e+002,  4.3100000e+002,  4.1300000e+002,  3.9700000e+002,  3.8200000e+002,  3.6600000e+002,  3.5200000e+002};
  //  co2_velocity[7] = std::vector<double>(tmp7, tmp7 + nt);

  //  double tmp8[] = {7.3700000e+002,  6.8700000e+002,  6.4100000e+002,  5.9600000e+002,  5.5200000e+002,  5.2100000e+002,  4.9400000e+002,  4.7000000e+002,  4.5000000e+002,  4.2900000e+002,  4.1200000e+002,  3.9700000e+002};
  //  co2_velocity[8] = std::vector<double>(tmp8, tmp8 + nt);

  //  double tmp9[] = { 8.1400000e+002,  7.6600000e+002,  7.1400000e+002,  6.7600000e+002,  6.3900000e+002,  6.1500000e+002,  5.9200000e+002,  5.6400000e+002,  5.4500000e+002,  5.2600000e+002,  5.0400000e+002,  4.8700000e+002};
  //  co2_velocity[9] = std::vector<double>(tmp9, tmp9 + nt);

  //} // local scope end co2 velocity

  { // local scope co2 density
    double tmp0[] = {1.8600000e-003,  1.8000000e-003,  1.7400000e-003,  1.6800000e-003,  1.6300000e-003,  1.5800000e-003,  1.5400000e-003,  1.4900000e-003,  1.4500000e-003,  1.4100000e-003,  1.3800000e-003,  1.3400000e-003};
    co2_density[0] = std::vector<double>(tmp0, tmp0 + nt);

    double tmp1[] = {1.9630000e-002,  1.8840000e-002,  1.8130000e-002,  1.7480000e-002,  1.6880000e-002,  1.6320000e-002,  1.5800000e-002,  1.5310000e-002,  1.4860000e-002,  1.4440000e-002,  1.4040000e-002,  1.3660000e-002};
    co2_density[1] = std::vector<double>(tmp1, tmp1 + nt);

    double tmp2[] = {1.2900000e-001,  9.3950000e-002,  8.7090000e-002,  8.1690000e-002,  7.7240000e-002,  7.3450000e-002,  7.0130000e-002,  6.7190000e-002,  6.4540000e-002,  6.2150000e-002,  6.0000000e-002,  5.7970000e-002};
    co2_density[2] = std::vector<double>(tmp2, tmp2 + nt);

    double tmp3[] = {8.3000000e-001,  6.8000000e-001,  2.0072000e-001,  1.8283000e-001,  1.6362000e-001,  1.4960000e-001,  1.3928000e-001,  1.3094000e-001,  1.2391000e-001,  1.1786000e-001,  1.1275000e-001,  1.0794000e-001};
    co2_density[3] = std::vector<double>(tmp3, tmp3 + nt);

    double tmp4[] = {9.1700000e-001,  8.0500000e-001,  6.8300000e-001,  4.4973000e-001,  3.2761000e-001,  2.6715000e-001,  2.3550000e-001,  2.1381000e-001,  1.9732000e-001,  1.8422000e-001,  1.7409000e-001,  1.6443000e-001};
    co2_density[4] = std::vector<double>(tmp4, tmp4 + nt);

    double tmp5[] = {9.3000000e-001,  8.6000000e-001,  7.8000000e-001,  6.9000000e-001,  5.7000000e-001,  4.8000000e-001,  3.9000000e-001,  3.3500000e-001,  3.0000000e-001,  2.7500000e-001,  2.5200000e-001,  2.3500000e-001};
    co2_density[5] = std::vector<double>(tmp5, tmp5 + nt);

    double tmp6[] = { 9.6000000e-001,  9.1000000e-001,  8.6000000e-001,  8.0000000e-001,  7.4000000e-001,  6.8000000e-001,  6.2000000e-001,  5.6000000e-001,  5.1000000e-001,  4.7000000e-001,  4.2000000e-001,  3.9000000e-001};
    co2_density[6] = std::vector<double>(tmp6, tmp6 + nt);

    double tmp7[] = {9.9000000e-001,  9.5000000e-001,  9.0000000e-001,  8.5000000e-001,  7.9000000e-001,  7.5000000e-001,  7.0000000e-001,  6.4000000e-001,  5.9000000e-001,  5.5000000e-001,  5.1000000e-001,  4.8000000e-001};
    co2_density[7] = std::vector<double>(tmp7, tmp7 + nt);

    double tmp8[] = {1.0080000e+000,  9.7000000e-001,  9.3000000e-001,  8.9000000e-001,  8.5000000e-001,  8.1000000e-001,  7.7000000e-001,  7.2000000e-001,  6.8000000e-001,  6.3000000e-001,  6.0000000e-001,  5.7000000e-001};
    co2_density[8] = std::vector<double>(tmp8, tmp8 + nt);

    double tmp9[] = {1.0400000e+000,  1.0000000e+000,  9.7000000e-001,  9.4000000e-001,  9.0000000e-001,  8.7000000e-001,  8.4000000e-001,  8.0000000e-001,  7.7000000e-001,  7.3000000e-001,  7.0000000e-001,  6.7000000e-001};
    co2_density[9] = std::vector<double>(tmp9, tmp9 + nt);

  } // local scope end co2 density

  { // local scope for co2 bulk
    double tmp0[] = {1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004,  1.3000000e-004 };
    co2_bulk[0] = std::vector<double>(tmp0, tmp0 + nt);

    double tmp1[] = {1.2900000e-003,  1.2800000e-003,  1.2800000e-003,  1.2800000e-003,  1.2800000e-003,  1.2700000e-003,  1.2700000e-003,  1.2700000e-003,  1.2600000e-003,  1.2600000e-003,  1.2600000e-003,  1.2600000e-003};
    co2_bulk[1] = std::vector<double>(tmp1, tmp1 + nt);

    double tmp2[] = {6.2400000e-003,  5.0100000e-003,  5.0400000e-003,  5.1300000e-003,  5.1800000e-003,  5.1400000e-003,  5.1700000e-003,  5.1000000e-003,  5.1000000e-003,  5.1100000e-003,  5.1100000e-003,  5.0700000e-003};
    co2_bulk[2] = std::vector<double>(tmp2, tmp2 + nt);

    double tmp3[] = {1.1922000e-001,  2.3270000e-002,  8.7300000e-003,  8.7900000e-003,  9.2300000e-003,  9.5100000e-003,  9.6500000e-003,  9.5500000e-003,  9.4100000e-003,  9.3700000e-003,  9.1600000e-003,  9.1000000e-003};
    co2_bulk[3] = std::vector<double>(tmp3, tmp3 + nt);

    double tmp4[] = {2.1393000e-001,  1.3598000e-001,  1.5990000e-002,  1.2820000e-002,  1.4690000e-002,  1.5200000e-002,  1.4600000e-002,  1.4120000e-002,  1.3860000e-002,  1.3630000e-002,  1.3460000e-002,  1.3410000e-002};
    co2_bulk[4] = std::vector<double>(tmp4, tmp4 + nt);

    double tmp5[] = {2.8543000e-001,  2.1157000e-001,  1.0449000e-001,  4.9930000e-002,  3.4780000e-002,  2.9520000e-002,  2.5360000e-002,  2.2650000e-002,  2.1710000e-002,  2.1250000e-002,  2.0470000e-002,  1.9900000e-002};
    co2_bulk[5] = std::vector<double>(tmp5, tmp5 + nt);

    double tmp6[] = {3.8954000e-001,  3.1249000e-001,  2.1672000e-001,  1.5558000e-001,  1.0798000e-001,  8.7150000e-002,  7.7260000e-002,  6.6650000e-002,  5.7920000e-002,  5.0560000e-002,  4.3010000e-002,  3.8450000e-002};
    co2_bulk[6] = std::vector<double>(tmp6, tmp6 + nt);

    double tmp7[] = {4.6861000e-001,  3.8912000e-001,  3.1435000e-001,  2.4603000e-001,  1.8813000e-001,  1.5870000e-001,  1.3003000e-001,  1.0916000e-001,  9.2990000e-002,  8.0260000e-002,  6.8320000e-002,  5.9470000e-002};
    co2_bulk[7] = std::vector<double>(tmp7, tmp7 + nt);

    double tmp8[] = {5.4751000e-001,  4.5781000e-001,  3.8212000e-001,  3.1614000e-001,  2.5900000e-001,  2.1987000e-001,  1.8791000e-001,  1.5905000e-001,  1.3770000e-001,  1.1595000e-001,  1.0185000e-001,  8.9840000e-002};
    co2_bulk[8] = std::vector<double>(tmp8, tmp8 + nt);

    double tmp9[] = {6.8910000e-001,  5.8676000e-001,  4.9450000e-001,  4.2956000e-001,  3.6749000e-001,  3.2906000e-001,  2.9439000e-001,  2.5448000e-001,  2.2871000e-001,  2.0197000e-001,  1.7781000e-001,  1.5890000e-001};
    co2_bulk[9] = std::vector<double>(tmp9, tmp9 + nt);

  } // local scope end


  if (ti < t.front() || ti > t.back() || pi < pmpa.front() || pi > pmpa.back())
    throw NRLib::Exception("Temperature or pressure outside valid range.");
  size_t i1 = 0;
  while (i1 < t.size() && ti >= t[i1])
    i1++;

  if (i1 >= 1)
    i1--;

  if (i1 > t.size() - 2)
    i1 = t.size() - 2;

  double t1 = (ti - t[i1])/(t[i1+1] - t[i1]);

  size_t j1 = 0;
  while (j1 < pmpa.size() && pi >= pmpa[j1])
    j1++;

  if (j1 >= 1)
    j1--;

  if (j1 > pmpa.size() - 2)
    j1 = pmpa.size() - 2;

  double t2 = (pi - pmpa[j1])/(pmpa[j1+1] - pmpa[j1]);

  double val_s1_1 = t1*co2_bulk[j1][i1+1] + (1.0 - t1)*co2_bulk[j1][i1];
  double val_s1_2 = t1*co2_bulk[j1+1][i1+1] + (1.0 - t1)*co2_bulk[j1+1][i1];
  bulk_modulus    = t2*val_s1_2 + (1-t2)*val_s1_1;

  val_s1_1 = t1*co2_density[j1][i1+1] + (1.0 - t1)*co2_density[j1][i1];
  val_s1_2 = t1*co2_density[j1+1][i1+1] + (1.0 - t1)*co2_density[j1+1][i1];
  density  = t2*val_s1_2 + (1-t2)*val_s1_1;

  /*val_s1_1 = t1*co2_velocity[j1][i1+1] + (1.0 - t1)*co2_velocity[j1][i1];
  val_s1_2 = t1*co2_velocity[j1+1][i1+1] + (1.0 - t1)*co2_velocity[j1+1][i1];
  velocity = t2*val_s1_2 + (1-t2)*val_s1_1;*/

}

double
DEMTools::CalcEffectiveElasticModuliUsingHill(double property1,
                                              double volumefraction1,
                                              double property2,
                                              double volumefraction2) {

  if (volumefraction2 < 0)
    volumefraction2 = 1.0 - volumefraction1;

  if (volumefraction1 + volumefraction2 > 1.0) {
    //NBNB fjellvoll add some warning
    volumefraction2 = 1 - volumefraction1;
  }

  return 0.5*(CalcEffectiveElasticModuliUsingReuss(property1,
                                                   volumefraction1,
                                                   property2,
                                                   volumefraction2) +
               CalcEffectiveElasticModuliUsingVoigt(property1,
                                                    volumefraction1,
                                                    property2,
                                                    volumefraction2));


}

double
DEMTools::CalcEffectiveElasticModuliUsingReuss(double property1,
                                               double volumefraction1,
                                               double property2,
                                               double volumefraction2) {

  if (volumefraction2 < 0)
    volumefraction2 = 1.0 - volumefraction1;

  if (volumefraction1 + volumefraction2 > 1.0) {
    //NBNB fjellvoll add some warning
    volumefraction2 = 1 - volumefraction1;
  }

  if (property1 == 0.0 || property2 == 0.0)
    throw NRLib::Exception("Invalid arguments:One of the properties are zero.");

  return 1.0/(volumefraction1/property1 + volumefraction2/property2);


}

double
DEMTools::CalcEffectiveElasticModuliUsingVoigt(double property1,
                                               double volumefraction1,
                                               double property2,
                                               double volumefraction2) {

  if (volumefraction2 < 0)
    volumefraction2 = 1.0 - volumefraction1;

  if (volumefraction1 + volumefraction2 > 1.0) {
    //NBNB fjellvoll add some warning
    volumefraction2 = 1 - volumefraction1;
  }

  return volumefraction1*property1 + volumefraction2*property2;

}

double
DEMTools::CalcEffectiveDensity(double rho1,
                               double volumefraction1,
                               double rho2,
                               double volumefraction2) {

  if (volumefraction2 < 0)
    volumefraction2 = 1.0 - volumefraction1;

  if (volumefraction1 + volumefraction2 > 1.0) {
    //NBNB fjellvoll add some warning
    volumefraction2 = 1 - volumefraction1;
  }

  return volumefraction1*rho1 + volumefraction2*rho2;

}


void
DEMTools::CalcEffectiveBulkAndShearModulus(const std::vector<double>&       bulk_modulus,
                                           const std::vector<double>&       shear_modulus,
                                           const std::vector<double>&       aspect_ratio,
                                           std::vector<double>&             concentration,
                                           double                           bulk_modulus_bg,
                                           double                           shear_modulus_bg,
                                           double&                    effective_bulk_modulus,
                                           double&                    effective_shear_modulus) {

  DEM dem(bulk_modulus,
          shear_modulus,
          aspect_ratio,
          concentration,
          bulk_modulus_bg,
          shear_modulus_bg);

  effective_bulk_modulus = effective_shear_modulus = 0;
  dem.CalcEffectiveModulus(effective_bulk_modulus, effective_shear_modulus);

}


double
DEMTools::CalcVelocityOfBrineFromTPS(double temperature,
                                     double pressure,
                                     double salinity) {

  return CalcVelocityOfWaterFromTP(temperature, pressure) +
         salinity*(1170.0 - 9.6*temperature + 0.055*(temperature*temperature) - 8.5E-5*(temperature*temperature*temperature) +
         2.6*pressure - 0.0029*temperature*pressure - 0.0476*(pressure*pressure)) +
         pow(salinity, 1.5)*(780.0 - 10.0*pressure + 0.16*(pressure*pressure)) - 1820*(salinity*salinity);

}

//calculate acoustic velocity of water
double
DEMTools::CalcVelocityOfWaterFromTP(double temperature,
                                    double pressure) {
  std::vector< std::vector<double> > omega(5, std::vector<double>(4, 0.0));

  omega[0][0] = 1402.85;
  omega[0][1] = 1.524;
  omega[0][2] = 3.437E-3;
  omega[0][3] = -1.197E-5;

  omega[1][0] = 4.871;
  omega[1][1] = -0.0111;
  omega[1][2] = 1.739E-4;
  omega[1][3] = -1.628E-6;

  omega[2][0] = -0.04783;
  omega[2][1] = 2.747E-4;
  omega[2][2] = -2.135E-6;
  omega[2][3] = 1.237E-8;

  omega[3][0] = 1.487E-4;
  omega[3][1] = -6.503E-7;
  omega[3][2] = -1.455E-8;
  omega[3][3] = 1.327E-10;

  omega[4][0] =  -2.197E-7;
  omega[4][1] = 7.987E-10;
  omega[4][2] = 5.230E-11;
  omega[4][3] = -4.614E-13;

  double ww = 0.0;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 4; j++)
      ww += omega[i][j]*pow(temperature, i)*pow(pressure, j);

  return ww;
}


double
DEMTools::CalcDensityOfWaterFromTP(double temperature,
                                   double pressure) {

  return  1.0 + 1E-6*(-80*temperature -3.3*(temperature*temperature) + 0.00175*(temperature*temperature*temperature) +
                      489*pressure - 2*temperature*pressure + 0.016*(temperature*temperature)*pressure -
                      1.3E-5*(temperature*temperature*temperature)*pressure -
                      0.333*(pressure*pressure) - 0.002*temperature*(pressure*pressure));

}

void
DEMTools::DebugTestCalcEffectiveModulus(double& effective_bulk_modulus,
                                        double& effective_shear_modulus,
                                        double& effective_density) {

  // Script for rock physics modelling using differential effective medium theory (DEM).
  // Specifications
  // Mineral properties
  // (According to Table A.4.1, p458-459, RP Handbook, Mavko et al. 2009).
  double quartz_k                = 37;   // gpa
  double quartz_mu               = 44;   // gpa
  double quartz_rho              = 2.65; // g/ccm

  double clay_k                  = 21;   // gpa
  double clay_mu                 = 7;    // gpa
  double clay_rho                = 2.6;  // g/ccm

  // Fluid properties
  // (Brine properties according to Batzle, M. and Wang, Z. 1992, and
  //  CO2 properties according to Z. Wang's compilation and measurements, RPH Tools)
  double porepressure            = 20;   // mpa
  double temperature             = 50;   // °c

  double brine_salinity          = 0.05; // 100*//
  double brine_k                 = CalcBulkModulusOfBrineFromTPS(temperature, porepressure, brine_salinity)*0.001 /*gpa*/;

  double brine_rho               = CalcDensityOfBrineFromTPS(temperature, porepressure, brine_salinity); // g/ccm

  double co2_k;
  double co2_rho;

  CalcCo2Prop(co2_k, co2_rho, temperature, porepressure);



  ////
  enum GeometryType {Spherical, Mixed};
  GeometryType my_geo_type = Mixed; // pore geometry, this enum chooses the one to use.
  // example 1: spherical
  double inclusiongeometry_spectrum          = 1.0;
  double inclusiongeometry_concentration     = 1.0;
  // example 2: mixed (oblate)
  double aspect_ratio2_tmp []          = {1, 0.5000, 0.1000, 0.0100, 1.0000e-003, 1.0000e-004};
  double concentration2_tmp []     = {0.6419, 0.3205, 0.0321, 0.0050, 5.0000e-004, 5.0000e-005};
  std::vector<double> aspect_ratio2(aspect_ratio2_tmp, aspect_ratio2_tmp + 6);
  std::vector<double> concentration2(concentration2_tmp, concentration2_tmp + 6);

  //// plf
  ////
  double porosity = 0.2;
  double lithology = 0.15;     // volume fraction of clay
  double saturation = 0.8;      // volume fraction of brine
  ////

  //// mixing
  // effective solid properties
  ////
  double solid_k                 = CalcEffectiveElasticModuliUsingHill(clay_k, lithology, quartz_k);
  double solid_mu                = CalcEffectiveElasticModuliUsingHill(clay_mu, lithology, quartz_mu);
  double solid_rho               = CalcEffectiveDensity(clay_rho, lithology, quartz_rho);


  ////
  // effective fluid properties
  double fluid_k                 = CalcEffectiveElasticModuliUsingReuss(brine_k, saturation, co2_k);     // homogeneous
  // fluid.k                 = geqhill(brine.k, saturation, co2.k);      // patchy
  double fluid_rho               = CalcEffectiveDensity(brine_rho, saturation, co2_rho);
  ////

  //// effective rock physics properties
  //
  ////
 effective_density                = CalcEffectiveDensity(fluid_rho, porosity, solid_rho);
  // without gassmann
  /*[rock.k rock.mu]        = geqdem(   [solid.k fluid.k*ones(1, nrofinclgeo)], ...
                                      [solid.mu 0*ones(1, nrofinclgeo)], ...
                                      inclusiongeometry.spectrum, ...
                                      porosity*inclusiongeometry.concentration);*/

 std::vector<double> bulk_modulus;
 std::vector<double> shear_modulus;
 std::vector<double> aspect_ratio;
 std::vector<double> concentration;

 if (my_geo_type == Spherical) {
   bulk_modulus =  std::vector<double>(1, fluid_k);
   shear_modulus = std::vector<double>(1, 0.0);
   aspect_ratio = std::vector<double>(1, inclusiongeometry_spectrum);
   concentration = std::vector<double>(1, inclusiongeometry_concentration*porosity);
 }
 else if (my_geo_type == Mixed) {
   bulk_modulus =  std::vector<double>(aspect_ratio2.size(), fluid_k);
   shear_modulus = std::vector<double>(aspect_ratio2.size(), 0.0);
   aspect_ratio  = aspect_ratio2;
   for (size_t i = 0; i < concentration2.size(); i++)
    concentration2[i] *= porosity;
   concentration = concentration2;
 }

  CalcEffectiveBulkAndShearModulus(bulk_modulus,
                                   shear_modulus,
                                   aspect_ratio,
                                   concentration,
                                   solid_k,
                                   solid_mu,
                                   effective_bulk_modulus,
                                   effective_shear_modulus);


}


void
DEMTools::DebugTestCalcEffectiveModulus2(double& effective_bulk_modulus,
                                         double& effective_shear_modulus,
                                         double& effective_density) {

  // Script for rock physics modelling using differential effective medium theory (DEM).
  // Specifications
  // Mineral properties
  // (According to Table A.4.1, p458-459, RP Handbook, Mavko et al. 2009).
  double quartz_k                = 37;   // gpa
  double quartz_mu               = 44;   // gpa
  double quartz_rho              = 2.65; // g/ccm

  Mineral quartz(quartz_k, quartz_mu, quartz_rho);

  double clay_k                  = 21;   // gpa
  double clay_mu                 = 7;    // gpa
  double clay_rho                = 2.6;  // g/ccm

  Mineral clay(clay_k, clay_mu, clay_rho);

  // Fluid properties
  // (Brine properties according to Batzle, M. and Wang, Z. 1992, and
  //  CO2 properties according to Z. Wang's compilation and measurements, RPH Tools)
  double porepressure            = 20;   // mpa
  double temperature             = 50;   // °c

  double brine_salinity          = 0.05; // 100*//

  Brine brine(brine_salinity);

  double brine_k;
  double brine_rho;

  brine.ComputeElasticParams(temperature, porepressure, brine_k, brine_rho);

  double co2_k;
  double co2_rho;

  CO2 co2;
  co2.ComputeElasticParams(temperature, porepressure, co2_k, co2_rho);

  ////
  enum GeometryType {Spherical, Mixed};
  GeometryType my_geo_type = Mixed; // pore geometry, this enum chooses the one to use.
  // example 1: spherical
  double inclusiongeometry_spectrum          = 1.0;
  double inclusiongeometry_concentration     = 1.0;
  // example 2: mixed (oblate)
  double aspect_ratio2_tmp []          = {1, 0.5000, 0.1000, 0.0100, 1.0000e-003, 1.0000e-004};
  double concentration2_tmp []     = {0.6419, 0.3205, 0.0321, 0.0050, 5.0000e-004, 5.0000e-005};
  std::vector<double> aspect_ratio2(aspect_ratio2_tmp, aspect_ratio2_tmp + 6);
  std::vector<double> concentration2(concentration2_tmp, concentration2_tmp + 6);

  //// plf
  ////
  double porosity = 0.2;
  double lithology = 0.15;     // volume fraction of clay
  double saturation = 0.8;      // volume fraction of brine
  ////

  //// mixing
  // effective solid properties
  ////
  std::vector<double> volume_fraction(2);

  volume_fraction[0] = lithology;
  volume_fraction[1] = 1.0 - volume_fraction[0];

  std::vector<Mineral*> mineral(2);
  mineral[0] = &clay;
  mineral[1] = &quartz;

  SolidMixed solidmixed(mineral, volume_fraction);

  // effective fluid properties
  //General problem: since a CO2 or Brine do not know k or rho we cannot use these classes itself to mix fluids! Is this really correct?
  //One solution would be to store the k and rho values inside the class or pointers to the values outside of class.

  std::vector<double> k(2);
  std::vector<double> rho(2);
  k[0] = brine_k;
  k[1] = co2_k;

  std::vector<double> volume_fraction2(2);
  volume_fraction2[0] = saturation;
  volume_fraction2[1] = 1.0 - volume_fraction2[0];

  rho[0] = brine_rho;
  rho[1] = co2_rho;


  FluidMixed fluid_mix(k, rho, volume_fraction2);
  ////


  //// effective rock physics properties
  // without gassmann
  /*[rock.k rock.mu]        = geqdem(   [solid.k fluid.k*ones(1, nrofinclgeo)], ...
                                      [solid.mu 0*ones(1, nrofinclgeo)], ...
                                      inclusiongeometry.spectrum, ...
                                      porosity*inclusiongeometry.concentration);*/

 std::vector<double> bulk_modulus;
 std::vector<double> shear_modulus;
 std::vector<double> aspect_ratio;
 std::vector<double> concentration;

 double fluid_k, dummy1=44.0, dummy2=5.0, dummy3=6.0;
 fluid_mix.ComputeElasticParams(dummy1, dummy2, fluid_k, dummy3);
 if (my_geo_type == Spherical) {
   bulk_modulus =  std::vector<double>(1, fluid_k);
   shear_modulus = std::vector<double>(1, 0.0);
   aspect_ratio = std::vector<double>(1, inclusiongeometry_spectrum);
   concentration = std::vector<double>(1, inclusiongeometry_concentration*porosity);
 }
 else if (my_geo_type == Mixed) {
   bulk_modulus =  std::vector<double>(aspect_ratio2.size(), fluid_k);
   shear_modulus = std::vector<double>(aspect_ratio2.size(), 0.0);
   aspect_ratio  = aspect_ratio2;
   for (size_t i = 0; i < concentration2.size(); i++)
    concentration2[i] *= porosity;
   concentration = concentration2;
 }

 RockInclusion rock_inclusion(solidmixed, fluid_mix, bulk_modulus, shear_modulus, aspect_ratio, concentration, porosity);

 rock_inclusion.GetSeimsmicParams(effective_bulk_modulus, effective_shear_modulus, effective_density);


}

