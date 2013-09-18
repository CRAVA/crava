#include "rplib/fluidco2.h"

#include "rplib/demmodelling.h"

#include "rplib/table_vp1.h"
#include "rplib/table_vp2.h"
#include "rplib/table_rho1.h"
#include "rplib/table_rho2.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/exception/exception.hpp"

#include <cassert>
#include <fstream>
#include <string>

FluidCO2::FluidCO2(double                      temp,
                   double                      pore_pressure,
                   const std::vector<double> & u)
: Fluid()
{
  u_ = u; //u contains independent samples used in quantiles of (temp,pore_pressure)
  ComputeElasticParams(temp, pore_pressure);

}

FluidCO2::FluidCO2(const FluidCO2 & rhs) : Fluid(rhs)
{
}

FluidCO2::~FluidCO2(){}


Fluid *
FluidCO2::Clone() const
{
  return new FluidCO2(*this);
}

void
FluidCO2::ComputeElasticParams(double temp, double pressure)
{
  static const NRLib::RegularSurface<double> surf_vp1   = ConstDataStoredAsSurface::CreateSurfaceVP1();
  static const NRLib::RegularSurface<double> surf_vp2   = ConstDataStoredAsSurface::CreateSurfaceVP2();
  static const NRLib::RegularSurface<double> surf_rho1  = ConstDataStoredAsSurface::CreateSurfaceRho1();
  static const NRLib::RegularSurface<double> surf_rho2  = ConstDataStoredAsSurface::CreateSurfaceRho2();

  static std::vector<double> p_func                     = GetPFunc();
  static double critical_temp                           = GetCritTemp();
  static double critical_pressure                       = GetCritPressure();

  // Sort input data into domains
  //Domain 1: Gas
  //Domain 2: Liquid and supercritical fluid

  double p_of_t             = p_func[0]*temp*temp*temp + p_func[1]*temp*temp + p_func[2]*temp + p_func[3];

  int scenario              = 1;

  if ((temp >= critical_temp && pressure >= critical_pressure) ||
      (temp < critical_temp && pressure >= p_of_t))
      scenario = 2;

  //bilinear interpolation, surfaces are multiplied with 100 in x, y and z-direction
  double vp = 0;
  double mu = 0;
  bool failed_getz = false;
  if (scenario == 1) {
    rho_      = surf_rho1.GetZ(100.0*pressure, 100.0*temp);
    if (surf_rho1.IsMissing(rho_))
      failed_getz = true;
    rho_     /= 100.0;

    vp        = surf_vp1.GetZ(100.0*pressure, 100.0*temp);
    if (surf_vp1.IsMissing(vp))
      failed_getz = true;
    vp       /= 100.0;
  }
  else {
    rho_      = surf_rho2.GetZ(100.0*pressure, 100.0*temp);
    if (surf_rho2.IsMissing(rho_))
      failed_getz = true;
    rho_     /= 100.0;

    vp        = surf_vp2.GetZ(100.0*pressure, 100.0*temp);
     if (surf_vp2.IsMissing(vp))
      failed_getz = true;
    vp       /= 100.0;
  }

  if (failed_getz)
    throw NRLib::Exception("CO2 Model: Interpolation failed.");

  //unit conversion from km/s -> m/s
  vp *= 1000.0;

  DEMTools::CalcElasticParamsFromSeismicParams(vp, 0.0, rho_, k_, mu);
}


void
FluidCO2::DebugTest() {
  //simple test function to test discontinuity line between gas and liquid between tripple and critical point

  //testing along discontinuity line
  double xTemp[]       = {-56.5570000000000, -55.1490000000000, -53.1490000000000, -51.1490000000000, -49.1490000000000, -47.1490000000000, -45.1490000000000, -43.1490000000000, -41.1490000000000, -39.1490000000000, -37.1490000000000, -35.1490000000000, -33.1490000000000, -31.1490000000000, -29.1490000000000, -27.1490000000000, -25.1490000000000, -23.1490000000000, -21.1490000000000 ,-19.1490000000000, -17.1490000000000, -15.1490000000000, -13.1490000000000, -11.1490000000000, -9.14900000000000, -7.14900000000000, -5.14900000000000, -3.14900000000000, -1.14900000000000, 0.850999999999999, 2.85100000000000, 4.85100000000000, 6.85100000000000, 8.85100000000000, 10.8510000000000, 12.8510000000000, 14.8510000000000, 16.8510000000000, 18.8510000000000, 20.8510000000000, 22.8510000000000, 24.8510000000000, 26.8510000000000, 27.8510000000000, 28.8510000000000, 29.8510000000000, 30.8510000000000, 30.9783000000000 };

  double yPres[]       = {0.517950000000000, 0.550410000000000, 0.599120000000000, 0.651010000000000, 0.706200000000000, 0.764830000000000, 0.827020000000000, 0.892900000000000 ,0.962610000000000 ,1.03620000000000, 1.11400000000000, 1.19600000000000 ,1.28240000000000 ,1.37330000000000 ,1.46890000000000 ,1.56920000000000, 1.67450000000000 ,1.78490000000000 ,1.90060000000000 ,2.02160000000000, 2.14820000000000, 2.28050000000000, 2.41870000000000, 2.56290000000000, 2.71330000000000, 2.87000000000000, 3.03330000000000, 3.20320000000000, 3.38010000000000, 3.56410000000000 ,3.75540000000000, 3.95410000000000, 4.16060000000000, 4.37510000000000, 4.59770000000000, 4.82880000000000 ,5.06870000000000, 5.31760000000000 ,5.57600000000000 ,5.84420000000000, 6.12260000000000, 6.41200000000000 ,6.71300000000000, 6.86820000000000, 7.02670000000000 ,7.18890000000000, 7.35540000000000, 7.37720000000000 };

  double epsilon = 0.50;
  double epsilon1 = 0.3;
  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i], yPres[i]);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] + epsilon1, yPres[i]);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] - epsilon1, yPres[i]);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i], yPres[i] - epsilon);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] + epsilon1, yPres[i] - epsilon);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] - epsilon1, yPres[i] -epsilon);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i], yPres[i] + epsilon);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] + epsilon1, yPres[i] + epsilon);

  for (int i = 0; i < 48; ++i)
    ComputeElasticParams(xTemp[i] - epsilon1, yPres[i] + epsilon);
}


