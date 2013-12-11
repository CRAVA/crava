#include "rplib/fluidco2.h"

#include "rplib/demmodelling.h"

#include "rplib/table_vp1.h"
#include "rplib/table_vp2.h"
#include "rplib/table_rho1.h"
#include "rplib/table_rho2.h"
#include "rplib/table_rho.h"
#include "rplib/table_vp.h"

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
  static const NRLib::RegularSurface<double> surf_vp    = ConstDataStoredAsSurface::CreateSurfaceVP();
  static const NRLib::RegularSurface<double> surf_rho   = ConstDataStoredAsSurface::CreateSurfaceRho();

  static const NRLib::RegularSurface<double> surf_vp1   = ConstDataStoredAsSurface::CreateSurfaceVP1();
  static const NRLib::RegularSurface<double> surf_vp2   = ConstDataStoredAsSurface::CreateSurfaceVP2();
  static const NRLib::RegularSurface<double> surf_rho1  = ConstDataStoredAsSurface::CreateSurfaceRho1();
  static const NRLib::RegularSurface<double> surf_rho2  = ConstDataStoredAsSurface::CreateSurfaceRho2();

  static std::vector<double> p_func                     = GetPFunc();
  static double critical_temp                           = GetCritTemp();
  static double critical_pressure                       = GetCritPressure();

  double scale = 100.0;
  double inv_scale = 1.0/scale;

  if (temp >= 1.0 && temp <= 100.0 && pressure >= 0.1 && pressure <= 100.0) { //very dense sampled table

    //bilinear interpolation, surfaces are multiplied with 100 in x, y and z-direction
    double vp = 0;
    double mu = 0;
    bool failed_getz = false;

    rho_  = surf_rho.GetZ(scale*pressure, scale*temp);
    if (surf_rho.IsMissing(rho_))
      failed_getz = true;

    rho_ *= inv_scale;

    vp    = surf_vp.GetZ(scale*pressure, scale*temp);
    if (surf_vp.IsMissing(vp))
      failed_getz = true;

    vp   *= inv_scale;

    if (failed_getz)
      throw NRLib::Exception("CO2 Model: Interpolation failed.");

    //unit conversion from km/s -> m/s
    vp *= 1000.0;
    DEMTools::CalcElasticParamsFromSeismicParams(vp, 0.0, rho_, k_, mu);

  }
  else {

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
      rho_      = surf_rho1.GetZ(scale*pressure, scale*temp);
      if (surf_rho1.IsMissing(rho_))
        failed_getz = true;
      rho_     *= inv_scale;

      vp        = surf_vp1.GetZ(scale*pressure, scale*temp);
      if (surf_vp1.IsMissing(vp))
        failed_getz = true;
      vp       *= inv_scale;
    }
    else {
      rho_      = surf_rho2.GetZ(scale*pressure, scale*temp);
      if (surf_rho2.IsMissing(rho_))
        failed_getz = true;
      rho_     *= inv_scale;

      vp        = surf_vp2.GetZ(scale*pressure, scale*temp);
       if (surf_vp2.IsMissing(vp))
        failed_getz = true;
      vp       *= inv_scale;
    }

    if (failed_getz)
      throw NRLib::Exception("CO2 Model: Interpolation failed.");

    //unit conversion from km/s -> m/s
    vp *= 1000.0;

    DEMTools::CalcElasticParamsFromSeismicParams(vp, 0.0, rho_, k_, mu);
  }
}


