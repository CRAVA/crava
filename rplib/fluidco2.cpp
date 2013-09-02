#include "rplib/fluidco2.h"

#include "rplib/demmodelling.h"

#include "rplib/surf_vp1.h"
#include "rplib/surf_vp2.h"
#include "rplib/surf_rho1.h"
#include "rplib/surf_rho2.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/exception/exception.hpp"

#include <cassert>
#include <fstream>
#include <string>



//Helper file for creating c++ file for storing surface
void WriteFileForStoringSurface(NRLib::RegularSurface<double>       & surf,
                                const std::string                   & first_part_of_path,
                                const std::string                   & out_file_name, //without extension and relative path
                                const std::string                   & name_of_function) {

  std::string full_file_name_cpp = first_part_of_path + out_file_name + ".cpp";
  std::ofstream ofile(full_file_name_cpp.c_str());

  if(!ofile)
    throw NRLib::Exception("Could not open file " + full_file_name_cpp + " for writing.");

  // Write first part
  ofile << "//This file is generated automatically with the function WriteHeaderFileForStoringSurface located in rplib/fluidco2.cpp." << "\n";
  ofile << "#include \"rplib/" << out_file_name << ".h\n";
  ofile << "#include \"nrlib/surface/regularsurface.hpp\"" << "\n";

  ofile << "namespace ConstDataStoredAsSurface {" << "\n";
  ofile << "const NRLib::RegularSurface<double> " << name_of_function << "() {" << "\n";

  //need grid2d
  ofile << "NRLib::Grid2D<double> grid2d(" << surf.GetNI() << "," << surf.GetNJ() << ", 0.0);" << "\n";

  int j = 0;
  NRLib::RegularSurface<double>::const_iterator i;
  for (i = surf.begin(); i != surf.end(); ++i) {
    ofile << "grid2d(" << j << ") =" <<  *i << "; ";

    if (j % 100 == 0 && j != 0)
      ofile << std::endl;
    j++;
  }

  ofile << "\nconst NRLib::RegularSurface<double> surf_final(" << surf.GetXMin() << "," << surf.GetYMin() << "," <<
           surf.GetLengthX() << "," << surf.GetLengthY() << "," << "grid2d);\n";
  ofile << "return surf_final;" << std::endl;

  ofile << "}" << std::endl; //end function
  ofile << "}" << std::endl; //end namespace
}

FluidCO2::FluidCO2(double                      temp,
                   double                      pore_pressure,
                   const std::vector<double> & u)
: Fluid()
{
  u_ = u; //u contains independent samples used in quantiles of (temp,pore_pressure)
  ComputeElasticParams(temp, pore_pressure);

  ////testing along discontinuity line
  //double xTemp[]       = {-56.5570000000000, -55.1490000000000, -53.1490000000000, -51.1490000000000, -49.1490000000000, -47.1490000000000, -45.1490000000000, -43.1490000000000, -41.1490000000000, -39.1490000000000, -37.1490000000000, -35.1490000000000, -33.1490000000000, -31.1490000000000, -29.1490000000000, -27.1490000000000, -25.1490000000000, -23.1490000000000, -21.1490000000000 ,-19.1490000000000, -17.1490000000000, -15.1490000000000, -13.1490000000000, -11.1490000000000, -9.14900000000000, -7.14900000000000, -5.14900000000000, -3.14900000000000, -1.14900000000000, 0.850999999999999, 2.85100000000000, 4.85100000000000, 6.85100000000000, 8.85100000000000, 10.8510000000000, 12.8510000000000, 14.8510000000000, 16.8510000000000, 18.8510000000000, 20.8510000000000, 22.8510000000000, 24.8510000000000, 26.8510000000000, 27.8510000000000, 28.8510000000000, 29.8510000000000, 30.8510000000000, 30.9783000000000 };

  //double yPres[]       = {0.517950000000000, 0.550410000000000, 0.599120000000000, 0.651010000000000, 0.706200000000000, 0.764830000000000, 0.827020000000000, 0.892900000000000 ,0.962610000000000 ,1.03620000000000, 1.11400000000000, 1.19600000000000 ,1.28240000000000 ,1.37330000000000 ,1.46890000000000 ,1.56920000000000, 1.67450000000000 ,1.78490000000000 ,1.90060000000000 ,2.02160000000000, 2.14820000000000, 2.28050000000000, 2.41870000000000, 2.56290000000000, 2.71330000000000, 2.87000000000000, 3.03330000000000, 3.20320000000000, 3.38010000000000, 3.56410000000000 ,3.75540000000000, 3.95410000000000, 4.16060000000000, 4.37510000000000, 4.59770000000000, 4.82880000000000 ,5.06870000000000, 5.31760000000000 ,5.57600000000000 ,5.84420000000000, 6.12260000000000, 6.41200000000000 ,6.71300000000000, 6.86820000000000, 7.02670000000000 ,7.18890000000000, 7.35540000000000, 7.37720000000000 };

  //double epsilon = 0.50;
  //double epsilon1 = 0.3;
  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i], yPres[i]);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] + epsilon1, yPres[i]);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] - epsilon1, yPres[i]);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i], yPres[i] - epsilon);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] + epsilon1, yPres[i] - epsilon);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] - epsilon1, yPres[i] -epsilon);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i], yPres[i] + epsilon);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] + epsilon1, yPres[i] + epsilon);

  //for (int i = 0; i < 48; ++i)
  //  ComputeElasticParams(xTemp[i] - epsilon1, yPres[i] + epsilon);

  //// other testing
  //ComputeElasticParams(17.0, 0.1);

  //ComputeElasticParams(17.0, 40.02);
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
  ////use this code part to generate new hardcoded surfaces if the surface changes
  //static NRLib::RegularSurface<double>  surf_vp1("K:/user/fjellvoll/GitHub/crava2/rplib/rox_surf/rox_surf_vp1_v2.storm");
  //static NRLib::RegularSurface<double>  surf_vp2("K:/user/fjellvoll/GitHub/crava2/rplib/rox_surf/rox_surf_vp2_v2.storm");
  //static NRLib::RegularSurface<double>  surf_rho1("K:/user/fjellvoll/GitHub/crava2/rplib/rox_surf/rox_surf_rho1_v2.storm");
  //static NRLib::RegularSurface<double>  surf_rho2("K:/user/fjellvoll/GitHub/crava2/rplib/rox_surf/rox_surf_rho2_v2.storm");

  //WriteFileForStoringSurface(surf_vp1, "K:/user/fjellvoll/GitHub/crava2/",  "rplib/rox_surf/surf_vp1", "CreateSurfaceVP1");
  //WriteFileForStoringSurface(surf_vp2, "K:/user/fjellvoll/GitHub/crava2/",  "rplib/rox_surf/surf_vp2", "CreateSurfaceVP2");
  //WriteFileForStoringSurface(surf_rho1, "K:/user/fjellvoll/GitHub/crava2/", "rplib/rox_surf/surf_rho1", "CreateSurfaceRho1);
  //WriteFileForStoringSurface(surf_rho2, "K:/user/fjellvoll/GitHub/crava2/", "rplib/rox_surf/surf_rho2", "CreateSurfaceRho2);


  static const NRLib::RegularSurface<double> surf_vp1   = ConstDataStoredAsSurface::CreateSurfaceVP1();
  static const NRLib::RegularSurface<double> surf_vp2   = ConstDataStoredAsSurface::CreateSurfaceVP2();
  static const NRLib::RegularSurface<double> surf_rho1  = ConstDataStoredAsSurface::CreateSurfaceRho1();
  static const NRLib::RegularSurface<double> surf_rho2  = ConstDataStoredAsSurface::CreateSurfaceRho2();

  static double critical_temp       = 30.9783;
  static double critical_pressure   = 7.3772;
  //static double flag_limit           = 0.5;

  static double p_func[]            = { 0.000003883701221,   0.000930453898625,   0.092759127015099,   3.480623887319762};

  // Sort input data into domains
  //Domain 1: Gas
  //Domain 2: Liquid and supercritical fluid

  double p_of_t             = p_func[0]*temp*temp*temp + p_func[1]*temp*temp + p_func[2]*temp + p_func[3];

  int scenario              = 1;
  //int    flag               = 0; //0: no warning,1: close to liquid-gas discontinuity, 2: close to critical point, 3: close to liquid-gas discontinuity and critical point

  if ((temp >= critical_temp && pressure >= critical_pressure) ||
      (temp < critical_temp && pressure >= p_of_t))
      scenario = 2;

  /* // for debugging
  if (temp <= critical_temp && std::abs(pressure - p_of_t) <= flag_limit)
    flag = 1;

  if ((std::abs(temp - critical_temp) <= flag_limit) &&
      (std::abs(pressure - critical_pressure) <= flag_limit))
    flag += 2;*/

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


