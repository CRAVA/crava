#ifndef CRAVA_DEFINITIONS_H
#define CRAVA_DEFINITIONS_H

#include <typeinfo>

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include "nrlib/segy/segy.hpp"

typedef  NRLib::RegularSurface<double>    Surface;
typedef  NRLib::Grid2D<double>            Grid2D;
typedef  NRLib::LogKit                    LogKit;
typedef  NRLib::StormContGrid             StormContGrid;
typedef  NRLib::SegY                      SegY;
typedef  NRLib::SegyGeometry              SegyGeometry;
typedef  NRLib::TraceHeaderFormat         TraceHeaderFormat;
typedef  NRLib::TextualHeader             TextualHeader;
/**
   Class to hold definitions and constants used by CRAVA.
*/

namespace Definitions
{
  inline static double AsciiIrapClassicUndefValue(void)       { return 9999900.0 ;}
  inline static double StormBinaryUndefValue(void)            { return    -999.0 ;}
  inline static float  SpatialFilterRegularisationValue(void) { return    0.0001f;}

  #ifndef RMISSING
  #define RMISSING -99999.000
  #endif

  #ifndef WELLMISSING
  #define WELLMISSING -999.00000
  #endif

  #define IMISSING -99999
}

#endif
