#ifndef RPLIB_FLUID_SURF_RHO_H
#define RPLIB_FLUID_SURF_RHO_H
#include "nrlib/surface/regularsurface.hpp"
namespace ConstDataStoredAsSurface {
const NRLib::RegularSurface<double> CreateSurfaceRho();
void  SetGridRhoPart1(NRLib::Grid2D<double>& grid2d);
}
#endif
