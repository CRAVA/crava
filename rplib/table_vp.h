#ifndef RPLIB_FLUID_SURF_VP_H
#define RPLIB_FLUID_SURF_VP_H
#include "nrlib/surface/regularsurface.hpp"
namespace ConstDataStoredAsSurface {
const NRLib::RegularSurface<double> CreateSurfaceVP();
void SetGridVPPart1(NRLib::Grid2D<double>& grid2d);
}
#endif
