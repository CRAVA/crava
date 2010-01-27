#include <fstream>
#include <math.h>

#include "nrlib/surface/surfaceio.hpp"

#include "lib/global_def.h"

#include "src/waveletfilter.h"

WaveletFilter::WaveletFilter(const std::string & /*fileName*/,
                             int               & /*errCode*/,
                             std::string       & /*errText*/)
 : hasHalpha_(false)
{
//  readFile(fileName, errCode, errText);
  createGrid();
}

WaveletFilter::~WaveletFilter(void)
{
}


bool WaveletFilter::readFile(const std::string & fileName,
                             int               & errCode,
                             std::string       & errText)
{
  std::vector<NRLib::RegularSurfaceRotated<double> > rot_surfaces = NRLib::ReadMultipleSgriSurf(fileName);
  if (rot_surfaces[0].GetAngle() != 0.0) {
    errText += "Grid for wavelet filter in file "+fileName+" is rotated. Must have rotation angle = 0.0.\n";
    errCode = 1;
  }
  else {
    alpha1_ = rot_surfaces[0].GetSurface();
    if (rot_surfaces.size() > 1) {
      Halpha_ = rot_surfaces[1].GetSurface();
      hasHalpha_ = true;
    }
  }

  return(true);
}

void WaveletFilter::createGrid()
{
  alpha1_ = NRLib::RegularSurface<double>(0,0,2*PI,0.5*PI,360,90,0.0);
  for (unsigned int i = 0; i<360; i++)
    for (unsigned int j = 0; j <= 45; j++)
      alpha1_(i,j) = 1.0;
}
