#include <fstream>
#include <math.h>

#include "nrlib/surface/surfaceio.hpp"

#include "lib/global_def.h"

#include "src/waveletfilter.h"

WaveletFilter::WaveletFilter(const std::string & /*fileName*/,
                             int               & /*errCode*/,
                             char              * /*errText*/)
 : hasHalpha_(false)
{
//  readFile(fileName, errCode, errText);
}

WaveletFilter::~WaveletFilter(void)
{
}


bool WaveletFilter::readFile(const std::string & fileName,
                             int               & errCode,
                             char              * errText)
{
  std::vector<NRLib::RegularSurfaceRotated<double> > rot_surfaces = NRLib::ReadMultipleSgriSurf(fileName);
  if (rot_surfaces[0].GetAngle() != 0.0) {
    sprintf(errText, "%sGrid for wavelet filter in file %s is rotated. Must have rotation angle = 0.0.\n", errText, fileName.c_str());
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

