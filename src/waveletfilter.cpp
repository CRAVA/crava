#include <fstream>
#include <math.h>

#include "nrlib/surface/surfaceio.hpp"

#include "lib/global_def.h"

#include "src/waveletfilter.h"

using namespace NRLib;

WaveletFilter::WaveletFilter(const std::string & fileName,
                             int               & errCode,
                             char              * errText)
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
  RegularSurfaceRotated<double> rot_surface = ReadSgriSurf(fileName);
  double angle = rot_surface.GetAngle();
  if (angle != 0.0) {
    sprintf(errText, "%sGrid for wavelet filter in file %s is rotated. Must have rotation angle = 0.0.\n", errText, fileName.c_str());
    errCode = 1;
  }
  else {
    alpha1_ = rot_surface.GetSurface();
  }

  return(true);
}

