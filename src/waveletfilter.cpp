#include <fstream>
#include <math.h>

#include "nrlib/surface/regularsurfacerotated.hpp"

#include "src/waveletfilter.h"


WaveletFilter::WaveletFilter(const std::string & fileName,
                             int               & errCode,
                             std::string       & errText)
 : hasHalpha_(false)
{
  readFile(fileName, errCode, errText);
//  createGrid();
}

WaveletFilter::~WaveletFilter(void)
{
}


void WaveletFilter::readFile(const std::string & fileName,
                             int               & errCode,
                             std::string       & errText)
{
  errCode = 0;
  try {
    std::vector<std::string> labels;
//    labels.push_back("arg1");  NBNB Odd
//    labels.push_back("arg2");
    std::vector<NRLib::RegularSurfaceRotated<float> > rot_surfaces = NRLib::ReadMultipleSgriSurf(fileName, labels);

    if(rot_surfaces[0].GetAngle() != 0.0) {
      errText += "Grid for wavelet filter in file "+fileName+" is rotated. Must have rotation angle = 0.0.\n";
      errCode = 1;
    }
    else {
      alpha1_ = rot_surfaces[0].ResampleSurface();
      if (rot_surfaces.size() > 1) {
        Halpha_ = rot_surfaces[1].ResampleSurface();
        hasHalpha_ = true;
      }
    }
  }
  catch (NRLib::Exception& e) {
    errText += e.what();
    errCode = 1;
  }
}

/*void WaveletFilter::createGrid()
{
  alpha1_ = NRLib::RegularSurface<float>(0,0,2*NRLib::Pi,0.5*NRLib::Pi,360,90,0.0);
  for (unsigned int i = 0; i<360; i++)
    for (unsigned int j = 0; j <= 45; j++)
      alpha1_(i,j) = 1.0;
}*/
const float
WaveletFilter::getHalpha(double phi, double psi)const
{
  if(hasHalpha_)
    return Halpha_.GetZ(phi*180/NRLib::Pi,psi*180/NRLib::Pi);
  else
    return 0.0f;
}

const float
WaveletFilter::getAlpha1(double phi, double psi) const 
{
  return alpha1_.GetZ(phi*180/NRLib::Pi,psi*180/NRLib::Pi);
}
  
const float
WaveletFilter::getBeta1(double phi, double psi)  const 
{
  return beta1_.GetZ(phi*180/NRLib::Pi,psi*180/NRLib::Pi)   ;
}

const float
WaveletFilter::getHbeta(double phi, double psi)  const 
{
  return Hbeta_.GetZ(phi*180/NRLib::Pi,psi*180/NRLib::Pi);
}
