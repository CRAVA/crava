#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lib/lib_misc.h"
#include "lib/global_def.h"

#include "src/vario.h"

Vario::Vario(float range1, float range2, float angle)
{
  type_   = new char[MAX_STRING];
  range1_ = range1;
  range2_ = range2;
  //
  // It seems (from txx and tyy below) that all angles are mathematical angles, 
  // that is, angles relative to the x-axis in counterclockwise direction.
  //
  angle_  = static_cast<float>(angle*M_PI/180.0); // conversion to radians, but no change of reference points.
}

Vario::~Vario()
{
  delete [] type_;
}

float Vario::findDist(float deltaX, float deltaY) const
{
  float dist; 
  float txx, txy, tyy;
  txx = static_cast<float>(cos(angle_)*cos(angle_)/(range1_*range1_) + sin(angle_)*sin(angle_)/(range2_*range2_));
  tyy = static_cast<float>(sin(angle_)*sin(angle_)/(range1_*range1_) + cos(angle_)*cos(angle_)/(range2_*range2_));
  txy = static_cast<float>(2*cos(angle_)*sin(angle_)*(1.0/(range1_*range1_)-1.0/(range2_*range2_)));
  dist =static_cast<float>(sqrt(txx*deltaX*deltaX+tyy*deltaY*deltaY+txy*deltaX*deltaY));
  return(dist);
}

void
Vario::convertRangesFromDegToRad()
{
  range1_ *= static_cast<float>(PI/180.0);
  range2_ *= static_cast<float>(PI/180.0);
}

void
Vario::rotateCounterClockwise(float rotAngle)
{
  angle_ += rotAngle;
}


SphericalVario::SphericalVario(float range1, float range2, float angle)
: Vario(range1, range2, angle)
{
  sprintf(type_,"Spherical");
}

SphericalVario::~SphericalVario()
{
}

float
SphericalVario::corr(float deltaX, float deltaY) const
{
  float dist = findDist(deltaX, deltaY);
  if(dist>=1.0)
    return(0.0);
  else
    return(float (1.0+dist*dist*dist/2.0-1.5*dist));
}


GenExpVario::GenExpVario(float pot, float range1, float range2, float angle)
: Vario(range1, range2, angle)
{
  sprintf(type_,"Generalised exponential");
  pot_ = pot;
}

GenExpVario::~GenExpVario()
{
}

float
GenExpVario::corr(float deltaX, float deltaY) const
{
  float dist;
  dist = findDist(deltaX, deltaY);
  return(float (exp(-3.0*pow(dist,pot_))));

}
