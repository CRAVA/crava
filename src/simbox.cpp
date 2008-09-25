#include <math.h>
#include <assert.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "nrlib/volume/volume.hpp"
#include "nrlib/surface/surfaceio.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "src/simbox.h"
#include "src/model.h"
#include "src/definitions.h"

Simbox::Simbox(void) 
  : Volume()
{
  status_      = EMPTY;
  topName_     = "";
  botName_     = "";
  inLine0_     = 0;
  crossLine0_  = 0;
  ilStep_      = 1;
  xlStep_      = 1;
  constThick_  = true;
  minRelThick_ = 1.0;
  dz_          = 0;
}

Simbox::Simbox(double x0, double y0, Surface * z0, double lx, 
               double ly, double lz, double rot, double dx, double dy, double dz) :
  Volume()
{
  status_      = BOXOK;
  topName_     = "";
  botName_     = "";
  SetDimensions(x0,y0,lx,ly);
  SetAngle(rot);
  
  Surface * z1 = new Surface(*z0);
  z1->Add(lz);
  SetSurfaces(z0,z1); //Automatically sets lz correct in this case.

  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  dz_          = dz;
  nx_          = int(0.5+lx/dx_);
  ny_          = int(0.5+ly/dy_);
  nz_          = int(0.5+lz/dz_);
  inLine0_     = 0;
  crossLine0_  = 0;
  constThick_  = true;
  minRelThick_ = 1.0;
}

Simbox::Simbox(const Simbox *simbox) : 
  Volume(*simbox)
{
  status_      = simbox->status_;
  cosrot_      = cos(GetAngle());
  sinrot_      = sin(GetAngle());
  dx_          = simbox->dx_;
  dy_          = simbox->dy_;
  dz_          = simbox->dz_;
  nx_          = simbox->nx_;
  ny_          = simbox->ny_;
  nz_          = simbox->nz_;
  inLine0_     = simbox->inLine0_;
  crossLine0_  = simbox->crossLine0_;
  constThick_  = simbox->constThick_;
  minRelThick_ = simbox->minRelThick_;
  topName_     = simbox->topName_;
  botName_     = simbox->botName_;
}   

Simbox::~Simbox()
{
}

int
Simbox::getIndex(double x, double y, double z) const

{
  int index = IMISSING;
  int i, j, k;
  getIndexes(x,y,z,i,j,k);
  if(k != IMISSING && j != IMISSING && i != IMISSING)
    index = int(i+j*nx_+k*nx_*ny_);
  return(index);
}

int
Simbox::getClosestZIndex(double x, double y, double z)
{
  int index = IMISSING;
  int i, j, k;
  getIndexesFull(x,y,z,i,j,k);
  if(i >=0 && i < nx_ && j >=0 && j < ny_)
  {
    if(k < 0)
      k = 0;
    else if(k >= nz_)
      k = nz_-1;
    index = i+j*nx_+k*nx_*ny_;
  }
  return(index);
}

void 
Simbox::getIndexes(double x, double y, double z, int & xInd, int & yInd, int & zInd) const
{
  xInd = IMISSING;
  yInd = IMISSING;
  zInd = IMISSING;
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx > 0 && rx < GetLX() && ry>0 && ry < GetLY())
  {
    double zBot, zTop = GetTopSurface().GetZ(x,y);
    if(GetTopSurface().IsMissing(zTop) == false)
    {
      zBot = GetBotSurface().GetZ(x,y);
      if(GetBotSurface().IsMissing(zBot) == false &&  z > zTop && z < zBot)
      {
        xInd = int(floor(rx/dx_));
        yInd = int(floor(ry/dy_));
        zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
        //LogKit::LogFormatted(LogKit::LOW,"rx,dx,xInd = %.4f %.4f %d   ry,dy,yInd = %.4f %.4f %d    %d\n",rx,dx_,xInd,ry,dy_,yInd,zInd);
      }
    }
  }
}

void 
Simbox::getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  xInd = int(floor(rx/dx_));
  yInd = int(floor(ry/dy_));
  zInd = IMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
      zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
  }
}

void
Simbox::getZInterpolation(double x, double y, double z, 
                          int & index1, int & index2, double & t) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  int xInd = int(floor(rx/dx_));
  int yInd = int(floor(ry/dy_));
  int zInd2, zInd1;
  index1 = IMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      zInd1 = static_cast<int>(floor((z-zTop)/dz)-0.5); //Find cell center above.
      if(zInd1 >=0 && zInd1 < nz_-1)
      {
        t = (z-zTop)/dz - 0.5 - static_cast<double>(zInd1);
        zInd2 = zInd1+1;
      }
      else 
      {
        t = 0;
        if(zInd1 < 0)
          zInd1 = 0;
        else
          zInd1 = nz_-1;
        zInd2 = zInd1;
      }
      index1 = xInd+yInd*nx_+zInd1*nx_*ny_;
      index2 = xInd+yInd*nx_+zInd2*nx_*ny_;
    }
  }
}

void  
Simbox::getCoord(int xInd, int yInd, int zInd, double &x, double &y, double &z) const
{
  double rx = (static_cast<double>(xInd) + 0.5)*dx_;
  double ry = (static_cast<double>(yInd) + 0.5)*dy_;
  x = rx*cosrot_-ry*sinrot_ + GetXMin();
  y = rx*sinrot_+ry*cosrot_ + GetYMin();
  z = RMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      z = zTop + (static_cast<double>(zInd) + 0.5)*dz;
    }
  }
}


void
Simbox::getMinMaxZ(double &minZ, double &maxZ) const
{
  minZ = GetTopSurface().Min();
  maxZ = GetBotSurface().Max();
}

int
Simbox::isInside(double x, double y) const
{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx < 0 || rx > GetLX() || ry<0 || ry > GetLY())
    return(0);
  else
    return(1);
}

int
Simbox::insideRectangle(const SegyGeometry *  geometry) const
{
  double xr   = geometry->getX0();
  double yr   = geometry->getY0();
  double rotr = geometry->getAngle();
  double lxr  = geometry->getlx();
  double lyr  = geometry->getly();
  double dxr  = geometry->getDx();
  double dyr  = geometry->getDy();

  // check that incoming rectangle is within simbox +-0.5 grid cells
  int allOk = 1;
  double cosrotr = cos(rotr);
  double sinrotr = sin(rotr);
  double x       = GetXMin();
  double y       = GetYMin();
  double rx      =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  double ry      = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;

  x  = GetXMin()+GetLX()*cosrot_;
  y  = GetYMin()+GetLX()*sinrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;
  
  x  = GetXMin()-GetLY()*sinrot_;
  y  = GetYMin()+GetLY()*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;
  
  x  = GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_;
  y  = GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.49*dx_ || rx > lxr+0.49*dx_ || ry<-0.49*dy_ || ry > lyr+0.49*dy_)
    allOk = 0;
  if(rotr<0)
    rotr+=2*PI;
  
  if (allOk==0) {
    double seisAzimuth = (-1)*rotr*(180/M_PI);
    double areaAzimuth = (-1)*GetAngle()*(180/M_PI);
    if (seisAzimuth < 0) seisAzimuth += 360.0;
    if (areaAzimuth < 0) areaAzimuth += 360.0;
    LogKit::LogFormatted(LogKit::LOW,"                        x0           y0            lx         ly         dx      dy     azimuth\n");
    LogKit::LogFormatted(LogKit::LOW,"--------------------------------------------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::LOW,"Model area:    %11.2f  %11.2f    %10.2f %10.2f    %7.2f %7.2f    %8.3f\n", 
                         GetXMin(), GetYMin(), GetLX(), GetLY(), dx_, dy_, areaAzimuth);
    LogKit::LogFormatted(LogKit::LOW,"Seismic area:  %11.2f  %11.2f    %10.2f %10.2f    %7.2f %7.2f    %8.3f\n", 
                         xr, yr, lxr, lyr, dxr, dyr, seisAzimuth);
    
    LogKit::LogFormatted(LogKit::HIGH,"\nCorner     XY Area                    XY Seismic\n");
    LogKit::LogFormatted(LogKit::HIGH,"-----------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::HIGH,"A %18.2f %11.2f    %11.2f %11.2f\n", GetXMin(),GetYMin(), xr,yr);
    LogKit::LogFormatted(LogKit::HIGH,"B %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()+GetLX()*cosrot_, GetYMin()+GetLX()*sinrot_,
                         xr+lxr*cosrotr, yr+lxr*sinrotr);
    LogKit::LogFormatted(LogKit::HIGH,"C %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()-GetLY()*sinrot_, GetYMin()+GetLY()*cosrot_,
                         xr -lyr*sinrotr, yr +lyr*cosrotr);
    LogKit::LogFormatted(LogKit::HIGH,"D %18.2f %11.2f    %11.2f %11.2f\n", 
                         GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_, GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_,
                         xr +lxr*cosrotr-lyr*sinrotr, yr +lxr*sinrotr+lyr*cosrotr);
    //
    // Calculate and write the largest possible AREA based on the (dx, dy, angle) given by user.
    //
    // Not implemented...
  }
  int error = 1 - allOk; 
  return error;
}


double
Simbox::getTop(double x, double y) const
{
  double zTop = GetTopSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(zTop))
    zTop = RMISSING;
  return(zTop);
}

double
Simbox::getBot(double x, double y) const
{
  double zBot = GetBotSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(zBot))
    zBot = RMISSING;
  return(zBot);
}

char *
Simbox::getStormHeader(int cubetype, int nx, int ny, int nz, bool flat, bool ascii) const
{
  if(flat == false)
    assert(topName_ != "");
  char * header = new char[500];
  if(ascii == false)
    sprintf(header,"storm_petro_binary\n");
  else
    sprintf(header,"storm_petro_ascii\n");

  sprintf(header,"%s0 %d %f\n",  header, cubetype, RMISSING);
  sprintf(header,"%sFFTGrid\n",header);
  if(flat == false)
    sprintf(header,"%s%f %f %f %f %s %s 0.0 0.0\n", header, GetXMin(), GetLX(), 
    GetYMin(), GetLY(), topName_.c_str(), botName_.c_str());
  else
    sprintf(header,"%s%f %f %f %f 0.0 %f 0.0 0.0\n", header, GetXMin(), GetLX(), 
    GetYMin(), GetLY(), GetLZ());

  sprintf(header,"%s%f %f\n\n", header, GetLZ(), GetAngle()*180/PI);
  sprintf(header,"%s%d %d %d\n", header, nx, ny, nz);
  return(header);
}

//NBNB Ragnar: Drep char * her ved bytte av logkit.
void
Simbox::writeTopBotGrids(const char * topname, const char * botname)
{
  char * tmpName;

  tmpName  = ModelSettings::makeFullFileName(topname);
  topName_ = NRLib2::RemovePath(std::string(tmpName));
  assert(typeid(GetTopSurface()) == typeid(Surface));
  const Surface & wtsurf = dynamic_cast<const Surface &>(GetTopSurface());
  NRLib2::WriteStormBinarySurf(wtsurf, std::string(tmpName));
  delete [] tmpName;

  tmpName  = ModelSettings::makeFullFileName(botname);
  botName_ = NRLib2::RemovePath(std::string(tmpName));
  assert(typeid(GetBotSurface()) == typeid(Surface));
  const Surface & wbsurf = dynamic_cast<const Surface &>(GetBotSurface());
  NRLib2::WriteStormBinarySurf(wbsurf, std::string(tmpName));
  delete [] tmpName;
}

int
Simbox::checkError(double lzLimit, char * errText)
{
  if(status_ == NODEPTH || status_ == EMPTY)
    status_ = EXTERNALERROR; //At this stage, lack of depth is an error

  if(status_ == EXTERNALERROR || status_ == INTERNALERROR)
    //Earlier internal errors are external for this purpose.
    return(EXTERNALERROR);

  if(status_ == NOAREA)
    return(BOXOK);

  if(dz_ < 0)
  {
    double z0, z1 = 0.0;
    double x, y, rx, ry = 0.5f*dy_;
    double lzCur, lzMin = double(1e+30);
    int i,j;
    for(j=0;j<ny_;j++)
    {
      rx = 0.5f*dx_;
      for(i=0;i<nx_;i++)
      {
        x = rx*cosrot_-ry*sinrot_ + GetXMin();
        y = rx*sinrot_+ry*cosrot_ + GetYMin();
        z0 = GetTopSurface().GetZ(x,y);
        z1 = GetBotSurface().GetZ(x,y);
        if(GetTopSurface().IsMissing(z0) == false && GetBotSurface().IsMissing(z1) == false )
        {
          lzCur = z1 - z0;
          if(lzCur < lzMin)
            lzMin = lzCur;
        }
        rx += dx_;
      }
      ry += dy_;
    }
  
    if(lzMin < 0.0)
    {
      status_ = INTERNALERROR;
      sprintf(errText,"At least parts of the Top surface is lower than the base surface. Are surfaces given in wrong order?\n");
    }
    else
    {
      double lzFac = lzMin/GetLZ();
      minRelThick_ = lzFac;
      if(lzFac < lzLimit) 
      {
        status_ = INTERNALERROR;
        sprintf(errText,"Error with top/bottom grids. Minimum thickness should be at least %f times maximum, is %f.\n", lzLimit, lzFac);
      }
      else 
      {
        dz_ = GetLZ()/static_cast<double>(nz_);
      }
    }
  }
  return(status_);
}


int
Simbox::setArea(const SegyGeometry * geometry, char * errText)
{
  double x0  = geometry->getX0();
  double y0  = geometry->getY0();
  double lx  = geometry->getlx();
  double ly  = geometry->getly();
  double rot = geometry->getAngle();
  double dx  = geometry->getDx();
  double dy  = geometry->getDy();

  ILxflag_ = geometry->getILxflag();
  try
  {
    SetDimensions(x0,y0,lx,ly);
  }
  catch (NRLib2::Exception & e)
  {
    sprintf(errText,"%s Could not set x0, y0, lx, and ly. %s",errText,e.what());
    return 1;
  }
  try
  {
    SetAngle(rot);
  }
  catch (NRLib2::Exception & e)
  {
    sprintf(errText,"%s Could not set rotation angle. %s",errText,e.what());
    return 1;
  }
  cosrot_ = cos(rot);
  sinrot_ = sin(rot);
  dx_     = dx;
  dy_     = dy;
  nx_     = int(0.5+lx/dx_);
  ny_     = int(0.5+ly/dy_);
  if(status_ == EMPTY)
    status_ = NODEPTH;
  else if(status_ == NOAREA)
    status_ = BOXOK;

  return 0;
}

void
Simbox::setSeisLines(int * lineParams)
{
  inLine0_    = lineParams[0];
  crossLine0_ = lineParams[1];
  ilStep_     = lineParams[2];
  xlStep_     = lineParams[3];
}

void
Simbox::setDepth(Surface * zref, double zShift, double lz, double dz)
{
  zref->Add(zShift);
  Surface * zBot = new Surface(*zref);
  zBot->Add(lz);
  SetSurfaces(zref,zBot);
  dz_ = dz;
  nz_ = int(0.5+lz/dz_);
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;
}

void
Simbox::setDepth(Surface * z0, Surface * z1, int nz)
{
  SetSurfaces(z0, z1);
  nz_ = nz;
  dz_ = -1;
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;

  constThick_ = false;
}

double
Simbox::getAvgRelThick(void) const
{
  double avgThick = 0.0f;
  for (int i = 0 ; i < nx_ ; i++) {
    for (int j = 0 ; j < ny_ ; j++) {
      avgThick += getRelThick(i, j);
    }
  }
  avgThick /= nx_*ny_; 
  return avgThick;
}

double
Simbox::getRelThick(int i, int j) const
{
  double rx = (static_cast<double>(i) + 0.5)*dx_;
  double ry = (static_cast<double>(j) + 0.5)*dy_;
  double x = rx*cosrot_-ry*sinrot_ + GetXMin();
  double y = rx*sinrot_+ry*cosrot_ + GetYMin();
  return(getRelThick(x, y));
}

double
Simbox::getRelThick(double x, double y) const
{
  double relThick = 1; //Default value to be used outside grid.
  double zTop = GetTopSurface().GetZ(x,y);
  double zBot = GetBotSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false && 
     GetBotSurface().IsMissing(zBot) == false)
    relThick = (zBot-zTop)/GetLZ();
  return(relThick);
}

void Simbox::findIJFromILXL(int IL, int XL, int &i, int &j)const
{

  if(ILxflag_==false)
  {
    i = (XL-crossLine0_)/xlStep_;
    j = (IL-inLine0_)/ilStep_;
  }
  else
  {
    i = (IL-inLine0_)/ilStep_;
    j = (XL-crossLine0_)/xlStep_;
  }
}
