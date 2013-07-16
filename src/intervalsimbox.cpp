/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#include <assert.h>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "nrlib/volume/volume.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/intervalsimbox.h"
#include "src/fftgrid.h"
#include "src/definitions.h"


IntervalSimbox::IntervalSimbox(double             x0,
                               double             y0,
                               const Surface    & z0,
                               double             lx,
                               double             ly,
                               double             lz,
                               double             rot,
                               double             dx,
                               double             dy,
                               double             dz){
  status_      = BOXOK;
  top_name_     = "";
  bot_name_     = "";
  SetDimensions(x0,y0,lx,ly);
  SetAngle(rot);

  Surface z1(z0);
  z1.Add(lz);
  SetSurfaces(z0,z1); //Automatically sets lz correct in this case.

  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  dz_          = dz;
  nx_          = int(0.5+lx/dx_);
  ny_          = int(0.5+ly/dy_);
  nz_          = int(0.5+lz/dz_);
  const_thick_  = true;
  min_rel_thick_ = 1.0;

  inline0_     = -0.5;
  crossline0_  = -0.5;
  xl_step_X_     =  cosrot_/dx_;
  xl_step_Y_     =  sinrot_/dx_;
  il_step_X_     = -sinrot_/dy_;
  il_step_Y_     =  cosrot_/dy_;
}



IntervalSimbox::~IntervalSimbox()
{
  delete top_correlation_surface_;
  delete bot_correlation_surface_;
}

int IntervalSimbox::GetIndex(double          x,
                             double          y,
                             double          z)   const{
  int index = IMISSING;
  int i, j, k;
  GetIndexes(x,y,z,i,j,k);
  if(k != IMISSING && j != IMISSING && i != IMISSING)
    index = int(i+j*nx_+k*nx_*ny_);
  return(index);
}

int
IntervalSimbox::GetClosestZIndex(double  x,
                                 double  y,
                                 double  z){
  int index = IMISSING;
  int i, j, k;
  GetIndexesFull(x,y,z,i,j,k);
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
IntervalSimbox::GetIndexes(double        x,
                           double        y,
                           double        z,
                           int         & x_ind,
                           int         & y_ind,
                           int         & z_ind) const{
  x_ind = IMISSING;
  y_ind = IMISSING;
  z_ind = IMISSING;
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx >= 0 && rx <= GetLX() && ry >= 0 && ry <= GetLY())
  {
    double z_bot, z_top = GetTopSurface().GetZ(x,y);
    if(GetTopSurface().IsMissing(z_top) == false)
    {
      z_bot = GetBotSurface().GetZ(x,y);
      if(GetBotSurface().IsMissing(z_bot) == false &&  z > z_top && z < z_bot)
      {
        x_ind = int(floor(rx/dx_));
        if(x_ind > nx_-1)
          x_ind = nx_-1;
        y_ind = int(floor(ry/dy_));
        if(y_ind > ny_-1)
          y_ind = ny_-1;
        z_ind = int(floor(static_cast<double>(nz_)*(z-z_top)/(z_bot-z_top)));
        //LogKit::LogFormatted(LogKit::Low,"rx,dx,xInd = %.4f %.4f %d   ry,dy,yInd = %.4f %.4f %d    %d\n",rx,dx_,xInd,ry,dy_,yInd,zInd);
      }
    }
  }
}

void
IntervalSimbox::GetIndexes(double        x,
                           double        y,
                           int         & x_ind,
                           int         & y_ind) const{
  x_ind = IMISSING;
  y_ind = IMISSING;
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx > 0 && rx < GetLX() && ry>0 && ry < GetLY())
  {
    x_ind = static_cast<int>(floor(rx/dx_));
    y_ind = static_cast<int>(floor(ry/dy_));
  }
}

void
IntervalSimbox::GetIndexesFull(double    x,
                               double    y,
                               double    z,
                               int     & x_ind,
                               int     & y_ind,
                               int     & z_ind) const{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  x_ind = int(floor(rx/dx_));
  y_ind = int(floor(ry/dy_));
  z_ind = IMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
      z_ind = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
  }
}

void
IntervalSimbox::GetInterpolationIndexes(double         x,
                                        double         y,
                                        double         z,
                                        double       & x_ind,
                                        double       & y_ind,
                                        double       & z_ind) const{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  x_ind = rx/dx_-0.5;
  y_ind = ry/dy_-0.5;
  z_ind = RMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
      z_ind = static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)-0.5;
  }
}

void
IntervalSimbox::GetZInterpolation(double       x,
                                  double       y,
                                  double       z,
                                  int        & index_1,
                                  int        & index_2,
                                  double     & t) const{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  int x_ind = int(floor(rx/dx_));
  int y_ind = int(floor(ry/dy_));
  int zInd2, zInd1;
  index_1 = IMISSING;
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
      index_1 = x_ind+y_ind*nx_+zInd1*nx_*ny_;
      index_2 = x_ind+y_ind*nx_+zInd2*nx_*ny_;
    }
  }
}

void  IntervalSimbox::GetCoord(int             x_ind,
                               int             y_ind,
                               int             z_ind,
                               double        & x,
                               double        & y,
                               double        & z) const{
  GetXYCoord(x_ind, y_ind, x, y);
  GetZCoord(z_ind, x, y, z);
}

void IntervalSimbox::GetXYCoord(int           x_ind,
                                int           y_ind,
                                double      & x,
                                double      & y) const{
  double rx = (static_cast<double>(x_ind) + 0.5)*dx_;
  double ry = (static_cast<double>(y_ind) + 0.5)*dy_;
  x = rx*cosrot_-ry*sinrot_ + GetXMin();
  y = rx*sinrot_+ry*cosrot_ + GetYMin();
}

void IntervalSimbox::GetZCoord(int            z_ind,
                               double         x,
                               double         y,
                               double       & z) const{
  z = RMISSING;
  double zBot, zTop = GetTopSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false)
  {
    zBot = GetBotSurface().GetZ(x,y);
    if(GetBotSurface().IsMissing(zBot) == false)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      z = zTop + (static_cast<double>(z_ind) + 0.5)*dz;
    }
  }
}

void IntervalSimbox::GetMinMaxZ(double      & min_z,
                                double      & max_z) const{
  min_z = GetZMin(nx_,ny_);
  max_z = GetZMax(nx_,ny_);
}

int IntervalSimbox::IsInside(double                x, 
                             double                y) const{
  double rx =  (x-GetXMin())*cosrot_ + (y-GetYMin())*sinrot_;
  double ry = -(x-GetXMin())*sinrot_ + (y-GetYMin())*cosrot_;
  if(rx < 0 || rx > GetLX() || ry<0 || ry > GetLY())
    return(0);
  else
    return(1);
}

int IntervalSimbox::InsideRectangle(const SegyGeometry *  geometry) const{
  double xr   = geometry->GetX0();
  double yr   = geometry->GetY0();
  double rotr = geometry->GetAngle();
  double lxr  = geometry->Getlx();
  double lyr  = geometry->Getly();
  double dxr  = geometry->GetDx();
  double dyr  = geometry->GetDy();

  // check that incoming rectangle is within IntervalSimbox +-0.5 grid cells
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
    rotr+=2*NRLib::Pi;

  if (allOk==0) {
    double seisAzimuth = (-1)*rotr*(180/NRLib::Pi);
    double areaAzimuth = (-1)*GetAngle()*(180/NRLib::Pi);
    if (seisAzimuth < 0) seisAzimuth += 360.0;
    if (areaAzimuth < 0) areaAzimuth += 360.0;
    LogKit::LogFormatted(LogKit::Low,"                        x0            y0           lx         ly     azimuth         dx      dy\n");
    LogKit::LogFormatted(LogKit::Low,"--------------------------------------------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::Low,"Model area:    %11.2f  %11.2f    %11.2f %11.2f    %8.3f    %7.2f %7.2f\n",
                         GetXMin(), GetYMin(), GetLX(), GetLY(), dx_, dy_, areaAzimuth);
    LogKit::LogFormatted(LogKit::Low,"Seismic area:  %11.2f  %11.2f    %10.2f %10.2f    %8.3f    %7.2f %7.2f\n",
                         xr, yr, lxr, lyr, dxr, dyr, seisAzimuth);

    LogKit::LogFormatted(LogKit::High,"\nCorner     XY Area                    XY Seismic\n");
    LogKit::LogFormatted(LogKit::High,"-----------------------------------------------------------\n");
    LogKit::LogFormatted(LogKit::High,"A %18.2f %11.2f    %11.2f %11.2f\n", GetXMin(),GetYMin(), xr,yr);
    LogKit::LogFormatted(LogKit::High,"B %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()+GetLX()*cosrot_, GetYMin()+GetLX()*sinrot_,
                         xr+lxr*cosrotr, yr+lxr*sinrotr);
    LogKit::LogFormatted(LogKit::High,"C %18.2f %11.2f    %11.2f %11.2f\n", GetXMin()-GetLY()*sinrot_, GetYMin()+GetLY()*cosrot_,
                         xr -lyr*sinrotr, yr +lyr*cosrotr);
    LogKit::LogFormatted(LogKit::High,"D %18.2f %11.2f    %11.2f %11.2f\n",
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

double IntervalSimbox::GetTop(int i, 
                              int j) const{
  double x, y;
  GetXYCoord(i,j,x,y);
  double zTop = GetTopSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(zTop))
    zTop = RMISSING;
  return(zTop);
}

double  IntervalSimbox::GetTop(double x, 
                               double y) const{
  double zTop = GetTopSurface().GetZ(x, y);
  if(GetTopSurface().IsMissing(zTop))
    zTop = RMISSING;
  return(zTop);
}

double  IntervalSimbox::GetBot(int  i, 
                               int  j) const{
  double x, y;
  GetXYCoord(i,j,x,y);
  double zBot = GetBotSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(zBot))
    zBot = RMISSING;
  return(zBot);
}

double  IntervalSimbox::GetBot(double   x, 
                               double   y) const
{
  double z_bot = GetBotSurface().GetZ(x, y);
  if(GetBotSurface().IsMissing(z_bot))
    z_bot = RMISSING;
  return(z_bot);
}

std::string IntervalSimbox::GetStormHeader(int  cubetype, 
                                           int  nx, 
                                           int  ny, 
                                           int  nz, 
                                           bool flat, 
                                           bool ascii) const{
  if(flat == false)
    assert(top_name_ != "");
  std::string header;
  if(ascii == false)
    header = "storm_petro_binary\n";
  else
    header = "storm_petro_ascii\n";

  header += "0 "+NRLib::ToString(cubetype) +" "+ NRLib::ToString(RMISSING,6)+"\n";
  header += "FFTGrid\n";
  if(flat == false)
    header += NRLib::ToString(GetXMin(),6) +" "+ NRLib::ToString(GetLX(),6) +" "+ NRLib::ToString(GetYMin(),6) +" "+ NRLib::ToString(GetLY(),6) +" "+ top_name_ +" "+ bot_name_ +" 0.0 0.0\n";
  else
    header += NRLib::ToString(GetXMin(),6) +" "+ NRLib::ToString(GetLX(),6) +" "+ NRLib::ToString(GetYMin(),6) +" "+ NRLib::ToString(GetLY(),6) +" 0.0 "+ NRLib::ToString(GetLZ(),6)+" 0.0 0.0\n";

  header += NRLib::ToString(GetLZ(),6) +" "+ NRLib::ToString(GetAngle()*180/NRLib::Pi,6)+"\n\n";
  header += NRLib::ToString(nx) +" "+ NRLib::ToString(ny) +" "+ NRLib::ToString(nz)+"\n";
  std::string str_header(header);

  return(str_header);
}

void  IntervalSimbox::WriteTopBotGrids(const std::string & top_name,
                                       const std::string & bot_name,
                                       const std::string & subdir,
                                       int                 output_format)
{
  assert(typeid(GetTopSurface()) == typeid(Surface));
  assert(typeid(GetBotSurface()) == typeid(Surface));

  const Surface & wtsurf = dynamic_cast<const Surface &>(GetTopSurface());
  const Surface & wbsurf = dynamic_cast<const Surface &>(GetBotSurface());

  IO::writeSurfaceToFile(wtsurf, top_name, subdir, output_format);
  IO::writeSurfaceToFile(wbsurf, bot_name, subdir, output_format);
}

void    IntervalSimbox::SetCorrelationSurfaces(Surface     * top_corr_surface,
                                               Surface     * bot_corr_surface){
  top_correlation_surface_ = top_corr_surface;
  bot_correlation_surface_ = bot_corr_surface;
}

void
IntervalSimbox::SetTopBotName(const std::string & topname,
                              const std::string & botname,
                              int                 outputFormat)
{
  std::string suffix;
  if ((outputFormat & IO::ASCII) > 0 && (outputFormat & IO::STORM) == 0)
    suffix = IO::SuffixAsciiIrapClassic();
  else
    suffix = IO::SuffixStormBinary();

  top_name_ = IO::getFilePrefix()+topname+suffix;
  bot_name_ = IO::getFilePrefix()+botname+suffix;
}

int IntervalSimbox::CalculateDz(double          lz_limit, 
                                std::string   & err_text){
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
      err_text += "At least parts of the top surface is lower than the base surface. Are surfaces given in wrong order?\n";
    }
    else
    {
      double lzFac = lzMin/GetLZ();
      min_rel_thick_ = lzFac;
      if(lzFac < lz_limit)
      {
        status_ = INTERNALERROR;
        err_text += "Error with top/bottom grids. Minimum thickness should be at least "+NRLib::ToString(lz_limit)+" times maximum, is "+NRLib::ToString(lzFac)+"\n";
      }
      else
      {
        dz_ = GetLZ()/static_cast<double>(nz_);
      }
    }
  }
  return(status_);
}


bool  IntervalSimbox::SetArea(const SegyGeometry  * geometry,
                              std::string         & err_text){
  double x0  = geometry->GetX0();
  double y0  = geometry->GetY0();
  double lx  = geometry->Getlx();
  double ly  = geometry->Getly();
  double rot = geometry->GetAngle();
  double dx  = geometry->GetDx();
  double dy  = geometry->GetDy();

  bool failed = false;

  try
  {
    SetDimensions(x0,y0,lx,ly);
  }
  catch (NRLib::Exception & e)
  {
    err_text += "Could not set x0, y0, lx, and ly.\n";
    err_text += e.what();
    return true; // Failed
  }
  try
  {
    SetAngle(rot);
  }
  catch (NRLib::Exception & e)
  {
    err_text += "Could not set rotation angle.\n";
    err_text += e.what();
    failed = true;
    return true; // Failed
  }
  cosrot_      = cos(rot);
  sinrot_      = sin(rot);
  dx_          = dx;
  dy_          = dy;
  nx_          = static_cast<int>(0.5+lx/dx_);
  ny_          = static_cast<int>(0.5+ly/dy_);

  // In case IL/XL information is not available, we fall back
  //  on the following base case values ...
  inline0_     = -0.5;
  crossline0_  = -0.5;
  il_step_X_     =  cosrot_/dx_;
  il_step_Y_     =  sinrot_/dx_;
  xl_step_X_     = -sinrot_/dy_;
  xl_step_Y_     =  cosrot_/dy_;

  if(status_ == EMPTY)
    status_ = NODEPTH;
  else if(status_ == NOAREA)
    status_ = BOXOK;

  return false; // OK
}

void
IntervalSimbox::SetDepth(const Surface  & z_ref,
                         double           z_shift,
                         double           lz,
                         double           dz,
                         bool             skip_check){
  Surface z_top(z_ref);
  z_top.Add(z_shift);
  Surface zBot(z_top);
  zBot.Add(lz);
  SetSurfaces(z_top,zBot,skip_check);
  dz_ = dz;
  nz_ = int(0.5+lz/dz_);
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;
}

void  IntervalSimbox::SetDepth(const Surface       & z0,
                               const Surface       & z1,
                               int                   nz,
                               bool                  skip_check = false){
  SetSurfaces(z0, z1, skip_check);
  nz_ = nz;
  dz_ = -1;
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;

  const_thick_ = false;
}

void  IntervalSimbox::SetILXL(const SegyGeometry * geometry){
  xl_step_X_ = geometry->GetXLStepX();
  xl_step_Y_ = geometry->GetXLStepY();
  il_step_X_ = geometry->GetILStepX();
  il_step_Y_ = geometry->GetILStepY();

  float x0 = static_cast<float>(GetXMin());
  float y0 = static_cast<float>(GetYMin());
  geometry->FindContILXL(x0, y0, inline0_, crossline0_); //Sets IL0 ,XL0
}

bool  IntervalSimbox::IsAligned(const SegyGeometry * geometry) const{
  double x,y;
  GetXYCoord(0, 0, x, y);
  int IL0, XL0;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL0, XL0);
  GetXYCoord(1, 0, x, y);
  int ILx1, XLx1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), ILx1, XLx1);
  GetXYCoord(0, 1, x, y);
  int ILy1, XLy1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), ILy1, XLy1);
  GetXYCoord(nx_-1, 0, x, y);
  int IL1, XL1;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL1, XL1);
  GetXYCoord(0, ny_-1, x, y);
  int IL2, XL2;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL2, XL2);
  GetXYCoord(nx_-1, ny_-1, x, y);
  int IL3, XL3;
  geometry->FindILXL(static_cast<float>(x), static_cast<float>(y), IL3, XL3);

  int XLdx = XLx1 - XL0;
  int XLdy = XLy1 - XL0;
  int ILdx = ILx1 - IL0;
  int ILdy = ILy1 - IL0;
  if(abs(XLdx*XLdy) > 0 || abs(ILdx*ILdy) > 0)
    return(false); //Moving along one axis lead to change in both il and xl

  int XLndx = XL1 - XL0;
  int XLndy = XL2 - XL0;
  int ILndx = IL1 - IL0;
  int ILndy = IL2 - IL0;
  if(XLndx != (nx_-1)*XLdx || XLndy != (ny_-1)*XLdy ||
     ILndx != (nx_-1)*ILdx || ILndy != (ny_-1)*ILdy)
    return(false); //XL or IL difference at corners not multiple of one-step difference

  if(XL3-XL0 != (nx_-1)*XLdx+(ny_-1)*XLdy || IL3-IL0 != (nx_-1)*ILdx+(ny_-1)*ILdy)
    return(false); //Check final corner, changes failed to match one-step.

  return(true);
}

double  IntervalSimbox::GetAvgRelThick(void) const{
  double avgThick = 0.0f;
  for (int i = 0 ; i < nx_ ; i++) {
    for (int j = 0 ; j < ny_ ; j++) {
      avgThick += GetRelThick(i, j);
    }
  }
  avgThick /= nx_*ny_;
  return avgThick;
}

double  IntervalSimbox::GetRelThick(int i, int j) const {
  double rx = (static_cast<double>(i) + 0.5)*dx_;
  double ry = (static_cast<double>(j) + 0.5)*dy_;
  double x = rx*cosrot_-ry*sinrot_ + GetXMin();
  double y = rx*sinrot_+ry*cosrot_ + GetYMin();
  return(GetRelThick(x, y));
}

double  IntervalSimbox::GetRelThick(double x, double y) const{
  double relThick = 1; //Default value to be used outside grid.
  double zTop = GetTopSurface().GetZ(x,y);
  double zBot = GetBotSurface().GetZ(x,y);
  if(GetTopSurface().IsMissing(zTop) == false &&
     GetBotSurface().IsMissing(zBot) == false)
    relThick = (zBot-zTop)/GetLZ();
  return(relThick);
}

void IntervalSimbox::GetMinAndMaxXY(double  & x_min,
                                    double  & x_max,
                                    double  & y_min,
                                    double  & y_max) const{
  x_min = std::min(GetXMin()+GetLX()*cosrot_, GetXMin());
  x_min = std::min(x_min,GetXMin()-GetLY()*sinrot_);
  x_min = std::min(x_min,GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_);

  x_max = std::max(GetXMin()+GetLX()*cosrot_, GetXMin());
  x_max = std::max(x_max,GetXMin()-GetLY()*sinrot_);
  x_max = std::max(x_max,GetXMin()+GetLX()*cosrot_-GetLY()*sinrot_);

  y_min = std::min(GetYMin()+GetLX()*sinrot_, GetYMin());
  y_min = std::min(y_min,GetYMin()+GetLY()*cosrot_);
  y_min = std::min(y_min,GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_);

  y_max = std::max(GetYMin(),GetYMin()+GetLX()*sinrot_);
  y_max = std::max(y_max,GetYMin()+GetLY()*cosrot_);
  y_max = std::max(y_max,GetYMin()+GetLX()*sinrot_+GetLY()*cosrot_);
}
