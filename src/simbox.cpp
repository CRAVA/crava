#include <math.h>
#include <assert.h>

#include "lib/irapgrid.h"
#include "lib/log.h"

#include "src/simbox.h"
#include "src/model.h"

Simbox::Simbox(void)
{
  status_      = EMPTY;
  topName_     = NULL;
  botName_     = NULL;
  z0Grid_      = NULL;
  z1Grid_      = NULL;
  inLine0_     = 0;
  crossLine0_  = 0;
  ilStep_      = 1;
  xlStep_      = 1;
  constThick_  = true;
  minRelThick_ = 1.0;
  dz_          = 0;
  lz_          = 0;
}

Simbox::Simbox(double x0, double y0, irapgrid * z0, double lx, double ly, double lz,
               double rot, double dx, double dy, double dz)
{
  status_      = BOXOK;
  topName_     = NULL;
  botName_     = NULL;
  x0_          = x0;
  y0_          = y0;
  z0Grid_      = z0;
  lx_          = lx;
  ly_          = ly;
  lz_          = lz;
  int error    = 0;
  z1Grid_      = irapgridArithmeticConstant(z0Grid_, lz_, 1, &error);
  rot_         = rot;
  cosrot_      = cos(rot_);
  sinrot_      = sin(rot_);
  dx_          = dx;
  dy_          = dy;
  dz_          = dz;
  nx_          = int(0.5+lx_/dx_);
  ny_          = int(0.5+ly_/dy_);
  nz_          = int(0.5+lz_/dz_);
  inLine0_     = 0;
  crossLine0_  = 0;
  constThick_  = true;
  minRelThick_ = 1.0;
}

Simbox::Simbox(const Simbox *simbox)
{
  status_      = simbox->status_;
  topName_     = NULL;
  botName_     = NULL;
  z0Grid_      = NULL;
  z1Grid_      = NULL;
  x0_          = simbox->x0_;
  y0_          = simbox->y0_;
  lx_          = simbox->lx_;
  ly_          = simbox->ly_;
  lz_          = simbox->lz_;
  rot_         = simbox->rot_;
  cosrot_      = simbox->cosrot_;
  sinrot_      = simbox->sinrot_;
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

  topName_ = new char [MAX_STRING]; 
  botName_ = new char [MAX_STRING]; 

  int failure = 0;
  if(simbox->z0Grid_ != NULL)
    z0Grid_ = irapgridCopyGrid(simbox->z0Grid_,&failure);
  if(simbox->z1Grid_ != NULL)
    z1Grid_ = irapgridCopyGrid(simbox->z1Grid_,&failure);
  if(simbox->topName_ != NULL)
    strcpy(topName_,simbox->topName_);
  if(simbox->botName_ != NULL )
    strcpy(botName_,simbox->botName_);
}   

Simbox::~Simbox()
{
  if(z0Grid_ != NULL)
  {
    freeIrapgrid(z0Grid_);
    free(z0Grid_);
  }
  if(z1Grid_ != NULL)
  {
    freeIrapgrid(z1Grid_);
    free(z1Grid_);
  }
  if(topName_ != NULL)
    delete [] topName_;
  if(botName_ != NULL)
    delete [] botName_;
}

int
Simbox::getIndex(double x, double y, double z) const

{
  int index = IMISSING;
  int i, j, k;
  getIndexes(x,y,z,i,j,k);
  if(k != IMISSING && j != IMISSING && i != IMISSING)
    index = int(i+j*floor(lx_/dx_+0.5)+k*floor((lx_/dx_)*(ly_/dy_)+0.5));
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
  double rx =  (x-x0_)*cosrot_ + (y-y0_)*sinrot_;
  double ry = -(x-x0_)*sinrot_ + (y-y0_)*cosrot_;
  if(rx > 0 && rx < lx_ && ry>0 && ry < ly_)
  {
    int error = 0;
    double zBot, zTop = irapgridGetValue(x, y, z0Grid_, &error);
    if(error == 0 && zTop != z0Grid_->missingcode)
    {
      zBot = irapgridGetValue(x, y, z1Grid_, &error);
      if(error == 0 && zBot != z1Grid_->missingcode &&  z > zTop && z < zBot)
      {
        xInd = int(floor(rx/dx_));
        yInd = int(floor(ry/dy_));
        zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
        //LogKit::writeLog("rx,dx,xInd = %.4f %.4f %d   ry,dy,yInd = %.4f %.4f %d    %d\n",rx,dx_,xInd,ry,dy_,yInd,zInd);
      }
    }
  }
}

void 
Simbox::getIndexesFull(double x, double y, double z, int & xInd, int & yInd, int & zInd) const
{
  double rx =  (x-x0_)*cosrot_ + (y-y0_)*sinrot_;
  double ry = -(x-x0_)*sinrot_ + (y-y0_)*cosrot_;
  xInd = int(floor(rx/dx_));
  yInd = int(floor(ry/dy_));
  zInd = IMISSING;
  int error = 0;
  double zBot, zTop = irapgridGetValue(x, y, z0Grid_, &error);
  if(error == 0 && zTop != z0Grid_->missingcode)
  {
    zBot = irapgridGetValue(x, y, z1Grid_, &error);
    if(error == 0 && zBot != z1Grid_->missingcode)
      zInd = int(floor(static_cast<double>(nz_)*(z-zTop)/(zBot-zTop)));
  }
}

void
Simbox::getZInterpolation(double x, double y, double z, 
                          int & index1, int & index2, double & t) const
{
  double rx =  (x-x0_)*cosrot_ + (y-y0_)*sinrot_;
  double ry = -(x-x0_)*sinrot_ + (y-y0_)*cosrot_;
  int xInd = int(floor(rx/dx_));
  int yInd = int(floor(ry/dy_));
  int zInd2, zInd1;
  index1 = IMISSING;
  int error = 0;
  double zBot, zTop = irapgridGetValue(x, y, z0Grid_, &error);
  if(error == 0 && zTop != z0Grid_->missingcode)
  {
    zBot = irapgridGetValue(x, y, z1Grid_, &error);
    if(error == 0 && zBot != z1Grid_->missingcode)
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
  x = rx*cosrot_-ry*sinrot_ + x0_;
  y = rx*sinrot_+ry*cosrot_ + y0_;
  z = RMISSING;
  int error = 0;
  double zBot, zTop = irapgridGetValue(x, y, z0Grid_, &error);
  if(error == 0 && zTop != z0Grid_->missingcode)
  {
    zBot = irapgridGetValue(x, y, z1Grid_, &error);
    if(error == 0 && zBot != z1Grid_->missingcode)
    {
      double dz = (zBot-zTop)/static_cast<double>(nz_);
      z = zTop + (static_cast<double>(zInd) + 0.5)*dz;
    }
  }
}


void
Simbox::getMinMaxZ(double &minZ, double &maxZ) const
{
  double min, max;
  irapgridMinMaxValue(z0Grid_, &min, &max);
  minZ = min;
  irapgridMinMaxValue(z1Grid_, &min, &max);
  maxZ = max;
}

int
Simbox::isInside(double x, double y) const
{
  double rx = (x-x0_)*cosrot_+(y-y0_)*sinrot_;
  double ry = -(x-x0_)*sinrot_ + (y-y0_)*cosrot_;
  if(rx < 0 || rx > lx_ || ry<0 || ry > ly_)
    return(0);
  else
    return(1);
}

int
Simbox::insideRectangle(double xr, double yr, double rotr, double lxr, double lyr) const
{
  int allOk = 1;
  double cosrotr = cos(rotr);
  double sinrotr = sin(rotr);
  double x = x0_;
  double y = y0_;
  double rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  double ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.01*dx_ || rx > lxr+0.01*dx_ || ry<-0.01*dy_ || ry > lyr+0.01*dy_)
    allOk = 0;

  x = x0_+lx_*cosrot_;
  y = y0_+lx_*sinrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.01*dx_ || rx > lxr+0.01*dx_ || ry<-0.01*dy_ || ry > lyr+0.01*dy_)
    allOk = 0;

  x = x0_-ly_*sinrot_;
  y = y0_+ly_*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.01*dx_ || rx > lxr+0.01*dx_ || ry<-0.01*dy_ || ry > lyr+0.01*dy_)
    allOk = 0;

  x = x0_+lx_*cosrot_-ly_*sinrot_;
  y = y0_+lx_*sinrot_+ly_*cosrot_;
  rx =  (x-xr)*cosrotr + (y-yr)*sinrotr;
  ry = -(x-xr)*sinrotr + (y-yr)*cosrotr;
  if(rx < -0.01*dx_ || rx > lxr+0.01*dx_ || ry<-0.01*dy_ || ry > lyr+0.01*dy_)
    allOk = 0;
  if(allOk == 0)
  {
    LogKit::writeLog("\n             X0         Y0              DeltaX       DeltaY    Angle\n");
    LogKit::writeLog("---------------------------------------------------------------------\n");
    LogKit::writeLog("Area:    %11.2f %11.2f   %11.2f %11.2f   %8.3f\n", x0_, y0_, lx_, ly_, (rot_*180)/PI);
    LogKit::writeLog("Seismic: %11.2f %11.2f   %11.2f %11.2f   %8.3f\n", xr, yr, lxr, lyr, (rotr*180/PI));
    LogKit::writeLog("\nCorner     XY Area                    XY Seismic\n");
    LogKit::writeLog("-----------------------------------------------------------\n");
    LogKit::writeLog("A %18.2f %11.2f    %11.2f %11.2f\n", x0_,y0_, xr,yr);
    LogKit::writeLog("B %18.2f %11.2f    %11.2f %11.2f\n", x0_+lx_*cosrot_, y0_+lx_*sinrot_,
      xr+lxr*cosrotr, yr+lxr*sinrotr);
    LogKit::writeLog("C %18.2f %11.2f    %11.2f %11.2f\n", x0_-ly_*sinrot_, y0_+ly_*cosrot_,
      xr -lyr*sinrotr, yr +lyr*cosrotr);
    LogKit::writeLog("D %18.2f %11.2f    %11.2f %11.2f\n", 
      x0_+lx_*cosrot_-ly_*sinrot_, y0_+lx_*sinrot_+ly_*cosrot_,
      xr +lxr*cosrotr-lyr*sinrotr, yr +lxr*sinrotr+lyr*cosrotr);
  }

  //
  // Calculate and write the largest possible AREA based on the (dx, dy, angle) given by user.
  //
  // Not implemented...

  return(allOk);
}


double
Simbox::getTop(double x, double y) const
{
  int error = 0;
  double zTop = irapgridGetValue(x, y, z0Grid_, &error);
  if(error != 0 || zTop == z0Grid_->missingcode)
    zTop = RMISSING;
  return(zTop);
}

double
Simbox::getBot(double x, double y) const
{
  int error = 0;
  double zBot = irapgridGetValue(x, y, z1Grid_, &error);
  if(error != 0 || zBot == z1Grid_->missingcode)
    zBot = RMISSING;
  return(zBot);
}

char *
Simbox::getStormHeader(int cubetype, int nx, int ny, int nz, bool flat, bool ascii) const
{
  if(flat == false)
    assert(topName_ != NULL);
  char * header = new char[500];
  if(ascii == false)
    sprintf(header,"storm_petro_binary\n");
  else
    sprintf(header,"storm_petro_ascii\n");

  sprintf(header,"%s0 %d %f\n",  header, cubetype, RMISSING);
  sprintf(header,"%sFFTGrid\n",header);
  if(flat == false)
    sprintf(header,"%s%f %f %f %f %s %s 0.0 0.0\n", header, x0_, lx_, 
    y0_, ly_, topName_, botName_);
  else
    sprintf(header,"%s%f %f %f %f 0.0 %f 0.0 0.0\n", header, x0_, lx_, 
    y0_, ly_, lz_);

  sprintf(header,"%s%f %f\n\n", header, lz_, rot_*180/PI);
  sprintf(header,"%s%d %d %d\n", header, nx, ny, nz);
  return(header);
}

void
Simbox::writeTopBotGrids(const char * topname, const char * botname)
{
  char dirsep = '/';
#ifdef _WINDOWS
  dirsep = '\\';
#endif

  char * tmpName = LogKit::makeFullFileName(topname);
  irapgridWriteBin(tmpName, z0Grid_);

  //Strip away path
  int i;
  for(i = strlen(tmpName)-1;i >=0; i--)
    if(tmpName[i] == dirsep)
      break;
  if(topName_ == NULL)
  {
    topName_ = new char[strlen(tmpName)-i+1];
    strcpy(topName_, &(tmpName[i+1]));
    //  botName_ = new char[strlen(topName_)+1];
    //  strcpy(botName_,topName_);
    //  i = strlen(botName_);
    //  botName_[i-12] = 'b';
    //  botName_[i-10] = 't';
  }
  delete [] tmpName;

  tmpName =LogKit::makeFullFileName(botname);
  if(botName_ ==NULL)
    {
      botName_ = new char[strlen(tmpName)-i+1];
    strcpy(botName_, &(tmpName[i+1]));
    }
  irapgridWriteBin(tmpName, z1Grid_);
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
        x = rx*cosrot_-ry*sinrot_ + x0_;
        y = rx*sinrot_+ry*cosrot_ + y0_;
        int error = 0;
        z0 = irapgridGetValue(x, y, z0Grid_, &error);
        if(error == 0)
          z1 = irapgridGetValue(x, y, z1Grid_, &error);
        if(error == 0 && z0 != z0Grid_->missingcode && z1 != z1Grid_->missingcode)
        {
          lzCur = z1 - z0;
          if(lzCur < lzMin)
            lzMin = lzCur;
          if(lzCur > lz_)
            lz_ = lzCur;
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
      double lzFac = lzMin/lz_;
      minRelThick_ = lzFac;
      if(lzFac < lzLimit) 
      {
        status_ = INTERNALERROR;
        sprintf(errText,"Error with top/bottom grids. Minimum thickness should be at least %f times maximum, is %f.\n", lzLimit, lzFac);
      }
      else 
      {
        dz_ = lz_/static_cast<double>(nz_);
      }
    }
  }
  return(status_);
}


void
Simbox::setArea(double x0, double y0, double lx, double ly, double rot, double dx, double dy)
{
  x0_     = x0;
  y0_     = y0;
  lx_     = lx;
  ly_     = ly;
  rot_    = rot;
  cosrot_ = cos(rot_);
  sinrot_ = sin(rot_);
  dx_     = dx;
  dy_     = dy;
  nx_     = int(0.5+lx_/dx_);
  ny_     = int(0.5+ly_/dy_);
  if(status_ == EMPTY)
    status_ = NODEPTH;
  else if(status_ == NOAREA)
    status_ = BOXOK;
//	LogKit::writeDebugLog(1,"%f %f %f %f %d %d\n",lx_,ly_,dx_,dy_,nx_,ny_);
}


void
Simbox::setDepth(irapgrid * zref, double zShift, double lz, double dz)
{
  int error = 0;
 if(z0Grid_ != NULL)
  {
    freeIrapgrid(z0Grid_);
    free(z0Grid_);
  }
  if(z1Grid_ != NULL)
  {
    freeIrapgrid(z1Grid_);
    free(z1Grid_);
  }
  z0Grid_ = irapgridArithmeticConstant(zref, zShift, 1, &error);
  z1Grid_ = irapgridArithmeticConstant(zref, zShift+lz, 1, &error);
  lz_ = lz;
  dz_ = dz;
  nz_ = int(0.5+lz_/dz_);
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;
}


void
Simbox::setDepth(irapgrid * z0, irapgrid * z1, int nz)
{
 
 if(z0Grid_ != NULL)
  {
    freeIrapgrid(z0Grid_);
    free(z0Grid_);
  }
  if(z1Grid_ != NULL)
  {
    freeIrapgrid(z1Grid_);
    free(z1Grid_);
  }
  z0Grid_ = z0;
  z1Grid_ = z1;
  nz_ = nz;
  dz_ = -1;
  if(status_ == EMPTY)
    status_ = NOAREA;
  else if(status_ == NODEPTH)
    status_ = BOXOK;

  constThick_ = false;
}


void
Simbox::setSeisLines(int il0, int xl0, int ilStep, int xlStep)
{
  inLine0_ = il0;
  crossLine0_ = xl0;
  ilStep_ = ilStep;
  xlStep_ = xlStep;
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
  double x = rx*cosrot_-ry*sinrot_ + x0_;
  double y = rx*sinrot_+ry*cosrot_ + y0_;
  return(getRelThick(x, y));
}

double
Simbox::getRelThick(double x, double y) const
{
  int error = 0;
  double zBot = 0.0f;
  double zTop = irapgridGetValue(x, y, z0Grid_, &error);
  double relThick = 1; //Default value to be used outside grid.
  if(error == 0)
    zBot = irapgridGetValue(x, y, z1Grid_, &error);
  if(error == 0 && zTop != z0Grid_->missingcode && 
    zBot != z1Grid_->missingcode)  
    relThick = (zBot-zTop)/lz_;
  return(relThick);
}
