#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "src/definitions.h"
#include "src/gridmapping.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/io.h"

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


GridMapping::GridMapping()
  : mapping_(NULL),
    simbox_(NULL),
    z0Grid_(NULL),
    z1Grid_(NULL),
    surfaceMode_(NONEGIVEN)
{
}

GridMapping::~GridMapping()
{
  if(simbox_!=NULL)
    delete simbox_;
  simbox_ = NULL;
  if(mapping_!=NULL)
    delete mapping_;
  mapping_ = NULL;
  delete z0Grid_;
  delete z1Grid_;
}

void 
GridMapping::makeTimeTimeMapping(const Simbox * timeCutSimbox)
{
  int nx   = timeCutSimbox->getnx();
  int ny   = timeCutSimbox->getny();
  int nz   = timeCutSimbox->getnz();
  mapping_ = new StormContGrid(*timeCutSimbox, nx, ny, nz); 
  simbox_  = new Simbox(timeCutSimbox);

  for(int i=0;i<nx;i++)
  {
    for(int j=0;j<ny;j++)
    {
      double x,y;
      timeCutSimbox->getXYCoord(i,j,x,y);
      double tTop   = timeCutSimbox->getTop(x,y);
      double tBase  = timeCutSimbox->getBot(x,y); 
      double deltaT = (tBase-tTop)/static_cast<double>(nz);    
      for(int k=0;k<nz;k++)      
        (*mapping_)(i,j,k) = static_cast<float>(tTop + static_cast<double>(k)*deltaT);
    }
  }
}

void 
GridMapping::makeTimeDepthMapping(FFTGrid      * velocity,
                                  const Simbox * timeSimbox)
{
  Simbox * depthSimbox = simbox_; // For readability

  int nx  = depthSimbox->getnx();
  int ny  = depthSimbox->getny();
  int nz  = depthSimbox->getnz();
  mapping_ = new StormContGrid(*depthSimbox, nx, ny, nz);
 // velocity->setAccessMode(FFTGrid::RANDOMACCESS);
  for(int i=0;i<nx;i++)
  {
    for(int j=0;j<ny;j++)
    {
      double x,y;
      depthSimbox->getXYCoord(i,j,x,y);
      double tTop   = timeSimbox->getTop(x,y);
      double tBase  = timeSimbox->getBot(x,y); 
      double zTop   = depthSimbox->getTop(x,y);
      double zBase  = depthSimbox->getBot(x,y); 
      double deltaT = (tBase-tTop)/(static_cast<double>(2000*timeSimbox->getnz()));
      double deltaZ = (zBase-zTop)/static_cast<double>(nz);    
      double sum    = 0.0;
      double sumz   = 0.0;    
      for(int k=0 ; k<timeSimbox->getnz() ; k++)
        sumz +=deltaT*velocity->getRealValue(i,j,k);
      double c = (zBase-zTop)/sumz;
      (*mapping_)(i,j,0) = static_cast<float>(tTop);
      int kk = 0;
      for(int k=1;k<nz;k++)
      {
        double z = k*deltaZ;
        while(sum<z && kk<nz)
        {
          kk++;
          sum+=deltaT*c*velocity->getRealValue(i,j,kk-1);
        }
        double v   = velocity->getRealValue(i,j,kk-1);
        double z0  = tTop;
        double dz1 = 2000.0*static_cast<double>(kk-1)*deltaT;
        double dz2 = 2000.0*(z - sum + deltaT*c*v)/(c*v);
        (*mapping_)(i,j,k) = static_cast<float>(z0 + dz1 + dz2);
      }
    }
  }
 // velocity->endAccess();
}


void
GridMapping::setMappingFromVelocity(FFTGrid * velocity, const Simbox * timeSimbox)
{
  if(simbox_!=NULL)  //Allow this to be called to override old mappings.
  {
    //z0Grid_, z1Grid_ points into simbox. Must copy the relevant before deleting.
    if(surfaceMode_ == TOPGIVEN) {
      z0Grid_ = new Surface(*z0Grid_);
      z1Grid_ = NULL;
    }
    else if(surfaceMode_ == BOTTOMGIVEN) {
      z0Grid_ = NULL;
      z1Grid_ = new Surface(*z1Grid_);
    }
    else {
      z0Grid_ = new Surface(*z0Grid_);
      z1Grid_ = new Surface(*z1Grid_);
    }
    delete simbox_;
    simbox_ = NULL;
  }
  if(mapping_!=NULL)
    delete mapping_;
  mapping_ = NULL;

  int format = velocity->getOutputFormat();
  bool failed = false;
  std::string errText("");
  calculateSurfaceFromVelocity(velocity, timeSimbox);
  setDepthSimbox(timeSimbox, timeSimbox->getnz(), format, failed, errText);
  makeTimeDepthMapping(velocity, timeSimbox);
  if (failed) {
    LogKit::LogFormatted(LogKit::Error,"\n%s\n",errText.c_str());
    exit(1);
  }
}



void 
GridMapping::calculateSurfaceFromVelocity(FFTGrid      * velocity, 
                                          const Simbox * timeSimbox)
{
  if(z0Grid_==NULL || z1Grid_==NULL)
  {
    Surface * isochore;
    if(z0Grid_==NULL)
      isochore = new Surface(*z1Grid_);
    else    
      isochore = new Surface(*z0Grid_);

    //
    // If a constant time surface has been used, it may have only four grid
    // nodes. To handle this situation we use the grid resolution whenever
    // this is larger than the surface resolution.
    //  
    int maxNx = std::max(timeSimbox->getnx(), static_cast<int>(isochore->GetNI()));
    int maxNy = std::max(timeSimbox->getny(), static_cast<int>(isochore->GetNJ())); 
    isochore->Resize(maxNx, maxNy);

    double dx = 0.5*isochore->GetDX();
    double dy = 0.5*isochore->GetDY();

    for(int j=0 ; j<static_cast<int>(isochore->GetNJ()) ; j++) 
    {
      for(int i=0 ; i<static_cast<int>(isochore->GetNI()) ; i++) 
      {
        double x, y;
        isochore->GetXY(i,j,x,y);   
        double tTop  = timeSimbox->getTop(x,y);
        double tBase = timeSimbox->getBot(x,y);
        int    nz    = timeSimbox->getnz();
        double dt    = (tBase - tTop)/(2000.0*static_cast<double>(nz));

        int ii, jj;
        timeSimbox->getIndexes(x,y,ii,jj);

        if(ii!=IMISSING && jj!=IMISSING)
        {
          double sum = 0.0;
          for(int k=0 ; k<nz ; k++)
            sum += velocity->getRealValue(ii,jj,k);
          (*isochore)(i,j) = sum*dt;
        }
        else
        {
          int i1, i2, i3, i4, j1, j2, j3, j4;
          timeSimbox->getIndexes(x+dx, y+dy, i1, j1);
          timeSimbox->getIndexes(x-dx, y-dy, i2, j2);
          timeSimbox->getIndexes(x+dx, y-dy, i3, j3);
          timeSimbox->getIndexes(x-dx, y+dy, i4, j4);

          int n = 0;
          if(i1!=IMISSING && j1!=IMISSING)
            n++;
          if(i2!=IMISSING && j2!=IMISSING)
            n++;
          if(i3!=IMISSING && j3!=IMISSING)
            n++;
          if(i4!=IMISSING && j4!=IMISSING)
            n++;

          if (n==0)
            isochore->SetMissing(i,j);
          else
          {
            double sum = 0.0;
            for(int k=0 ; k<nz ; k++) { 
              if(i1!=IMISSING && j1!=IMISSING)
                sum += velocity->getRealValue(i1,j1,k);
              if(i2!=IMISSING && j2!=IMISSING)
                sum += velocity->getRealValue(i2,j2,k);
              if(i3!=IMISSING && j3!=IMISSING)
                sum += velocity->getRealValue(i3,j3,k);
              if(i4!=IMISSING && j4!=IMISSING)
                sum += velocity->getRealValue(i4,j4,k);
            }
            (*isochore)(i,j) = sum*dt/static_cast<double>(n);
          }
        }
      }
    }

    if(z0Grid_==NULL) 
    {
      z0Grid_ = new Surface(*isochore);
      z0Grid_->Multiply(-1.0);
      z0Grid_->AddNonConform(z1Grid_);
    }
    else 
    {
      z1Grid_ = new Surface(*isochore);
      z1Grid_->AddNonConform(z0Grid_);
    }
    delete isochore;
  }
}

void 
GridMapping::setDepthSurfaces(const std::vector<std::string> & surfFile, 
                              bool                           & failed, 
                              std::string                    & errText)
{
  if(surfFile[0]=="" && surfFile[1]=="")
  {
    errText += "Both top and base depth surfaces are missing.\n";
    failed = 1;
  }
  if(surfFile[0] != "")
  {
    try {
      Surface tmpSurf(surfFile[0]);
      z0Grid_ = new Surface(tmpSurf);
      surfaceMode_ = TOPGIVEN;
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = 1;
    }
  }
  if(surfFile[1] != "")
  {
    try {
      Surface tmpSurf(surfFile[1]);
      z1Grid_ = new Surface(tmpSurf);
      if(surfaceMode_ == TOPGIVEN)
        surfaceMode_ = BOTHGIVEN;
      else
        surfaceMode_ = BOTTOMGIVEN;
    }
    catch (NRLib::Exception & e) {
      errText += e.what();
      failed = 1;
    }
  }
}

void GridMapping::setDepthSimbox(const Simbox * timeSimbox,
                                 int            nz,
                                 int            outputFormat,
                                 bool         & failed,
                                 std::string  & errText)
{
  simbox_ = new Simbox(timeSimbox);
  simbox_->setDepth(*z0Grid_, *z1Grid_, nz);

  std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixDepth();
  std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixDepth();
  simbox_->writeTopBotGrids(topSurf, baseSurf, IO::PathToInversionResults(), outputFormat);

  double dummyLzLimit = 0.0; // The other LzLimit is only for inversion, not depth conversion
  int status = simbox_->calculateDz(dummyLzLimit,errText);
  if(status == Simbox::INTERNALERROR)
  {
    errText += "A problem was encountered for depth output grid.\n";
    failed = true;
  }
  double zmin, zmax;
  simbox_->getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::Low,"\nDepth output interval:\n");
  LogKit::LogFormatted(LogKit::Low,"  True vertical depth   avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       zmin+simbox_->getlz()*simbox_->getAvgRelThick()*0.5,
                       zmin, zmax); 
  LogKit::LogFormatted(LogKit::Low,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
                       simbox_->getlz()*simbox_->getAvgRelThick(),
                       simbox_->getlz()*simbox_->getMinRelThick(),
                       simbox_->getlz());
  LogKit::LogFormatted(LogKit::Low,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
                       simbox_->getdz()*simbox_->getAvgRelThick(),
                       simbox_->getdz(),
                       simbox_->getdz()*simbox_->getMinRelThick());
}
