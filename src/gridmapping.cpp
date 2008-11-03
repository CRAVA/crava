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
#include "src/modelfile.h"
#include "lib/global_def.h"

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


GridMapping::GridMapping()
  : mapping_(NULL),
    simbox_(NULL),
    z0Grid_(NULL),
    z1Grid_(NULL)
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
}

void 
GridMapping::calculateSurfaceFromVelocity(FFTGrid      * velocity, 
                                          const Simbox * timeSimbox)
{
  if(z0Grid_==NULL || z1Grid_==NULL)
  {
  Surface * isochore;
  if(z0Grid_==NULL)
  {
    isochore = new Surface(*z1Grid_);
    isochore->Resize(timeSimbox->getnx(), timeSimbox->getny());
  }
  else    
  {
    isochore = new Surface(*z0Grid_);
    isochore->Resize(timeSimbox->getnx(), timeSimbox->getny());
  }
   
  for(int j=0 ; j<isochore->GetNJ() ; j++) 
  {
    for(int i=0 ; i<isochore->GetNI() ; i++) 
    {
      double x, y;
      isochore->GetXY(i,j,x,y);   
      int    nz    = timeSimbox->getnz();
      double tTop  = timeSimbox->getTop(x,y);
      double tBase = timeSimbox->getBot(x,y);
      double dt    = (tBase - tTop)/(2000.0*static_cast<double>(nz));
      
      int ii, jj, i1, i2, i3, i4, j1, j2, j3, j4;
      timeSimbox->getIndexes(x,y,ii,jj);
      double dx = 0.5*isochore->GetDX();
      double dy = 0.5*isochore->GetDY();
      timeSimbox->getIndexes(x+dx, y+dy, i1, j1);
      timeSimbox->getIndexes(x-dx, y-dy, i2, j2);
      timeSimbox->getIndexes(x+dx, y-dy, i3, j3);
      timeSimbox->getIndexes(x-dx, y+dy, i4, j4);
      double w1 = 0.5;
      double w2 = 0.125;
      int n = 0;
      if(ii==IMISSING || jj==IMISSING)
        n = 4;
      if(i1==IMISSING || j1==IMISSING)
        n++;
      if(i2==IMISSING || j2==IMISSING)
        n++;
      if(i3==IMISSING || j3==IMISSING)
        n++;
      if(i4==IMISSING || j4==IMISSING)
        n++;

      if(ii!=IMISSING &&jj!=IMISSING)
      {
        w1 = w1 +n*0.125*4.0/(8-n);
        w2 = w2 +n*0.125/(8-n);
      }
      else
      {
        w2 = w2+n*0.125/(8-n);
      }
    
     if (n==8)
        isochore->SetMissing(i,j);
     else
      {
        double sum = 0.0;
        for(int k=0 ; k<nz ; k++) { 
          if(ii!=IMISSING && jj!=IMISSING)
            sum += w1*velocity->getRealValue(ii,jj,k)*dt;
          if(i1!=IMISSING && j1!=IMISSING)
            sum += w2*velocity->getRealValue(i1,j1,k)*dt;
          if(i2!=IMISSING && j2!=IMISSING)
            sum += w2*velocity->getRealValue(i2,j2,k)*dt;
          if(i3!=IMISSING && j3!=IMISSING)
            sum += w2*velocity->getRealValue(i3,j3,k)*dt;
          if(i4!=IMISSING && j4!=IMISSING)
            sum += w2*velocity->getRealValue(i4,j4,k)*dt;
        }
        if(z0Grid_==NULL)
        {
          (*isochore)(i,j) = -sum+z1Grid_->GetZ(x,y,false);
        }
        else
          (*isochore)(i,j) = sum+z0Grid_->GetZ(x,y,false);
      }
    }
  }
  if(z0Grid_==NULL) 
  {
    z0Grid_ = new Surface(*isochore);
  }
  else 
  {
    z1Grid_ = new Surface(*isochore);
  }

  delete isochore;
  }
}

void 
GridMapping::setDepthSurfaces(char ** surfFile, 
                              bool  & failed, 
                              char  * errText)
{
  if(surfFile[0]==NULL && surfFile[1]==NULL)
  {
    sprintf(errText,"%s Both top and base depth surfaces are missing.", errText);
    failed = 1;
  }
  if(surfFile[0] != NULL)
  {
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(surfFile[0]);
      z0Grid_ = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
      sprintf(errText,"%s%s",errText,e.what());
      failed = 1;
    }
  }
  if(surfFile[1] != NULL)
  {
    try {
      Surface tmpSurf = NRLib2::ReadStormSurf(surfFile[1]);
      z1Grid_ = new Surface(tmpSurf);
    }
    catch (NRLib2::Exception & e) {
    sprintf(errText,"%s%s",errText,e.what());
    failed = 1;
    }
  }
}

void GridMapping::setDepthSimbox(const Simbox * timeSimbox,
                                 int            nz,
                                 int            outputFormat,
                                 bool         & failed,
                                 char         * errText)
{
  simbox_ = new Simbox(timeSimbox);
  simbox_->setDepth(z0Grid_, z1Grid_, nz);
  simbox_->writeTopBotGrids("Surface_Top_Depth", "Surface_Base_Depth", outputFormat);

  double dummyLzLimit = 0.0; // The other LzLimit is only for inversion, not depth conversion
  int error = simbox_->checkError(dummyLzLimit,errText);
  if(error == Simbox::INTERNALERROR)
    {
    sprintf(errText,"%sA problems was encountered for depth output grid\n", errText);
    failed = true;
  }
  double zmin, zmax;
  simbox_->getMinMaxZ(zmin,zmax);
  LogKit::LogFormatted(LogKit::LOW,"\nDepth output interval:\n");
  LogKit::LogFormatted(LogKit::LOW,"  True vertical depth   avg / min / max    : %7.1f /%7.1f /%7.1f\n",
                       zmin+simbox_->getlz()*simbox_->getAvgRelThick()*0.5,
                       zmin, zmax); 
  LogKit::LogFormatted(LogKit::LOW,"  Interval thickness    avg / min / max    : %7.1f /%7.1f /%7.1f\n", 
                       simbox_->getlz()*simbox_->getAvgRelThick(),
                       simbox_->getlz()*simbox_->getMinRelThick(),
                       simbox_->getlz());
  LogKit::LogFormatted(LogKit::LOW,"  Sampling density      avg / min / max    : %7.2f /%7.2f /%7.2f\n", 
                       simbox_->getdz()*simbox_->getAvgRelThick(),
                       simbox_->getdz(),
                       simbox_->getdz()*simbox_->getMinRelThick());
}
