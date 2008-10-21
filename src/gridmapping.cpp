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
#include "src/modelsettings.h"
#include "lib/global_def.h"

#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"


GridMapping::GridMapping(const Simbox  * simbox, 
                         ModelSettings * modelSettings, 
                         char          * errText,
                         int             format, 
                         FFTGrid       * velocity)
  : mapping_(NULL),
    simbox_(NULL),
    z0Grid_(NULL),
    z1Grid_(NULL),
    format_(format),
    surfmissing_(0)
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
  makeTimeDepthMapping(mapping_,
                       timeSimbox, 
                       simbox_, 
                       velocity);
}

void 
GridMapping::makeTimeDepthMapping(StormContGrid *& mapping,
                                  const Simbox   * timeSimbox,
                                  const Simbox   * depthSimbox,
                                  FFTGrid        * velocity)
{
  int nx  = depthSimbox->getnx();
  int ny  = depthSimbox->getny();
  int nz  = depthSimbox->getnz();
  mapping = new StormContGrid(*depthSimbox, nx, ny, nz);

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
      (*mapping)(i,j,0) = static_cast<float>(tTop);
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
        (*mapping)(i,j,k) = static_cast<float>(z0 + dz1 + dz2);
      }
    }
  }
}

void 
GridMapping::calculateSurfaceFromVelocity(FFTGrid       * velocity, 
                                          const Simbox  * timeSimbox, 
                                          ModelSettings * modelSettings, 
                                          bool          & failed, 
                                          char          * errText)
{
  if(surfmissing_>0)
  {
    // calculate new surface
    int nx = timeSimbox->getnx();
    int ny = timeSimbox->getny();
    int nz = timeSimbox->getnz();
    double * values = new double[nx*ny];

    for(int i=0 ; i<nx ; i++)
    {
      for(int j=0 ; j<ny ; j++)
      {
        double x, y;
        timeSimbox->getXYCoord(i,j,x,y);
        double sum = 0.0;
        double dt  = (timeSimbox->getBot(x,y)-timeSimbox->getTop(x,y))/(2000.0*static_cast<double>(nz));
        for(int k=0 ; k<nz ; k++)          
          sum+=velocity->getRealValue(i,j,k)*dt;
        if(surfmissing_==1)
          values[i*nx+j] = z1Grid_->GetZ(x,y) - sum;
        else if(surfmissing_==2)
          values[i*nx+j] = z0Grid_->GetZ(x,y) + sum;
      }
    }
    if(surfmissing_==1)
      z0Grid_ = new Surface(z1Grid_->GetXMin(), 
                            z1Grid_->GetYMin(), 
                            z1Grid_->GetLengthX(), 
                            z1Grid_->GetLengthY(), 
                            z1Grid_->GetNI(),
                            z1Grid_->GetNJ(), 
                            (*values));
    else
      z1Grid_ = new Surface(z0Grid_->GetXMin(), 
                            z0Grid_->GetYMin(), 
                            z0Grid_->GetLengthX(), 
                            z0Grid_->GetLengthY(), 
                            z0Grid_->GetNI(),
                            z0Grid_->GetNJ(), 
                            (*values));
    setSimbox(modelSettings, failed, errText,nz);
  }
  surfmissing_ = 0;
}

void 
GridMapping::setDepthSurfaces(ModelFile     * modelFile, 
                              ModelSettings * modelSettings, 
                              bool          & failed, 
                              char          * errText, 
                              int             nz)
{
  surfmissing_ = 0;
  char **surfFile = modelFile->getDepthSurfFile();
  z0Grid_ = NULL;
  if(surfFile[0]==0 && surfFile[1]==0)
  {
    sprintf(errText,"%s Both top and base depth surfaces are missing.", errText);
    failed = 1;
  }

  if(surfFile[0]!=NULL)
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
  else
    surfmissing_ = 1;

  z1Grid_ = NULL;
  if(surfFile[1]!=NULL)
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
  else
    surfmissing_ = 2;

  if(surfmissing_==0)
  {
    setSimbox(modelSettings,failed, errText, nz);
  } 
}

void GridMapping::setSimbox(ModelSettings *modelSettings, bool &failed, char *errText, int nz)
{
  simbox_ = new Simbox();
  simbox_->setDepth(z0Grid_, z1Grid_, nz);
 
  const char * topname = "topdepth.storm";
  const char * botname = "botdepth.storm";
  simbox_->writeTopBotGrids(topname, botname);
  int error = 0;
  const SegyGeometry * areaParams = modelSettings->getAreaParameters(); 
  if (areaParams != NULL)
  {
    simbox_->setArea(areaParams, errText);
    if(error ==1)
    {
      sprintf(errText," %s Problems with definition of depth simbox.",errText);
      failed = 1;
    }
    else
    {
      error = simbox_->checkError(modelSettings->getLzLimit(),errText);
      if(error == Simbox::INTERNALERROR)
      {
        sprintf(errText," %s A problems was encountered for depth output grid\n", errText);
        failed = true;
      }
    }
  }
  else
  {
    LogKit::LogFormatted(LogKit::ERROR," No area available for simbox\n");
    failed = true;
  }
  if (simbox_->status() == Simbox::BOXOK)
  {
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
}
