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


GridMapping::GridMapping(const Simbox *simbox, ModelFile *modelFile, ModelSettings *modelSettings, bool depthmode, bool&failed,
                         char *errText,int format, FFTGrid *velocity)
{
  depthmode_ = depthmode;
  format_ = format;
  failed = false;
  if(depthmode_==1)
    simbox_ = NULL;
  else
    simbox_ = new Simbox(simbox);

  int nz = simbox->getnz();
  if(depthmode==1)
    setSurfaces(modelFile, modelSettings, failed, errText, nz); // If two surfaces are read from file, the simbox is completed in this routine.
  
  if(velocity!=NULL && surfmissing_>0 && depthmode_ ==1)
    calculateSurfaceFromVelocity(velocity, simbox, modelSettings, failed, errText); // simbox is set in this routine
 
  mapping_ = NULL;
  if(velocity!=NULL || depthmode_==0)
    makeMapping(velocity,simbox);
}
GridMapping::~GridMapping()
{
  if(simbox_!=NULL)
    delete simbox_;
  simbox_ = 0;
  if(mapping_!=NULL)
    delete mapping_;
  mapping_ = 0;
}

void GridMapping::makeMapping(FFTGrid *velocity, const Simbox *timeSimbox)
{
  int nx = simbox_->getnx();
  int ny = simbox_->getny();
  int nz = simbox_->getnz();
  mapping_ = new StormContGrid(*simbox_, nx, ny, nz);
  
  int i,j,k;
  if(velocity==NULL) // time-time mapping
  {
    double deltaz, x, y;
    
    for(i=0;i<nx;i++)
    {
      x = simbox_->getx0()+i*simbox_->getdx();
      for(j=0;j<ny;j++)
      {
        y = simbox_->gety0()+j*simbox_->getdy();
        deltaz = (simbox_->getBot(x,y)-simbox_->getTop(x,y))/nz;    
        for(k=0;k<nz;k++)      
          (*mapping_)(i,j,k) = static_cast<float> (simbox_->getTop(x,y)+k*deltaz);
      }
    }

  }
  else // time depth mapping
  {
    int kk;
    double c, sumz;
    double deltat,deltaz,sum,x,y,z;
    float value;
    for(i=0;i<nx;i++)
    {
      x = simbox_->getx0()+i*simbox_->getdx();
      for(j=0;j<ny;j++)
      {
        y = simbox_->gety0()+j*simbox_->getdy();
        deltat = (timeSimbox->getBot(x,y)-timeSimbox->getTop(x,y))/(2000*timeSimbox->getnz());
        deltaz = (simbox_->getBot(x,y)-simbox_->getTop(x,y))/nz;    
        sumz = 0.0;    
        for(k=0;k<timeSimbox->getnz();k++)
          sumz +=deltat*velocity->getRealValue(i,j,k);
        c = (simbox_->getBot(x,y)-simbox_->getTop(x,y))/sumz;
        (*mapping_)(i,j,0) = static_cast<float> (timeSimbox->getTop(x,y));
        kk = 0;
        sum = 0;
        for(k=1;k<nz;k++)
        {
          z = k*deltaz;
          while(sum<z && kk<nz)
          {
            kk++;
            sum+=deltat*c*velocity->getRealValue(i,j,kk-1);
          }
          value = float(timeSimbox->getTop(x,y)+2000*(kk-1)*deltat+2000*(z-sum+deltat*c*velocity->getRealValue(i,j,kk-1))/(c*velocity->getRealValue(i,j,kk-1)));
          (*mapping_)(i,j,k) = value;
        }
      }
    }

  }

}

void GridMapping::calculateSurfaceFromVelocity(FFTGrid *velocity, const Simbox *simbox, ModelSettings *modelSettings, bool &failed, char *errText)
{
    
  if(surfmissing_>0 && velocity!=NULL)
  {
    // calculate new surface
    int nz = simbox->getnz();
    int nx,ny;
    nx = simbox->getnx();
    ny = simbox->getny();
    double * values = new double[nx*ny];
    double dt, sum, x, y;
    int i,j,k;
    for(i=0;i<nx;i++)
    {
      x = simbox->getx0()+i*simbox->getdx();
      for(j=0;j<ny;j++)
      {
        y = simbox->gety0()+j*simbox->getdy();
        sum = 0;
        dt = (simbox->getBot(x,y)-simbox->getTop(x,y))/(2000.0*simbox->getnz());
        for(k=0;k<simbox->getnz();k++)          
          sum+=velocity->getRealValue(i,j,k)*dt;
        if(surfmissing_==2)
            values[i*nx+j] = z0Grid_->GetZ(x,y)+ sum;
        else if(surfmissing_==1)
             values[i*nx+j] = z1Grid_->GetZ(x,y)-sum;
      }
    }
    
   
    if(surfmissing_==2)
    {
  
        z1Grid_ = new Surface(z0Grid_->GetXMin(), z0Grid_->GetYMin(), 
                                     z0Grid_->GetLengthX(), z0Grid_->GetLengthY(), 
                                     z0Grid_->GetNI(),z0Grid_->GetNJ(), (*values));
    }
    else
    {
        z0Grid_ = new Surface(z1Grid_->GetXMin(), z1Grid_->GetYMin(), 
                                     z1Grid_->GetLengthX(), z1Grid_->GetLengthY(), 
                                     z1Grid_->GetNI(),z1Grid_->GetNJ(), (*values));
    }
    setSimbox(modelSettings, failed, errText,nz);
  }
  surfmissing_ = 0;
  
  
}

void GridMapping::setSurfaces(ModelFile *modelFile, ModelSettings *modelSettings, bool &failed, char *errText, int nz)
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
                             zmin+simbox_->getlz()*simbox_->getAvgRelThick(),
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
