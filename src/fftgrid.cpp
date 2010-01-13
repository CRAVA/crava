#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <string>

#include "lib/global_def.h"
#include "nrlib/segy/segy.hpp"
#include "lib/lib_misc.h"
#include "lib/random.h"
#include "lib/timekit.hpp"

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/exception/exception.hpp"

#include "src/fftgrid.h"
#include "src/corr.h"
#include "src/simbox.h"
#include "src/model.h"
#include "src/timings.h"
#include "src/definitions.h"
#include "src/gridmapping.h"
#include "src/io.h"


FFTGrid::FFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp)
{
  cubetype_       = CTMISSING;
  theta_          = RMISSING;
  nx_             = nx;
  ny_             = ny;
  nz_             = nz;
  nxp_            = nxp;
  nyp_            = nyp;
  nzp_            = nzp;
  scale_          = 1.0;

  cnxp_           = nxp_/2+1;
  rnxp_           = 2*(cnxp_);   

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = 0; 
  counterForSet_  = 0;
  istransformed_  = false;
  rvalue_         = NULL;
  add_            = true;
 // nGrids_        += 1;
  // index= i+rnxp_*j+k*rnxp_*nyp_;
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 

 
}

FFTGrid::FFTGrid(FFTGrid  * fftGrid)
{
  cubetype_       = fftGrid->cubetype_;
  theta_          = fftGrid->theta_;
  nx_             = fftGrid->nx_;
  ny_             = fftGrid->ny_;
  nz_             = fftGrid->nz_;
  nxp_            = fftGrid->nxp_;
  nyp_            = fftGrid->nyp_;
  nzp_            = fftGrid->nzp_;
  scale_          = fftGrid->scale_;

  cnxp_           = nxp_/2+1;
  rnxp_           = 2*(cnxp_);   

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = fftGrid->getCounterForGet(); 
  counterForSet_  = fftGrid->getCounterForSet();
  istransformed_  = false;
  add_ = fftGrid->add_;
  
  createRealGrid(add_);
  for(int k=0;k<nzp_;k++)
    for(int j=0;j<nyp_;j++)
      for(int i=0;i<rnxp_;i++)   
      {
        float value = fftGrid->getNextReal();
        setNextReal(value);
      }// k,j,i

 
}

FFTGrid::~FFTGrid()
{
if (rvalue_!=NULL)
{
  if(add_==true)
    nGrids_ = nGrids_ - 1;
  fftw_free(rvalue_);
  LogKit::LogFormatted(LogKit::DEBUGLOW,"\nFFTGrid Destructor: nGrids_ = %d",nGrids_);
}
}

void
FFTGrid::fillInFromSegY(SegY* segy, Simbox *simbox)
{
  assert(cubetype_  !=  CTMISSING);

  createRealGrid();

  int i,j,k,refi,refj,refk;
  float distx,disty,distz,mult;
  double x,y,z;
  float* meanvalue = NULL;
  bool  isParameter=false;

  if(cubetype_==PARAMETER)  isParameter=true;

  int outMode = SegY::MISSING;
  if(cubetype_ == DATA)
    outMode = SegY::ZERO;
  else if(cubetype_ == PARAMETER)
    outMode = SegY::CLOSEST;

  fftw_real value  = 0.0;
  float val1,val2;
  if(isParameter)
  { 
    meanvalue= static_cast<float*>(fftw_malloc(sizeof(float)*nyp_*nxp_));

    for(j=0;j<nyp_;j++)
      for(i=0;i<nxp_;i++)   
      {
        refi = getXSimboxIndex(i);
        refj = getYSimboxIndex(j);
        refk = 0;
        simbox->getCoord(refi,refj,refk,x,y,z);
        val1 = segy->GetValue(x,y,z, outMode);
        refk = nz_-1;
        simbox->getCoord(refi,refj,refk,x,y,z);
        val2 = segy->GetValue(x,y,z, outMode);
        meanvalue[i+j*nxp_] = static_cast<float>((val1+val2)/2.0);
      }
  } // endif

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  LogKit::LogFormatted(LogKit::LOW,"\nResampling seismic data into %dx%dx%d grid:",nxp_,nyp_,nzp_);
  setAccessMode(WRITE);

  int monitorSize = MAXIM(1,int(nyp_*nzp_*0.02));
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |  ");
  printf("\n  ^");

  for( k = 0; k < nzp_; k++)
  {
    for( j = 0; j < nyp_; j++)
    {
      for( i = 0; i < rnxp_; i++)   
      {
        if(i<nxp_)
        {
          refi   = getXSimboxIndex(i);
          refj   = getYSimboxIndex(j);
          refk   = getZSimboxIndex(k);
          simbox->getCoord(refi,refj,refk,x,y,z);
          distx  = getDistToBoundary(i,nx_,nxp_);
          disty  = getDistToBoundary(j,ny_,nyp_);
          distz  = getDistToBoundary(k,nz_,nzp_);         
          mult   = float(pow(MAXIM(1.0-distx*distx-disty*disty-distz*distz,0.0),3));
          value=segy->GetValue(x,y,z, outMode );
          if(isParameter)
            value=float ( ((mult*value+(1.0-mult)*meanvalue[i+j*nxp_])) );
          else
            value*=mult;
        }
        else
          value=RMISSING;        

        setNextReal(value);

      } //for k,j,i
      if ((nyp_*k+j+1)%monitorSize == 0) 
      { 
        printf("^");
        fflush(stdout);
      }
    }
  }
  LogKit::LogFormatted(LogKit::LOW,"\n");
  endAccess();
  if(isParameter) fftw_free(meanvalue);

  Timings::setTimeResamplingSeismic(wall,cpu);
}

void
FFTGrid::fillInFromStorm(Simbox            * actSimBox,
                         StormContGrid     * grid, 
                         const std::string & parName, 
                         bool                scale)
{
  assert(cubetype_ != CTMISSING);
  createRealGrid();
  float scalevert, scalehor;
  if(scale==false)
  {
    scalevert = 1.0;
    scalehor = 1.0;
  }
  else //from sgri file
  {
    LogKit::LogFormatted(LogKit::LOW,"Sgri file read. Rescaling z axis from s to ms, x and y from km to m. \n");
    scalevert = 1000.0;
    scalehor  = 1000.0;
  }
  int i,j,k,refi,refj,refk;
  float distx,disty,distz,mult;
  float* meanvalue;

  //assert(cubetype_==PARAMETER);

  fftw_real value  = 0.0;
  
  double x,y,z;
  double val;

  meanvalue= static_cast<float*>(fftw_malloc(sizeof(float)*nyp_*nxp_));

  for(j=0;j<nyp_;j++)
    for(i=0;i<nxp_;i++)   
    {
      refi   = getXSimboxIndex(i);
      refj   = getYSimboxIndex(j);
      actSimBox->getCoord(refi, refj, 0, x, y, z);
      val=grid->GetValueZInterpolated(x*scalehor,y*scalehor,z*scalevert);
      actSimBox->getCoord(refi, refj, nz_-1, x, y, z);
 
        if(val != RMISSING)
   
          val = (grid->GetValueZInterpolated(x*scalehor,y*scalehor,z*scalevert)+val)/2.0;
        else     
          val = grid->GetValueZInterpolated(x*scalehor,y*scalehor,z*scalevert);
      meanvalue[i+j*nxp_] = static_cast<float>(val);
    }

  LogKit::LogFormatted(LogKit::LOW,"\nResampling %s into %dx%dx%d grid:\n",parName.c_str(),nxp_,nyp_,nzp_);
  setAccessMode(WRITE);
  
  int monitorSize = MAXIM(1,int(nyp_*nzp_*0.02));
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |  ");
  printf("\n  ^");
  
  for( k = 0; k < nzp_; k++)
  {
    for( j = 0; j < nyp_; j++)
    {
      for( i = 0; i < rnxp_; i++)   
      {
        refi   = getXSimboxIndex(i);
        refj   = getYSimboxIndex(j);
        refk   = getZSimboxIndex(k);
        distx  = getDistToBoundary(i,nx_,nxp_);
        disty  = getDistToBoundary(j,ny_,nyp_);
        distz  = getDistToBoundary(k,nz_,nzp_);         
        mult   = float(pow(MAXIM(1.0-distx*distx-disty*disty-distz*distz,0.0),3));
        
        if(i<nxp_)  // computes the index reference from the cube puts it in value
        {
          actSimBox->getCoord(refi, refj, refk, x, y, z);
          value = static_cast<fftw_real>(grid->GetValueZInterpolated(x*scalehor,y*scalehor,z*scalevert));
          if(value==RMISSING)
            value = 0;
          value=static_cast<float>( ((mult*value+(1.0-mult)*meanvalue[i+j*nxp_])) );
        }
        else
          value=RMISSING;        
        
        setNextReal(value);
      } //for k,j,i
      if ((nyp_*k+j+1)%monitorSize == 0) 
      { 
        printf("^");
        fflush(stdout);
      }
    }
  }
  LogKit::LogFormatted(LogKit::LOW,"\n");
  endAccess();
  fftw_free(meanvalue);
}

void
FFTGrid::fillInFromRealFFTGrid(FFTGrid& fftGrid)
{
  assert(cubetype_  !=  CTMISSING);
  createRealGrid();
  int i,j,k,refi,refj,refk;
  float distx,disty,distz,mult;
  float* meanvalue;

  assert(cubetype_==PARAMETER);

  fftw_real value  = 0.0;

  meanvalue= static_cast<float*>(fftw_malloc(sizeof(float)*nyp_*nxp_));

  for(j=0;j<nyp_;j++)
    for(i=0;i<nxp_;i++)   
    {
      refi   = getXSimboxIndex(i);
      refj   = getYSimboxIndex(j);
      float val = fftGrid.getRealValue(refi, refj, 0);
      float val2 = fftGrid.getRealValue(refi, refj, nz_-1);
      meanvalue[i+j*nxp_] = (val == RMISSING || val2 == RMISSING ? 0.f : (val+val2)*0.5f);
    }

    setAccessMode(WRITE);
    for( k = 0; k < nzp_; k++)
    {
      for( j = 0; j < nyp_; j++)
        for( i = 0; i < rnxp_; i++)   
        {
          refi   = getXSimboxIndex(i);
          refj   = getYSimboxIndex(j);
          refk   = k<nz_ ? k : k<(nz_+nzp_)/2 ? nz_-1 : 0;
          distx  = getDistToBoundary(i,nx_,nxp_);
          disty  = getDistToBoundary(j,ny_,nyp_);
          distz  = getDistToBoundary(k,nz_,nzp_);         
          mult   = float(pow(MAXIM(1.0-distx*distx-disty*disty-distz*distz,0.0),3));

          if(i<nxp_)  // computes the index reference from the cube puts it in value
          {
            float te = fftGrid.getRealValue(refi, refj, refk);
            value = te; 
            value=float ( ((mult*value+(1.0-mult)*meanvalue[i+j*nxp_])) );
            te = 0;
          }
          else
            value=RMISSING;        

          setNextReal(value);
        } //for k,j,i
    }

    endAccess();
    fftw_free(meanvalue);

}


void
FFTGrid::fillInConstant(float value)
{  
  createRealGrid();
  int i,j,k;
  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++) 
      {
        if(i<nxp_)
          setNextReal(value);
        else
          setNextReal(RMISSING);
      }  

      endAccess();
}

void
FFTGrid::fillInTest(float value1,float value2)
{  
  createRealGrid();
  int i,j,k;
  float value; 
  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++) 
      {
        if(k > nz_/2 && k < (nz_+nzp_)/2) value = value2; else value = value1;
        if(i<nxp_)
          setNextReal(value);
        else
          setNextReal(RMISSING);
      }  
      endAccess();
}

void
FFTGrid::fillInFromArray(float *value) //NB works only for padding size up to nxp=2*nx, nyp=2*ny, nzp=2*nz
{
  createRealGrid();
  cubetype_ = PARAMETER;
  int i,j,k, ii,jj,kk, iii, jjj, kkk;
  setAccessMode(WRITE);
  kkk = 1;
  for( k = 0; k < nzp_; k++)
  {
    jjj = 1;
    if(k<nz_)
      kk = k;
    else
    {
      kk = k-kkk;
      kkk++;
    }

    for( j = 0; j < nyp_; j++)
    {
      iii = 1;
      if(j<ny_)
        jj = j;
      else
      {
        jj = j-jjj;
        jjj++;
      }
      for( i = 0; i < rnxp_; i++) 
      {
        if(i<nx_)
          ii = i;
        else
        {
          ii = i-iii;
          iii++;
        }
        if(i<nxp_)
          setNextReal(value[ii+jj*nx_+kk*nx_*ny_]);
        else
          setNextReal(RMISSING);
      }  
    }
  }

  endAccess();
}


fftw_real*
FFTGrid::fillInParamCorr(Corr* corr,int minIntFq, float gradI, float gradJ)
{
  assert(corr->getnx() == nxp_);
  assert(corr->getny() == nyp_);
  assert(istransformed_== false);

  int i,j,k,baseK,cycleI,cycleJ;
  float value,subK;
  fftw_real* circCorrT;
  circCorrT = reinterpret_cast<fftw_real*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)));
  computeCircCorrT(corr,circCorrT);
  makeCircCorrTPosDef(circCorrT, minIntFq);

  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++)   
      {    
        if(i < nxp_)  // computes the index reference from the cube puts it in value
        {
          cycleI = i-nxp_;
          if(i < -cycleI)
            cycleI = i;
          cycleJ = j-nyp_;
          if(j < -cycleJ)
            cycleJ = j;

          subK  =  k+cycleI*gradI+cycleJ*gradJ; //Subtract to counter rotation.
          baseK =  int(floor(subK));
          subK  -= baseK;
          while(baseK < 0)
            baseK += nzp_;       //Use cyclicity
          while(baseK >= nzp_)
            baseK -= nzp_;       //Use cyclicity
          value = (1-subK)*circCorrT[baseK];
          if(baseK != nzp_-1)
            value += subK*circCorrT[baseK+1];
          else
            value += subK*circCorrT[0];
          value *= float( (*(corr->getPriorCorrXY()))(i+nxp_*j));
        }
        else
          value = RMISSING;        
        
        setNextReal(value);
      } //for k,j,i
  
  endAccess();
  return circCorrT;//fftw_free(circCorrT);
}


void
FFTGrid::fillInErrCorr(Corr* parCorr, float gradI, float gradJ)
{
  // Note:  this contain the latteral correlation and the multiplyers for the 
  // time correlation. The time correlation is further adjusted by the wavelet
  // and the derivative of the wavelet.
  // the angular correlation is given by a functional Expression elsewhere

  assert(parCorr->getnx() == nxp_);
  assert(parCorr->getny() == nyp_);
  assert(istransformed_== false);

  int i,j,k,baseK,cycleI,cycleJ;
  float value,subK;
  int range = 1;
  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++)   
      {    
        cycleI = i-nxp_;
        if(i < -cycleI)
          cycleI = i;
        cycleJ = j-nyp_;
        if(j < -cycleJ)
          cycleJ = j;

        subK  =  k+cycleI*gradI+cycleJ*gradJ;
        if(i < nxp_) {
          if(fabs(subK) < range*1.0f || fabs(subK-nzp_) < range*1.0f) {
            baseK =  int(subK);
            subK  -= baseK;
            while(baseK < -range)
              baseK += nzp_;       //Use cyclicity
            while(baseK >= range)
              baseK -= nzp_;       //Use cyclicity

            if(baseK < 0) {
              baseK = -1-baseK;
              subK  = 1-subK;
            }
            value = (1-subK)* float( (*(parCorr->getPriorCorrXY()))(i+nxp_*j) );
          }
          else
            value = 0;
        }
        else
          value = RMISSING;        

        setNextReal(value);
      } //for k,j,i
      endAccess();
}


void
FFTGrid::fillInComplexNoise(RandomGen * ranGen)
{
  assert(ranGen);
  istransformed_ = true;
  int i;
  cubetype_=PARAMETER;
  float std = float(1/sqrt(2.0));
  for(i=0;i<csize_;i++)
  {
    //if(xind == 0 || xind == nx-1 && nx is even)
    if(((i % cnxp_) == 0) || (((i % cnxp_) == cnxp_-1) && ((i % 2) == 1)))
    {
      int xshift = i % cnxp_;
      int jkind  = (i-xshift)/cnxp_;
      int jind   = jkind % nyp_;       //Index j along y-direction
      int kind   = int(jkind/nyp_);    //Index k along z-direction
      int jccind, kccind, jkccind;     //Indexes for complex conjugated
      if(jind == 0)
        jccind = 0;
      else
        jccind = nyp_-jind;
      if(kind == 0)
        kccind = 0;
      else
        kccind = nzp_-kind;
      jkccind = jccind+kccind*nyp_;
      if(jkccind == jkind)             //Number is its own cc, i. e. real
      {
        cvalue_[i].re = float(ranGen->rnorm01());
        cvalue_[i].im = 0;
      }
      else if(jkccind > jkind)         //Have not simulated cc yet.
      {
        cvalue_[i].re = float(std*ranGen->rnorm01());
        cvalue_[i].im = float(std*ranGen->rnorm01());
      }
      else                             //Look up cc value
      {
        int cci = jkccind*cnxp_+xshift;
        cvalue_[i].re = cvalue_[cci].re;
        cvalue_[i].im = -cvalue_[cci].im;
      }
    }
    else
    {
      cvalue_[i].re = float(std*ranGen->rnorm01());
      cvalue_[i].im = float(std*ranGen->rnorm01());
    }
  }
}

void
FFTGrid::createRealGrid(bool add)
{
  //  long int timestart, timeend;
  //  time(&timestart);
  istransformed_=false;
  rvalue_    = static_cast<fftw_real*>(fftw_malloc(rsize_ * sizeof(fftw_real)));
  cvalue_    = reinterpret_cast<fftw_complex*>(rvalue_); //
  counterForGet_  = 0; 
  counterForSet_  = 0;
  add_ = add;
  if(add==true)
    nGrids_        += 1;
 // LogKit::LogFormatted(LogKit::ERROR,"\nFFTGrid createRealGrid : nGrids = %d    maxGrids = %d\n",nGrids_,maxAllowedGrids_);
  if (nGrids_ > maxAllowedGrids_) {
    std::string text;
    text += "\n\nERROR in FFTGrid createRealGrid. You have allocated to many FFTGrids. The fix";
    text += "\nis to increase the nGrids variable calculated in Model::checkAvailableMemory().\n";
    text += "\nDo you REALLY need to allocate more grids?\n";
    text += "\nAre there no grids that can be released?\n";
    LogKit::LogFormatted(LogKit::ERROR, text);
    exit(1);
  }    
  maxAllocatedGrids_ = std::max(nGrids_, maxAllocatedGrids_);

  //  time(&timeend);
  //  LogKit::LogFormatted(LogKit::LOW,"\nReal grid created in %ld seconds.\n",timeend-timestart);
}

void
FFTGrid::createComplexGrid()
{
  //  long int timestart, timeend;
  //  time(&timestart);
  istransformed_  = true;
  rvalue_         = static_cast<fftw_real*>(fftw_malloc(rsize_ * sizeof(fftw_real)));
  cvalue_         = reinterpret_cast<fftw_complex*>(rvalue_); //
  counterForGet_  = 0; 
  counterForSet_  = 0;
  nGrids_        += 1;
 // LogKit::LogFormatted(LogKit::ERROR,"\nFFTGrid createComplexGrid : nGrids = %d    maxGrids = %d\n",nGrids_,maxAllowedGrids_);
  if (nGrids_ > maxAllowedGrids_) {
    std::string text;
    text += "\n\nERROR in FFTGrid createComplexGrid. You have allocated to many FFTGrids. The fix";
    text += "\nis to increase the nGrids variable calculated in Model::checkAvailableMemory().\n";
    text += "\nDo you REALLY need to allocate more grids?\n";
    text += "\nAre there no grids that can be released?\n";
    LogKit::LogFormatted(LogKit::ERROR, text);
    exit(1);
  }    
  maxAllocatedGrids_ = std::max(nGrids_, maxAllocatedGrids_);
  //  time(&timeend);
  //  LogKit::LogFormatted(LogKit::LOW,"\nComplex grid created in %ld seconds.\n",timeend-timestart);
}

int 
FFTGrid::getFillNumber(int i, int n, int np )
{
  //  for the series                 i = 0,1,2,3,4,5,6,7 
  //  GetFillNumber(i, 5 , 8)  returns   0,1,2,3,4,4,1,0 (cut middle, i.e 3,2)
  //  GetFillNumber(i, 4 , 8)  returns   0,1,2,3,3,2,1,0 (copy)
  //  GetFillNumber(i, 3 , 8)  returns   0,1,2,2,1,1,1,0 (drag middle out, i.e. 1)

  int refi     =  0; 
  int BeloWnp, AbovEn;

  if (i< np)  
  {
    if (i<n)
      // then it is in the main cube
      refi  =  i;
    else
    { 
      // Get cyclic extention 
      BeloWnp  = np-i-1;   
      AbovEn   = i-n+1;
      if( AbovEn < BeloWnp )
      { 
        // Then the index is closer to end than start.
        refi=MAXIM(n-AbovEn,n/2);
      }else{
        // The it is closer to  start than the end
        refi=MINIM(BeloWnp,n/2);
      }//endif
    }//endif
  }//endif
  else
  {
    // This happens when the index is larger than the padding size
    // this happens in some cases because rnxp_ is larger than nxp_
    // and the x cycle is of length rnxp_ 
    refi=IMISSING;
  }//endif 
  return(refi);
}

int  
FFTGrid::getZSimboxIndex(int k)
{
  int refk;

  if(k < (nz_+nzp_)/2)
    refk=k;
  else
    refk=k-nzp_;

  return refk;
}


float 
FFTGrid::getDistToBoundary(int i, int n, int np )
{
  //  for the series                 i = 0,1,2,3,4,5,6,7 
  //  GetFillNumber(i, 5 , 8)  returns   0,0,0,0,0,p,r,p  p is between 0 and 1, r is larger than 1
  //  GetFillNumber(i, 4 , 8)  returns   0,0,0,0,p,r,r,p  p is between 0 and 1, r's are larger than 1
  //  GetFillNumber(i, 3 , 8)  returns   0,0,0,p,r,r,r,p  p is between 0 and 1, r's are larger than 1

  float dist     =    0.0;
  float taperlength = 0.0;
  int   BeloWnp, AbovEn;

  if (i< np)  
  {
    if (i<n)
      // then it is in the main cube
      dist  =  0.0;
    else
    { 
      taperlength= static_cast<float>((MINIM(n,np-n)/2.1)) ;// taper goes to zero  at taperlength
      BeloWnp  = np-i;   
      AbovEn   = i-(n-1);
      if( AbovEn < BeloWnp )
      { 
        // Then the index is closer to end than start.
        dist = static_cast<float>(AbovEn/taperlength);
      }
      else
      {
        // The it is closer to  start than the end (or identical to)
        dist = static_cast<float>(BeloWnp/taperlength); 
      }//endif
    }//endif
  }//endif
  else
  {
    // This happens when the index is larger than the padding size
    // this happens in some cases because rnxp_ is larger than nxp_
    // and the x cycle is of length rnxp_ 
    dist=RMISSING;
  }//endif 
  return(dist);
}


fftw_complex 
FFTGrid::getNextComplex() 
{
  assert(istransformed_==true);
  assert(counterForGet_ < csize_);
  counterForGet_  +=  1; 
  if(counterForGet_ == csize_)
  {
    counterForGet_=0;
    return(cvalue_[csize_ - 1]);
  }
  else 
    return(cvalue_[counterForGet_ - 1] );
}



float 
FFTGrid::getNextReal()
{
  assert(istransformed_ == false);
  assert(counterForGet_ < rsize_);
  counterForGet_  +=  1;
  float r;
  if(counterForGet_==rsize_)  
  {
    counterForGet_=0; 
    r = static_cast<float>(rvalue_[rsize_-1]);
  } 
  else  
    r = static_cast<float>(rvalue_[counterForGet_-1]);  
  return r;
} 

float        
FFTGrid::getRealValue(int i, int j, int k, bool extSimbox)
{ 
  // when index is in simbox (or the extended simbox if extSimbox is true) it returns the grid value  
  // else it returns RMISSING
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 
  float value;

  // assert(istransformed_==false);

  bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
    ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));


  if( inSimbox && notMissing )
  { // if index in simbox
    int index=i+rnxp_*j+k*rnxp_*nyp_;
    value = static_cast<float>(rvalue_[index]); 
  }
  else
  {
    value = RMISSING;
  }

  return( value );
}
float *
FFTGrid::getRealTrace(int i, int j)
{
  float *value = new float[nz_];
  for(int k=0;k<nz_;k++)
    value[k] = FFTGrid::getRealValue(i,j,k);

  return value;

}

std::vector<float>   
FFTGrid::getRealTrace2(int i, int j)
{
  std::vector<float> value;
  for(int k = 0; k < nz_; k++)
    value.push_back(getRealValue(i,j,k));
  return value;
}

float        
FFTGrid::getRealValueCyclic(int i, int j, int k)
{
  float value;
  if(i<0) 
    i = nxp_+i;
  if(j<0)
    j = nyp_+j;
  if(k<0)
    k = nzp_+k;
  
  if(i<nxp_ && j<nyp_ && k<nzp_)
  {
    int index=i+rnxp_*j+k*rnxp_*nyp_;
    value = static_cast<float>(rvalue_[index]); 
  }
  else
  {
    value = RMISSING;
  }

  return( value );

}

float        
FFTGrid::getRealValueInterpolated(int i, int j, float kindex, bool extSimbox)
{ 
  // when index is in simbox (or the extended simbox if extSimbox is true) it returns the grid value  
  // else it returns RMISSING
  // i index in x direction 
  // j index in y direction 
  // k index in z direction, float, should interpolate

  float value, val1, val2;

  int k1 = int(floor(kindex));
  val1 = getRealValue(i,j,k1,extSimbox);
  if(val1==RMISSING)
    return(RMISSING);
  int k2 = k1+1;
  val2 = getRealValue(i,j,k2,extSimbox);
  if(val2==RMISSING)
    return(val1);
  value = float(1.0-(kindex-k1))*val1+float(kindex-k1)*val2;
  return(value);

}

fftw_complex        
FFTGrid::getComplexValue(int i, int j, int k, bool extSimbox) const
{ 
  // when index is in simbox (or the extended simbox if extSimbox is true) it returns the grid value  
  // else it returns RMISSING
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 
  fftw_complex value;

  assert(istransformed_==true);

  bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
    ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));


  if( inSimbox && notMissing )
  { // if index in simbox
    int index=i + j*cnxp_ + k*cnxp_*nyp_;
    value = fftw_complex (cvalue_[index]); 
  }
  else
  {
    value.re = RMISSING;
    value.im = RMISSING;
  }

  return( value );
}


fftw_complex       
FFTGrid::getFirstComplexValue()
{  
  assert(istransformed_);  
  fftw_complex value;

  setAccessMode(READ);
  counterForGet_=0;

  value = getNextComplex();  

  counterForGet_=0; 
  endAccess();

  return( value );
}

int 
FFTGrid::setNextComplex(fftw_complex value)
{
  assert(istransformed_==true);
  assert(counterForSet_ < csize_);
  counterForSet_  +=  1;
  if(counterForSet_==csize_)
  {
    counterForSet_=0;
    cvalue_[csize_-1]=value;
  }
  else 
    cvalue_[counterForSet_-1]=value;
  return(0);  
}

int   
FFTGrid::setNextReal(float  value)
{    
  assert(istransformed_== false);
  assert(counterForSet_ < rsize_);
  counterForSet_  +=  1; 
  if(counterForSet_==rsize_)
  {
    counterForSet_=0;
    rvalue_[rsize_-1] = static_cast<fftw_real>(value);
  }
  else 
    rvalue_[counterForSet_-1] = static_cast<fftw_real>(value);
  return(0);  
}


int   
FFTGrid::setRealValue(int i, int j ,int k, float  value, bool extSimbox)
{    
  assert(istransformed_== false);

  bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
    ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));

  if( inSimbox && notMissing )
  { // if index in simbox
    int index=i+rnxp_*j+k*rnxp_*nyp_;
    rvalue_[index] = value; 
    return( 0 );
  }
  else
    return(1);
}

int FFTGrid::setRealTrace(int i, int j, float *value)
{
  int notok;
  for(int k=0;k<nz_;k++)
  {
    notok = setRealValue(i,j,k,value[k]);
    if(notok==1)
      return(1);
  }
  return(0);


}
int   
FFTGrid::setComplexValue(int i, int j ,int k, fftw_complex value, bool extSimbox)
{    
  assert(istransformed_== true);

  bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
    ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));

  if( inSimbox && notMissing )
  { // if index in simbox
    int index=i + j*cnxp_ + k*cnxp_*nyp_;
    cvalue_[index] = value; 
    return( 0 );
  }
  else
    return(1);
}

int
FFTGrid::square()
{
  int i;

  if(istransformed_==true)
  {
    for(i = 0;i < csize_; i++)
    {
      if ( cvalue_[i].re == RMISSING || cvalue_[i].im == RMISSING)
      {
        cvalue_[i].re = RMISSING;
        cvalue_[i].im = RMISSING;
      }
      else
      {
        cvalue_[i].re = cvalue_[i].re * cvalue_[i].re + cvalue_[i].im * cvalue_[i].im;
        cvalue_[i].im = 0.0;
      }
    } // i
  }
  else
  {
    for(i = 0;i < rsize_; i++)
    {
      if( rvalue_[i]== RMISSING)
      {
        rvalue_[i] = RMISSING;
      }
      else
      {
        rvalue_[i] = float( rvalue_[i]*rvalue_[i] );
      }

    }// i
  }

  return(0);  
}

int
FFTGrid::expTransf()
{
  assert(istransformed_==false);
  int i;
  for(i = 0;i < rsize_; i++)
  {
    if( rvalue_[i]== RMISSING)
    {
      rvalue_[i] = RMISSING;
    }
    else
    {
      rvalue_[i] = float( exp(rvalue_[i]) );
    }

  }// i
  return(0);  
}

int
FFTGrid::logTransf()
{
  assert(istransformed_==false);
  int i;
  for(i = 0;i < rsize_; i++)
  {
    if( rvalue_[i]== RMISSING ||  rvalue_[i] <= 0.0 )
    {
      rvalue_[i] = 0;
    }
    else
    {
      rvalue_[i] = float( log(rvalue_[i]) );
    }

  }// i
  return(0);  
}

int
FFTGrid::collapseAndAdd(float * grid)
{
  assert(istransformed_==false);
  int   i,j;
  float value;


  for(j = 0; j < nyp_; j++)
    for(i=0;i<nxp_;i++)   
    { 
      value = rvalue_[i + j*rnxp_ ];
      grid[i + j*nxp_] += value ;  
    }
    return(0);  
}

void
FFTGrid::fftInPlace()
{ 
  // uses norm preserving transform for parameter and data
  // in case of correlation and cross correlation it  
  // scale  by 1/N on the inverse such that it maps between 
  // the correlation function and eigen values of the corresponding circular matrix 

  time_t timestart, timeend;
  time(&timestart);

  assert(istransformed_==false);   
  assert(cubetype_!= CTMISSING);

  if( cubetype_!= COVARIANCE )  
    FFTGrid::multiplyByScalar(static_cast<float>(1.0/( sqrt(static_cast<float>(nxp_*nyp_*nzp_)) ) ) );

  int flag;
  rfftwnd_plan plan;  
  flag = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan= rfftw3d_create_plan(nzp_,nyp_,nxp_,FFTW_REAL_TO_COMPLEX,flag);
  rfftwnd_one_real_to_complex(plan,rvalue_,cvalue_);
  fftwnd_destroy_plan(plan);
  istransformed_=true;
  time(&timeend);
  LogKit::LogFormatted(LogKit::DEBUGLOW,"\nFFT of grid type %d finished after %ld seconds \n",cubetype_, timeend-timestart);  
}

void
FFTGrid::invFFTInPlace()
{  
  // uses norm preserving transform for parameter and data
  // in case of correlation and cross correlation it  
  // scale  by 1/N on the inverse such that it maps between 
  // the correlation function and eigen values of the corresponding circular matrix 

  time_t timestart, timeend;
  time(&timestart);

  assert(istransformed_==true);
  assert(cubetype_!= CTMISSING);

  float scale;
  int flag;
  rfftwnd_plan plan; 
  if(cubetype_==COVARIANCE)
    scale=float( 1.0/(nxp_*nyp_*nzp_));
  else
    scale=float( 1.0/sqrt(float(nxp_*nyp_*nzp_)));

  flag = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan= rfftw3d_create_plan(nzp_,nyp_,nxp_,FFTW_COMPLEX_TO_REAL,flag);
  rfftwnd_one_complex_to_real(plan,cvalue_,rvalue_);
  fftwnd_destroy_plan(plan);
  istransformed_=false; 
  FFTGrid::multiplyByScalar( scale);
  time(&timeend);
  LogKit::LogFormatted(LogKit::DEBUGLOW,"\nInverse FFT of grid type %d finished after %ld seconds \n",cubetype_, timeend-timestart);  
}

void
FFTGrid::realAbs()
{
  assert(istransformed_==true);
  int i;
  for(i=0;i<csize_;i++)
  {
    cvalue_[i].re = float (  sqrt( cvalue_[i].re * cvalue_[i].re ) );
    cvalue_[i].im = 0.0; 
  }
}

void
FFTGrid::add(FFTGrid* fftGrid)
{
  assert(nxp_==fftGrid->getNxp());
  if(istransformed_==true)
  {
    int i;
    for(i=0;i<csize_;i++)
    {
      cvalue_[i].re += fftGrid->cvalue_[i].re;
      cvalue_[i].im += fftGrid->cvalue_[i].im; 
    }
  }

  if(istransformed_==false)
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] += fftGrid->rvalue_[i];
    }
  }
}
void
FFTGrid::subtract(FFTGrid* fftGrid)
{
  assert(nxp_==fftGrid->getNxp());
  if(istransformed_==true)
  {
    int i;
    for(i=0;i<csize_;i++)
    {
      cvalue_[i].re -= fftGrid->cvalue_[i].re;
      cvalue_[i].im -= fftGrid->cvalue_[i].im; 
    }
  }

  if(istransformed_==false)
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] -= fftGrid->rvalue_[i];
    }
  }
}
void
FFTGrid::changeSign()
{
  if(istransformed_==true)
  {
    int i;
    for(i=0;i<csize_;i++)
    {
      cvalue_[i].re = -cvalue_[i].re;
      cvalue_[i].im = -cvalue_[i].im; 
    }
  }

  if(istransformed_==false)
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] = -rvalue_[i];
    }
  }
}
void
FFTGrid::multiply(FFTGrid* fftGrid)
{
  assert(nxp_==fftGrid->getNxp());
  if(istransformed_==true)
  {
    int i;
    for(i=0;i<csize_;i++)
    {
      fftw_complex tmp = cvalue_[i];
      cvalue_[i].re = fftGrid->cvalue_[i].re*tmp.re - fftGrid->cvalue_[i].im*tmp.im;
      cvalue_[i].im = fftGrid->cvalue_[i].im*tmp.re + fftGrid->cvalue_[i].re*tmp.im;
    }
  }

  if(istransformed_==false)
  {
    int i;
    for(i=0;i < rsize_;i++)
    {
      rvalue_[i] *= fftGrid->rvalue_[i];
    }
  }
}

void 
FFTGrid::multiplyByScalar(float scalar)
{
  assert(istransformed_==false);
  int i;
  for(i=0;i<rsize_;i++)
  {
    rvalue_[i]*=scalar;
  }
}

void
FFTGrid::computeCircCorrT(Corr* corr,fftw_real* CircCorrT)
{
  int k,n,refk;
  float dummy;
  const float* CorrT;

  CorrT = corr->getPriorCorrT(n,dummy);

  assert(CorrT[0] != 0);

  for(k = 0 ; k < 2*(nzp_/2+1) ; k++ )
  {
    if(k < nzp_)
    {
      if( k < nzp_/2+1)
        refk = k;
      else
        refk = nzp_ - k;
      if(refk < n)
        CircCorrT[k] = CorrT[refk];
      else
        CircCorrT[k] = 0.0;
    }
    else
    {
      CircCorrT[k]=RMISSING;
    }
  }
}


void
FFTGrid::makeCircCorrTPosDef(fftw_real* CircCorrT,int minIntFq)
{
  int k;
  fftw_complex* fftCircCorrT;
  fftCircCorrT = fft1DzInPlace(CircCorrT);
  for(k = 0; k < nzp_/2+1; k++)
  {
    if(k <= minIntFq)
      fftCircCorrT[k].re = 0.0 ;
    else
      fftCircCorrT[k].re = float(sqrt(fftCircCorrT[k].re * fftCircCorrT[k].re + 
      fftCircCorrT[k].im * fftCircCorrT[k].im ));
    fftCircCorrT[k].im = 0.0;
  }

  CircCorrT   = invFFT1DzInPlace(fftCircCorrT);
  //
  // NBNB-PAL: If the number of layers is too small CircCorrT[0] = 0. How 
  //           do we avoid this, or how do we flag the problem?
  //
  float scale;
  if (CircCorrT[0] > 1.e-5f) // NBNB-PAL: Temporary solution for above mentioned problem
    scale = float( 1.0/CircCorrT[0] );
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"\nERROR: The circular temporal correlation (CircCorrT) is undefined. You\n");
    LogKit::LogFormatted(LogKit::LOW,"       probably need to increase the number of layers...\n\nAborting\n");
    exit(1);
  }    
  for(k = 0; k < nzp_; k++)
  {
    CircCorrT[k] *= scale;
  }
}


fftw_complex* 
FFTGrid::fft1DzInPlace(fftw_real*  in)
{  
  // in is over vritten by out 
  // not norm preservingtransform ifft(fft(funk))=N*funk

  int flag;

  rfftwnd_plan plan;  
  fftw_complex* out;
  out = reinterpret_cast<fftw_complex*>(in);

  flag    = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan    = rfftwnd_create_plan(1, &nzp_ ,FFTW_REAL_TO_COMPLEX,flag);
  rfftwnd_one_real_to_complex(plan,in ,out);
  fftwnd_destroy_plan(plan);

  return out;
}

fftw_real*  
FFTGrid::invFFT1DzInPlace(fftw_complex* in)
{
  // in is over vritten by out 
  // not norm preserving transform  ifft(fft(funk))=N*funk

  int flag;
  rfftwnd_plan plan; 
  fftw_real*  out;
  out = reinterpret_cast<fftw_real*>(in); 

  flag = FFTW_ESTIMATE | FFTW_IN_PLACE;
  plan= rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);
  rfftwnd_one_complex_to_real(plan,in,out);
  fftwnd_destroy_plan(plan);
  return out;
}

bool
FFTGrid::consistentSize(int nx,int ny, int nz, int nxp, int nyp, int nzp)
{
  bool consistent = (nx==nx_);
  if( consistent ) consistent = (ny==ny_);
  if( consistent ) consistent = (nz==nz_);
  if( consistent ) consistent = (nxp==nxp_);
  if( consistent ) consistent = (nyp==nyp_);
  if( consistent ) consistent = (nzp==nzp_);
  return consistent;
}


void 
FFTGrid::writeFile(const std::string & fName, 
                   const std::string & subDir, 
                   const Simbox      * simbox, 
                   const std::string   label, 
                   const float         z0, 
                   GridMapping       * depthMap, 
                   GridMapping       * timeMap,
                   const TraceHeaderFormat & thf)
{
  std::string fileName = IO::makeFullFileName(subDir, fName);      
  
  if(formatFlag_ > 0) //Output format specified.
  {
    bool expTrans = (label == "exptrans"); // use sgriLabel to active exponential transformations
    if((domainFlag_ & IO::TIMEDOMAIN) > 0) {
      if(timeMap == NULL) { //No resampling of storm 
        if(expTrans) {
          if((formatFlag_ & IO::STORM) > 0) 
            FFTGrid::writeStormFile(fileName, simbox, true, false);
          if((formatFlag_ & IO::ASCII) > 0)
            FFTGrid::writeStormFile(fileName, simbox, true, true);
        }
        else {
          if((formatFlag_ & IO::STORM) > 0) 
            FFTGrid::writeStormFile(fileName, simbox, false, false);
          if((formatFlag_ & IO::ASCII) > 0)
            FFTGrid::writeStormFile(fileName, simbox, false, true);
        }
      }
      else {
        FFTGrid::writeResampledStormCube(timeMap, fileName, simbox, formatFlag_, expTrans);
      }

      //SEGY, SGRI, and CRAVA are never resampled in time.
      if((formatFlag_ & IO::SEGY) >0)
        FFTGrid::writeSegyFile(fileName, simbox, z0, thf);
      if((formatFlag_ & IO::SGRI) >0)
        FFTGrid::writeSgriFile(fileName, simbox, label);
      if((formatFlag_ & IO::CRAVA) >0)
        FFTGrid::writeCravaFile(fileName, simbox);
    }

    if(depthMap != NULL && (domainFlag_ & IO::DEPTHDOMAIN) > 0) { //Writing in depth. Currently, only stormfiles are written in depth.
      std::string depthName = fileName+"_Depth";
      if(depthMap->getMapping() == NULL) {
        if(depthMap->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::WARNING,
                               "Depth interval lacking when trying to write %s. Write cancelled.\n",depthName.c_str());
          return;
        }
        if(expTrans) {
          if((formatFlag_ & IO::STORM) > 0) 
            FFTGrid::writeStormFile(depthName, depthMap->getSimbox(), true, false);
          if((formatFlag_ & IO::ASCII) > 0)
            FFTGrid::writeStormFile(depthName, depthMap->getSimbox(), true, true);
        }
        else {
          if((formatFlag_ & IO::STORM) >0) 
            FFTGrid::writeStormFile(depthName, depthMap->getSimbox(), false, false);
          if((formatFlag_ & IO::ASCII) >0)
            FFTGrid::writeStormFile(depthName, depthMap->getSimbox(), false, true);
        }

        if((formatFlag_ & IO::SEGY) >0)
          makeDepthCubeForSegy(depthMap->getSimbox(),depthName);
      }
      else
      {
        if(depthMap->getSimbox() == NULL) {
          LogKit::LogFormatted(LogKit::WARNING,
                               "Depth mapping incomplete when trying to write %s. Write cancelled.\n",depthName.c_str());
          return;
        }
        // Writes also segy in depth if required
        FFTGrid::writeResampledStormCube(depthMap, depthName, simbox, formatFlag_, expTrans);
      }
    }
  }
}

void
FFTGrid::writeStormFile(const std::string & fileName, const Simbox * simbox, bool expTrans, bool ascii, bool padding, bool flat)
{
  int nx, ny, nz;
  if(padding == true)
  {
    nx = nxp_;
    ny = nyp_;
    nz = nzp_;
  }
  else
  {
    nx = nx_;
    ny = ny_;
    nz = nz_;
  }

  std::string gfName;
  std::string header = simbox->getStormHeader(cubetype_, nx, ny, nz, flat, ascii);
  int i, j, k;
  float value;

  if(ascii == false) {
    gfName = fileName + IO::SuffixStormBinary();
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM binary file "+gfName+"...");
    std::ofstream binFile;
    NRLib::OpenWrite(binFile, gfName, std::ios::out | std::ios::binary);
    binFile << header;
    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++)
        for(i=0;i<nx;i++)
        {
          if (expTrans)
            value = exp(getRealValue(i,j,k,true));
          else
            value = getRealValue(i,j,k,true);          
          NRLib::WriteBinaryFloat(binFile, value);
        }
    //NRLib::WriteBinaryInt(binFile, 0);
    binFile << "0\n";
    binFile.close();
  }
  else {
    gfName = fileName + IO::SuffixGeneralData();
    std::ofstream file;
    NRLib::OpenWrite(file, gfName);
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM ascii file "+gfName+"...");
    file << header;
    file << std::fixed << std::setw(11) << std::setprecision(6);
    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++) {
        for(i=0;i<nx;i++) {
          value = getRealValue(i,j,k,true);
          if (expTrans)
            file << exp(value) << " ";
          else
            file << value << " ";
        }
        file << "\n";
      }
    file << "0\n";
  }

  LogKit::LogFormatted(LogKit::LOW,"done\n");
}


int
FFTGrid::writeSegyFile(const std::string & fileName, const Simbox * simbox, float z0, const TraceHeaderFormat &thf)
{
  //  long int timestart, timeend;
  //  time(&timestart);

  std::string gfName = fileName + IO::SuffixSegy();
  //SegY * segy = new SegY(gfName, simbox);
  TextualHeader header = TextualHeader::standardHeader();
  float dz = float(floor(simbox->getdz()+0.5));
  if(dz==0)
    dz = 1.0;
  int segynz;
  double zmin,zmax;
  simbox->getMinMaxZ(zmin,zmax);
  segynz = int(ceil((zmax-z0)/dz));
  SegY *segy = new SegY(gfName.c_str(),z0,segynz,dz,header, thf);
  SegyGeometry * geometry = new SegyGeometry(simbox->getx0(), simbox->gety0(), simbox->getdx(), simbox->getdy(), 
                                             simbox->getnx(), simbox->getny(),simbox->getIL0(), simbox->getXL0(), 
                                             simbox->getILStepX(), simbox->getILStepY(), 
                                             simbox->getXLStepX(), simbox->getXLStepY(),  
                                             simbox->getAngle());
  segy->SetGeometry(geometry);
  delete geometry; //Call above takes a copy.
  LogKit::LogFormatted(LogKit::LOW,"\nWriting SEGY file "+gfName+"...");

  int i,j,k;
  double x,y,z;
  std::vector<float> trace(segynz);//Maximum amount of data needed.
  for(j=0;j<simbox->getny();j++)
  {
    for(i=0;i<simbox->getnx();i++)
    {
      simbox->getCoord(i, j, 0, x, y, z);
      z = simbox->getTop(x,y);

      if(z == RMISSING || z == WELLMISSING)
      {
        //printf("Missing trace.\n");
        for(k=0;k<nz_;k++)
          trace[k] = 0;
      }
      else
      {
        double gdz       = simbox->getdz()*simbox->getRelThick(i,j);
        int    firstData = static_cast<int>(floor(0.5+(z-z0)/dz));
        int    endData   = static_cast<int>(floor(0.5+((z-z0)+nz_*gdz)/dz));

        if(endData > segynz)
        {
          printf("Internal warning: SEGY-grid too small (%d, %d needed). Truncating data.\n", nz_, endData); 
          endData = segynz;
        }
        for(k=0;k<firstData;k++)
          trace[k] = 0;
        for(;k<endData;k++)
          trace[k] = this->getRegularZInterpolatedRealValue(i,j,z0,dz,k,z,gdz);
        for(;k<segynz;k++)
          trace[k] = 0;

        float xx = float(x);
        float yy = float(y);
        segy->StoreTrace(xx, yy, trace, NULL);
      }
    }
  }

  segy->WriteAllTracesToFile();

  delete segy; //Closes file.
  // delete [] value;
  LogKit::LogFormatted(LogKit::LOW,"done\n");
  //  time(&timeend);
  //printf("\n Write SEGY was performed in %ld seconds.\n",timeend-timestart);
  return(0);
}


void
FFTGrid::writeResampledStormCube(GridMapping       * gridmapping, 
                                 const std::string & fileName, 
                                 const Simbox      * simbox,
                                 const int           format,
                                 bool                expTrans)
{
  // simbox is related to the cube we resample from. gridmapping contains simbox for the cube we resample to.
 
  float time, kindex;
  StormContGrid *mapping = gridmapping->getMapping();
  StormContGrid *outgrid = new StormContGrid(*mapping);
 
  double x,y;
  int nz = mapping->GetNK();
  for(int i=0;i<nx_;i++)
  {
    for(int j=0;j<ny_;j++)
    {
      simbox->getXYCoord(i,j,x,y);
      for(int k=0;k<nz;k++)
      {
        time = (*mapping)(i,j,k);
        kindex = float((time - static_cast<float>(simbox->getTop(x,y)))/simbox->getdz());
        float value = getRealValueInterpolated(i,j,kindex);
        if (expTrans)
          (*outgrid)(i,j,k) = exp(value);
        else
          (*outgrid)(i,j,k) = value;
      }
    }
  }

  std::string gfName;
  std::string header;
  if ((format & IO::ASCII) > 0) // ASCII
  {
    gfName = fileName + IO::SuffixGeneralData();
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER,nx_,ny_,nz, 0, 1);
    outgrid->SetFormat(StormContGrid::STORM_ASCII);
    outgrid->WriteToFile(gfName,header);
  }

  if ((format & IO::STORM) > 0)
  {
    gfName =  fileName + IO::SuffixStormBinary();
    header = gridmapping->getSimbox()->getStormHeader(FFTGrid::PARAMETER,nx_,ny_,nz, 0, 0);
    outgrid->SetFormat(StormContGrid::STORM_BINARY);
    outgrid->WriteToFile(gfName,header);
  }
  if((formatFlag_ & IO::SEGY) > 0)
  {
    gfName =  fileName + IO::SuffixSegy();
    writeSegyFromStorm(outgrid,gfName);
  }
  delete outgrid;
}

int
FFTGrid::writeSgriFile(const std::string & fileName, const Simbox *simbox, const std::string label)
{
  double vertScale = 0.001;
  double horScale  = 0.001;
  std::string fName = fileName + IO::SuffixSgriHeader();

  std::ofstream headerFile;
  NRLib::OpenWrite(headerFile, fName);

  LogKit::LogFormatted(LogKit::LOW,"\nWriting SGRI header file "+fName+"...");

  headerFile << "NORSAR General Grid Format v1.0\n";
  headerFile << "3\n";
  headerFile << "X (km)\n";
  headerFile << "Y (km)\n";
  headerFile << "Z (s)\n";
  headerFile << "FFT-grid\n";
  headerFile << "1\n";
  headerFile << label << std::endl;
  headerFile << "1 1 1\n";

  double zMin, zMax;
  simbox->getMinMaxZ(zMin, zMax);
/*  float dz = static_cast<float> (floor(simbox->getdz()+0.5)); //To have the same sampling as in SegY
  if (dz == 0.0)
    dz = 1.0; */
  float dz = static_cast<float> (simbox->getdz());
  int nz = static_cast<int> (ceil((zMax - zMin)/dz));
  int ny = simbox->getny();
  int nx = simbox->getnx();
  headerFile << nx << " " << ny << " " << nz << std::endl;
  headerFile << std::setprecision(10);
  headerFile << simbox->getdx()*horScale << " " << simbox->getdy()*horScale << " " << dz*vertScale << std::endl;
  double x0 = simbox->getx0() + 0.5 * simbox->getdx();
  double y0 = simbox->gety0() + 0.5 * simbox->getdy();
  double z0 = zMin + 0.5 * dz;
  headerFile << x0*horScale << " " << y0*horScale << " " << z0*vertScale << std::endl;
  headerFile << simbox->getAngle() << " 0\n";
  headerFile << RMISSING << std::endl;
  
  fName = fileName + IO::SuffixSgri();
  headerFile << fName << std::endl;
  headerFile << "0\n";

  std::ofstream binFile;
  NRLib::OpenWrite(binFile, fName, std::ios::out | std::ios::binary);

  int i,j,k;
  double x, y, z, zTop, zBot;
  float value;
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      simbox->getXYCoord(i, j, x, y);
      zTop = simbox->getTop(x, y);
      zBot = simbox->getBot(x, y);
      double scale = 1.0/(simbox->getdz()*simbox->getRelThick(i,j));
      for (k=0; k<nz; k++) {
        z = zMin + k*dz;
        if (z < zTop || z > zBot)
          value = RMISSING;
        else {
          int simboxK = static_cast<int> ((z - zTop)*scale + 0.5);
          value = getRealValue(i,j,simboxK);
        }
#ifndef BIGENDIAN
        NRLib::WriteBinaryFloat(binFile, value);
#else
        NRLib::WriteBinaryFloat(binFile, value, END_LITTLE_ENDIAN);
#endif
      }
    }
  }
  return(0);
}


void
FFTGrid::writeCravaFile(const std::string & fileName, const Simbox * simbox)
{
  try {
    std::ofstream binFile;
    std::string fName = fileName + IO::SuffixCrava();
    NRLib::OpenWrite(binFile, fName, std::ios::out | std::ios::binary);

    std::string fileType = "crava_fftgrid_binary";
    binFile << fileType << "\n";

    NRLib::WriteBinaryDouble(binFile, simbox->getx0());
    NRLib::WriteBinaryDouble(binFile, simbox->gety0());
    NRLib::WriteBinaryDouble(binFile, simbox->getdx());
    NRLib::WriteBinaryDouble(binFile, simbox->getdy());
    NRLib::WriteBinaryInt(binFile, simbox->getnx());
    NRLib::WriteBinaryInt(binFile, simbox->getny());
    NRLib::WriteBinaryDouble(binFile, simbox->getIL0());
    NRLib::WriteBinaryDouble(binFile, simbox->getXL0());
    NRLib::WriteBinaryDouble(binFile, simbox->getILStepX());
    NRLib::WriteBinaryDouble(binFile, simbox->getILStepY());
    NRLib::WriteBinaryDouble(binFile, simbox->getXLStepX());
    NRLib::WriteBinaryDouble(binFile, simbox->getXLStepY());
    NRLib::WriteBinaryDouble(binFile, simbox->getAngle());
    NRLib::WriteBinaryInt(binFile, rnxp_);
    NRLib::WriteBinaryInt(binFile, nyp_);
    NRLib::WriteBinaryInt(binFile, nzp_);
    int i;
    for(i=0;i<rsize_;i++)
      NRLib::WriteBinaryFloat(binFile, rvalue_[i]);

    binFile.close();
  }
  catch (NRLib::Exception & e) {
    std::string message = "Error: "+std::string(e.what())+"\n";
    LogKit::LogMessage(LogKit::ERROR, message);
  }
}


void
FFTGrid::readCravaFile(const std::string & fileName, std::string & errText) 
{
  std::string error;
  try {
    std::ifstream binFile;
    NRLib::OpenRead(binFile, fileName, std::ios::in | std::ios::binary);

    std::string fileType;
    getline(binFile,fileType);

    double dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryInt(binFile);
    dummy = NRLib::ReadBinaryInt(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    dummy = NRLib::ReadBinaryDouble(binFile);
    int nx = NRLib::ReadBinaryInt(binFile);
    int ny = NRLib::ReadBinaryInt(binFile);
    int nz = NRLib::ReadBinaryInt(binFile);

    if(nx != rnxp_ || ny != nyp_ || nz != nzp_) {
      std::string message = "Grid dimension is wrong for file read from the direct crava format '"+
        fileName+"'.";
      binFile.close();
      throw(NRLib::Exception(message));
    }
    createRealGrid();
    int i;
    for(i=0;i<rsize_;i++)
      rvalue_[i] = NRLib::ReadBinaryFloat(binFile);
    
    binFile.close();
  }
  catch (NRLib::Exception & e) {
    error = std::string("Error: ") + e.what() + "\n";
  }
  errText += error;
}


float
FFTGrid::getRegularZInterpolatedRealValue(int i, int j, double z0Reg, 
                                          double dzReg, int kReg, 
                                          double z0Grid, double dzGrid)
{
  float z     = static_cast<float> (z0Reg+dzReg*kReg);
  float t     = static_cast<float> ((z-z0Grid)/dzGrid);
  int   kGrid = int(t);

  t -= kGrid;

  if(kGrid<0) {
    if(kGrid<-1)
      return(RMISSING);
    else {
      kGrid = 0;
      t = 0;
    }
  }
  else if(kGrid>nz_-2) {
    if(kGrid > nz_-1)
      return(RMISSING);
    else {
      kGrid = nz_-2;
      t = 1;
    }
  }

  float v1 = getRealValue(i,j,kGrid);
  float v2 = getRealValue(i,j,kGrid+1);
  float v  = (1-t)*v1+t*v2;
  return(v);
}

void
FFTGrid::writeAsciiFile(const std::string & fileName)
{
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  LogKit::LogFormatted(LogKit::LOW,"\nWriting ASCII file "+fileName+"...");
  int i, j, k;
  float value;
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)
      {
        value = getNextReal(); 
        file << value << "\n";
      }
  file.close();
  LogKit::LogFormatted(LogKit::LOW,"done.\n");
}

void
FFTGrid::writeAsciiRaw(const std::string & fileName)
{
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  LogKit::LogFormatted(LogKit::LOW,"\nWriting ASCII file "+fileName+"...");
  int i, j, k;
  float value;
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)
      {
        value = getNextReal(); 
        file << i << " " << j << " " << k << " " << value;
      }
  file.close();
  LogKit::LogFormatted(LogKit::LOW,"done.\n");
}


void
FFTGrid::interpolateSeismic(float energyTreshold)
{
  assert(cubetype_ == DATA);
  int i, j, k, index = 0;
  short int * flags = new short int[nx_*ny_];
  int curFlag, flag = 0; //Flag rules: bit 0 = this trace bad, bit 1 = any prev. bad
  int imin = nx_;
  int imax = 0;
  int jmin = ny_;
  int jmax = 0;
  float value, energy, totalEnergy = 0;;
  float * energyMap = new float[nx_*ny_];
  for(j=0;j<ny_;j++)
    for(i=0;i<nx_;i++)
    {
      energy = 0;
      for(k=0;k<nz_;k++) {
        value = getRealValue(i, j, k);
        energy += value*value;
      }
      energyMap[index] = energy;
      totalEnergy += energy;
      index++;
    }

    index = 0;
    float energyLimit = energyTreshold*totalEnergy/float(nx_*ny_);
    int nInter = 0;    //#traces interpolated.
    int nInter0 = 0;   //#traces interpolated where there was no response at all.
    for(j=0;j<ny_;j++)
      for(i=0;i<nx_;i++)
      {
        curFlag = 0;
        if(energyMap[index] <= energyLimit) {//Values in this trace are bogus, interpolate.
          curFlag = 1;
          nInter++;
          if(energyMap[index] == 0.0f)
            nInter0++;
        }
        flags[index] = short(flag+curFlag);
        if(curFlag == 1)
          flag = 2;
        else
        {
          if(i < imin)
            imin = i;
          if(i > imax)
            imax = i;
          if(j < jmin)
            jmin = j;
          if(j > jmax)
            jmax = j;
        }
        index++;
      }

      LogKit::LogFormatted(LogKit::LOW,"\nInterpolated %d of %d traces (%d with zero response).",
        nInter, nx_*ny_, nInter0);
      int curIndex = 0;
      for(j=0;j<ny_;j++)
        for(i=0;i<nx_;i++)
        {
          if((flags[curIndex] % 2) == 1)
          {
            if((interpolateTrace(curIndex, flags, i, j) % 2) == 0)
            {
              index = curIndex-1;
              while(index >= 0 && flags[index] > 1)
              {
                if((flags[index] % 2) == 1)
                  interpolateTrace(index, flags, (index % nx_), index/nx_);
                index--;
              }
              if(index < 0)
                index = 0;
              assert(flags[index] < 2);
              index++;
              while(flags[index-1] == 0 && index <= curIndex)
              {
                if(flags[index] > 1)
                  flags[index] -= 2;
                index++;
              }
            }
          }
          curIndex++;
        }
        extrapolateSeismic(imin, imax, jmin, jmax);
        delete [] energyMap;
        delete [] flags;
}


int
FFTGrid::interpolateTrace(int index, short int * flags, int i, int j)
{
  int  k, nt = 0;
  bool left = (i > 0 && (flags[index-1] % 2) == 0);
  bool right = (i < nx_-1 && (flags[index+1] % 2) == 0);
  bool up = (j > 0 && (flags[index-nx_] % 2) == 0);
  bool down = (j < ny_-1 && (flags[index+nx_] % 2) == 0);
  if(left == true) nt++;
  if(right == true) nt++;
  if(up == true) nt++;
  if(down == true) nt++;
  if(nt > 1)
  {
    float * mean = new float[nzp_];
    for(k=0;k<nzp_;k++)
    {
      if(left == true)
        mean[k] = getRealValue(i-1, j, k,true);
      else
        mean[k] = 0;
    }
    if(right == true)
      for(k=0;k<nzp_;k++)
        mean[k] += getRealValue(i+1,j,k,true);
    if(up == true)
      for(k=0;k<nzp_;k++)
        mean[k] += getRealValue(i,j-1,k,true);
    if(down == true)
      for(k=0;k<nzp_;k++)
        mean[k] += getRealValue(i,j+1,k,true);
    for(k=0;k<nzp_;k++)
      setRealValue(i, j, k, mean[k]/float(nt),true);
    delete [] mean;
    assert((flags[index] % 2) == 1);
    flags[index]--;
  }
  return(flags[index]);
}


void
FFTGrid::extrapolateSeismic(int imin, int imax, int jmin, int jmax)
{
  int i, j, k, refi, refj;
  float value, distx, disty, distz, mult;
  for(j=0;j<nyp_;j++)
  {
    refj = getYSimboxIndex(j);
    if(refj < jmin)
      refj = jmin;
    else if(refj > jmax)
      refj = jmax;
    for(i=0;i<nxp_;i++)
    {
      if(i < imin || i > imax || j < jmin || j > jmax)
      {
        refi = getXSimboxIndex(i);
        if(refi < imin)
          refi = imin;
        else if(refi > imax)
          refi = imax;
        for(k=0;k<nzp_;k++)
        {
          value = getRealValue(refi, refj, k, true);
          distx  = getDistToBoundary(i,nx_,nxp_);
          disty  = getDistToBoundary(j,ny_,nyp_);
          distz  = getDistToBoundary(k,nz_,nzp_);         
          mult   = float(pow(MAXIM(1.0-distx*distx-disty*disty-distz*distz,0.0),3));
          setRealValue(i, j, k, mult*value, true);
        }
      }
    }
  }
}




void
FFTGrid::checkNaN()
{
  /*
  #ifndef UNIX_ELLER_NOE_ANNET
  int i;
  if(istransformed_)
  {
  for(i=0;i<csize_;i++)
  if(_isnan(cvalue_[i].re) ||
  _isnan(cvalue_[i].im))
  break;
  if(i == csize_)
  i = -1;
  }
  else
  {
  for(i=0;i<rsize_;i++)
  if(_isnan(rvalue_[i]))
  break;
  if(i == rsize_)
  i = -1;
  }
  assert(i == -1);
  #endif
  */
}

void FFTGrid::writeSegyFromStorm(StormContGrid *data, std::string fileName)
{

  int i,k,j;
  TextualHeader header = TextualHeader::standardHeader();
  int nx = data->GetNI();
  int ny = data->GetNJ();
  const SegyGeometry *geometry = new SegyGeometry(data->GetXMin(),data->GetYMin(),data->GetDX(),data->GetDY(),
               nx,ny,data->GetAngle());
  float dz = float(floor((data->GetLZ()/data->GetNK())));
  //int nz = int(data->GetZMax()/dz);
  //float z0 = float(data->GetZMin());
  float z0 = 0.0;
  int nz = int(ceil((data->GetZMax())/dz));
  SegY segyout(fileName,0,nz,dz,header);
  segyout.SetGeometry(geometry);

  std::vector<float> datavec;
  datavec.resize(nz);
  float x, y, xt, yt, z;
  for(j=0;j<ny;j++)
    for(i=0;i<nx;i++)
    {
      xt = float((i+0.5)*geometry->GetDx());
      yt = float((j+0.5)*geometry->GetDy());
      x = float(geometry->GetX0()+xt*geometry->GetCosRot()-yt*geometry->GetSinRot());
      y = float(geometry->GetY0()+yt*geometry->GetCosRot()+xt*geometry->GetSinRot());
  
      double zbot= data->GetBotSurface().GetZ(x,y);
      double ztop = data->GetTopSurface().GetZ(x,y);
      int    firstData = static_cast<int>(floor((ztop)/dz));
      int    endData   = static_cast<int>(floor((zbot)/dz));

      if(endData > nz)
      {
        printf("Internal warning: SEGY-grid too small (%d, %d needed). Truncating data.\n", nz, endData); 
        endData = nz;
      }
      for(k=0;k<firstData;k++)
      {
        datavec[k] = 0.0;
      }

      for(k=firstData;k<endData;k++)
      {
         z = z0+k*dz;
        datavec[k] = float(data->GetValueZInterpolated(x,y,z));
      }
      for(k=endData;k<nz;k++)
      {
        datavec[k] = 0.0;
      }
      segyout.StoreTrace(x,y,datavec,NULL);
    }
  
  segyout.WriteAllTracesToFile();

}

void FFTGrid::makeDepthCubeForSegy(Simbox *simbox,const std::string & fileName)
{
  StormContGrid * stormcube = new StormContGrid((*simbox),nx_,ny_,nz_);
  
  int i,j,k;
  for(k=0;k<nz_;k++)
    for(j=0;j<ny_;j++)
      for(i=0;i<nx_;i++)
        (*stormcube)(i,j,k) = getRealValue(i,j,k,true);
  
  std::string fullFileName = fileName + IO::SuffixSegy();
  FFTGrid::writeSegyFromStorm(stormcube, fullFileName);
}

int FFTGrid::formatFlag_        = 0;
int FFTGrid::domainFlag_        = IO::TIMEDOMAIN;
int FFTGrid::maxAllowedGrids_   = 1;   // One grid is allocated and deallocated before memory check.
int FFTGrid::maxAllocatedGrids_ = 0; 
int FFTGrid::nGrids_            = 0;
