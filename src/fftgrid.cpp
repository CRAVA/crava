#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "lib/global_def.h"
#include "lib/segy.h"
#include "lib/lib_misc.h"
#include "lib/random.h"

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "nrlib/iotools/logkit.hpp"

#include "fftgrid.h"
#include "corr.h"
#include "simbox.h"
#include "model.h"

using namespace NRLib2;

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
  rnxp_	          = 2*(cnxp_);	   

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = 0; 
  counterForSet_  = 0;
  istransformed_  = false;	
  rvalue_ = NULL;
  // index= i+rnxp_*j+k*rnxp_*nyp_;
  // i index in x direction 
  // j index in y direction 
  // k index in z direction 
}

FFTGrid::FFTGrid(FFTGrid  * fftGrid)
{
  float value;
  int   i,j,k;


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
  rnxp_	          = 2*(cnxp_);	   

  csize_          = cnxp_*nyp_*nzp_;
  rsize_          = rnxp_*nyp_*nzp_;
  counterForGet_  = fftGrid->getCounterForGet(); 
  counterForSet_  = fftGrid->getCounterForSet();
  istransformed_  = false;		
  createRealGrid();
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)   
      {
        value=fftGrid->getNextReal();
        setNextReal(value);
      }// k,j,i
}



FFTGrid::~FFTGrid()
{
  if (rvalue_!=NULL)
    fftw_free(rvalue_);
}

void
FFTGrid::fillInFromSegY(SegY* segy)
{
  assert(cubetype_  !=  CTMISSING);

  createRealGrid();

  int i,j,k,refi,refj,refk;
  float distx,disty,distz,mult;
  float* meanvalue = NULL;
  bool  isParameter=false;

  if(cubetype_==PARAMETER)  isParameter=true;

  int outMode = SegY::MISSING;
  if(cubetype_ == DATA)
    outMode = SegY::ZERO;
  else if(cubetype_ == PARAMETER)
    outMode = SegY::CLOSEST;

  fftw_real value  = 0.0;

  if(isParameter)
  { 
    meanvalue= (float*)  fftw_malloc(sizeof(float)*nyp_*nxp_);

    for(j=0;j<nyp_;j++)
      for(i=0;i<nxp_;i++)   
      {
        refi = getXSimboxIndex(i);
        refj = getYSimboxIndex(j);
        meanvalue[i+j*nxp_] = float ((segy->getValue(refi,refj,0, outMode) + segy->getValue(refi,refj,nz_-1, outMode) )/2.0);
      }
  } // endif

  LogKit::LogFormatted(LogKit::LOW,"\nResampling seismic data into %dx%dx%d grid:\n",nxp_,nyp_,nzp_);
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
          value=segy->getValue(refi,refj,refk, outMode );
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
}

void
FFTGrid::fillInFromStorm(Simbox * tmpSimBox, Simbox * actSimBox,
                         float * grid, const char * parName)
{
  assert(cubetype_ != CTMISSING);
  createRealGrid();
  int i,j,k,refi,refj,refk;
  float distx,disty,distz,mult;
  float* meanvalue;

  assert(cubetype_==PARAMETER);

  fftw_real value  = 0.0;

  int index, index1, index2;
  double x,y,z, t;
  float val;

  meanvalue= (float*)  fftw_malloc(sizeof(float)*nyp_*nxp_);

  for(j=0;j<nyp_;j++)
    for(i=0;i<nxp_;i++)   
    {
      refi   = getXSimboxIndex(i);
      refj   = getYSimboxIndex(j);
      actSimBox->getCoord(refi, refj, 0, x, y, z);
      index = tmpSimBox->getClosestZIndex(x,y,z);
      if(index != IMISSING)
        val = grid[index];
      else
        val = 0;
      actSimBox->getCoord(refi, refj, nz_-1, x, y, z);
      index = tmpSimBox->getClosestZIndex(x,y,z);
      if(index != IMISSING)
        if(val != 0)
          val = float((grid[index]+val)/2.0);
        else
          val = grid[index];
      meanvalue[i+j*nxp_] = val;
    }

  LogKit::LogFormatted(LogKit::LOW,"Resampling %s into %dx%dx%d grid:\n",parName,nxp_,nyp_,nzp_);
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
          
          /*          index = tmpSimBox->getClosestZIndex(x,y,z);
                      if(index != RMISSING)
                      value=grid[index];
          */
          tmpSimBox->getZInterpolation(x, y, z, index1, index2, t);
          if(index1 != IMISSING)
            value = static_cast<float>((1-t)*grid[index1]+t*grid[index2]);
          else
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

  meanvalue= (float*)  fftw_malloc(sizeof(float)*nyp_*nxp_);

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
FFTGrid::fillInParamCorr(Corr* corr,int minIntFq)
{
  assert(corr->getnx() == nxp_);
  assert(corr->getny() == nyp_);
  assert(istransformed_== false);
  cubetype_ = COVARIANCE;

  createRealGrid();

  int i,j,k;
  float value;
  fftw_real* circCorrT;
  //  float * corrXY;
  circCorrT = (fftw_real*) fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real));
  computeCircCorrT(corr,circCorrT);
  makeCircCorrTPosDef(circCorrT, minIntFq);
  //makeCorrXYPosDef();
  // setAccessMode(RANDOMACCESS);
  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++)   
      {    
        if(i < nxp_)  // computes the index reference from the cube puts it in value
        {
          value = float( (*(corr->getCorrXY()))(i+nxp_*j)) * circCorrT[k];
        }
        else
          value = RMISSING;        

        setNextReal(value);
      } //for k,j,i

      endAccess();
      return circCorrT;//fftw_free(circCorrT);
}

void
FFTGrid::fillInErrCorr(Corr* parCorr)
{
  // Note:  this contain the latteral correlation and the multiplyers for the 
  // time correlation. The time correlation is further adjusted by the wavelet
  // and the derivative of the wavelet.
  // the angular correlation is given by a functional Expression elsewhere

  assert(parCorr->getnx() == nxp_);
  assert(parCorr->getny() == nyp_);
  assert(istransformed_== false);
  cubetype_ = COVARIANCE;

  createRealGrid();

  int i,j,k;
  float value,mult;
  setAccessMode(WRITE);
  for( k = 0; k < nzp_; k++)
    for( j = 0; j < nyp_; j++)
      for( i = 0; i < rnxp_; i++)   
      {    
        //if( k == 0) mult = 1.0; else  if( k == 1 || k == nzp_ - 1 ) mult = float( -0.33333 ); else mult = 0.0;

        if( k == 0) 
          mult = 1.0; 
        else
          mult = 0.0;

        if(i < nxp_)  // computes the index reference from the cube puts it in value
        {
          value = float( (*(parCorr->getCorrXY()))(i+nxp_*j) )*mult ;
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
FFTGrid::createRealGrid()
{
  //  long int timestart, timeend;
  //  time(&timestart);
  istransformed_=false;
  rvalue_    = (fftw_real*) fftw_malloc(rsize_ * sizeof(fftw_real));
  cvalue_    = (fftw_complex*) rvalue_; //
  counterForGet_  = 0; 
  counterForSet_  = 0;
  //  time(&timeend);
  //  LogKit::LogFormatted(LogKit::LOW,"\nReal grid created in %ld seconds.\n",timeend-timestart);
}

void
FFTGrid::createComplexGrid()
{
  //  long int timestart, timeend;
  //  time(&timestart);
  istransformed_  = true;
  rvalue_         = (fftw_real*) fftw_malloc(rsize_ * sizeof(fftw_real));
  cvalue_         = (fftw_complex*) rvalue_; //
  counterForGet_  = 0; 
  counterForSet_  = 0;
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
      taperlength= (float) (MINIM(n,np-n)/2.1) ;// taper goes to zero  at taperlength
      BeloWnp  = np-i;   
      AbovEn   = i-(n-1);
      if( AbovEn < BeloWnp )
      { 
        // Then the index is closer to end than start.
        dist= (float) (AbovEn/taperlength);
      }
      else
      {
        // The it is closer to  start than the end (or identical to)
        dist=(float)  (BeloWnp/taperlength); 
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
    r = (float)   rvalue_[rsize_-1];
  } 
  else  
    r = (float) rvalue_[counterForGet_-1] ;  
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

  assert(istransformed_==false);

  bool  inSimbox   = (extSimbox ? ( (i < nxp_) && (j < nyp_) && (k < nzp_)):
  ((i < nx_) && (j < ny_) && (k < nz_)));
  bool  notMissing = ( (i > -1) && (j > -1) && (k > -1));


  if( inSimbox && notMissing )
  { // if index in simbox
    int index=i+rnxp_*j+k*rnxp_*nyp_;
    value = float (rvalue_[index]); 
  }
  else
  {
    value = RMISSING;
  }

  return( value );
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
    rvalue_[rsize_-1] = (fftw_real) value;
  }
  else 
    rvalue_[counterForSet_-1] = (fftw_real) value;
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
    FFTGrid::multiplyByScalar(float (1.0/( sqrt((float)(nxp_*nyp_*nzp_)) ) ) );

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

  CorrT =  corr->getCorrT(n,dummy);

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
  out = (fftw_complex*) in;

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
  out = (fftw_real*) in; 

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
FFTGrid::writeFile(const char * fileName, const Simbox * simbox, bool writeSegy)
{
  if(formatFlag_ != NONE)
  {
    if((formatFlag_ & STORMFORMAT) == STORMFORMAT)
      writeStormFile(fileName, simbox);
    if((formatFlag_ & SEGYFORMAT) == SEGYFORMAT && writeSegy==1)
      writeSegyFile(fileName, simbox);
    if((formatFlag_ & STORMASCIIFORMAT) == STORMASCIIFORMAT)
      writeStormFile(fileName, simbox, true);
  }
}

void
FFTGrid::writeStormFile(const char * fileName, const Simbox * simbox, bool ascii, bool padding, bool flat)
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

  char * gfName;
  char * header = simbox->getStormHeader(cubetype_, nx, ny, nz, flat, ascii);
  FILE * file;
  int i, j, k;
  float value;

  if(ascii == false) {
    gfName = ModelSettings::makeFullFileName(fileName, ".storm");
    file = fopen(gfName,"wb");
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM binary file %s...",gfName);
    fwrite(header, 1, strlen(header), file);
    char * output = (char * ) &value;
    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++)
        for(i=0;i<nx;i++)
        {
          value = getRealValue(i,j,k,true);
  #ifndef UNIX
          fwrite(&(output[3]),1,1,file);
          fwrite(&(output[2]),1,1,file);
          fwrite(&(output[1]),1,1,file);
          fwrite(&(output[0]),1,1,file);
  #else
          fwrite(output, 1, 4, file);
  #endif
        }
    output[0] = '0';
    output[1] = '\n';
    fwrite(output, 2, 1, file);
  }
  else {
    gfName = ModelSettings::makeFullFileName(fileName, ".txt");
    file = fopen(gfName,"w");
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM ascii file %s...",gfName);
    fprintf(file,"%s",header);
    for(k=0;k<nz;k++)
      for(j=0;j<ny;j++) {
        for(i=0;i<nx;i++) {
          value = getRealValue(i,j,k,true);
          fprintf(file,"%f ", value);
        }
        fprintf(file,"\n");
      }
    fprintf(file,"0\n");
  }
  delete [] gfName;

  fclose(file);
  LogKit::LogFormatted(LogKit::LOW,"done\n");
  //  time(&timeend);
  //printf("\n Write Storm was performed in %ld seconds.\n",timeend-timestart);
  delete [] header;
}


int
FFTGrid::writeSegyFile(const char * fileName, const Simbox * simbox)
{
  //  long int timestart, timeend;
  //  time(&timestart);

  char * gfName = ModelSettings::makeFullFileName(fileName, ".segy");
  SegY * segy = new SegY(gfName, simbox);
  LogKit::LogFormatted(LogKit::LOW,"\nWriting SEGY file %s...",gfName);
  //  LogKit::LogFormatted(LogKit::LOW,"%d x %d traces.\n",nx_, ny_);
  delete [] gfName;
  char errMsg[MAX_STRING];
  if(segy->checkError(errMsg) > 0)
  {
    LogKit::LogFormatted(LogKit::LOW,"failed\n");
    return(1);
  }
  int i, j, jStart, jEnd, dj;
  if(simbox->getILStep() > 0)
  {
    jStart = 0;
    jEnd = ny_;
    dj = 1;
  }
  else
  {
    jStart = ny_-1;
    jEnd = -1;
    dj = -1;
  }

  float * value = new float[nz_];
  for(j=jStart;j != jEnd;j+= dj)
    for(i=0;i<nx_;i++)
    {
      segy->writeTrace(i, j, this);


      /* Old code, assumes dz equal everywhere
      //      LogKit::LogFormatted(LogKit::LOW,"*** Trace %d %d of %d %d\n",i, j, nx_, ny_); 
      for(k=0;k<nz_;k++)
      value[k] = getRealValue(i,j,k);
      segy->writeTrace(i, j, value);
      */
    }
    delete segy; //Closes file.
    delete [] value;
    LogKit::LogFormatted(LogKit::LOW,"done\n");
    //  time(&timeend);
    //printf("\n Write SEGY was performed in %ld seconds.\n",timeend-timestart);
    return(0);
}

void
FFTGrid::writeAsciiFile(char * fileName)
{
  char * gfName = ModelSettings::makeFullFileName(fileName, ".ascii");
  FILE *file = fopen(gfName,"wb");
  LogKit::LogFormatted(LogKit::LOW,"\nWriting ASCII file %s...",gfName);
  delete [] gfName;
  int i, j, k;
  float value;
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)
      {
        value = getNextReal(); 
        fprintf(file,"%f\n",value);
      }
      fclose(file);
      LogKit::LogFormatted(LogKit::LOW,"done.\n");
}

void
FFTGrid::writeAsciiRaw(char * fileName)
{
  char * gfName = ModelSettings::makeFullFileName(fileName, ".ascii");
  FILE *file = fopen(gfName,"wb");
  LogKit::LogFormatted(LogKit::LOW,"\nWriting ASCII file %s...",fileName);
  int i, j, k;
  float value;
  for(k=0;k<nzp_;k++)
    for(j=0;j<nyp_;j++)
      for(i=0;i<rnxp_;i++)
      {
        value = getNextReal(); 
        fprintf(file,"%d %d %d %f\n",i, j, k, value);
      }
      fclose(file);
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

      LogKit::LogFormatted(LogKit::LOW,"Interpolated %d of %d traces (%d with zero response).\n\n",
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

int FFTGrid::formatFlag_ = FFTGrid::NOOUTPUT;
