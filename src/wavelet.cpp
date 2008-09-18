#include <string.h>
#include <assert.h>
#include <math.h>

#include "wavelet.h"
#include "wavelet1D.h"
#include "wavelet3D.h"

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/sgri.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "src/modelsettings.h"
#include "src/model.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/definitions.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet::Wavelet(int dim)
  : dim_(dim)
{
}

Wavelet::Wavelet(ModelSettings * modelSettings, int dim)
  : dim_(dim)
{
	maxShift_       = modelSettings->getMaxWaveletShift();
	minRelativeAmp_ = modelSettings->getMinRelWaveletAmp();

	isReal_     = true;
	inFFTorder_ = false;
	scale_=1; 
	gridNI_=0;   
	gridNJ_=0;
	shiftGrid_=NULL;  
	gainGrid_=NULL; 
}

Wavelet::Wavelet(Wavelet * wavelet, int dim)
  : dim_(dim)
{
  gridNI_    = 0;   
  gridNJ_    = 0;
  shiftGrid_ = NULL;  
  gainGrid_  = NULL; 
  if(! wavelet->getIsReal() ) wavelet->invFFT1DInPlace();
  theta_     = wavelet->theta_;
  dz_        = wavelet->dz_;
  nz_        = wavelet->nz_;
  nzp_       = wavelet->nzp_;
  cz_        = wavelet->cz_;
  inFFTorder_= wavelet->inFFTorder_;
  isReal_    = wavelet->isReal_;  
  norm_      = wavelet->norm_;
}

// Makes error correlations 
//
//
/*
Wavelet::Wavelet(Wavelet * wavelet, int difftype, int dim)
  : dim_(dim)
{
  gridNI_    = 0;   
  gridNJ_    = 0;
  shiftGrid_ = NULL;  
  gainGrid_  = NULL; 
  if(! wavelet->getIsReal() ) wavelet->invFFT1DInPlace();
  theta_     = wavelet->theta_;
  dz_        = wavelet->dz_;
  nz_        = wavelet->nz_;
  nzp_       = wavelet->nzp_;
  cz_        = wavelet->cz_;
  inFFTorder_= wavelet->inFFTorder_;
  isReal_    = wavelet->isReal_;  
  norm_      = wavelet->norm_;
}

Wavelet::Wavelet(int difftype, int nz, int nzp, int dim)
  : dim_(dim)
{
  gridNI_     = 0;   
  gridNJ_     = 0;
  shiftGrid_  = NULL;  
  gainGrid_   = NULL; 
  theta_      = RMISSING;
  dz_         = RMISSING;
  nz_         = nz;
  nzp_        = nzp;
  cz_         = 0;
  inFFTorder_ = true;
  errCode_    = 0;
}
*/

Wavelet::~Wavelet()
{
  if(shiftGrid_!=NULL)
    delete [] shiftGrid_;
  if(gainGrid_!=NULL)
    delete [] gainGrid_; 
}

Wavelet*  Wavelet::getLocalWavelet(int i, int j)
{
  Wavelet * localWavelet;
  if (dim_ == 1) {
    localWavelet = new Wavelet1D(this);
    float  ts    = this->getLocalTimeshift(i,j);
    float gain   = this->getLocalGainFactor(i,j);
    localWavelet->shiftAndScale(ts,gain);
  }
  else
    localWavelet = new Wavelet3D(this);

  return localWavelet;
}


float  Wavelet::getLocalTimeshift(int i, int j) const
{
  float shift = 0.0f;

  int ind;
  if(shiftGrid_ != NULL && i < gridNI_ &&  i>=0 && j < gridNJ_ &&  j>=0)
  {
    ind = gridNI_*j + i;
    if(shiftGrid_[ind] != WELLMISSING)
      shift = shiftGrid_[ind];
  }


  return shift;
}

float  Wavelet::getLocalGainFactor(int i, int j) const
{
  float gain = 1.0f;
  int ind;
  if(gainGrid_ != NULL && i < gridNI_ &&  i>=0 && j < gridNJ_ &&  j>=0)
  {
    ind = gridNI_*j + i;
    if(gainGrid_[ind] != WELLMISSING)
      gain = gainGrid_[ind];
  }

  return(gain);
}

void
Wavelet::scale(float scale)
{
  if (scale != 1.0f)
    LogKit::LogFormatted(LogKit::LOW,"  Scaling wavelet with gain factor         : %.3e\n",scale);
  scale_ = scale;
}

void 
Wavelet::setShiftGrid(Surface * grid, Simbox * simbox)
{
  gridNI_ = simbox->getnx();
  gridNJ_ = simbox->getny();
  if(shiftGrid_ != NULL)
    delete [] shiftGrid_;
  shiftGrid_ = new float[gridNI_*gridNJ_];
  for(int j=0;j<gridNJ_;j++)
    for(int i=0;i<gridNI_;i++)
    {
      double x, y, z;
      simbox->getCoord(i, j, 0, x, y, z);
      shiftGrid_[i+gridNI_*j] = static_cast<float>(grid->GetZ(x,y));
    }
}

void
Wavelet::printVecToFile(char* fileName,fftw_real* vec, int nzp) const
{
  if( ModelSettings::getDebugLevel() > 0) { 
      char * fName = ModelSettings::makeFullFileName(fileName, ".dat");
      FILE *file = fopen(fName,"w");
      int i;
      for(i=0;i<nzp;i++)
        fprintf(file,"%f\n",vec[i]);

      fclose(file);
      delete [] fName;
    }  
}

void           
Wavelet::fftInv(fftw_complex* cAmp,fftw_real* rAmp,int nt)
{
  rfftwnd_plan p2 = rfftwnd_create_plan(1, &nt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_complex_to_real(p2, cAmp, rAmp);
  fftwnd_destroy_plan(p2);
  double sf = 1.0/double(nt);
  for(int i=0;i<nt;i++)
    rAmp[i]*=fftw_real(sf);
}

void
Wavelet::fft(fftw_real* rAmp,fftw_complex* cAmp,int nt)
{
  rfftwnd_plan p1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_one_real_to_complex(p1, rAmp, cAmp);
  fftwnd_destroy_plan(p1);
}

void 
Wavelet::shiftReal(float shift, fftw_real* rAmp,int nt)
{
  fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);
  fft(rAmp,cAmp, nt);
  int cnzp= nt/2+1;
  float expo;
  fftw_complex tmp,mult;
  for(int i=0;i<cnzp;i++)
  {
    tmp     = cAmp[i];
    expo    = static_cast<float>(-2.0*shift*PI*static_cast<float>(i)/static_cast<float>(nt));
    mult.re = cos(expo);
    mult.im = sin(expo);
    cAmp[i].re = tmp.re*mult.re-tmp.im*mult.im;
    cAmp[i].im = tmp.re*mult.im+tmp.im*mult.re;
  }
  fftInv(cAmp,rAmp, nt);
}

/*void 
Wavelet::shiftReal(int shift, fftw_real* rAmp,int nt)
{
  float* tmp=new float[nt];
  int i,index;
  for(i=0;i<nt;i++)
    tmp[i]=rAmp[i];
  
  for(i=0;i<nt;i++)
  {
    index = i-shift;
    if(index<0)
      index+=nt;
    if(index>=nt)
      index-=nt;

    rAmp[i] = tmp[index];
  }
  delete [] tmp;
}
*/
void
Wavelet::fillInCpp(float* alpha,float* beta,float* rho,int start,int length,fftw_real* cpp_r,int nzp)
{
  int i;

  for(i=0;i<nzp;i++)
    cpp_r[i]=0;

  for(i=start;i < start+length-1;i++)
  {
    float ei1 = computeElasticImpedance(alpha[i],beta[i],rho[i]);
    float ei2 = computeElasticImpedance(alpha[i+1],beta[i+1],rho[i+1]);
    cpp_r[i] =  ei2-ei1;
  } 
}


void
Wavelet::fillInSeismic(float* seisData,int start, int length,fftw_real* seis_r,int nzp) const
{ 
  int i;
  for(i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for(i=start; i<start+length; i++)
  {
    seis_r[i] = seisData[i];
  }
/*
  int lTregion = 3;
  int* modify  = getIndexPrior(start,lTregion,nzp);
  int* conditionto = getIndexPost(start-1,lTregion,nzp);
  //NBNB Odd: interpolate endpoints?
*/

}

void
Wavelet::estimateCor(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* ccor_1_2_c,int cnzp) const
{
  for(int i=0;i<cnzp;i++)
  {
    ccor_1_2_c[i].re = var1_c[i].re*var2_c[i].re+var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }

}

void
Wavelet::convolve(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* out_c,int cnzp) const
{
  for(int i=0;i<cnzp;i++)
  {
    out_c[i].re = var1_c[i].re*var2_c[i].re+var1_c[i].im*var2_c[i].im; 
    out_c[i].im = var1_c[i].im*var2_c[i].re - var1_c[i].re*var2_c[i].im;
  }

}

float
Wavelet::computeElasticImpedance(float vp, float vs, float rho) const
{
  // vp, vs, rho are logtransformed
  float angImp;

  angImp = float(coeff_[0]*vp+coeff_[1]*vs+coeff_[2]*rho );
  
  return(angImp); 
}


void
Wavelet::findContiniousPartOfData(bool* hasData,int nz,int &start, int &length) const
{ 
  int i;
  int lPice=0;
  int lengthMaxPice=-1;
  int startLongestPice=0;
  bool previousHadData = false;

  for(i = 0; i < nz ;i++)
  {
    if(hasData[i])
    {
      if(! previousHadData)
        lPice=1;
      else
        lPice++;
      previousHadData = true;
    }
    else
    {
      if(previousHadData)
      {
        if(lengthMaxPice < lPice)
        {
          lengthMaxPice  = lPice;
          startLongestPice = i-lPice;
        }
      }
      previousHadData=false;
    }
  }

  if(previousHadData)
  {
    if(lengthMaxPice < lPice)
    {
      lengthMaxPice  = lPice;
      startLongestPice = i-lPice;
    }
  }

  start  = startLongestPice;
  length = lengthMaxPice; 
}

/*
int* 
Wavelet::getIndexPrior(int start,int nInd,int nzp)
{
   //returns indexbefore start np cyclic;
  int* index= new int[nInd];
  int i,tmp;
  for(i=0;i<nInd;i++)
  {
    tmp=start-i-1;
    if(tmp<0)
      tmp = tmp+nzp;
    if(tmp>=nzp)
      tmp = tmp-nzp;
    index[i]= tmp;
  }
  return index;
}

int* 
Wavelet::getIndexPost(int start,int nInd,int nzp)
{
  //returns index after start np cyclic;
  int* index = new int[nInd];
  int i,tmp;
  for(i=0;i<nInd;i++)
  {
    tmp=start+i+1;
    if(tmp<0)
      tmp = tmp+nzp;
    if(tmp>=nzp)
      tmp = tmp-nzp;
    index[i]= tmp;
  }
  return index;
}
*/

int Wavelet::getWaveletLengthI()
{
  bool trans=false;
  if(isReal_==false)
  {
    invFFT1DInPlace();
    trans=true;
  }
  
  float maxAmp =  fabs(getRAmp(0)); // gets max amp 
  for(int i=1;i <nzp_;i++)
    if(fabs(getRAmp(i)) > maxAmp)
      maxAmp = fabs(getRAmp(i));

  float minAmp= maxAmp*minRelativeAmp_; // minimum relevant amplitude

  int wLength=nzp_;

  for(int i=nzp_/2;i>0;i--)
  {
    if(fabs(getRAmp(i)) >minAmp)
    {
      wLength= (i*2+1);// adds both sides 
      break;
    }
    if(fabs(getRAmp(nzp_-i)) > minAmp)
    {
      wLength= (2*i+1);// adds both sides 
      break;
    }
  }
  wLength =MINIM(wLength,2*(nzp_/2)-1); // allways even number
  if(trans==true)
    fft1DInPlace();

  return wLength;
}

float
Wavelet::getWaveletLengthF()
{
  return dz_*float( getWaveletLengthI() );
}

float         
Wavelet::getNoiseStandardDeviation(Simbox * simbox, FFTGrid * seisCube, WellData ** wells, int nWells, char *errText, int &error)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"\n  Estimating noise from seismic data and (nonfiltered) blocked wells");
  float errStd  = 0.0f;
  float dataVar = 0.0f;
  // initialization
  scale_=1; 
  gridNI_=0;   
  gridNJ_=0;
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  
  float * dz = new float[nWells];
  int i, k, w;

  int   nz   = simbox->getnz();
  float dz0  = static_cast<float>(simbox->getdz());
  int   nzp  = seisCube->getNzp();
  int   cnzp = (nzp/2+1);
  int   rnzp = 2*cnzp;

  float * alpha    = new float[nz];
  float * beta     = new float[nz];
  float * rho      = new float[nz];
  float * seisData = new float[nz];
  bool  * hasData  = new bool[nz];

  //Noise estimation
  fftw_real**       cpp_r = new fftw_real*[nWells];
  fftw_complex**    cpp_c = reinterpret_cast<fftw_complex**>(cpp_r);
  
  fftw_real**       seis_r = new fftw_real*[nWells];
  fftw_complex**    seis_c = reinterpret_cast<fftw_complex**>(seis_r); 

  fftw_real**       synt_r = new fftw_real*[nWells];
  fftw_complex**    synt_c = reinterpret_cast<fftw_complex**>(synt_r);

  fftw_real**       cor_seis_synt_r = new fftw_real*[nWells];
  fftw_complex**    cor_seis_synt_c = reinterpret_cast<fftw_complex**>(cor_seis_synt_r); 

  //fftw_real**       err_r = new fftw_real*[nWells];
//  fftw_complex**    err_c = (fftw_complex** ) err_r;  //Useful for debug

  fftw_real**      wavelet_r = new fftw_real*[nWells];
  fftw_complex**   wavelet_c = reinterpret_cast<fftw_complex**>(wavelet_r);
  float* shiftWell   = new float[nWells];
  float* errVarWell  = new float[nWells];
  float* dataVarWell = new float[nWells];
  int*   nActiveData = new int[nWells];

  int maxBlocks = 0;
  //float wShift  = getWaveletShift();
  for(i=0;i<nWells;i++)
  {
    cpp_r[i]           = new fftw_real[rnzp];
    seis_r[i]          = new fftw_real[rnzp];
    cor_seis_synt_r[i] = new fftw_real[rnzp];
    wavelet_r[i]       = new fftw_real[rnzp];
    synt_r[i]          = new fftw_real[rnzp];
    nActiveData[i]     = 0;
    errVarWell[i]      = 0.0f;
    dataVarWell[i]     = 0.0f;
    const int * ipos   = wells[i]->getBlockedLogsPropThick()->getIpos();
    const int * jpos   = wells[i]->getBlockedLogsPropThick()->getJpos();
    dz[i]              = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])*simbox->getdz());
    int nBlocks        = wells[i]->getBlockedLogsPropThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (w = 0 ; w < nWells ; w++) 
  {
    if (wells[w]->getUseForWaveletEstimation())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsPropThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bl->getVerticalTrend(bl->getAlpha(), alpha);
      bl->getVerticalTrend(bl->getBeta(), beta);
      bl->getVerticalTrend(bl->getRho(), rho);
      bl->getVerticalTrend(seisLog, seisData);

      for (k = 0 ; k < nz ; k++) {
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      }

      int start,length;
      findContiniousPartOfData(hasData,nz,start,length);

      if(length*dz0 > waveletLength_) // must have enough data
      {
        fillInCpp(alpha,beta,rho,start,length,cpp_r[w],nzp);  // fills in reflection coefficients
        fft(cpp_r[w],cpp_c[w],nzp);
        fillInnWavelet(wavelet_r[w],nzp,dz[w]); // fills inn wavelet
        //flipVec(wavelet_r[w],nzp);
        fft(wavelet_r[w],wavelet_c[w],nzp);
        convolve(cpp_c[w],wavelet_c[w],synt_c[w],cnzp);
        fillInSeismic(seisData,start, length,seis_r[w],nzp);
        fft(seis_r[w],seis_c[w],nzp);
        estimateCor(synt_c[w],seis_c[w],cor_seis_synt_c[w],cnzp);
        fftInv(cor_seis_synt_c[w],cor_seis_synt_r[w],nzp);
        float shift=findBulkShift(cor_seis_synt_r[w],dz[w], nzp);
        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w]=shift;
        fftInv(synt_c[w],synt_r[w],nzp);
        shiftReal(-shift/dz[w],synt_r[w],nzp);
        fillInSeismic(seisData,start, length,seis_r[w],nzp);
        if(ModelSettings::getDebugLevel() > 0)
        {
          char* fileName = new char[MAX_STRING];
          sprintf(fileName,"seismic_Well_%d_%d",w,int(floor(theta_/PI*180.0+0.5)));
          printVecToFile(fileName,seis_r[w], nzp);
          sprintf(fileName,"synthetic_seismic_Well_%d_%d",w,int(floor(theta_/PI*180.0+0.5)));
          printVecToFile(fileName,synt_r[w], nzp);
          sprintf(fileName,"wellTime_Well_%d",w);
          char * fName = ModelSettings::makeFullFileName(fileName, ".dat");
          FILE *file = fopen(fName,"wb");
          //fprint(file,"%f %f  %f \n",,dz[w],);
          fclose(file);
          delete [] fileName;
          delete fName;
        }
        for(i=start;i<start+length;i++)
        { 
          float err=(seis_r[w][i] - synt_r[w][i]);
          errVarWell[w]+=err*err;
          dataVarWell[w]+=seis_r[w][i] *seis_r[w][i] ;
        }
        nActiveData[w]=length;
      }
      else
      {
        LogKit::LogFormatted(LogKit::LOW,"\n  Not using vertical well %s for error estimation (length=%.1fms  required length=%.1fms).",
                         wells[w]->getWellname(),length*dz0,waveletLength_);
      }
    }
  }

  float * scaleOptWell    = new float[nWells];
  float * errWellOptScale = new float[nWells];
  float * errWell         = new float[nWells];

  float errOptScale;
  float optScale = findOptimalWaveletScale(synt_r,seis_r,nWells,nzp,dataVarWell,
                                           errOptScale,errWell,scaleOptWell,errWellOptScale);
  delete [] seisLog;
  delete [] dz;
  delete [] hasData;

  int nData=0;
  for(i=0;i<nWells;i++)
  {
    nData   += nActiveData[i];
    errStd  += errVarWell[i];
    dataVar += dataVarWell[i];
    if(nActiveData[i]>0)
    {    
      errVarWell[i]  /= nActiveData[i];
      dataVarWell[i] /= nActiveData[i];
    }
  }

  if (nData == 0)
  {
    sprintf(errText, "%s Cannot estimate signal-to-noise ratio. No legal well data available.\n", errText);
    error += 1;
  }
  dataVar /= float(nData);
  errStd  /= float(nData);
  errStd   = sqrt(errStd);
  
  LogKit::LogFormatted(LogKit::MEDIUM,"\n  Reporting errors (as standard deviations) estimated in different ways:\n\n");

  if (readtype_ == ESTIMATE)
  {
    LogKit::LogFormatted(LogKit::LOW,"\n");
    LogKit::LogFormatted(LogKit::LOW,"                                     SeisData       OptimalGlobal      OptimalLocal\n");
    LogKit::LogFormatted(LogKit::LOW,"  Well                  shift[ms]     StdDev         Gain   S/N         Gain   S/N \n");
    LogKit::LogFormatted(LogKit::LOW,"  ----------------------------------------------------------------------------------\n");
    for(i=0;i<nWells;i++)
    {
      if(nActiveData[i]>0) {
        float SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
        float SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);
        LogKit::LogFormatted(LogKit::LOW,"  %-20s   %6.2f     %9.2e      %6.2f %6.2f      %6.2f %6.2f\n", 
                             wells[i]->getWellname(),shiftWell[i],sqrt(dataVarWell[i]),
                             optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
      }
      else
        LogKit::LogFormatted(LogKit::LOW,"  %-20s      -            -             -      -           -      -\n",
                             wells[i]->getWellname()); 
    }
  }
  else
  {
    LogKit::LogFormatted(LogKit::LOW,"\n");
    LogKit::LogFormatted(LogKit::LOW,"                                     SeisData        ActuallyUsed       OptimalGlobal      OptimalLocal\n");
    LogKit::LogFormatted(LogKit::LOW,"  Well                  shift[ms]     StdDev          Gain   S/N         Gain   S/N         Gain   S/N \n");
    LogKit::LogFormatted(LogKit::LOW,"  ------------------------------------------------------------------------------------------------------\n");
    for(i=0;i<nWells;i++)
    {
      if(nActiveData[i]>0) {
        float SNActuallyUsed  = dataVarWell[i]/errVarWell[i];
        float SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
        float SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);
        LogKit::LogFormatted(LogKit::LOW,"  %-20s   %6.2f     %9.2e        1.00 %7.2f     %6.2f %7.2f     %6.2f %7.2f\n", 
                             wells[i]->getWellname(),shiftWell[i],sqrt(dataVarWell[i]),
                             SNActuallyUsed,optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
      }
      else
        LogKit::LogFormatted(LogKit::LOW,"  %-20s      -            -             -      -           -      -           -      -   \n",
                             wells[i]->getWellname()); 
    }
  }  

  delete [] shiftWell;
  delete [] errVarWell;
  delete [] dataVarWell;
  delete [] nActiveData;
  delete [] scaleOptWell;
  delete [] errWellOptScale;
  delete [] errWell;

  float empSNRatio = dataVar/(errStd*errStd);
  LogKit::LogFormatted(LogKit::LOW,"\n  Signal to noise ratio used for this angle stack is: %6.2f\n", empSNRatio);

  if (empSNRatio < 1.1f) 
  {
    sprintf(errText, "%sThe empirical signal-to-noise ratio Var(seismic data)/Var(noise) is %.2f. Ratios\n",errText, empSNRatio);
    sprintf(errText, "%s smaller than 1 are illegal and CRAVA has to stop. CRAVA was for some reason not able\n", errText);
    sprintf(errText,"%s to estimate this ratio reliably, and you must give it as input to the model file.\n", errText);
    sprintf(errText,"%s \nNOTE: If the wavelet was estimated by CRAVA the solution may be to remove one or more wells\n", errText);
    sprintf(errText,"%s from the wavelet estimation (compare shifts and SN-ratios for different wells).\n\n", errText);
    error += 1;
  }
 
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  for(i=0;i<nWells;i++)
  {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    //delete [] err_r[i] ;
    delete [] synt_r[i] ;
    delete [] wavelet_r[i];
    delete [] cor_seis_synt_r[i];
  }
  delete [] cpp_r;
  delete [] seis_r;
  //delete [] err_r;
  delete [] synt_r;
  delete [] wavelet_r;
  delete [] cor_seis_synt_r;
      
  return errStd;
}

float          
Wavelet::findOptimalWaveletScale(fftw_real ** synt_seis_r,
                                 fftw_real ** seis_r,
                                 int          nWells,
                                 int          nzp,
                                 float      * wellWeight,
                                 float      & err, // NBNB-PAL: Det er uheldig å returnere err slik
                                 float      * errWell,
                                 float      * scaleOptWell,
                                 float      * errWellOptScale) const
{
  float   optScale   = 1.0f;
  float   scaleLimit = 3.0f;
  int     nScales    = 51; // should be odd to include 1.00
  float * scales     = new float[nScales];
  float * error      = new float[nScales];

  for(int i=0;i<nScales;i++)
  {
    scales[i] = exp(-log(scaleLimit)+i*2*(log(scaleLimit))/(nScales-1));
    error[i]  = 0.0f;
  }

  int    * counter  = new int[nWells];
  float  * seisNorm = new float[nWells];
  float ** resNorm  = new float*[nWells];

  for(int i=0;i<nWells;i++)
  {
    resNorm[i]  = new float[nScales];
    seisNorm[i] = 0.0f;
  }

  float minSeisAmp = 1e-7;
  int totCount=0;

  for(int i=0;i<nWells;i++)
  {
    counter[i]=0;
    if(wellWeight[i]>0)
    {
      // Count number of layers with seismic data
      for(int k=0;k<nzp;k++)
        if(fabs(seis_r[i][k]) > minSeisAmp)
          counter[i]++;
      totCount+=counter[i];

      for(int j=0;j<nScales;j++)
      {
        resNorm[i][j]=0.0;
        for(int k=0;k<nzp;k++)
        {
          if(fabs(seis_r[i][k]) > minSeisAmp)
          {
            seisNorm[i]   += seis_r[i][k] * seis_r[i][k];
            float      foo = scales[j]*synt_seis_r[i][k] - seis_r[i][k];
            resNorm[i][j] += foo*foo;
          }
          error[j]+=resNorm[i][j];
        }
      }
    }//if
  }

  int   optInd=0; 
  float optValue=error[0];
  for(int i=1;i<nScales;i++)
    if(error[i]<optValue)
    {
      optValue=error[i];
      optInd=i;
    }

    delete [] error;

    err = sqrt(optValue/static_cast<float>(totCount));
    optScale = scales[optInd];

    for(int i=0;i<nWells;i++)
    {
      if(counter[i]>0)
        errWell[i] = sqrt(resNorm[i][optInd]/counter[i]);
      else
        errWell[i] = 0.0f;
    }

    for(int i=0;i<nWells;i++)
    {
      if(wellWeight[i]>0)
      {
        optValue = resNorm[i][0];
        optInd=0;
        for(int j=1;j<nScales;j++)
          if(resNorm[i][j]<optValue)
          {
            optValue=resNorm[i][j];
            optInd=j;
          }
          scaleOptWell[i]    = scales[optInd];
          errWellOptScale[i] = sqrt(optValue/float(counter[i]));
      }
      else
      {
        scaleOptWell[i]    = 0.0f;
        errWellOptScale[i] = 0.0f;
      }
    }

    for (int i=0; i<nWells; i++)
      delete [] resNorm[i];
    delete [] resNorm;
    delete [] seisNorm;
    delete [] counter;
    delete [] scales;

    return optScale;
}

float
Wavelet::findBulkShift(fftw_real* vec_r,float dz,int nzp)
{
  
  
  float shift=0.0f;
  
  float sum=0;
  int i,polarity;
  // if the sum from -maxShift_ to maxShift_ ms is 
  // positive then polarity is positive   
  for(i=0;i<ceil(maxShift_/dz);i++)//zero included
    sum+=vec_r[i];
  for(i=0;i<floor(maxShift_/dz);i++)
    sum+=vec_r[nzp-i-1];
    
  polarity=-1;
  if(sum > 0)
    polarity=1;

  // gets optimal shift
  float maxValue;
  float shiftF;
  int shiftI;
  float f1,f2,f3;

 
  maxValue = 0.0f;
  shiftI=0;

  for(i=0;i<ceil(maxShift_/dz);i++)
  {
    if(vec_r[i]*polarity > maxValue)
    {
      maxValue = vec_r[i]*polarity;
      shiftI = i;
    }
  }
  for(i=0;i<floor(maxShift_/dz);i++)
  {
    if(vec_r[nzp-1-i]*polarity > maxValue)
    {
      maxValue = vec_r[nzp-1-i]*polarity;
      shiftI = -1-i;
    }
  }
  if(shiftI < 0)
  {
    if(vec_r[nzp+shiftI-1]*polarity < maxValue) //then local max
    {
      f1 = vec_r[nzp+shiftI-1];
      f2 = vec_r[nzp+shiftI];
      int ind3;
      if(shiftI==-1)
        ind3 = 0;
      else
        ind3=nzp+shiftI+1;
      f3 = vec_r[ind3];
      float x0=(f1-f3)/(2*(f1+f3-2*f2));
      shiftF=float(shiftI)+x0;
    }
    else  // do as good as we can
      shiftF=float(shiftI);
  }
  else //positive or zero shift
  {
    if(vec_r[shiftI+1]*polarity < maxValue) //then local max
    {
      f3 = vec_r[shiftI+1];
      f2 = vec_r[shiftI];
      int ind1;
      if(shiftI==0)
        ind1 = nzp-1;
      else
        ind1=shiftI-1;
      f1 = vec_r[ind1];
      float x0=(f1-f3)/(2*(f1+f3-2*f2));
      shiftF=shiftI+x0;
    }
    else  // do as good as we can
      shiftF=float(shiftI);
  }
  shift = shiftF*dz;
  shiftReal(-shiftF,vec_r,nzp);// for testing
 
  return shift;
}

void
Wavelet::fillInnWavelet(fftw_real* wavelet_r,int nzp,float dz)
{
  //Note: Assumes that dz_ > dz
  // NBNB OddK 
  fftw_real previous1 = getRAmp(0);
  fftw_real previous2 = getRAmp(0);
  fftw_real current1 = 0.0f;
  fftw_real current2 = 0.0f;
  wavelet_r[0]        = previous1;
  //wavelet_r[nzp/2]    = 0;
  //wavelet_r[nzp/2+1]    = 0;
  int counterForWavelet   = 1;

  float w;
  int i;
  for(i=1;i<= nzp/2;i++)
  {
    if(counterForWavelet*dz_ < i*dz)
    {
      counterForWavelet++;
      previous1=current1;
      previous2=current2;
    }
    current1=getRAmp(counterForWavelet);
    current2=getRAmp(nzp-counterForWavelet);
    w = (counterForWavelet*dz_-i*dz) /dz_;
    wavelet_r[i]     = (1-w)*current1+ w*previous1;
    wavelet_r[nzp-i] = (1-w)*current2+ w*previous2;
  }
}


void 
Wavelet::setGainGrid(Surface * grid, Simbox * simbox)
{
  double sum = 0.0;
  int nData = 0;
  gridNI_ = simbox->getnx();
  gridNJ_ = simbox->getny();
  if(gainGrid_ != NULL)
    delete [] gainGrid_;
  gainGrid_ = new float[gridNI_*gridNJ_];
  for(int j=0;j<gridNJ_;j++)
    for(int i=0;i<gridNI_;i++)
    {
      double x, y, z;
      simbox->getCoord(i, j, 0, x, y, z);
      double value = grid->GetZ(x,y);
      gainGrid_[i+gridNI_*j] = static_cast<float>(value);
      if(value != WELLMISSING)
      {
        sum += log(value);
        nData++;
      }
    }
    float invGeoMean = float(exp(-sum/static_cast<double>(nData)));
    for(int j=0;j<gridNJ_;j++)
      for(int i=0;i<gridNI_;i++)
      {
        if(gainGrid_[i+gridNI_*j] != WELLMISSING)
          gainGrid_[i+gridNI_*j] *= invGeoMean;
      }
      float geoMean = 1/invGeoMean;
      norm_ *= geoMean;
      bool fftflag = isReal_;
      if(fftflag == true)
        invFFT1DInPlace();
      for(int i=0;i<nzp_;i++) {
        float rAmp = getRAmp(i);
        rAmp *= geoMean;
        setRAmp(rAmp,i);
        //        rAmp_[i] *= geoMean;
      }
        if(fftflag == true)
          fft1DInPlace();
}


/*void           
Wavelet::flipVec(fftw_real* vec, int n)
{
  int i;
  fftw_real* tmp= new fftw_real[n];
  for(i=1;i<n;i++)
  {
    tmp[i]=vec[i];
  }
  
  for(i=1;i<n;i++)
  {
    vec[i]=tmp[n-i];
  }
  delete [] tmp;
}
*/
