#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"

#include "nrlib/surface/regularsurface.hpp"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/wavelet3D.h"
#include "src/modelsettings.h"
#include "src/model.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/definitions.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/vario.h"
#include "src/krigingdata2d.h"
#include "src/kriging2d.h"
#include "src/covgrid2d.h"
#include "src/io.h"

Wavelet::Wavelet(int dim)
  : shiftGrid_(NULL),
    gainGrid_(NULL),
    dim_(dim)
{
}

Wavelet::Wavelet(int dim, float * reflCoef)
  : shiftGrid_(NULL),
    gainGrid_(NULL),
    dim_(dim)
{
  coeff_[0] = reflCoef[0];
  coeff_[1] = reflCoef[1];
  coeff_[2] = reflCoef[2];
}

Wavelet::Wavelet(ModelSettings * modelSettings, int dim, float * reflCoef)
  : inFFTorder_(false),
    isReal_(true),
    maxShift_(modelSettings->getMaxWaveletShift()),
    minRelativeAmp_(modelSettings->getMinRelWaveletAmp()),
    scale_(1),
    shiftGrid_(NULL),
    gainGrid_(NULL),
    dim_(dim)
{
  coeff_[0]       = reflCoef[0];
  coeff_[1]       = reflCoef[1];
  coeff_[2]       = reflCoef[2];
}

Wavelet::Wavelet(Wavelet * wavelet, int dim)
  : shiftGrid_(NULL),
    gainGrid_(NULL),
    dim_(dim)
{
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

Wavelet::~Wavelet()
{
  if(shiftGrid_!=NULL)
    delete shiftGrid_;
  if(gainGrid_!=NULL)
    delete gainGrid_; 
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

  if(shiftGrid_ != NULL && i < static_cast<int>(shiftGrid_->GetNI()) && i>=0 && 
     j < static_cast<int>(shiftGrid_->GetNJ()) &&  j>=0)
  {
   // ind = gridNI_*j + i;
    if((*shiftGrid_)(i,j) != WELLMISSING)
      shift = float((*shiftGrid_)(i,j));
  }


  return shift;
}

float  Wavelet::getLocalGainFactor(int i, int j) const
{
  float gain = 1.0f;
  
  if(gainGrid_ != NULL && i < static_cast<int>(gainGrid_->GetNI()) &&  i>=0 && 
    j < static_cast<int>(gainGrid_->GetNJ()) &&  j>=0)
  {
    //ind = gridNI_*j + i;
    if((*gainGrid_)(i,j) != WELLMISSING)
      gain = float((*gainGrid_)(i,j));
  }

  return(gain);
}

void
Wavelet::scale(float scale)
{
  if (scale != 1.0f)
    LogKit::LogFormatted(LogKit::LOW,"  Scaling wavelet with factor         : %.3e\n",scale);
  scale_ = scale;
}

void 
Wavelet::setShiftGrid(Grid2D *grid)
{
  
  if(shiftGrid_ != NULL)
    delete [] shiftGrid_;
  shiftGrid_ = grid;
  
}

void
Wavelet::printVecToFile(const std::string & fileName, fftw_real* vec, int nzp) const
{
  if( ModelSettings::getDebugLevel() > 0) { 
    std::string fName = fileName + IO::SuffixGeneralData();
    fName = IO::makeFullFileName(IO::PathToWavelets(), fName);
    std::ofstream file;
    NRLib::OpenWrite(file,fName);
    for(int i=0;i<nzp;i++)
      file << vec[i] << "\n";
    file.close();
  }  
}

void 
Wavelet::shiftReal(float shift, fftw_real* rAmp,int nt)
{
  fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);
  Utils::fft(rAmp,cAmp, nt);
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
  Utils::fftInv(cAmp,rAmp, nt);
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
  wLength =MINIM(wLength,2*((nzp_+1)/2) - 1); // always odd number
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
Wavelet::calculateSNRatioAndLocalWavelet(Simbox        * simbox, 
                                         FFTGrid       * seisCube, 
                                         WellData     ** wells, 
                                         Grid2D       *& shift, 
                                         Grid2D       *& gain, 
                                         ModelSettings * modelSettings,
                                         char          * errText, 
                                         int           & error, 
                                         Grid2D       *& noiseScaled, 
                                         int             number, 
                                         float           globalScale)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"\n  Estimating noise from seismic data and (nonfiltered) blocked wells");
  float errStd  = 0.0f;
  float dataVar = 0.0f;
  // initialization
 // scale_=1; 
  
  shiftGrid_=NULL;  
  gainGrid_=NULL; 
  
  Vario  * localWaveletVario = modelSettings->getLocalWaveletVario();
  int      nWells            = modelSettings->getNumberOfWells();
 // bool     useLocalWavelet   = modelSettings->getUseLocalWavelet();
 // int      outputFormat      = modelSettings->getOutputFormatFlag();
 // int      otherOutput       = modelSettings->getOtherOutputFlag();
  bool     doEstimateLocalShift = modelSettings->getEstimateLocalShift(number);
  bool     doEstimateLocalScale = modelSettings->getEstimateLocalScale(number);
  bool     doEstimateLocalNoise = modelSettings->getEstimateLocalNoise(number);
  bool     doEstimateGlobalScale = modelSettings->getEstimateGlobalWaveletScale(number);
  bool     doEstimateSNRatio = modelSettings->getEstimateSNRatio(number);
// bool     estimationMode = modelSettings->getEstimationMode();

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
  fftw_real    ** cpp_r = new fftw_real*[nWells];
  fftw_complex ** cpp_c = reinterpret_cast<fftw_complex**>(cpp_r);
  
  fftw_real    ** seis_r = new fftw_real*[nWells];
  fftw_complex ** seis_c = reinterpret_cast<fftw_complex**>(seis_r); 

  fftw_real    ** synt_r = new fftw_real*[nWells];
  fftw_complex ** synt_c = reinterpret_cast<fftw_complex**>(synt_r);

  fftw_real    ** cor_seis_synt_r = new fftw_real*[nWells];
  fftw_complex ** cor_seis_synt_c = reinterpret_cast<fftw_complex**>(cor_seis_synt_r); 

  //fftw_real**       err_r = new fftw_real*[nWells];
  //fftw_complex**    err_c = (fftw_complex** ) err_r;  //Useful for debug

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
    const int * ipos   = wells[i]->getBlockedLogsOrigThick()->getIpos();
    const int * jpos   = wells[i]->getBlockedLogsOrigThick()->getJpos();
    dz[i]              = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])*simbox->getdz());
    int nBlocks        = wells[i]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
 // seisCube->setAccessMode(FFTGrid::RANDOMACCESS);
  for (w = 0 ; w < nWells ; w++) 
  {
    if (wells[w]->getUseForWaveletEstimation())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
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
      bl->findContiniousPartOfData(hasData,nz,start,length);

      if(length*dz0 > waveletLength_) // must have enough data
      {
        bl->fillInCpp(coeff_,start,length,cpp_r[w],nzp);  // fills in reflection coefficients
        Utils::fft(cpp_r[w],cpp_c[w],nzp);
        fillInnWavelet(wavelet_r[w],nzp,dz[w]); // fills inn wavelet
        //flipVec(wavelet_r[w],nzp);
        Utils::fft(wavelet_r[w],wavelet_c[w],nzp);
        convolve(cpp_c[w],wavelet_c[w],synt_c[w],cnzp);
        bl->fillInSeismic(seisData,start,length,seis_r[w],nzp);
        Utils::fft(seis_r[w],seis_c[w],nzp);
        bl->estimateCor(synt_c[w],seis_c[w],cor_seis_synt_c[w],cnzp);
        Utils::fftInv(cor_seis_synt_c[w],cor_seis_synt_r[w],nzp);
        //Estimate shift. Do not run if shift given, use given shift.        
        float shift=findBulkShift(cor_seis_synt_r[w],dz[w], nzp);
        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w]=shift;
        Utils::fftInv(synt_c[w],synt_r[w],nzp);
        shiftReal(-shift/dz[w],synt_r[w],nzp);
        bl->fillInSeismic(seisData,start,length,seis_r[w],nzp);
        if(ModelSettings::getDebugLevel() > 0)
        {
          std::string angle = NRLib::ToString(theta_/(M_PI*180.0),1);
          std::string fileName;
          fileName = "seismic_Well_" + NRLib::ToString(w) + "_" + angle;
          printVecToFile(fileName,seis_r[w], nzp);
          fileName = "synthetic_seismic_Well_" + NRLib::ToString(w) + "_" + angle;
          printVecToFile(fileName,synt_r[w], nzp);
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
  for(i=0;i<nWells;i++)
    scaleOptWell[i] = -1.0;
  float * errWellOptScale = new float[nWells];
  float * errWell         = new float[nWells];
  bool writelog = false;
  float errOptScale = 1.0;
  //Estimate global scale, local scale and error.
  //If global scale given, do not use return value. Do kriging with global scale as mean.
  //If local scale given, run separate routine to find local noise if wanted.
  float optScale;
  if(doEstimateLocalScale==true || doEstimateGlobalScale==true)
  {
    optScale = findOptimalWaveletScale(synt_r,seis_r,nWells,nzp,dataVarWell,
                                           errOptScale,errWell,scaleOptWell,errWellOptScale);
    
    writelog = true;
    if(doEstimateGlobalScale==false)
      optScale = globalScale;
    else
    {
      scale(optScale);
      for(i=0;i<nWells;i++)
        scaleOptWell[i]/=optScale;
    }
  }
  // local scale given means gain !=NULL
  else if(doEstimateLocalNoise==true && gain!=NULL)
  {
    optScale = globalScale; // Global scale must be given if local scale is given
    findLocalNoiseWithGainGiven(synt_r,seis_r,nWells,nzp,dataVarWell, errOptScale, errWell, errWellOptScale, scaleOptWell,gain, wells, simbox);
    writelog = true;
  }
  else
  {
    optScale = globalScale;
    // only for loging
    findOptimalWaveletScale(synt_r,seis_r,nWells,nzp,dataVarWell,
                                           errOptScale,errWell,scaleOptWell,errWellOptScale);
  }
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

   // if (readtype_ == ESTIMATE)
  //  {
      LogKit::LogFormatted(LogKit::LOW,"\n");
      LogKit::LogFormatted(LogKit::LOW,"                                     SeisData       OptimalGlobal      OptimalLocal\n");
      LogKit::LogFormatted(LogKit::LOW,"  Well                  shift[ms]     StdDev         Gain   S/N         Gain   S/N \n");
      LogKit::LogFormatted(LogKit::LOW,"  ----------------------------------------------------------------------------------\n");
      for(i=0;i<nWells;i++)
      {
        if(nActiveData[i]>0)
         {
           float SNOptimalGlobal, SNOptimalLocal;    
            SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
            SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);   
            LogKit::LogFormatted(LogKit::LOW,"  %-20s   %6.2f     %9.2e      %6.2f %6.2f      %6.2f %6.2f\n", 
            wells[i]->getWellname(),shiftWell[i],sqrt(dataVarWell[i]),
            optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
        }
        else
          LogKit::LogFormatted(LogKit::LOW,"  %-20s      -            -             -      -           -      -\n",
          wells[i]->getWellname()); 
      }
      for(i=0;i<nWells;i++)
      {
        if((scaleOptWell[i]>=3.0 || scaleOptWell[i]<=0.3334) && nActiveData[i]>0)
        {
          LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: The well %s has a optimal local gain value indicating that this well should not be used for wavelet estimation\n",
          wells[i]->getWellname());
        }
      }
 //   }
  /*  else //commented out because Actually used is wrong.
    {
      LogKit::LogFormatted(LogKit::LOW,"\n");
      LogKit::LogFormatted(LogKit::LOW,"                                     SeisData        ActuallyUsed       OptimalGlobal      OptimalLocal\n");
      LogKit::LogFormatted(LogKit::LOW,"  Well                  shift[ms]     StdDev          Gain   S/N         Gain   S/N         Gain   S/N \n");
      LogKit::LogFormatted(LogKit::LOW,"  ------------------------------------------------------------------------------------------------------\n");
      for(i=0;i<nWells;i++)
      {
        if(nActiveData[i]>0) {
          float SNActuallyUsed  = dataVarWell[i]/errVarWell[i];
          float SNOptimalGlobal, SNOptimalLocal;
          SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
          SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);
       
          LogKit::LogFormatted(LogKit::LOW,"  %-20s   %6.2f     %9.2e        1.00 %7.2f     %6.2f %7.2f     %6.2f %7.2f\n", 
            wells[i]->getWellname(),shiftWell[i],sqrt(dataVarWell[i]),
            SNActuallyUsed,optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
        }
        else
          LogKit::LogFormatted(LogKit::LOW,"  %-20s      -            -             -      -           -      -           -      -   \n",
          wells[i]->getWellname()); 
      }
    }*/
//  }

 // if(useLocalWavelet && (shift==NULL || gain==NULL))
 // {
 //   estimateLocalWavelets(shift, gain, shiftWell, scaleOptWell,
 //                         localWaveletVario, nActiveData, simbox, 
 //                         wells, nWells, outputFormat, otherOutput);
  if(doEstimateLocalScale==true)
  {
// Estimate global noise with local waveletscale
    dataVar = 0.0;
    errStd = 0.0;
    for(i=0;i<nWells;i++)
    {
      dataVar+=(dataVarWell[i]*nActiveData[i]);
      errStd+=(errWellOptScale[i]*errWellOptScale[i]*nActiveData[i]);
    }
    dataVar/=nData;
    errStd/=nData;
    errStd = sqrt(errStd);
  }
  else if(doEstimateGlobalScale==true)
   errStd = errOptScale;
 //   
 // }
  

  if(doEstimateLocalShift || doEstimateLocalScale || doEstimateLocalNoise) 
  {
    //
    // Pretabulate correlations
    //
    const CovGrid2D cov(localWaveletVario, 
                        simbox->getnx(),
                        simbox->getny(),
                        simbox->getdx(), 
                        simbox->getdy());
    
    if (ModelSettings::getDebugLevel() > 0) {
      std::string baseName = std::string("Local_Wavelet_Correlation") + IO::SuffixAsciiIrapClassic();
      std::string fileName = IO::makeFullFileName(IO::PathToWavelets(), baseName);
      cov.writeToFile(fileName);
    }

    if(doEstimateLocalShift)
      estimateLocalShift(cov, shift, shiftWell, nActiveData, simbox,wells, nWells);
    
    if(doEstimateLocalScale)
      estimateLocalGain(cov, gain, scaleOptWell, 1.0, nActiveData, simbox,wells, nWells);

    if(doEstimateLocalNoise)
    {
      float errStdLN;
      if(doEstimateSNRatio==true)
        errStdLN = errStd;
      else //SNRatio given in model file
        errStdLN = sqrt(dataVar/modelSettings->getSNRatio(number));
      if(gain==NULL && doEstimateLocalScale==false && doEstimateGlobalScale==false) // No local wavelet scale
      {
        for(i=0;i<nWells;i++)
        {
          errVarWell[i] = sqrt(errVarWell[i]);
        }
        estimateLocalNoise(cov, noiseScaled, errStdLN, errVarWell, nActiveData, simbox,wells, nWells); 
      }
      else if(doEstimateGlobalScale==true && doEstimateLocalScale==false) // global wavelet scale
        estimateLocalNoise(cov, noiseScaled, errStdLN,errWell, nActiveData, simbox,wells, nWells); 
      else
        estimateLocalNoise(cov, noiseScaled, errStdLN, errWellOptScale, nActiveData, simbox,wells, nWells); 
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
  if(doEstimateSNRatio==true)
    LogKit::LogFormatted(LogKit::LOW,"\n  Signal to noise ratio used for this angle stack is: %6.2f\n", empSNRatio);

  if (empSNRatio < 1.1f) 
  {
    LogKit::LogFormatted(LogKit::WARNING,"\nERROR: The empirical signal-to-noise ratio Var(data)/Var(noise) is %.2f. Ratios smaller",empSNRatio);
    LogKit::LogFormatted(LogKit::WARNING,"\n       than 1.1 are not acceptable. The signal-to-noise ratio was not reliably estimated");
    LogKit::LogFormatted(LogKit::WARNING,"\n       and you must give it as input in the model file.\n");
    LogKit::LogFormatted(LogKit::WARNING,"\n       If the wavelet was estimated by CRAVA the solution may be to remove one or more wells");
    LogKit::LogFormatted(LogKit::WARNING,"\n       from the wavelet estimation (compare shifts and SN-ratios for different wells).\n");

    sprintf(errText, "%sInvalid signal-to-noise ratio obtained for the angle-gather of %.1f degrees.\n",
            errText, static_cast<float>(180.0/M_PI)*seisCube->getTheta());
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
 // seisCube->endAccess();    
  return empSNRatio;
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

  float minSeisAmp = static_cast<float> (1e-7);
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
        }
        error[j]+=resNorm[i][j];
        
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

void Wavelet::
findLocalNoiseWithGainGiven(fftw_real ** synt_seis_r,
                            fftw_real ** seis_r,
                            int nWells,
                            int nzp,
                            float * wellWeight,
                            float & err,
                            float * errWell,
                            float      * scaleOptWell,
                            float * errWellOptScale, 
                            Grid2D * gain, 
                            WellData **wells, Simbox *simbox) const
{
  double *scale = new double[nWells];   
  float error = 0.0; 

  int    * counter  = new int[nWells];
  float  * seisNorm = new float[nWells];
  float * resNorm  = new float[nWells];

  for(int i=0;i<nWells;i++)
  {
    resNorm[i]  = 0.0f;
    seisNorm[i] = 0.0f;
  }

  float minSeisAmp = static_cast<float> (1e-7);
  int totCount=0;
  const double *x, *y;
  int nData;
  for(int i=0;i<nWells;i++)
  {
    x = wells[i]->getXpos(nData);
    y = wells[i]->getYpos(nData);
    int ix, iy;
    simbox->getIndexes(x[0],y[0],ix,iy);
    //scale[i] = gain->GetZ(x[0],y[0]);
    scale[i] = (*gain)(ix,iy);
    counter[i]=0;
    if(wellWeight[i]>0)
    {
      // Count number of layers with seismic data
      for(int k=0;k<nzp;k++)
        if(fabs(seis_r[i][k]) > minSeisAmp)
          counter[i]++;
      totCount+=counter[i];
      for(int k=0;k<nzp;k++)
      {
        if(fabs(seis_r[i][k]) > minSeisAmp)
        {
          seisNorm[i]   += seis_r[i][k] * seis_r[i][k];
          float      foo = float(scale[i]*synt_seis_r[i][k] - seis_r[i][k]);
          resNorm[i] += foo*foo;
        }
        error+=resNorm[i];
      }
    }//if
  }

  float optValue=error;
  err = sqrt(optValue/static_cast<float>(totCount));
  for(int i=0;i<nWells;i++)
  {
    if(counter[i]>0)
      errWell[i] = sqrt(resNorm[i]/counter[i]);
    else
      errWell[i] = 0.0f;
  }

  for(int i=0;i<nWells;i++)
  {
    if(wellWeight[i]>0)
    {
      optValue = resNorm[i];

      scaleOptWell[i]    = float(scale[i]);
      errWellOptScale[i] = sqrt(optValue/float(counter[i]));
    }
    else
    {
      scaleOptWell[i]    = 0.0f;
      errWellOptScale[i] = 0.0f;
    }
  }

  delete [] resNorm;
  delete [] seisNorm;
  delete [] counter;
  delete [] scale;
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
Wavelet::setGainGrid(Grid2D * grid)
{
  double sum = 0.0;
  int nData = 0;
  
  if(gainGrid_ != NULL)
    delete [] gainGrid_;
  gainGrid_ = grid;
  for(int j=0;j<static_cast<int>(gainGrid_->GetNJ());j++)
    for(int i=0;i<static_cast<int>(gainGrid_->GetNI());i++)
    {
      sum+=log((*gainGrid_)(i,j));
      nData++;
    }

    float invGeoMean = float(exp(-sum/static_cast<double>(nData)));
    for(int j=0;j<static_cast<int>(gainGrid_->GetNJ());j++)
      for(int i=0;i<static_cast<int>(gainGrid_->GetNI());i++)
      {
        if((*gainGrid_)(i,j) != WELLMISSING)
          (*gainGrid_)(i,j) *= invGeoMean;
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

void
Wavelet::estimateLocalShift(const CovGrid2D  & cov,
                            Grid2D          *& shift,
                            float            * shiftWell,
                            int              * nActiveData,
                            Simbox           * simbox,
                            WellData        ** wells,
                            int                nWells)
{
  //
  // NBNB-PAL: Since slightly deviated wells are accepted, we should
  // eventually make gain- and shift-cubes rather than single maps.
  //

  //
  // Collect data for kriging
  //
  KrigingData2D shiftData;
  
  for(int i=0;i<nWells;i++)
  {
    if(nActiveData[i]>0) 
    {
      //
      // Coordinates for data point must be chosed from blocked 
      // logs and not from wells
      //
      BlockedLogs * bl = wells[i]->getBlockedLogsOrigThick();
      const double * xPos = bl->getXpos(); 
      const double * yPos = bl->getYpos();
      int xInd, yInd;
      simbox->getIndexes(xPos[0],yPos[0],xInd,yInd);        
      shiftData.addData(xInd,yInd,shiftWell[i]);
    }
  }
  shiftData.findMeanValues();
 
  //
  // Perform kriging
  //
  if(shift==NULL) {
    shift = new Grid2D(simbox->getnx(), 
                       simbox->getny(), 
                       0.0f);
    Kriging2D::krigSurface(*shift, shiftData, cov);
  }
 
}

void
Wavelet::estimateLocalGain(const CovGrid2D  & cov,
                           Grid2D          *& gain,
                           float            * scaleOptWell,
                           float              globalScale, 
                           int              * nActiveData,
                           Simbox           * simbox,
                           WellData        ** wells,
                           int                nWells)
{
  //
  // Collect data for kriging
  //
  KrigingData2D gainData;
  
  for(int i=0;i<nWells;i++)
  {
    if(nActiveData[i]>0) 
    {
      //
      // Coordinates for data point must be chosed from blocked 
      // logs and not from wells
      //
      BlockedLogs * bl = wells[i]->getBlockedLogsOrigThick();
      const double * xPos = bl->getXpos(); 
      const double * yPos = bl->getYpos();
      int xInd, yInd;
      simbox->getIndexes(xPos[0],yPos[0],xInd,yInd);        
      gainData.addData(xInd,yInd,scaleOptWell[i]);
    }
  }
  gainData.findMeanValues();
  
  //
  // Perform kriging
  //
  if(gain==NULL) {
    gain = new Grid2D(simbox->getnx(), 
                       simbox->getny(), 
                       globalScale);
    Kriging2D::krigSurface(*gain, gainData, cov);
  }
}

// Estimate local scaled noise
void
Wavelet::estimateLocalNoise(const CovGrid2D  & cov,
                            Grid2D          *& noiseScaled,
                            float              globalNoise,
                            float            * errWellOptScale,
                            int              * nActiveData,
                            Simbox           * simbox,
                            WellData        ** wells,
                            int                nWells)
{
  //
  // Collect data for kriging
  //
  KrigingData2D noiseData;
  
  for(int i=0;i<nWells;i++)
  {
    if(nActiveData[i]>0) 
    {
      //
      // Coordinates for data point must be chosed from blocked 
      // logs and not from wells
      //
      BlockedLogs * bl = wells[i]->getBlockedLogsOrigThick();
      const double * xPos = bl->getXpos(); 
      const double * yPos = bl->getYpos();
      int xInd, yInd;
      simbox->getIndexes(xPos[0],yPos[0],xInd,yInd);        
      noiseData.addData(xInd,yInd,errWellOptScale[i]/globalNoise);
    }
  }
  noiseData.findMeanValues();
  
  //
  // Perform kriging
  //
  if(noiseScaled==NULL) {
    noiseScaled = new Grid2D(simbox->getnx(), 
                             simbox->getny(), 
                             1.0);
    Kriging2D::krigSurface(*noiseScaled, noiseData, cov);
  }
}

float Wavelet::findGlobalScaleForGivenWavelet(ModelSettings *modelSettings, 
                                                Simbox *simbox,
                                                FFTGrid        * seisCube, 
                                                WellData ** wells)
{
  int i,j,k;
  int nWells        = modelSettings->getNumberOfWells();
  int     nz            = simbox->getnz();                     
  int     nzp           = seisCube->getNzp();
  int     rnzp          = 2*(nzp/2+1);
  fftw_real    ** cpp_r = new fftw_real*[nWells];
  fftw_real    ** seis_r = new fftw_real*[nWells];
  float *seisData = new float[nz];
  bool *hasData = new bool[nz];
  int maxBlocks  = 0;
  for(i=0;i<nWells;i++)
  {
    cpp_r[i]      = new fftw_real[rnzp];
    seis_r[i]     = new fftw_real[rnzp];
    for(j=0;j<rnzp;j++)
    {
      cpp_r[i][j] = 0;
      seis_r[i][j] = 0;
    }
    int nBlocks        = wells[i]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];
  float maxSeis = 0.0;
  float maxCpp = 0.0;
  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) 
  {
    if (wells[w]->getUseForWaveletEstimation())
    {
      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
    //  seisCube->setAccessMode(FFTGrid::RANDOMACCESS);
      bl->getBlockedGrid(seisCube,seisLog);
    //  seisCube->endAccess();
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      
      bl->getVerticalTrend(seisLog, seisData);

      for (k = 0 ; k < nz ; k++) {
        hasData[k] = seisData[k] != RMISSING;
      }
      int start,length;
      bl->findContiniousPartOfData(hasData,nz,start,length);
      
      bl->fillInCpp(coeff_,start,length,cpp_r[w],nzp); 
      bl->fillInSeismic(seisData,start,length,seis_r[w],nzp);

      for(i=0;i<nzp;i++)
      {
        if(cpp_r[w][i]>maxCpp)
          maxCpp = cpp_r[w][i];
        if(seis_r[w][i]>maxSeis)
          maxSeis = seis_r[w][i];
      }

    }
  }
  delete [] seisLog;
  delete [] seisData;
  delete [] hasData;
  for(i=0;i<nWells;i++)
  {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
  }
  delete [] cpp_r;
  delete [] seis_r;

  float kk = maxSeis/maxCpp;


  float maxwavelet = 0.0;
  for(int i=0;i<nzp_;i++)
   if(getRAmp(i)> maxwavelet)
     maxwavelet = getRAmp(i);

  float scale1 = kk/maxwavelet;
 
  return scale1;
}


