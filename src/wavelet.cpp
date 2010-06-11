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
#include "src/io.h"

Wavelet::Wavelet(int dim)
  : dim_(dim),
    shiftGrid_(NULL),
    gainGrid_(NULL)
{
}

Wavelet::Wavelet(int       dim, 
                 Wavelet * wavelet)
  : theta_(wavelet->getTheta()),
    dz_(wavelet->getDz()),
    nz_(wavelet->getNz()),
    nzp_(wavelet->getNzp()),
    cz_(wavelet->getCz()),
    inFFTorder_(wavelet->getInFFTOrder()),
    norm_(wavelet->getNorm()),
    waveletLength_(wavelet->getWaveletLength()),
    dim_(dim),
    shiftGrid_(NULL),
    gainGrid_(NULL)
{
  if(! wavelet->getIsReal()) wavelet->invFFT1DInPlace();
  isReal_ = wavelet->getIsReal(); //NBNB-Frode: Always true?
  coeff_[0] = wavelet->coeff_[0];
  coeff_[1] = wavelet->coeff_[1];
  coeff_[2] = wavelet->coeff_[2];

  cnzp_ = nzp_/2+1;
  rnzp_ = 2*cnzp_;
  rAmp_ = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);

  if(isReal_)
    for(int i = 0; i < rnzp_; i++) {
      rAmp_[i] = wavelet->getRAmp(i);  
    }
  else //NBNB-Frode: This will never happen?
    for(int i = 0; i < cnzp_; i++) {
      cAmp_[i].re = wavelet->getCAmp(i).re; 
      cAmp_[i].im = wavelet->getCAmp(i).im; 
    }
}


Wavelet::Wavelet(const std::string & fileName, 
                 int                 fileFormat, 
                 ModelSettings     * modelSettings, 
                 float             * reflCoef,
                 float               theta,
                 int                 dim,
                 int               & errCode, 
                 std::string       & errText)
  : theta_(theta),
    inFFTorder_(false),
    isReal_(true),
    dim_(dim),
    scale_(1),
    shiftGrid_(NULL),
    gainGrid_(NULL)
{
  coeff_[0]       = reflCoef[0];
  coeff_[1]       = reflCoef[1];
  coeff_[2]       = reflCoef[2];
  switch (fileFormat) {
  case JASON: 
    WaveletReadJason(fileName, errCode, errText);
    break;
  case NORSAR:
    WaveletReadNorsar(fileName, errCode, errText);
    break;
  }
  formats_       = modelSettings->getWaveletFormatFlag();
  waveletLength_ = findWaveletLength(modelSettings->getMinRelWaveletAmp());
  LogKit::LogFormatted(LogKit::LOW,"\n  Estimated wavelet length:  %.1fms.\n",waveletLength_);

  if(errCode == 0) {
    for(int i=0; i < rnzp_ ;i++) {  
      if(i < nzp_)
        rAmp_[i]*=scale_;
      else
        rAmp_[i]=RMISSING;
    }//end for i
  }
}

Wavelet::Wavelet(Wavelet * wavelet, 
                 int       difftype)
  : theta_(wavelet->getTheta()),
    dz_(wavelet->getDz()),
    nz_(wavelet->getNz()),
    nzp_(wavelet->getNzp()),
    cz_(wavelet->getCz()),
    inFFTorder_(wavelet->getInFFTOrder()),
    norm_(wavelet->getNorm()),
    dim_(wavelet->getDim()),
    shiftGrid_(NULL),
    gainGrid_(NULL)
{
  if(! wavelet->getIsReal() ) wavelet->invFFT1DInPlace();
  isReal_    = wavelet->isReal_;  //NBNB-Frode: Always true?
  norm_      = wavelet->norm_;

  cnzp_ = nzp_/2+1;
  rnzp_ = 2*cnzp_;
  rAmp_ = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);
  int i;

  if(difftype != FOURIER) {
//    double norm2 = 0.0;
    for( i = 0; i < rnzp_; i++) {
      if(i < nzp_) {
        switch(difftype) {
        case FIRSTORDERFORWARDDIFF:
          if(i == nzp_-1 )
            rAmp_[i] = wavelet->getRAmp(0) - wavelet->getRAmp(i);
          else
            rAmp_[i] = wavelet->getRAmp(i+1)-wavelet->getRAmp(i);      
          break;
        case FIRSTORDERBACKWARDDIFF:
          if(i == 0 )
            rAmp_[i] = wavelet->getRAmp(0) - wavelet->getRAmp(nzp_-1);
          else
            rAmp_[i] = wavelet->getRAmp(i)-wavelet->getRAmp(i-1);      
          break;
        case FIRSTORDERCENTRALDIFF:
          if(i == 0 )
            rAmp_[i] = static_cast<fftw_real> (0.5 * (wavelet->getRAmp(i+1) - wavelet->getRAmp(nzp_-1)));
          else {
            if(i == nzp_-1 )
              rAmp_[i] = static_cast<fftw_real> (0.5 * (wavelet->getRAmp(0) - wavelet->getRAmp(i-1)));
            else
              rAmp_[i] = static_cast<fftw_real> (0.5 * (wavelet->getRAmp(i+1) - wavelet->getRAmp(i-1)));
          }      
          break;
        }
 //       norm2 += static_cast<double> (rAmp_[i]*rAmp_[i]);
      }
      else
        rAmp_[i] = RMISSING;
    }
//    norm_= findNorm(); //static_cast<float> (sqrt(norm2));
  }
  else {
    fftw_complex  iValue;
    for(i=0; i < nzp_; i++ )
      rAmp_[i] = wavelet->getRAmp(i);
    fft1DInPlace();
    for(i=0; i < cnzp_; i++ ) {
      iValue  =  cAmp_[i];
      cAmp_[i].re = static_cast<float> (- iValue.im * 2.0 * M_PI * i/(static_cast<float>(nzp_)));
      cAmp_[i].im = float(iValue.re * 2 * M_PI * i/(static_cast<float>(nzp_)));
    }
    invFFT1DInPlace();
//    norm_ = findNorm();
  }
  norm_ = findNorm();
}

Wavelet::Wavelet(int difftype, 
                 int nz, 
                 int nzp)
  : theta_(RMISSING),
    dz_(RMISSING),
    nz_(nz),
    nzp_(nzp),
    cz_(0),
    inFFTorder_(true),
    dim_(1),
    shiftGrid_(NULL),
    gainGrid_(NULL)
{
  cnzp_       = nzp_/2+1;
  rnzp_       = 2*cnzp_;
  rAmp_       = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  
  cAmp_       = reinterpret_cast<fftw_complex*>(rAmp_);
  norm_       = RMISSING;
  int i;

  if(difftype != FOURIER) {
    isReal_    = true;
    for(i=0;i < rnzp_; i++) {
      if( i < nzp_)
        rAmp_[i] = 0.0;
      else
        rAmp_[i] = RMISSING;
    }

    switch(difftype) {
    case FIRSTORDERFORWARDDIFF:
      rAmp_[0] = -1.0; 
      rAmp_[nzp_-1] = 1.0; 
      norm_    = sqrt(2.0f);
      break;
    case FIRSTORDERBACKWARDDIFF:
      rAmp_[0]      = 1.0;
      rAmp_[nzp_-1] =  -1.0;      
      norm_    = sqrt(2.0f); 
      break;
    case FIRSTORDERCENTRALDIFF:
      rAmp_[1]      = 0.5;
      rAmp_[nzp_-1] = -0.5;
      norm_    = sqrt(0.5f); 
      break;
    }
  }
  else {
    isReal_    = false;
    for(i=0;i < cnzp_; i++) {
      cAmp_[i].re = 0.0f;
      cAmp_[i].im = static_cast<float>( 2.0 * M_PI * i/(static_cast<float>(nzp_)) );   
    }

    invFFT1DInPlace();
    norm_ = findNorm();
  }       
}


Wavelet::~Wavelet()
{
  if(shiftGrid_!=NULL)
    delete shiftGrid_;
  if(gainGrid_!=NULL)
    delete gainGrid_;
  fftw_free(rAmp_);
}

void Wavelet::fft1DInPlace()
{
  // use the operator version of the fourier transform
  if(isReal_) {
    int flag; 
    rfftwnd_plan plan;  
    flag    = FFTW_ESTIMATE | FFTW_IN_PLACE;
    plan    = rfftwnd_create_plan(1, &nzp_ ,FFTW_REAL_TO_COMPLEX,flag);
    //
    // NBNB-PAL: The call rfftwnd_on_real_to_complex is causing UMRs in Purify.
    //
    rfftwnd_one_real_to_complex(plan,rAmp_,cAmp_);
    fftwnd_destroy_plan(plan);
    isReal_ = false;
  }
}

void Wavelet::invFFT1DInPlace()
{
  // use the operator version of the fourier transform
  if(!isReal_) {
    int flag;
    rfftwnd_plan plan; 

    flag = FFTW_ESTIMATE | FFTW_IN_PLACE;
    plan= rfftwnd_create_plan(1,&nzp_,FFTW_COMPLEX_TO_REAL,flag);
    rfftwnd_one_complex_to_real(plan,cAmp_,rAmp_);
    fftwnd_destroy_plan(plan);
    isReal_=true;
    double scale= static_cast<double>(1.0/static_cast<double>(nzp_));
    for(int i=0; i < nzp_; i++)
      rAmp_[i] = static_cast<fftw_real>(rAmp_[i]*scale);  
  }
}

fftw_real  
Wavelet::getRAmp(int k)
{
  fftw_real value;

  if(isReal_) {
    if(k < nzp_)
      value = rAmp_[k];
    else
      value = 0.0;
  }
  else {
    invFFT1DInPlace();
    if(k < nzp_)
      value = rAmp_[k];
    else
      value = 0.0;

    fft1DInPlace();
  }
  return value;
}


fftw_complex   
Wavelet::getCAmp(int k) const
{
  assert(!isReal_);
  fftw_complex  value;

  if(k < cnzp_) {
    value.re =  cAmp_[k].re;
    value.im =  cAmp_[k].im;
  }
  else {
    int refk =  nzp_-k;
    value.re =  cAmp_[refk].re;
    value.im =  - cAmp_[refk].im;
  }
  return value;
}


fftw_complex   
Wavelet::getCAmp(int    k, 
                 float  scale) const
{
  ///////////////////////////////////////////////////////////
  //
  // Get the  fourier transform of the streched wavelet.
  // scale is  in [0 1] wavelet
  // scale = 1 this is identical to getCAmp(int k)
  // Note we do not use the normal scale relation:  
  //  
  // FFT( w(s*t) ) = fw(w/s) * ( 1/s )   with  fw( w ) = FFT( w(t) )
  //
  // We return:  fw(w/s) 
  ///////////////////////////////////////////////////////////// 

  assert(!isReal_);

  fftw_complex  value;
  float omega,dOmega;  
  int   omU, omL;
  if( k < cnzp_) {
    omega  = float(k) / scale;
    if(omega >= cnzp_) {
      value.re =  0.0f;
      value.im =  0.0f;
    }
    else {
      omL        = int(floor( omega ));
      omU        = int( floor( omega )) + 1;
      if(omU >= cnzp_) 
        omU -= 1;    

      dOmega     = omega - float(omL);
      value.re =  cAmp_[omL].re * ( 1.0f - dOmega ) + cAmp_[omU].re * dOmega;
      value.im =  (cAmp_[omL].im * ( 1.0f - dOmega ) + cAmp_[omU].im * dOmega);
    }
  }
  else {
    int refk =  nzp_-k;
    omega = float(refk)/scale;

    if(omega >= cnzp_) {
      value.re =  0.0f;
      value.im =  0.0f;
    }
    else {
      omL        = int(floor( omega ));
      omU        = int( floor( omega )) + 1;
      if(omU >= cnzp_) 
        omU -= 1;
      dOmega     = omega - float(omL);
      value.re =  cAmp_[omL].re * ( 1.0f - dOmega ) + cAmp_[omU].re * dOmega;
      value.im =  (- cAmp_[omL].im * ( 1.0f - dOmega ) - cAmp_[omU].im * dOmega);
    }
  }
  return value;
}

void           
Wavelet::setRAmp(float  value, 
                 int    k)
{
  rAmp_[k] = value;
}

Wavelet*  
Wavelet::getLocalWavelet(int i, 
                         int j)
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

void
Wavelet::resample(float dz, 
                  int   nz, 
                  int   nzp) 
{
  //LogKit::LogFormatted(LogKit::LOW,"  Resampling wavelet\n");
  assert(isReal_);
  assert(!inFFTorder_);
  
  int cnzp  =  nzp/2 + 1;
  int rnzp  =  2*cnzp;

  fftw_real * wlet  = static_cast<fftw_real *>(fftw_malloc( sizeof(fftw_real)*rnzp ));

  float z;
  for(int k=0; k < rnzp; k++) {
    if(k < nzp) {
      if(k < nzp/2+1)
        z = static_cast<float>( dz*k );
      else 
        z = static_cast<float>( dz*(k-nzp) );
      wlet[k] = getWaveletValue(z, rAmp_ , cz_, nz_, dz_);
    }
    else
      wlet[k] =RMISSING;
  }
  fftw_free( rAmp_);

  double norm2 = 0.0; 
  for(int k=0; k < nzp; k++) 
    norm2 += static_cast<double> (wlet[k]*wlet[k]);
  norm_       = static_cast<float>(sqrt( norm2));

  rAmp_       = static_cast<fftw_real *>(wlet); // rAmp_ is not allocated 
  cAmp_       = reinterpret_cast<fftw_complex*>(rAmp_);
  nzp_        = nzp;
  rnzp_       = rnzp;
  cnzp_       = cnzp;
  cz_         = 0;
  nz_         = nz;
  dz_         = dz;
  inFFTorder_ = true;

  if( ModelSettings::getDebugLevel() > 0 ) {
    std::string fileName = "resampled_wavelet_";
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
  }
}

void 
Wavelet::multiplyRAmpByConstant(float c)
{
  for(int i=0; i < rnzp_ ;i++) {  
    if(i < nzp_)
      rAmp_[i]*=c;
    else
      rAmp_[i]=RMISSING;
  }//end for i
}

void
Wavelet::scale(float scale)
{
  if (scale != 1.0f) {
    LogKit::LogFormatted(LogKit::LOW,"  Scaling wavelet with factor         : %.3e\n",scale);
    for(int i=0; i < rnzp_ ; i++)
      if(rAmp_[i] != RMISSING)
        rAmp_[i]=rAmp_[i]*scale;
  }
  scale_ = scale;
}

void
Wavelet::printToFile(const std::string & fileName, 
                     bool                overrideDebug)
{
  if(overrideDebug == true || ModelSettings::getDebugLevel() > 0) {
    std::string fName = IO::makeFullFileName(IO::PathToWavelets(), fileName);
    std::ofstream file;
    NRLib::OpenWrite(file, fName);
    for(int i=0;i<nzp_;i++)
      file << getRAmp(i) << "\n";
    file.close();
  }
}

void
Wavelet::writeWaveletToFile(const std::string & fileName, 
                            float               approxDzIn)
{
  float approxDz = std::min(approxDzIn,static_cast<float>(floor(dz_*10)/10));
  approxDz = std::min(approxDzIn,dz_);
  
  //Trick: Written wavelet may be shorter than the actual.
  //This gives inconsistency if a wavelet is read and written.
  //Make consistent by truncating wavelet to writing range before interpolation.
  int     activeCells = int(floor(waveletLength_/2/dz_));
  float * remember    = new float[nzp_];
  for(int i=activeCells+1;i<=nz_;i++) {
    remember[i]        = rAmp_[i];
    rAmp_[i]           = 0;
    remember[nzp_-i] = rAmp_[nzp_-i];
    rAmp_[nzp_-i]    = 0;
  }
  float          T            = nzp_*dz_;
  int            nzpNew       = int(ceil(T/approxDz - 0.5));  
  float          dznew        = T/float(nzpNew);
  int            cnzpNew      = (nzpNew/2)+1;
  
  fftw_real    * waveletNew_r =  new fftw_real[2*cnzpNew];
  fftw_complex * waveletNew_c =  reinterpret_cast<fftw_complex*>(waveletNew_r);

  fft1DInPlace();

  double         multiplyer = static_cast<double>(nzpNew)/static_cast<double>(nzp_);
  
  for(int i=0;i<cnzpNew;i++) {
    if(i < cnzp_) {
      waveletNew_c[i].re = static_cast<fftw_real>(cAmp_[i].re*multiplyer);
      waveletNew_c[i].im = static_cast<fftw_real>(cAmp_[i].im*multiplyer);
      if((i==(cnzp_-1)) & (2*((cnzp_-1)/2) != cnzp_-1)) //boundary effect in fft domain
        waveletNew_c[i].re*=0.5;
    }
    else { 
      waveletNew_c[i].re = 0;
      waveletNew_c[i].im = 0;
    }
  }
  invFFT1DInPlace();
  for(int i=activeCells+1;i<=nz_;i++) {
    rAmp_[i] = remember[i];
    rAmp_[nzp_-i] = remember[nzp_-i];
  }
  delete [] remember;
  
  
  Utils::fftInv(waveletNew_c,waveletNew_r,nzpNew );// note might be n^2 algorithm for some nzpNew
  
  int wLength = int(floor(waveletLength_/dznew+0.5));
  int halfLength = wLength/2; // integer division
  wLength =  halfLength*2+1;// allways odd
  if( wLength>nzpNew) {  
    wLength=2*(nzpNew/2)-1;// allways odd
    halfLength=wLength/2;
  }
  
  float shift = -dznew*halfLength;
  
  std::string fName;
  std::ofstream file;
  
  if ((formats_ & IO::JASONWAVELET)>0)
    {
    fName = std::string(fileName) + NRLib::ToString(theta_*(180/M_PI), 1) + "deg" + IO::SuffixJasonWavelet();
    fName = IO::makeFullFileName(IO::PathToWavelets(), fName);
    
    LogKit::LogFormatted(LogKit::MEDIUM,"  Writing Wavelet to file \'"+fName+"\'\n");
    
    NRLib::OpenWrite(file, fName);
    
    file << "\"* Export format using Comma Separated Values\"\n"
         << "\"* Wavelet written from CRAVA\"\n"
         << "\"* Generated \"\n"
         << "\"*\"\n"
         << "\"* File format: \"\n"
         << "\"* - N lines starting with * are comment (such as this line)\"\n"
         << "\"* - 1 line with four fields (data type, data unit, depth type, depth unit) \"\n"
         << "\"* - 1 line with start time  \"\n"
         << "\"* - 1 line with sample interval \"\n"
         << "\"* - 1 line with number of data lines \"\n"
         << "\"* - N lines with trace data \"\n"
         << "\"* Data values are represented as floating point numbers,\"\n"
         << "\"* \"\n"
         << "\"wavelet\",\"none\",\"time\",\"ms\"\n"
         << std::fixed 
         << std::setprecision(0)
         << shift   << "\n"
         << std::setprecision(2)
         << dznew   << "\n"
         << wLength << "\n";
    
    for(int i=halfLength ; i > 0 ; i--)
      file << std::setprecision(6) << waveletNew_r[nzpNew-i] << "\n";
    for(int i=0;i<=halfLength;i++)
      file << std::setprecision(6) << waveletNew_r[i] << "\n";
    file.close();
  }

  if ((formats_ & IO::NORSARWAVELET)>0)  
    {
    //Writing wavelet in swav-format
    fName = std::string(fileName) + NRLib::ToString(theta_*(180/M_PI), 1) + "deg" + IO::SuffixNorsarWavelet();
    fName = IO::makeFullFileName(IO::PathToWavelets(), fName);
  
    int   cznew = static_cast<int>(floor((fabs(shift/dznew))+0.5));
    float t0    = static_cast<float>(-dznew * cznew);

    LogKit::LogFormatted(LogKit::MEDIUM,"  Writing Wavelet to file \'"+fName+"\'\n");

    NRLib::OpenWrite(file, fName);

    file << "pulse file-3\n"
         << "* Written by CRAVA\n"
         << "* time - millisec - sample value legend\n"
         << "Time 1 Amplitude\n" 
         << "* Min, main and max frequency/wavenumber\n"
         << "-999.99 -999.99 -999.99\n"
         << "* Number of samples - millisec sampling - time at first sample\n"
         << std::fixed
         << std::setprecision(2)
         << wLength
         << " "
         << dznew
         << " "
         << std::setprecision(0)
         << t0  
         << "\n";

    for (int i = halfLength ; i > 0 ; i--)
      file << std::setprecision(6) << waveletNew_r[nzpNew-i] << "\n";
    for (int i = 0; i <= halfLength; i++)
      file << std::setprecision(6) << waveletNew_r[i] << "\n";
    file.close();
  }

  delete [] waveletNew_r;
}

void 
Wavelet::setShiftGrid(Grid2D *grid)
{
  if(shiftGrid_ != NULL)
    delete [] shiftGrid_;
  shiftGrid_ = grid; 
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
    for(int i=0;i<static_cast<int>(gainGrid_->GetNI());i++) {
      sum+=log((*gainGrid_)(i,j));
      nData++;
    }

    float invGeoMean = float(exp(-sum/static_cast<double>(nData)));
    for(int j=0;j<static_cast<int>(gainGrid_->GetNJ());j++)
      for(int i=0;i<static_cast<int>(gainGrid_->GetNI());i++) {
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
    }
    if(fftflag == true)
      fft1DInPlace();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////PROTECTED MEMBER METHODS////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
void
Wavelet::shiftAndScale(float shift,
                       float gain)
{
  int k;

  fftw_complex  ampMultiplier,tmp;

  if(isReal_) 
    fft1DInPlace();

  float iShift=shift/dz_;

  for(k=0;k < cnzp_; k++) {
    ampMultiplier.re = float(gain*cos(2.0*(M_PI*(iShift)*k)/float(nzp_)));
    ampMultiplier.im = float(gain*sin(-2.0*(M_PI*(iShift)*k)/float(nzp_)));

    tmp.re = ampMultiplier.re*cAmp_[k].re - ampMultiplier.im*cAmp_[k].im;
    tmp.im = ampMultiplier.im*cAmp_[k].re + ampMultiplier.re*cAmp_[k].im;

    cAmp_[k] =tmp;
  }
}

float
Wavelet::findWaveletLength(float minRelativeAmp)
{
  bool trans=false;
  if(isReal_==false) {
    invFFT1DInPlace();
    trans=true;
  }
  
  float maxAmp =  fabs(getRAmp(0)); // gets max amp 
  for(int i=1;i <nzp_;i++)
    if(fabs(getRAmp(i)) > maxAmp)
      maxAmp = fabs(getRAmp(i));

  float minAmp= maxAmp*minRelativeAmp; // minimum relevant amplitude

  int wLength=nzp_;

  for(int i=nzp_/2;i>0;i--) {
    if(fabs(getRAmp(i)) >minAmp) {
      wLength= (i*2+1);// adds both sides 
      break;
    }
    if(fabs(getRAmp(nzp_-i)) > minAmp) {
      wLength= (2*i+1);// adds both sides 
      break;
    }
  }
  wLength = std::min(wLength,2*((nzp_+1)/2) - 1); // always odd number
  if(trans==true)
    fft1DInPlace();

  return (dz_*static_cast<float>(wLength));
}

float
Wavelet::findNorm() const
{
  double norm2=0.0;
  for(int i=0; i < nzp_; i++ )
    norm2 += static_cast<double> (rAmp_[i]*rAmp_[i]);
  float norm = static_cast<float>(sqrt(norm2));
  return norm;
}

fftw_real* 
Wavelet::averageWavelets(const std::vector<std::vector<fftw_real> > & wavelet_r,
                         int                                          nWells,
                         int                                          nzp,
                         const std::vector<float>                   & wellWeight,
                         const std::vector<float>                   & dz,
                         float                                        dzOut) const
{
  // assumes dz[w] < dzOut for all w
  fftw_real* wave= static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  
  for(int i=0; i<nzp; i++)
    wave[i] = 0.0; // initialize

  std::vector<float> weight(nWells);// weight is length of data interval in well
  float sum = 0.0f;
  for(int w=0; w<nWells; w++)
    sum += wellWeight[w];

  for(int w=0; w<nWells; w++)
    weight[w] = wellWeight[w]/sum;
  
  float ww;

  for(int w=0;w<nWells;w++) {
    if(wellWeight[w] > 0.0) {
      wave[0] += weight[w]* wavelet_r[w][0]; // central
      for(int i=1;i <nzp/2-1;i++) { // positive time
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];  // fraction of distance to ind*dz[w]
        if(ind<nzp/2)
            ww= (1-t)* wavelet_r[w][ind]+t*(wavelet_r[w][ind+1]);
        else
          ww=0;
        wave[i]    += weight[w]*ww;
      }
      for(int i=1;i < nzp/2-1;i++) {// negative time Note dz[w] <= dzOut for all w
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];
        if(ind<nzp/2)
          ww    = (1-t)* wavelet_r[w][nzp-ind]+t*(wavelet_r[w][nzp-(ind+1)]);
        else
          ww=0;
        wave[nzp-i] += weight[w]*ww;
      }
    }
  }
  
  std::string fileName;
  fileName = "wavelet_"+NRLib::ToString(int(floor(theta_/M_PI*180+0.5)))+"_fftOrder_noshift";
  printVecToFile(fileName,wave,nzp_);

  return wave;
}

void
Wavelet::fillInnWavelet(fftw_real* wavelet_r,
                        int nzp,
                        float dz)
{
  //Note: Assumes that dz_ > dz
  // NBNB OddK 
  fftw_real previous1 = getRAmp(0);
  fftw_real previous2 = getRAmp(0);
  fftw_real current1 = 0.0f;
  fftw_real current2 = 0.0f;
  wavelet_r[0]        = previous1;

  int counterForWavelet   = 1;
  float w;
  for(int i=1;i<= nzp/2;i++) {
    if(counterForWavelet*dz_ < i*dz) {
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

float
Wavelet::findBulkShift(fftw_real* vec_r,
                       float dz,
                       int nzp,
                       float maxShift)
{
  float shift=0.0f;
  float sum=0;
  int i,polarity;
  // if the sum from -maxShift to maxShift ms is 
  // positive then polarity is positive   
  for(i=0;i<ceil(maxShift/dz);i++)//zero included
    sum+=vec_r[i];
  for(i=0;i<floor(maxShift/dz);i++)
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

  for(i=0;i<ceil(maxShift/dz);i++) {
    if(vec_r[i]*polarity > maxValue) {
      maxValue = vec_r[i]*polarity;
      shiftI = i;
    }
  }
  for(i=0;i<floor(maxShift/dz);i++) {
    if(vec_r[nzp-1-i]*polarity > maxValue) {
      maxValue = vec_r[nzp-1-i]*polarity;
      shiftI = -1-i;
    }
  }
  if(shiftI < 0) {
    if(vec_r[nzp+shiftI-1]*polarity < maxValue) { //then local max
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
  else { //positive or zero shift
    if(vec_r[shiftI+1]*polarity < maxValue) { //then local max
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
Wavelet::shiftReal(float        shift, 
                   fftw_real  * rAmp,
                   int          nt)
{
  fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);
  Utils::fft(rAmp,cAmp, nt);
  int cnzp= nt/2+1;
  float expo;
  fftw_complex tmp,mult;
  for(int i=0;i<cnzp;i++) {
    tmp     = cAmp[i];
    expo    = static_cast<float>(-2.0*shift*M_PI*static_cast<float>(i)/static_cast<float>(nt));
    mult.re = cos(expo);
    mult.im = sin(expo);
    cAmp[i].re = tmp.re*mult.re-tmp.im*mult.im;
    cAmp[i].im = tmp.re*mult.im+tmp.im*mult.re;
  }
  Utils::fftInv(cAmp,rAmp, nt);
}

void
Wavelet::printVecToFile(const std::string & fileName, 
                        fftw_real         * vec, 
                        int                 nzp) const
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

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////PRIVATE MEMBER METHODS//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

float 
Wavelet::getWaveletValue(float   z, 
                         float * Wavelet, 
                         int     center, 
                         int     nz, 
                         float   dz)
{
  // returns the value of the vavelet in the location z. Wavelet have the length nz 
  // and the center value is Wavelet[center]
  // uses kriging with ricker 20Hz wavelet as correlation function.
  float value;
  int    k,l;
  int*   ind=new int[6];// iL1,iL2,iL3,iR1,iR2,iR3;
  double* val=new double[6];//vL1,vL2,vL3,vR1,vR2,vR3;

  ind[2]= int( floor( (z/dz) ) );
  for(k=0;k<6;k++)
    ind[k]=  ind[2]+k-2;

  for(k=0;k<6;k++)
    val[k]=  getArrayValueOrZero(ind[k]+center , Wavelet,  nz); 

  double** Cov = new double*[6];
  double*  cov = new double[6];
  double   nu = 20;
  double   deltaT;

  for(k=0;k<6;k++)
  { 
    Cov[k] = new double[6];
    deltaT = (dz*ind[k]-z)*0.001;
    cov[k] = (1-2*nu*nu*M_PI*M_PI*(deltaT)*(deltaT))*exp(-nu*nu*M_PI*M_PI*(deltaT)*(deltaT));
    for(l=0;l<6;l++)
    {
      deltaT =(dz*ind[k]-dz*ind[l])*0.001;
      Cov[k][l] = (1-2*nu*nu*M_PI*M_PI*deltaT * deltaT )*exp(-nu*nu*M_PI*M_PI*deltaT*deltaT);
    }
  }
  //OK not very intellegent implementation since chol is done for each time step.
  lib_matrCholR(6,  Cov);
  lib_matrAxeqbR(6, Cov, cov); // cov contains now the kriging weigths; 

  value = 0.0;
  for(k=0;k<6;k++)
  {
    value+= float(val[k]*cov[k]);
    delete [] Cov[k];
  }
  delete [] Cov;
  delete [] cov;
  delete [] val;
  delete [] ind;
  return value;
}



float 
Wavelet::getArrayValueOrZero(int     i,
                             float * Wavelet, 
                             int     nz) const
{
  float value;

  if(i > -1 && i < nz)
    value = Wavelet[i];
  else
    value = 0.0; 
  return value;
}


float  
Wavelet::getLocalTimeshift(int i, 
                           int j) const
{
  float shift = 0.0f;

  if(shiftGrid_ != NULL && i < static_cast<int>(shiftGrid_->GetNI()) && i>=0 && 
     j < static_cast<int>(shiftGrid_->GetNJ()) &&  j>=0) {
    if((*shiftGrid_)(i,j) != WELLMISSING)
      shift = float((*shiftGrid_)(i,j));
  }

  return shift;
}

float  
Wavelet::getLocalGainFactor(int i, 
                            int j) const
{
  float gain = 1.0f;
  
  if(gainGrid_ != NULL && i < static_cast<int>(gainGrid_->GetNI()) &&  i>=0 && 
    j < static_cast<int>(gainGrid_->GetNJ()) &&  j>=0) {
    if((*gainGrid_)(i,j) != WELLMISSING)
      gain = float((*gainGrid_)(i,j));
  }

  return(gain);
}

void
Wavelet::WaveletReadJason(const std::string & fileName, 
                          int               & errCode, 
                          std::string       & errText)
{
  std::ifstream file;
  NRLib::OpenRead(file,fileName);
  std::string dummyStr;
  
  bool lineIsComment = true; 
  int  line          = 0;
  int  thisLine      = 0;

  while( lineIsComment == true) {
    if(NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    }

    NRLib::ReadNextToken(file,dummyStr,line);
    if (line == thisLine)
      NRLib::DiscardRestOfLine(file,line,false);
    thisLine = line;        
    if((dummyStr[0]!='*') &  (dummyStr[0]!='"')) {
      lineIsComment = false;
    }
  }

  float shift = NRLib::ParseType<float>(dummyStr);

  if (NRLib::CheckEndOfFile(file))  {
    errText += "Error: End of file "+fileName+" premature.\n";
    errCode=1; 
    return;
  } 
  NRLib::ReadNextToken(file,dummyStr,line);
  if (line == thisLine)
    NRLib::DiscardRestOfLine(file,line,false);
  thisLine = line;

  dz_ = NRLib::ParseType<float>(dummyStr);

  if (NRLib::CheckEndOfFile(file)) {
    errText += "Error: End of file "+fileName+" premature.\n";
    errCode=1; 
    return;
  }
  NRLib::ReadNextToken(file,dummyStr,line); 
  if (line == thisLine)
    NRLib::DiscardRestOfLine(file,line,false);
  thisLine = line;
 
  nz_ = NRLib::ParseType<int>(dummyStr);
  cz_   =  static_cast<int>(floor((fabs(shift/dz_))+0.5));
  nzp_  = nz_;
  cnzp_ = nzp_/2+1;
  rnzp_ = 2*cnzp_; 
  rAmp_ = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp_));
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);
  norm_ = RMISSING;

  for(int i=0; i<nz_;i++) {
    if (NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    } 
    NRLib::ReadNextToken(file,dummyStr,line);

    rAmp_[i] = static_cast<fftw_real>(NRLib::ParseType<float>(dummyStr));
  }
  file.close();
}

void
Wavelet::WaveletReadNorsar(const std::string & fileName, 
                           int               & errCode, 
                           std::string       & errText)
{
  std::ifstream file;
  NRLib::OpenRead(file,fileName);
  std::string dummyStr;
  
  bool  lineIsComment = true; 
  int   line          = 0;
  int   timeUnits;
  float milliSec;
  float samplingInt;
  float firstSample;

  NRLib::DiscardRestOfLine(file, line, false);

  // Domain
  while (lineIsComment == true) {
    if(NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    }

    NRLib::ReadNextToken(file,dummyStr,line);
    if (dummyStr[0] != '*') 
      lineIsComment = false;
    else
      NRLib::DiscardRestOfLine(file,line,false);
  }

  if (NRLib::Uppercase(dummyStr) != "TIME") {
    errText += "Error: The NORSAR wavelet in file "+fileName+" needs to be given in Time domain.\n";
    errCode  = 1;
    return;
  }

  timeUnits = NRLib::ReadNext<int>(file,line);
  milliSec = static_cast<float>(std::pow(10.0, 3*(1-timeUnits)));

  NRLib::DiscardRestOfLine(file,line,false); //Contains text for sample values not used by CRAVA

  // Frequency info
  lineIsComment = true;
  while( lineIsComment == true) {
    if(NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    }

    NRLib::ReadNextToken(file,dummyStr,line);
    if (dummyStr[0] != '*') 
      lineIsComment = false;
    else
      NRLib::DiscardRestOfLine(file,line,false);
  }
 
  NRLib::DiscardRestOfLine(file,line,false); //Line contains frequency information not used by CRAVA

  // Sample info
  lineIsComment = true; 
  while( lineIsComment == true) {
    if(NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    }

    NRLib::ReadNextToken(file,dummyStr,line);
    if (dummyStr[0] != '*') 
      lineIsComment = false;
    else
      NRLib::DiscardRestOfLine(file,line,false);
  }

  nz_         = NRLib::ParseType<int>(dummyStr);
  samplingInt = NRLib::ReadNext<float>(file,line);
  firstSample = NRLib::ReadNext<float>(file,line);
  
  dz_   = samplingInt * milliSec;
  cz_   = static_cast<int>((-firstSample / dz_) + 0.5);
  nzp_  = nz_;
  cnzp_ = nzp_/2+1;
  rnzp_ = 2*cnzp_; 
  rAmp_ = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnzp_));
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);
  norm_ = RMISSING;

  // Pulse samples
  for(int i=0; i<nz_;i++) {
    if (NRLib::CheckEndOfFile(file)) {
      errText += "Error: End of file "+fileName+" premature.\n";
      errCode=1; 
      return;
    } 
    NRLib::ReadNextToken(file,dummyStr,line);

    rAmp_[i] = static_cast<fftw_real>(NRLib::ParseType<float>(dummyStr));
  }
  file.close();
}
