#include "wavelet3D.h"

#include <string.h>
#include <assert.h>
#include <math.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_misc.h"
#include "lib/lib_matr.h"
#include "lib/sgri.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/modelsettings.h"
#include "src/model.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/definitions.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet3D::Wavelet3D(const std::string & fileName, 
                     ModelSettings     * modelSettings, 
                     Simbox            * simBox, 
                     float               theta, 
                     float             * reflCoef,  
                     int               & errCode, 
                     char              * errText)
: Wavelet(modelSettings, 3, reflCoef)
{
  nx_ = simBox->getnx();
  ny_ = simBox->getny();
  nz_ = simBox->getnz();
  dx_ = static_cast<float>(simBox->getdx());
  dy_ = static_cast<float>(simBox->getdy());
  dz_ = static_cast<float>(simBox->getdz());
  theta_ = theta;
  norm_ = RMISSING;

  double xPadFac = modelSettings->getXPadFac();
  nxp_   =  findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(nx_)*(1.0+xPadFac) )) );
  double yPadFac = modelSettings->getYPadFac();
  nyp_   =  findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(ny_)*(1.0+yPadFac) )) );
  double zPadFac = modelSettings->getZPadFac();
  nzp_   =  findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(nz_)*(1.0+zPadFac) )) );
  
  ampCube_ = FFTGrid(nx_, ny_, nz_, nxp_, nyp_, nzp_);
  ampCube_.createRealGrid();
  ampCube_.setType(FFTGrid::COVARIANCE);
  ampCube_.setAccessMode(FFTGrid::RANDOMACCESS);

  inFFTorder_ = true;
  readtype_ = Wavelet::SGRI;
  Sgri *sgri = new Sgri(fileName.c_str(), errText, errCode);
  
  if (errCode == 0) {
    float xLim = 0.5f * dx_ * nxp_;
    float yLim = 0.5f * dy_ * nyp_;
    float zLim = 0.5f * dz_ * nzp_;
    int axisOk = sgri->sizeOk(xLim, yLim, zLim);
    if (axisOk == 1) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in x-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 2) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in y-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 3) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in x- and y-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 4) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in z-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 5) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in x- and z-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 6) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in y- and z-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    else if (axisOk == 7) {
      sprintf(errText,"%s3-D wavelet read from file %s has too big size in x-, y- and z-direction compared to padded fft-grid.\n",errText,fileName.c_str()); 
      errCode=1;
    }
    if (axisOk != 0) {
      sprintf(errText,"%s\nToo big wavelet grid read from file  %s.\n",errText, fileName.c_str());
      errCode=1;
    }
  }

  if (errCode == 0) {
    int i,j,k;
    float x, y, z, value;
    for (k=0; k<=nzp_/2; k++) {
      z = k * dz_;
      for (j=0; j<=nyp_/2; j++) {
        y = j * dy_;
        for (i=0; i<=nxp_/2; i++) {
          x = i * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
        for (i=(nxp_/2)+1; i<nxp_; i++) {
          x = (i-nxp_) * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
      }
      for (j=(nyp_/2)+1; j<nyp_; j++) {
        y = (j-nyp_) * dy_;
        for (i=0; i<=nxp_/2; i++) {
          x = i * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
        for (i=(nxp_/2)+1; i<nxp_; i++) {
          x = (i-nxp_) * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        } 
      }
    }
    for (k=(nzp_/2)+1; k<nzp_; k++) {
      z = (k-nzp_) * dz_;
      for (j=0; j<=nyp_/2; j++) {
        y = j * dy_;
        for (i=0; i<=nxp_/2; i++) {
          x = i * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
        for (i=(nxp_/2)+1; i<nxp_; i++) {
          x = (i-nxp_) * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
      }
      for (j=(nyp_/2)+1; j<nyp_; j++) {
        y = (j-nyp_) * dy_;
        for (i=0; i<=nxp_/2; i++) {
          x = i * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        }
        for (i=(nxp_/2)+1; i<nxp_; i++) {
          x = (i-nxp_) * dx_;
          value = sgri->getWaveletValue(x, y, z);
          ampCube_.setRealValue(i, j, k, value, true);
        } 
      }
    }
    waveletLength_ = getWaveletLengthF();
    
    //For debugging purposes, the shifted FFTGrid is written to file to compare to input sgri
    FFTGrid *shiftAmp = new FFTGrid(nx_, ny_, nz_, nx_, ny_, nz_);
    shiftAmp->fillInConstant(0.0);
    shiftAmp->setType(FFTGrid::DATA);
    shiftAmp->setAccessMode(FFTGrid::RANDOMACCESS);

    shiftFFTGrid(shiftAmp);
    char fName[200];
    int thetaDeg = int ( theta/PI*180 + 0.5 );
    sprintf(fName, "WL_as_shiftedFFTGrid_%d", thetaDeg);
    std::string sgriLabel("3D Wavelet as shifted FFT-grid for incidence angle ");
    sgriLabel += NRLib::ToString(thetaDeg);
    shiftAmp->writeFile(fName, simBox, sgriLabel);
    delete shiftAmp;
    //End for debugging purposes
  }

  delete sgri;
}

Wavelet3D::Wavelet3D(Wavelet * wavelet, int difftype)
  : Wavelet(wavelet, 3)
{
  assert(wavelet->getDim() == 3);
  assert(wavelet->getIsReal());
  nx_ = wavelet->getNx();
  ny_ = wavelet->getNy();
  nxp_ = wavelet->getNxp();
  nyp_ = wavelet->getNyp();
  dx_ = wavelet->getDx();
  dy_ = wavelet->getDy();
  ampCube_ = FFTGrid(nx_, ny_, nz_, nxp_, nyp_, nzp_);
  ampCube_.createRealGrid();
  ampCube_.setType(FFTGrid::COVARIANCE);
  ampCube_.setAccessMode(FFTGrid::RANDOMACCESS);
  float rValue;
  for (int i=0; i<nxp_; i++) {
    for (int j=0; j<nyp_; j++) {
      for (int k=0; k<nzp_; k++) {
        if (difftype == FIRSTORDERFORWARDDIFF) {
          if (k == nzp_-1) 
            rValue = wavelet->getRAmp(0,j,i) - wavelet->getRAmp(k,j,i);
          else
            rValue = wavelet->getRAmp(k+1,j,i) - wavelet->getRAmp(k,j,i);
        }
        else { //(difftype == FIRSTORDERBACKWARDDIFF)
          if (k == 0)
            rValue = wavelet->getRAmp(k,j,i) - wavelet->getRAmp(nzp_-1,j,i);
          else
            rValue = wavelet->getRAmp(k,j,i) -wavelet->getRAmp(k-1,j,i);
        }
        ampCube_.setRealValue(i,j,k,rValue,true);
      }
    }
  }
}

Wavelet3D::Wavelet3D(Wavelet * wavelet)
  : Wavelet(wavelet, 3)
{
  assert(wavelet->getDim() == 3);
  assert(wavelet->getIsReal());
  nx_ = wavelet->getNx();
  ny_ = wavelet->getNy();
  nxp_ = wavelet->getNxp();
  nyp_ = wavelet->getNyp();
  dx_ = wavelet->getDx();
  dy_ = wavelet->getDy();
  ampCube_ = FFTGrid(nx_, ny_, nz_, nxp_, nyp_, nzp_);
  ampCube_.createRealGrid();
  ampCube_.setType(FFTGrid::COVARIANCE);
  ampCube_.setAccessMode(FFTGrid::RANDOMACCESS);
  for (int k=0; k<nzp_; k++) {
    for (int j=0; j<nyp_; j++) {
      for (int i=0; i<nxp_; i++) {
        float rvalue = wavelet->getRAmp(k,j,i);
        ampCube_.setRealValue(i,j,k,rvalue,true);
      }
    }
  }
}

void           
Wavelet3D::fft1DInPlace()
{
  ampCube_.fftInPlace();
  isReal_ = false;
}

void
Wavelet3D::invFFT1DInPlace()
{
  ampCube_.invFFTInPlace();
  isReal_ = true;
}

bool           
Wavelet3D::consistentSize(int nzp, int nyp, int nxp) const 
{ 
  bool ok = true;
  if (nzp!=nzp_) {
    printf("nzp=%d  nzp_wavelet3D=%d\n",nzp,nzp_);
    ok = false;
  }
  if (nyp != nyp_) {
    printf("nyp=%d  nyp_wavelet3D=%d\n",nyp,nyp_);
    ok = false;
  }
  if (nxp != nxp_) {
    printf("nxp=%d  nxp_wavelet3D=%d\n",nxp,nxp_);
    ok = false;
  }
  return (ok);
}

fftw_complex   
Wavelet3D::getCAmp(int k, int j, int i) const
{
  return(ampCube_.getComplexValue(i,j,k,true));
}

fftw_real      
Wavelet3D::getRAmp(int k, int j, int i)
{
  return(ampCube_.getRealValue(i,j,k,true));
}

fftw_complex   
Wavelet3D::getCAmp(int k, float, int j, int i) const
{
  return(ampCube_.getComplexValue(i,j,k,true));
}

void           
Wavelet3D::setRAmp(float value, int k, int j, int i)
{
  ampCube_.setRealValue(i,j,k,value);
}

void
Wavelet3D::setCAmp(fftw_complex value, int k, int j, int i)
{
  ampCube_.setComplexValue(i, j ,k, value, true);
}

void
Wavelet3D::scale(float scale)
{
  Wavelet::scale(scale);

  float rAmp;
  for(int i=0; i < nxp_ ; i++) {
    for (int j=0; j < nyp_; j++) {
      for (int k=0; k < nzp_; k++) {
        rAmp = getRAmp(k,j,i);
        if (rAmp != RMISSING) {
          rAmp *= scale;
          setRAmp(rAmp,k,j,i);
        }
      }
    }
  }
}

void
Wavelet3D::multiplyByR(float p)
{
  assert(!getIsReal());
  float scale = static_cast<float> (1.0/(nxp_*nyp_*nzp_));
  int cnxp = nxp_/2+1;
  int i, j, k;
  float kx, ky, kz;
  float radius;
  fftw_complex cValue;

  for (k=0; k<=nzp_/2; k++) {
    kz = static_cast<float> (k);
    for (j=0; j<=nyp_/2; j++) {
      ky = static_cast<float> (j);
      for (i=0; i<cnxp; i++) {
        kx = static_cast<float> (i);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        cValue = getCAmp(k,j,i);
        cValue.re *= scale * pow(radius, p);
        cValue.im *= scale * pow(radius, p);
        setCAmp(cValue,k,j,i);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      ky = static_cast<float> (nyp_-j);
      for (i=0; i<cnxp; i++) {
        kx = static_cast<float> (i);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        cValue = getCAmp(k,j,i);
        cValue.re *= scale * pow(radius, p);
        cValue.im *= scale * pow(radius, p);
        setCAmp(cValue,k,j,i);
      }
    }
  }
  for (k=(nzp_/2)+1; k<nzp_; k++) {
    kz = static_cast<float> (nzp_-k);
    for (j=0; j<=nyp_/2; j++) {
      ky = static_cast<float> (j);
      for (i=0; i<cnxp; i++) {
        kx = static_cast<float> (i);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        cValue = getCAmp(k,j,i);
        cValue.re *= scale * pow(radius, p);
        cValue.im *= scale * pow(radius, p);
        setCAmp(cValue,k,j,i);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      ky = static_cast<float> (nyp_-j);
      for (i=0; i<cnxp; i++) {
        kx = static_cast<float> (i);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        cValue = getCAmp(k,j,i);
        cValue.re *= scale * pow(radius, p);
        cValue.im *= scale * pow(radius, p);
        setCAmp(cValue,k,j,i);
      }
    }
  }
}



void
Wavelet3D::printToFile(std::string fileName, bool overrideDebug) 
{
  if(overrideDebug == true || ModelSettings::getDebugLevel() > 0) {
    std::string fName = ModelSettings::makeFullFileName(fileName+".txt");
    FILE *file = fopen(fName.c_str(),"w");
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM ascii file "+fName+"...");
    int i,j,k;
    for(k=0;k<nzp_;k++)
      for(j=0;j<nyp_;j++) {
        for(i=0;i<nxp_;i++) {
          fprintf(file,"%f ", getRAmp(k,j,i));
        }
        fprintf(file,"\n");
      }
    fprintf(file,"0\n");
    fclose(file);
  }
}

void
Wavelet3D::writeWaveletToFile(char* fileName, float, Simbox *simbox)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"  Writing 3D-Wavelet to file. \n");
  ampCube_.writeFile(fileName, simbox, "", false);
}

void
Wavelet3D::shiftFFTGrid(FFTGrid *shiftAmp)
{
  int cx = nx_/2;
  int cy = ny_/2;
  int cz = nz_/2;
  int xOdd=0 , yOdd=0, zOdd=0;
  if (nx_/2 != (nx_+1)/2)
    xOdd=1;
  if (ny_/2 != (ny_+1)/2)
    yOdd=1;
  if (nz_/2 != (nz_+1)/2)
    zOdd=1;

  int i,j,k;
  float value;
  for (i=0; i<cx+xOdd; i++) {
    for (j=0; j<cy+yOdd; j++) {
      for (k=0; k<cz+zOdd; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i+cx,j+cy,k+cz,value);
      }
      for (k=nzp_-cz; k<nzp_; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i+cx,j+cy,k-(nzp_-cz),value);
      }
    }
    for (j=nyp_-cy; j<nyp_; j++) {
      for (k=0; k<cz; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i+cx,j-(nyp_-cy),k+cz,value);
      }
      for (k=nzp_-cz; k<nzp_; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i+cx,j-(nyp_-cy),k-(nzp_-cz),value);
      }
    }
  }
  for (i=nxp_-cx; i<nxp_; i++) {
    for (j=0; j<cy; j++) {
      for (k=0; k<cz; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i-(nxp_-cx),j+cy,k+cz,value);
      }
      for (k=nzp_-cz; k<nzp_; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i-(nxp_-cx),j+cy,k-(nzp_-cz),value);
      }
    }
    for (j=nyp_-cy; j<nyp_; j++) {
      for (k=0; k<cz; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i-(nxp_-cx),j-(nyp_-cy),k+cz,value);
      }
      for (k=nzp_-cz; k<nzp_; k++) {
        value = ampCube_.getRealValue(i,j,k,true);
        shiftAmp->setRealValue(i-(nxp_-cx),j-(nyp_-cy),k-(nzp_-cz),value);
      }
    }
  }
}
