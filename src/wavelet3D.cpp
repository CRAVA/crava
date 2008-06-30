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
#include "lib/log.h"
#include "lib/sgri.h"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

Wavelet3D::Wavelet3D(char * fileName, ModelSettings * modelSettings, Simbox *simBox, float theta, int dim)
  :Wavelet(modelSettings, dim)
{
  nx_ = simBox->getnx();
  ny_ = simBox->getny();
  nz_ = simBox->getnz();
  dx_ = static_cast <float> (simBox->getdx());
  dy_ = static_cast <float> (simBox->getdy());
  dz_ = static_cast <float> (simBox->getdz());
  theta_ = theta;

  float xPad = modelSettings->getXpad();
  nxp_   =  findClosestFactorableNumber( (int) ceil( nx_*(1.0f+xPad) ) );
  float yPad = modelSettings->getYpad();
  nyp_   =  findClosestFactorableNumber( (int) ceil( ny_*(1.0f+yPad) ) );
  float zPad = modelSettings->getZpad();
  nzp_   =  findClosestFactorableNumber( (int) ceil( nz_*(1.0f+zPad) ) );

  ampCube_ = FFTGrid(nx_, ny_, nz_, nxp_, nyp_, nzp_);
  ampCube_.createRealGrid();
  ampCube_.setType(FFTGrid::COVARIANCE);
  ampCube_.setAccessMode(FFTGrid::RANDOMACCESS);
  ampCube_.setOutputFormat(modelSettings->getFormatFlag());

  readtype_ = Wavelet::SGRI;
	Sgri *sgri = new Sgri(fileName, errText_, errCode_);
  if (errCode_ != 0) {
    LogKit::writeLog("Error when reading sgri from file %s.\n",fileName);
    LogKit::writeLog("%s \n", errText_);
    exit(1);
  }

  float xLim = 0.5f * dx_ * nxp_;
  float yLim = 0.5f * dy_ * nyp_;
  float zLim = 0.5f * dz_ * nzp_;
  if(!sgri->sizeOk(xLim, yLim, zLim)) {
    if (errCode_ == 0)  
      sprintf(errText_,"3-D wavelet read from file %s has too big size compared to padded fft-grid.\n", fileName);
    else  
      sprintf(errText_,"%s3-D wavelet read from file %s has too big size compared to padded fft-grid.\n",errText_,fileName); 
    errCode_=1; 
    LogKit::writeLog("\nToo big wavelet grid read from file  %s.\n",fileName);
    LogKit::writeLog("%s \n", errText_);
    exit(1);
  }

  float sf = 100.0;
  int i,j,k;
  float x, y, z, value;
  for (k=0; k<=nzp_/2; k++) {
    z = k * dz_;
    for (j=0; j<=nyp_/2; j++) {
      y = j * dy_;
      for (i=0; i<=nxp_/2; i++) {
        x = i * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
      for (i=(nxp_/2)+1; i<nxp_; i++) {
        x = (i-nxp_) * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      y = (j-nyp_) * dy_;
      for (i=0; i<=nxp_/2; i++) {
        x = i * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
      for (i=(nxp_/2)+1; i<nxp_; i++) {
        x = (i-nxp_) * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
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
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
      for (i=(nxp_/2)+1; i<nxp_; i++) {
        x = (i-nxp_) * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      y = (j-nyp_) * dy_;
      for (i=0; i<=nxp_/2; i++) {
        x = i * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      }
      for (i=(nxp_/2)+1; i<nxp_; i++) {
        x = (i-nxp_) * dx_;
        value = sf * sgri->getWaveletValue(x, y, z);
        ampCube_.setRealValue(i, j, k, value, true);
      } 
    }
  }
  
  FFTGrid *shiftAmp = new FFTGrid(nx_, ny_, nz_, nx_, ny_, nz_);
  shiftAmp->fillInConstant(0.0);
  shiftAmp->setType(FFTGrid::DATA);
  shiftAmp->setAccessMode(FFTGrid::RANDOMACCESS);
  shiftAmp->setOutputFormat(modelSettings->getFormatFlag());

  shiftFFTGrid(shiftAmp);
  shiftAmp->writeFile("3DWavelet", simBox);
  delete shiftAmp;

  delete sgri;
}

Wavelet3D::Wavelet3D(Wavelet * wavelet, int dim)
  : Wavelet(wavelet, dim)
{
}

void           
Wavelet3D::fft1DInPlace()
{
  ampCube_.fftInPlace();
}

void
Wavelet3D::invFFT1DInPlace()
{
  ampCube_.invFFTInPlace();
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
Wavelet3D::printToFile(char* fileName, bool overrideDebug) 
{
  if(overrideDebug == true || LogKit::getDebugLevel() > 0) {
    char * fName = LogKit::makeFullFileName(fileName, ".txt");
    FILE *file = fopen(fName,"w");
    LogKit::writeLog("\nWriting STORM ascii file %s...",fName);
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
      delete [] fName;
  }
}

void
Wavelet3D::writeWaveletToFile(char* fileName, float, Simbox *simbox)
{
  ampCube_.writeFile(fileName, simbox, false);
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
   
