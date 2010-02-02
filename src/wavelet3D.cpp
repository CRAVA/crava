#include <iostream>
#include <fstream>

#include <string.h>
#include <assert.h>
#include <math.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/global_def.h"
#include "lib/lib_matr.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/surface/surfaceio.hpp"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/definitions.h"
#include "src/wavelet3D.h"
#include "src/wavelet1D.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/model.h"
#include "src/io.h"
#include "src/waveletfilter.h"

Wavelet3D::Wavelet3D(const std::string            & filterFile,
                     const std::vector<Surface *> & estimInterval,
                     const NRLib::Grid2D<float>   & refTimeGradX,
                     const NRLib::Grid2D<float>   & refTimeGradY,
                     FFTGrid                      * seisCube,
                     ModelSettings                * modelSettings,
                     WellData                    ** wells,
                     Simbox                       * simBox,
                     float                        * reflCoef,
                     int                            angle_index,
                     float                          theta,
                     int                          & errCode,
                     std::string                  & errText)
  : Wavelet(3),
    filter_(filterFile, errCode, errText)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"  Estimating 3D wavelet pulse from seismic data and (nonfiltered) blocked wells\n");
  maxShift_             = modelSettings->getMaxWaveletShift();
  minRelativeAmp_       = modelSettings->getMinRelWaveletAmp();  theta_ = theta;
  norm_ = RMISSING;
  cz_         = 0;
  inFFTorder_ = true;
  coeff_[0] = reflCoef[0];
  coeff_[1] = reflCoef[1];
  coeff_[2] = reflCoef[2];

  nz_               = simBox->getnz();
  float dx          = static_cast<float>(simBox->getdx());
  float dy          = static_cast<float>(simBox->getdy());
  dz_               = static_cast<float>(simBox->getdz());
  nzp_              = seisCube->getNzp();

  unsigned int nWells = modelSettings->getNumberOfWells();
  float v0            = modelSettings->getAverageVelocity();
  float stretch       = modelSettings->getStretchFactor(angle_index);
  bool hasHalpha      = filter_.hasHalpha();
  
  int nhalfWl         = static_cast<int> (0.5 * modelSettings->getWaveletTaperingL() / dz_);
  int nWl             = 2 * nhalfWl + 1;
  std::vector<std::vector<float> > wlestWell(nWells, std::vector<float>(nWl, 0.0));

  int nActiveWells = 0;
  for (unsigned int w=0; w<nWells; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::MEDIUM, "  Well :  %s\n", wells[w]->getWellname().c_str());

      BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick();  
 
      std::vector<double> tGradX;//(bl->getNumberOfBlocks(),0.0);
      std::vector<double> tGradY;//(bl->getNumberOfBlocks(),0.0);

//    NBNB-Frode: Fyll gradX og gradY for alle nzp_ ved hjelp av kall til
      bl->findSeismicGradient(seisCube,
                              simBox,
                              tGradX,
                              tGradY);
      const int * iPos = bl->getIpos();
      const int * jPos = bl->getJpos();
      
      unsigned int nBlocks = bl->getNumberOfBlocks();
      float * zGradX  = new float[nBlocks];
      float * zGradY  = new float[nBlocks];
      float * t0GradX = new float[nBlocks];
      float * t0GradY = new float[nBlocks];
      for (unsigned int b = 0; b<nBlocks; b++) {
        t0GradX[b]    = refTimeGradX(iPos[b], jPos[b]);
        zGradX[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradX[b]) - t0GradX[b]); //0.001f is due to ms/s conversion
        t0GradY[b]    = refTimeGradY(iPos[b], jPos[b]);
        zGradY[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradY[b]) - t0GradY[b]);
      }

      float * az = new float[nz_]; 
      bl->getVerticalTrend(zGradX, az);
      delete [] zGradX;
      float * bz = new float[nz_];
      bl->getVerticalTrend(zGradY, bz);
      delete [] zGradY;
      float * at0 = new float[nz_]; 
      bl->getVerticalTrend(t0GradX, at0);
      delete [] t0GradX;
      float * bt0 = new float[nz_];
      bl->getVerticalTrend(t0GradY, bt0);
      delete [] t0GradY;

      bool  *hasWellData = new bool[nz_];
      findLayersWithData(estimInterval, 
                         bl, 
                         seisCube, 
                         az, 
                         bz, 
                         hasWellData);

      int start, length;
      bl->findContiniousPartOfData(hasWellData, 
                                   nz_, 
                                   start, 
                                   length);
      if (length > nWl) {
        nActiveWells++;
        float *cpp         = new float[nzp_];
        bl->fillInCpp(coeff_,start,length,cpp,nzp_);

        std::vector<float> cppAdj(length, 0.0);
        std::vector<float> Halpha;
        if (hasHalpha)
          Halpha.resize(length,0.0);

        for (int i=start; i < start+length-1; i++) {
          double phi    = findPhi(az[i], bz[i]);
          double r      = sqrt(az[i]*az[i] + bz[i]*bz[i] + 1);
          double psi    = acos(1.0 / r);
          float alpha1  = filter_.getAlpha1(phi, psi);
          float f       = static_cast<float> (1.0 * r / stretch);
          cppAdj[i-start]     = cpp[i] * alpha1 / f;
          if (hasHalpha)
            Halpha[i-start]   = filter_.getHalpha(phi, psi);
        }
        delete [] cpp;

        std::ofstream cppFile;
        NRLib::OpenWrite(cppFile, "cppAdj.dat");
        for (int i=0; i<length; i++)
          cppFile << std::setprecision(4) << std::setw(10) << cppAdj[i] << std::endl;
        cppFile.close();

        float *zLog = new float[nBlocks];
        for (unsigned int b=0; b<nBlocks; b++) {
          double zTop     = simBox->getTop(iPos[b], jPos[b]);
          zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(iPos[b], jPos[b]) * dz_);
        }
        float * zPosWell = new float[nz_];
        bl->getVerticalTrend(zLog, zPosWell);
        delete zLog;

        int nTracesX    = static_cast<int> (modelSettings->getEstRangeX(angle_index) / dx);
        int nTracesY    = static_cast<int> (modelSettings->getEstRangeY(angle_index) / dy);

        std::vector<std::vector<float> > gMat;

        std::vector<float> dVec;
        int nPoints = 0;
        for (int xTr = -nTracesX; xTr <= nTracesX; xTr++) {
          for (int yTr = -nTracesY; yTr <= nTracesY; yTr++) {
            float * seisLog = new float[nBlocks];
            bl->getBlockedGrid(seisCube, seisLog, xTr, yTr);
            float *seisData = new float[nz_];
            bl->getVerticalTrend(seisLog, seisData);
            delete [] seisLog;
            float *zLog = new float[nBlocks];
            for (unsigned int b=0; b<nBlocks; b++) {
              int xIndex      = iPos[b] + xTr;
              int yIndex      = jPos[b] + yTr;
              double zTop     = simBox->getTop(xIndex, yIndex);
              zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(xIndex, yIndex) * dz_);
            }
            float * zPosTrace = new float[nz_];
            bl->getVerticalTrend(zLog, zPosTrace);
            delete zLog;
            for (int t=start; t < start+length; t++) {
              if (seisData[t] != RMISSING) {
                dVec.push_back(seisData[t]);
                float **lambda = new float *[length];
                for (int i=0; i<length; i++) {
                  lambda[i]     = new float[nWl];
                  for (int j=0; j<nWl; j++)
                    lambda[i][j] = 0.0;
                }
                for (int tau = start; tau < start+length; tau++) {
                  //Hva gjør vi hvis zData[t] er RMISSING. Kan det skje?
                  float at = at0[tau] + 2.0f*az[tau]/v0;
                  float bt = bt0[tau] + 2.0f*bz[tau]/v0;
                  float u = static_cast<float> (zPosTrace[t] - zPosWell[tau] - at*(xTr*dx) - bt*(yTr*dy));
//                  float u = static_cast<float> (zPosTrace[t] - zPosWell[tau] - az[tau]*(xTr*dx_) - bz[tau]*(yTr*dy_));
                  if (hasHalpha) {
                    for (int i=0; i<nWl; i++) {
                      float v = u - static_cast<float>((i - nhalfWl)*dz_);
                      float h = Halpha[tau-start];
                      lambda[tau-start][i] = static_cast<float> (h / (PI *(h*h + v*v)));
                    }
                  }
                  else {
                    int tLow  = static_cast<int> (floor(u / dz_));
                    int tHigh = tLow + 1;
                    float lambdaValue = u - static_cast<float> (tLow * dz_);
                    if (u >= 0.0 && tLow <= nhalfWl) { 
                      lambda[tau-start][tLow]   = 1 -lambdaValue;
                      lambda[tau-start][tHigh]  = lambdaValue; 
                    }
                    else if (u < 0.0 && -tHigh <= nhalfWl) {
                      if (tHigh < 0)
                        lambda[tau-start][tHigh+nWl] = lambdaValue;
                      else
                        lambda[tau-start][0] = lambdaValue;
                      if (tLow >= -nhalfWl)
                        lambda[tau-start][tLow+nWl] = 1 - lambdaValue;
                    }
                  }
                }
                std::vector<float> gVec(nWl);
                for (int i=0; i<nWl; i++) {
                  gVec[i] = 0.0f;
                  for (int j=0; j<length; j++)
                    gVec[i] += cppAdj[j] * lambda[j][i];
                }
                gMat.push_back(gVec);
                nPoints++;
                for (int i=0; i<length; i++)
                  delete [] lambda[i];
                delete [] lambda;
              } //if (seisData[t] != RMISSING)
            }
          }
        }

        std::ofstream gMatFile;
        NRLib::OpenWrite(gMatFile, "gMat.dat");
        for (int i=0; i<nPoints; i++) {
          gMatFile << "gMat[" << std::setw(2) << i << "]: ";
          for (int j=0; j<nWl; j++)
            gMatFile << std::setprecision(4) << std::setw(10) << gMat[i][j];
          gMatFile << std::endl;
        }
        gMatFile.close();

        double **gTrg = new double *[nWl];
        for (int i=0; i<nWl; i++) {
          gTrg[i]     = new double[nWl];
          for (int j=0; j<nWl; j++) {
            gTrg[i][j] = 0.0;
            for (int k=0; k<nPoints; k++)
              gTrg[i][j] += gMat[k][i]*gMat[k][j];
          }
        }

        std::ofstream gTrgMatFile;
        NRLib::OpenWrite(gTrgMatFile, "gTrgMat.dat");
        for (int i=0; i<nWl; i++) {
          gTrgMatFile << "gTrg[" << std::setw(2) << i << "]: ";
          for (int j=0; j<nWl; j++)
            gTrgMatFile << std::setprecision(4) << std::setw(10) << gTrg[i][j];
          gTrgMatFile << std::endl;
        }
        gTrgMatFile.close();

        std::ofstream dFile;
        NRLib::OpenWrite(dFile, "data.dat");
        for (int i=0; i<nPoints; i++)
          dFile << std::setprecision(4) << std::setw(10) << dVec[i] << std::endl;
        dFile.close();

        double *gTrd = new double[nWl];
        for (int i=0; i<nWl; i++) {
          gTrd[i] = 0.0;
          for (int j=0; j<nPoints; j++)
            gTrd[i] += gMat[j][i] * dVec[j];
        }

        lib_matrCholR(nWl, gTrg);
        lib_matrAxeqbR(nWl, gTrg, gTrd);
        for (int i=0; i<nWl; i++)
          wlestWell[w][i] = static_cast<float> (gTrd[i]);

        cnzp_ = nzp_/2+1;
        rnzp_ = 2*cnzp_;
        rAmp_ = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  
        cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);

        int nWl = wlestWell.size();
        int nHalfWl = static_cast<int> (nWl/2);
        rAmp_[0] = static_cast<fftw_real> (wlestWell[w][0]);
        for (int i=1; i<=nHalfWl; i++) {
          rAmp_[i]      = static_cast<fftw_real> (wlestWell[w][i]);
          rAmp_[nzp_-i] = static_cast<fftw_real> (wlestWell[w][nWl-i]);
        }
        for (int i=nHalfWl+1; i<nzp_-nHalfWl; i++)
          rAmp_[i] = 0.0;
        for (int i=nzp_; i<rnzp_; i++)
          rAmp_[i] = RMISSING;

        waveletLength_ = getWaveletLengthF();

        delete [] gTrd;

        for (int i=0; i<nWl; i++)
          delete [] gTrg[i];
        delete [] gTrg;

        delete [] az;
        delete [] bz;

        delete [] hasWellData;

      }
      else {
        LogKit::LogFormatted(LogKit::MEDIUM,"     No enough data for 3D wavelet estimation in well %s\n", wells[w]->getWellname().c_str());
      }
    }
  }

  std::vector<float> wlest(nWl,0.0);
  for (int j=0; j<nWl; j++) {
    for(unsigned int w=0; w<nWells; w++)
      wlest[j] += wlestWell[w][j];
    wlest[j] = wlest[j] / static_cast<float>(nActiveWells);
  }

  std::ofstream wlFile;
  NRLib::OpenWrite(wlFile, "wlest.dat");
  for (int i=0; i<nWl; i++)
    wlFile << std::setprecision(4) << std::setw(10) << wlest[i] << std::endl;
  wlFile.close();
}

Wavelet3D::Wavelet3D(const std::string & fileName, 
            int                 fileFormat, 
            ModelSettings     * modelSettings, 
            float             * reflCoef,
            float               theta,
            int               & errCode, 
            std::string       & errText,
            const std::string & filterFile)
  : Wavelet(fileName, fileFormat, modelSettings, reflCoef, theta, 3, errCode, errText),
    filter_(filterFile, errCode, errText)
{
}

Wavelet3D::Wavelet3D(Wavelet * wavelet)
  : Wavelet(3, wavelet)
{
  assert(wavelet->getDim() == 3);
  assert(wavelet->getIsReal());

}

Wavelet3D::~Wavelet3D()
{
}

void
Wavelet3D::findLayersWithData(const std::vector<Surface *> & estimInterval,
                              BlockedLogs                  * bl,
                              FFTGrid                      * seisCube,
                              float                        * az,
                              float                        * bz,
                              bool                         * hasWellData) const
{
  float *seisLog     = new float[bl->getNumberOfBlocks()];
  bl->getBlockedGrid(seisCube, seisLog);
  float *seisData    = new float[nz_];
  bl->getVerticalTrend(seisLog, seisData);
  delete [] seisLog;

  float *alpha       = new float[nz_];
  float *beta        = new float[nz_];
  float *rho         = new float[nz_];
  bl->getVerticalTrend(bl->getAlpha(), alpha);
  bl->getVerticalTrend(bl->getBeta(), beta);
  bl->getVerticalTrend(bl->getRho(), rho);

  for (int k=0; k<nz_; k++) 
    hasWellData[k] = (alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING && az[k] != RMISSING && bz[k] != RMISSING && seisData[k] != RMISSING);
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;

  //Check that data are within wavelet estimation interval
  if (estimInterval.size() > 0) {
    const double *xPos = bl->getXpos();
    const double *yPos = bl->getYpos();
    const double *zPos = bl->getZpos();
    for (int k=0; k<nz_; k++) {
      const double zTop  = estimInterval[0]->GetZ(xPos[k],yPos[k]);
      const double zBase = estimInterval[1]->GetZ(xPos[k],yPos[k]);
      if ((zPos[k]-0.5*dz_) < zTop || (zPos[k]+0.5*dz_) > zBase)
        hasWellData[k] = false;
    }
  }
}

double 
Wavelet3D::findPhi(float a, float b) const
//Return value should be between 0 and 2*PI
{
  double phi;
  double epsilon = 0.000001;
  if (a > epsilon && b >= 0.0) //1. quadrant 
    phi = atan(b/a);
  else if (a > epsilon && b < 0.0) //4. quadrant
    phi = 2*PI + atan(b/a);
  else if (a < - epsilon && b >= 0.0) //2. quadrant
    phi = PI + atan(b/a);
  else if (a < - epsilon && b < 0.0) //3. quadrant
    phi = PI + atan(b/a);
  else if (b  >= 0.0) //kx very small
    phi = 0.5 * PI;
  else //kx very small
    phi = 1.5 * PI;

  return(phi);
}

double 
Wavelet3D::findPsi(float radius, float kz) const
//Return value should be between 0 and 0.5*PI
{
  double epsilon = 0.000001;
  double psi = 0.0;
  if (kz < 0.0)
    kz = -kz;
  if (radius > epsilon)
    psi = acos(kz/radius);

  return(psi);
}

fftw_complex 
Wavelet3D::findWLvalue(Wavelet1D       * wavelet1d,
                                    float             omega) const
{
  int lowindex = static_cast<int> (omega / wavelet1d->getDz());
  fftw_complex c_low, c_high;
  if (lowindex >= wavelet1d->getNz()) {
    c_low.re = 0.0; 
    c_low.im = 0.0;
    c_high.re = 0.0;
    c_high.im = 0.0;
  }
  else if (lowindex == wavelet1d->getNz()-1) {
    c_high.re = 0.0;
    c_high.im = 0.0;
    c_low = wavelet1d->getCAmp(lowindex);
  }
  else {
    c_low = wavelet1d->getCAmp(lowindex);
    c_high = wavelet1d->getCAmp(lowindex + 1);
  }
  float fac = omega - lowindex * wavelet1d->getDz();
  fftw_complex cValue;
  cValue.re = (1-fac) * c_low.re + fac * c_high.re;
  cValue.im = (1-fac) * c_low.im + fac * c_high.im;

  return cValue;
}
