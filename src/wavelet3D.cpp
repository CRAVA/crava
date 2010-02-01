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
  : Wavelet(modelSettings, 3, reflCoef),
    filter_(filterFile, errCode, errText)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"  Estimating 3D wavelet pulse from seismic data and (nonfiltered) blocked wells\n");
  readtype_     = ESTIMATE;
  theta_ = theta;
  norm_ = RMISSING;

  nx_             = simBox->getnx();
  ny_             = simBox->getny();
  nz_             = simBox->getnz();
  dx_             = static_cast<float>(simBox->getdx());
  dy_             = static_cast<float>(simBox->getdy());
  dz_             = static_cast<float>(simBox->getdz());
  nxp_            = seisCube->getNxp();
  nyp_            = seisCube->getNyp();
  nzp_            = seisCube->getNzp();

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

        int nTracesX    = static_cast<int> (modelSettings->getEstRangeX(angle_index) / dx_);
        int nTracesY    = static_cast<int> (modelSettings->getEstRangeY(angle_index) / dy_);

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
                  float u = static_cast<float> (zPosTrace[t] - zPosWell[tau] - at*(xTr*dx_) - bt*(yTr*dy_));
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

  ampCube_ = FFTGrid(0, 0, 0, 0, 0, 0);
  ampCube_.createComplexGrid();
  wavelet1D_ = new Wavelet1D(static_cast<Wavelet *> (this), modelSettings, reflCoef, wlest);
}


Wavelet3D::Wavelet3D(Wavelet1D           * wavelet1D,
                     const std::string   & filterFile,
                     ModelSettings       * modelSettings,
                     int                   angle_index,
                     Simbox              * simBox,
                     float                 theta,
                     int                 & errCode,
                     std::string         & errText)
  : Wavelet(3),
    wavelet1D_(wavelet1D),
    filter_(filterFile, errCode, errText)
{
  float v0 = modelSettings->getAverageVelocity();
  nx_ = simBox->getnx();
  ny_ = simBox->getny();
  nz_ = simBox->getnz();
  dx_ = static_cast<float>(simBox->getdx());
  dy_ = static_cast<float>(simBox->getdy());
  dz_ = static_cast<float>(simBox->getdz() * 0.5f * v0 * 0.001f);

  double xPadFac = modelSettings->getXPadFac();
  nxp_   =  FFTGrid::findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(nx_)*(1.0+xPadFac) )) );
  double yPadFac = modelSettings->getYPadFac();
  nyp_   =  FFTGrid::findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(ny_)*(1.0+yPadFac) )) );
  double zPadFac = modelSettings->getZPadFac();
  nzp_   =  FFTGrid::findClosestFactorableNumber( static_cast<int>(ceil( static_cast<double>(nz_)*(1.0+zPadFac) )) );
  
  theta_ = theta;
  norm_ = RMISSING;

  wavelet1D_->resample(static_cast<float>(simBox->getdz()), 
                      nz_, 
                      static_cast<float> (zPadFac), 
                      theta);
  wavelet1D_->fft1DInPlace();
    
  ampCube_ = FFTGrid(nx_, ny_, nz_, nxp_, nyp_, nzp_);
  ampCube_.createComplexGrid();
  ampCube_.setType(FFTGrid::COVARIANCE);
  ampCube_.setAccessMode(FFTGrid::RANDOMACCESS);

  int i, j, k;
  float kx, ky, kz;
  float radius, alpha1, hAlpha, alpha2, omega;
  double phi, psi;
  fftw_complex cValue;
  
  float stretch = modelSettings->getStretchFactor(angle_index);
  float minus2pi = static_cast<float> (-2.0 * PI);

  for (k=0; k<=nzp_/2; k++) {
    kz = static_cast<float> (k / dz_);
    for (j=0; j<=nyp_/2; j++) {
      ky = static_cast<float> (j / dy_);
      for (i=0; i<=nxp_/2; i++) {
        kx = static_cast<float> (i / dx_);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        phi = findPhi(kx, ky);
        psi = findPsi(radius, kz);
        alpha1 = static_cast<float> (filter_.getAlpha1(phi, psi));
        hAlpha = static_cast<float> (filter_.getHalpha(phi, psi));
        omega = (0.5f * v0 * radius) / stretch;
        cValue = findWLvalue(wavelet1D_, omega);
        alpha2 = exp(minus2pi * omega * hAlpha);
        cValue.re *= static_cast<fftw_real> (alpha1 * alpha2);
        setCAmp(cValue,k,j,i);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      ky = static_cast<float> ((j-nyp_) / dy_);
      for (i=0; i<=nxp_/2; i++) {
        kx = static_cast<float> (i / dx_);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        phi = findPhi(kx, ky);
        psi = findPsi(radius, kz);
        alpha1 = static_cast<float> (filter_.getAlpha1(phi, psi));
        hAlpha = static_cast<float> (filter_.getHalpha(phi, psi));
        omega = (0.5f * v0 * radius) / stretch;
        cValue = findWLvalue(wavelet1D_, omega);
        alpha2 = exp(minus2pi * omega * hAlpha);
        cValue.re *= static_cast<fftw_real> (alpha1 * alpha2);
        setCAmp(cValue,k,j,i);
      }
    }
  }
  for (k=(nzp_/2)+1; k<nzp_; k++) {
    kz = static_cast<float> ((nzp_-k) / dz_);
    for (j=0; j<=nyp_/2; j++) {
      ky = static_cast<float> (j / dy_);
      for (i=0; i<=nxp_/2; i++) {
        kx = static_cast<float> (i / dx_);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        phi = findPhi(kx, ky);
        psi = findPsi(radius, kz);
        alpha1 = static_cast<float> (filter_.getAlpha1(phi, psi));
        hAlpha = static_cast<float> (filter_.getHalpha(phi, psi));
        omega = (0.5f * v0 * radius) / stretch;
        cValue = findWLvalue(wavelet1D_, omega);
        alpha2 = exp(minus2pi * omega * hAlpha);
        cValue.re *= static_cast<fftw_real> (alpha1 * alpha2);
        setCAmp(cValue,k,j,i);
      }
    }
    for (j=(nyp_/2)+1; j<nyp_; j++) {
      ky = static_cast<float> ((nyp_-j) / dy_);
      for (i=0; i<=nxp_/2; i++) {
        kx = static_cast<float> (i / dx_);
        radius = sqrt(kx*kx + ky*ky + kz*kz);
        phi = findPhi(kx, ky);
        psi = findPsi(radius, kz);
        alpha1 = static_cast<float> (filter_.getAlpha1(phi, psi));
        hAlpha = static_cast<float> (filter_.getHalpha(phi, psi));
        omega = (0.5f * v0 * radius) / stretch;
        cValue = findWLvalue(wavelet1D_, omega);
        alpha2 = exp(minus2pi * omega * hAlpha);
        cValue.re *= static_cast<fftw_real> (alpha1 * alpha2);
        setCAmp(cValue,k,j,i);
      }
    }
  }
  invFFT1DInPlace();
  wavelet1D_->invFFT1DInPlace();
/*
  FFTGrid *shiftAmp = new FFTGrid(nx_, ny_, nz_, nx_, ny_, nz_);
  shiftAmp->fillInConstant(0.0);
  shiftAmp->setType(FFTGrid::DATA);
  shiftAmp->setAccessMode(FFTGrid::RANDOMACCESS);

  shiftFFTGrid(shiftAmp);

  std::ofstream headerFile;
  NRLib::OpenWrite(headerFile, "WL_as_shiftedFFTGrid.Sgrh");

  headerFile << "NORSAR General Grid Format v1.0\n";
  headerFile << "3\n";
  headerFile << "x (km)\n";
  headerFile << "y (km)\n";
  headerFile << "z (km)\n";
  headerFile << "FFT-grid\n";
  headerFile << "1\n";
  headerFile << "3D-wavelet" << std::endl;
  headerFile << "1 1 1\n";
  headerFile << nx_ << " " << ny_ << " " << nz_ << std::endl;
  headerFile << std::setprecision(10);
  headerFile << dx_*0.001 << " " << dy_*0.001 << " " << dz_*0.001 << std::endl;
  double x0 = 0.001 * (simBox->getx0() + 0.5 * dx_ * nx_);
  double y0 = 0.001 * (simBox->gety0() + 0.5 * dy_ * ny_);
  double z0 = 0.001 * (simBox->GetZMin(0,0) + 0.5 * dz_ * nz_);
//  headerFile << "0.0 0.0 0.0\n";
  headerFile << x0 << " " << y0 << " " << z0 << std::endl;
  headerFile << "0.0\n";
  headerFile << RMISSING << std::endl;
  headerFile << "WL_as_shiftedFFTGrid.Sgri\n";
  headerFile.close();

  std::ofstream binFile;
  NRLib::OpenWrite(binFile, "WL_as_shiftedFFTGrid.Sgri" , std::ios::out | std::ios::binary);
  float value;
  for (k=0; k<nz_; k++)
    for (j=0; j<ny_; j++)
      for (i=0; i<nx_; i++) {
          value = shiftAmp->getRealValue(i,j,k);
#ifndef BIGENDIAN
        NRLib::WriteBinaryFloat(binFile, value);
#else
        NRLib::WriteBinaryFloat(binFile, value, END_LITTLE_ENDIAN);
#endif
      }
  binFile.close();
  */
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
  nx_         = wavelet->getNx();
  ny_         = wavelet->getNy();
  nxp_        = wavelet->getNxp();
  nyp_        = wavelet->getNyp();
  dx_         = wavelet->getDx();
  dy_         = wavelet->getDy();

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

Wavelet3D::~Wavelet3D()
{
  if (wavelet1D_)
    delete wavelet1D_;
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
Wavelet3D::printToFile(const std::string & fileName, bool overrideDebug) 
{
  if(overrideDebug == true || ModelSettings::getDebugLevel() > 0) {
    std::string fName = fileName + IO::SuffixGeneralData();
    fName = IO::makeFullFileName(IO::PathToWavelets(), fName);
    std::ofstream file;
    NRLib::OpenWrite(file, fName);
    LogKit::LogFormatted(LogKit::LOW,"\nWriting STORM ascii file "+fName+"...");
    for(int k=0;k<nzp_;k++)
      for(int j=0;j<nyp_;j++) {
        for(int i=0;i<nxp_;i++) {
          file << getRAmp(k,j,i);
        }
        file << "\n";
      }
    file << "0\n";
    file.close();
  }
}

void
Wavelet3D::writeWaveletToFile(const std::string & fileName, float approxDzIn)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"  Writing 3D-Wavelet to file. \n");
  wavelet1D_->writeWaveletToFile(fileName, approxDzIn);
//  ampCube_.writeFile(fileName, IO::PathToWavelets(), simbox, "", false);
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
