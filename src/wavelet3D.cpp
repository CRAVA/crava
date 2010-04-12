#include <iostream>
#include <fstream>

#include <string.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

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

Wavelet3D::Wavelet3D(const std::string                          & filterFile,
                     const std::vector<Surface *>               & estimInterval,
                     const NRLib::Grid2D<float>                 & refTimeGradX,
                     const NRLib::Grid2D<float>                 & refTimeGradY,
                     const std::vector<std::vector<double> >    & tGradX,
                     const std::vector<std::vector<double> >    & tGradY,
                     FFTGrid                                    * seisCube,
                     ModelSettings                              * modelSettings,
                     WellData                                  ** wells,
                     Simbox                                     * simBox,
                     float                                      * reflCoef,
                     int                                          angle_index,
                     int                                        & errCode,
                     std::string                                & errText)
  : Wavelet(3),
    filter_(filterFile, errCode, errText)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"  Estimating 3D wavelet pulse from seismic data and (nonfiltered) blocked wells\n");

  theta_      = seisCube->getTheta();;
  norm_       = RMISSING;
  cz_         = 0;
  inFFTorder_ = true;
  isReal_     = true;
  coeff_[0]   = reflCoef[0];
  coeff_[1]   = reflCoef[1];
  coeff_[2]   = reflCoef[2];
  formats_    = modelSettings->getWaveletFormatFlag();

  nz_         = simBox->getnz();
  float dx    = static_cast<float>(simBox->getdx());
  float dy    = static_cast<float>(simBox->getdy());
  dz_         = static_cast<float>(simBox->getdz());
  nzp_        = seisCube->getNzp();
  cnzp_       = nzp_/2+1;
  rnzp_       = 2*cnzp_;

  unsigned int nWells    = modelSettings->getNumberOfWells();
  float        v0        = modelSettings->getAverageVelocity();
//  float        stretch   = modelSettings->getStretchFactor(angle_index);
   
  int nhalfWl         = static_cast<int> (0.5 * modelSettings->getWaveletTaperingL() / dz_);
  int nWl             = 2 * nhalfWl + 1;
  std::vector<std::vector<fftw_real> > wellWavelets(nWells/*, std::vector<float>(rnzp_, 0.0)*/);
  std::vector<float>                   wellWeight(nWells, 0.0);  
  std::vector<float>                   dzWell(nWells, 0.0);
  for (unsigned int w=0; w<nWells; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::MEDIUM, "  Well :  %s\n", wells[w]->getWellname().c_str());

      BlockedLogs *bl    = wells[w]->getBlockedLogsOrigThick();  
      const std::vector<int> iPos = bl->getIposVector();
      const std::vector<int> jPos = bl->getJposVector();
       
      std::vector<float> az(nz_); 
      std::vector<float> bz(nz_);
      std::vector<float> at0(nz_);
      std::vector<float> bt0(nz_);
      unsigned int nBlocks = bl->getNumberOfBlocks();
      calculateGradients(bl,
                         iPos,
                         jPos,
                         refTimeGradX,
                         refTimeGradY,
                         tGradX[w],
                         tGradY[w],
                         v0,
                         az,
                         bz,
                         at0,
                         bt0);


      std::vector<bool> hasWellData(nz_);
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

      dzWell[w]          = static_cast<float>(simBox->getRelThick(iPos[0],jPos[0]) * dz_);
      if (length > nWl) {
        std::vector<float> Halpha;
        std::vector<fftw_real> cppAdj = adjustCpp(bl,
                                                  az,
                                                  bz,
                                                  Halpha,
                                                  start,
                                                  length);

        printVecToFile("cpp_dipadjust", &cppAdj[0], length);

        std::vector<float> zLog(nBlocks);
        for (unsigned int b=0; b<nBlocks; b++) {
          float zTop     = static_cast<float> (simBox->getTop(iPos[b], jPos[b]));
//          zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(iPos[b], jPos[b]) * dz_);
          zLog[b]         = static_cast<float> (zTop + b * dzWell[w]);
        }
        std::vector<float> zPosWell(nz_);
        bl->getVerticalTrend(&zLog[0], &zPosWell[0]);

        int nTracesX    = static_cast<int> (modelSettings->getEstRangeX(angle_index) / dx);
        int nTracesY    = static_cast<int> (modelSettings->getEstRangeY(angle_index) / dy);

        std::vector<std::vector<float> > gMat;
        std::vector<float> dVec;
        int nPoints = 0;
        for (int xTr = -nTracesX; xTr <= nTracesX; xTr++) {
          for (int yTr = -nTracesY; yTr <= nTracesY; yTr++) {
            std::vector<float> seisLog(nBlocks);
            bl->getBlockedGrid(&seisCube[0], &seisLog[0], xTr, yTr);
            std::vector<float> seisData(nz_);
            bl->getVerticalTrend(&seisLog[0], &seisData[0]);
            for (unsigned int b=0; b<nBlocks; b++) {
              //int xIndex      = iPos[b] + xTr; //NBNB Frode: Gir warning fordi den ikkje blir brukt. Slett om du ikkje skal bruka.
              //int yIndex      = jPos[b] + yTr; //NBNB Frode: Gir warning fordi den ikkje blir brukt. Slett om du ikkje skal bruka.
              float zTop     = static_cast<float> (simBox->getTop(iPos[b], jPos[b]));
//              zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(xIndex, yIndex) * dz_);
              zLog[b]         = static_cast<float> (zTop + b * dzWell[w]);
            }
            std::vector<float> zPosTrace(nz_);
            bl->getVerticalTrend(&zLog[0], &zPosTrace[0]);
            for (int t=start; t < start+length; t++) {
              if (seisData[t] != RMISSING) {
                dVec.push_back(seisData[t]);
                std::vector<std::vector<float> > lambda(length, std::vector<float>(nWl,0.0)); //NBNB-Frode: Bør denne deklareres utenfor?
                for (int tau = start; tau < start+length; tau++) {
                  //Hva gjør vi hvis zData[t] er RMISSING. Kan det skje?
                  float at = at0[tau] + 2.0f*az[tau]/v0;
                  float bt = bt0[tau] + 2.0f*bz[tau]/v0;
                  float u = static_cast<float> (zPosTrace[t] - zPosWell[tau] - at*xTr*dx - bt*yTr*dy);
                  if (filter_.hasHalpha()) {
                    for (int i=0; i<nWl; i++) {
                      float v = u - static_cast<float>((i - nhalfWl)*dzWell[w]);
                      float h = Halpha[tau-start];
                      lambda[tau-start][i] = static_cast<float> (h / (M_PI *(h*h + v*v)));
                    }
                  }
                  else {
                    int tLow  = static_cast<int> (floor(u / dzWell[w]));
                    int tHigh = tLow + 1;
                    float lambdaValue = (u/dzWell[w]) - static_cast<float> (tLow);
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
                    } // else if
                  } // else
                } // for (tau=start...start+length)
                std::vector<float> gVec(nWl);
                for (int i=0; i<nWl; i++) {
                  gVec[i] = 0.0f;
                  for (int j=0; j<length; j++)
                    gVec[i] += cppAdj[j] * lambda[j][i];
                }
                gMat.push_back(gVec);
                nPoints++;
              } //if (seisData[t] != RMISSING)
            }
          }
        }

        printMatToFile("design_matrix", gMat, nPoints, nWl);
        printVecToFile("seismic", &dVec[0], nPoints);
        wellWavelets[w] = calculateWellWavelet(gMat,
                                               dVec,
                                               nWl,
                                               nhalfWl,
                                               nPoints);
        printVecToFile("well_wl", &wellWavelets[w][0], nzp_);
 
        wellWeight[w] = calculateWellWeight(nWl,
                                            nPoints,
                                            gMat,
                                            wellWavelets[w],
                                            dVec);
      } //if (length > nWl)
      else {
        LogKit::LogFormatted(LogKit::MEDIUM,"     No enough data for 3D wavelet estimation in well %s\n", wells[w]->getWellname().c_str());
      }
    } // if(wells->getUseForEstimation)
  } // for (w=0...nWells) 

  rAmp_ = averageWavelets(wellWavelets, nWells, nzp_, wellWeight, dzWell, dz_);
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);
  waveletLength_ = findWaveletLength(modelSettings->getMinRelWaveletAmp());
  LogKit::LogFormatted(LogKit::LOW,"  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( ModelSettings::getDebugLevel() > 0 )
    writeWaveletToFile("estimated_wavelet_", 1.0f);

  norm_ = findNorm();

  fftw_real * trueAmp = rAmp_;
  rAmp_               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
  cAmp_               = reinterpret_cast<fftw_complex *>(rAmp_);

  for(unsigned int w=0; w<nWells; w++) {
    if (wells[w]->getUseForWaveletEstimation() && 
       ((modelSettings->getWaveletOutputFlag() & IO::WELL_WAVELETS)>0 || modelSettings->getEstimationMode()))
    {
      for(int i=0; i<nzp_; i++)
        rAmp_[i] = wellWavelets[w][i];
      std::string wellname(wells[w]->getWellname());
      NRLib::Substitute(wellname,"/","_");
      NRLib::Substitute(wellname," ","_");
      std::string fileName = IO::PrefixWellWavelet() + wellname + "_"; 
      writeWaveletToFile(fileName, 1.0f);
    }
  }
  fftw_free(rAmp_);
  rAmp_ = trueAmp;
  cAmp_ = reinterpret_cast<fftw_complex *>(rAmp_);
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
}

Wavelet3D::~Wavelet3D()
{
}

float         
Wavelet3D::calculateSNRatioAndLocalWavelet(Simbox        * /*simbox*/, 
                                           FFTGrid       * /*seisCube*/, 
                                           WellData     ** /*wells*/, 
                                           ModelSettings * /*modelSettings*/,
                                           std::string   & errText, 
                                           int           & error, 
                                           int             /*number*/,
                                           Grid2D       *& /*noiseScaled*/, 
                                           Grid2D       *& /*shift*/, 
                                           Grid2D       *& /*gain*/)
{

 //Not possible to estimate signal-to-noise ratio for 3D wavelets
  errText += "Estimation of signal-to-noise ratio is not possible for 3D wavelets.\n";
  errText += "The s/n ratio must be specified in the model file\n";
  error++;
  return (1.0f);
}

void
Wavelet3D::findLayersWithData(const std::vector<Surface *> & estimInterval,
                              BlockedLogs                  * bl,
                              FFTGrid                      * seisCube,
                              const std::vector<float>     & az,
                              const std::vector<float>     & bz,
                              std::vector<bool>            & hasWellData) const
{
  std::vector<float> seisLog(bl->getNumberOfBlocks());
  bl->getBlockedGrid(seisCube, &seisLog[0]);
  std::vector<float> seisData(nz_);
  bl->getVerticalTrend(&seisLog[0], &seisData[0]);

  std::vector<float> alpha(nz_);
  std::vector<float> beta(nz_);
  std::vector<float> rho(nz_);
  bl->getVerticalTrend(bl->getAlpha(), &alpha[0]);
  bl->getVerticalTrend(bl->getBeta(), &beta[0]);
  bl->getVerticalTrend(bl->getRho(), &rho[0]);

  for (int k=0; k<nz_; k++) 
    hasWellData[k] = (alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING && az[k] != RMISSING && bz[k] != RMISSING && seisData[k] != RMISSING);

  //Check that data are within wavelet estimation interval
  if (estimInterval.size() > 0) {
    const std::vector<double> xPos = bl->getXposVector();
    const std::vector<double> yPos = bl->getYposVector();
    const std::vector<double> zPos = bl->getZposVector();
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
//Return value should be between 0 and 360
{
  double phi;
  double epsilon = 0.001;
  if (a > epsilon && b >= 0.0) //1. quadrant 
    phi = atan(b/a);
  else if (a > epsilon && b < 0.0) //4. quadrant
    phi = 2*M_PI + atan(b/a);
  else if (a < - epsilon && b >= 0.0) //2. quadrant
    phi = M_PI + atan(b/a);
  else if (a < - epsilon && b < 0.0) //3. quadrant
    phi = M_PI + atan(b/a);
  else if (b  >= 0.0) //kx very small
    phi = 0.5 * M_PI;
  else //kx very small
    phi = 1.5 * M_PI;
  phi = phi * 180 / M_PI;
  return(phi);
}

double 
Wavelet3D::findPsi(float r) const
//Return value should be between 0 and 90
{
  double psi    = acos(1.0 / r);
  psi = psi * 180 / M_PI;
  return(psi);
}

fftw_complex 
Wavelet3D::findWLvalue(float omega) const
{
  int lowindex = static_cast<int> (omega / dz_);
  fftw_complex c_low, c_high;
  if (lowindex >= nz_) {
    c_low.re = 0.0; 
    c_low.im = 0.0;
    c_high.re = 0.0;
    c_high.im = 0.0;
  }
  else if (lowindex == nz_-1) {
    c_high.re = 0.0;
    c_high.im = 0.0;
    c_low = getCAmp(lowindex);
  }
  else {
    c_low = getCAmp(lowindex);
    c_high = getCAmp(lowindex + 1);
  }
  float fac = omega - lowindex * dz_;
  fftw_complex cValue;
  cValue.re = (1-fac) * c_low.re + fac * c_high.re;
  cValue.im = (1-fac) * c_low.im + fac * c_high.im;

  return cValue;
}

std::vector<fftw_real>
Wavelet3D::adjustCpp(BlockedLogs              * bl,
                     const std::vector<float> & az,
                     const std::vector<float> & bz,
                     std::vector<float>       & Halpha,
                     int start,
                     int length) const
{
  std::vector<fftw_real> cpp(nzp_);
  bl->fillInCpp(coeff_, start, length, &cpp[0], nzp_);
  printVecToFile("cpp", &cpp[0], length);

  std::vector<fftw_real> cppAdj(length, 0.0);

  if (filter_.hasHalpha())
    Halpha.resize(length,0.0);
  for (int i=start; i < start+length-1; i++) {
    double phi    = findPhi(az[i], bz[i]);
    float r       = sqrt(az[i]*az[i] + bz[i]*bz[i] + 1);
    double psi    = findPsi(r);
    float alpha1  = filter_.getAlpha1(phi, psi);
    float stretch = std::cos(theta_);
    cppAdj[i-start]     = static_cast<fftw_real> (cpp[i] * alpha1 * stretch / r);
//    cppAdj[i-start]     = static_cast<fftw_real> (cpp[i]);
    if (filter_.hasHalpha())
      Halpha[i-start]   = filter_.getHalpha(phi, psi);
  }
  return cppAdj;
}


void
Wavelet3D::calculateGradients(BlockedLogs                * bl,
                              const std::vector<int>     & iPos,
                              const std::vector<int>     & jPos,
                              const NRLib::Grid2D<float> & refTimeGradX,
                              const NRLib::Grid2D<float> & refTimeGradY,
                              const std::vector<double>  & tGradX,
                              const std::vector<double>  & tGradY,
                              float                        v0,
                              std::vector<float>         & az,
                              std::vector<float>         & bz,
                              std::vector<float>         & at0,
                              std::vector<float>         & bt0) const
{
  unsigned int nBlocks = bl->getNumberOfBlocks();
  std::vector<float> zGradX(nBlocks);
  std::vector<float> zGradY(nBlocks);
  std::vector<float> t0GradX(nBlocks);
  std::vector<float> t0GradY(nBlocks);
  for (unsigned int b = 0; b<nBlocks; b++) {
    t0GradX[b]    = refTimeGradX(iPos[b], jPos[b]);
    zGradX[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradX[b]) - t0GradX[b]); //0.001f is due to ms/s conversion
    t0GradY[b]    = refTimeGradY(iPos[b], jPos[b]);
    zGradY[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradY[b]) - t0GradY[b]);
  }

  bl->getVerticalTrend(&zGradX[0], &az[0]);
  bl->getVerticalTrend(&zGradY[0], &bz[0]);
  bl->getVerticalTrend(&t0GradX[0], &at0[0]);
  bl->getVerticalTrend(&t0GradY[0], &bt0[0]);
}


std::vector<fftw_real>
Wavelet3D::calculateWellWavelet(const std::vector<std::vector<float> > & gMat,
                                const std::vector<float>               & dVec,
                                int                                      nWl,
                                int                                      nhalfWl,
                                int                                      nPoints) const
{
  std::vector<fftw_real> wellWavelet(rnzp_, 0.0);

  double **gTrg = new double *[nWl];
  for (int i=0; i<nWl; i++) {
    gTrg[i]     = new double[nWl];
    for (int j=0; j<nWl; j++) {
      gTrg[i][j] = 0.0;
      for (int k=0; k<nPoints; k++)
        gTrg[i][j] += gMat[k][i]*gMat[k][j];
    }
  }

  double *gTrd = new double[nWl];
  for (int i=0; i<nWl; i++) {
    gTrd[i] = 0.0;
    for (int j=0; j<nPoints; j++)
      gTrd[i] += gMat[j][i] * dVec[j];
  }

  lib_matrCholR(nWl, gTrg);
  lib_matrAxeqbR(nWl, gTrg, gTrd);
  for (int i=0; i<nWl; i++)
    delete [] gTrg[i];
  delete [] gTrg;

  wellWavelet[0] = static_cast<fftw_real> (gTrd[0]);
  for (int i=1; i<nhalfWl; i++) {
    wellWavelet[i]      = static_cast<fftw_real> (gTrd[i]);
    wellWavelet[nzp_-i] = static_cast<fftw_real> (gTrd[nWl-i]);
  }
  for (int i=nhalfWl+1; i<nzp_-nhalfWl; i++)
    wellWavelet[i] = 0.0f;
  for (int i=nzp_; i<rnzp_; i++)
    wellWavelet[i] = RMISSING;
  delete [] gTrd;

  return wellWavelet;
}

float
Wavelet3D::calculateWellWeight(int nWl,
                               int nPoints,
                               const std::vector<std::vector<float> > & gMat,
                               const std::vector<float>               & wellWavelet,
                               const std::vector<float>               & dVec) const
{
  double s2 = 0.0;
  for (int i=0; i<nPoints; i++) {
    double prod = 0.0;
    for (int j=0; j<nWl; j++)
      prod += static_cast<double> (gMat[i][j]*wellWavelet[j]);
    double residual = prod - dVec[i];
    s2 += residual * residual;
  }
  float weight = static_cast<float> (1/s2);
  return weight;
}

void
Wavelet3D::printMatToFile(const std::string                       & fileName, 
                          const std::vector<std::vector<float> >  & mat, 
                          int                                       n,
                          int                                       m) const
{
  if( ModelSettings::getDebugLevel() > 0) { 
    std::string fName = fileName + IO::SuffixGeneralData();
    fName = IO::makeFullFileName(IO::PathToWavelets(), fName);
    std::ofstream file;
    NRLib::OpenWrite(file,fName);
    for(int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++)
        file << mat[i][j] << " ";
      file << "\n";
    }
    file.close();
  }  
}
