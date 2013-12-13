/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <fstream>

#include <string.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "nrlib/flens/nrlib_flens.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/regularsurface.hpp"

#include "src/modelsettings.h"
//#include "src/blockedlogs.h"
#include "src/definitions.h"
#include "src/wavelet3D.h"
#include "src/wavelet1D.h"
#include "src/wavelet.h"
//#include "src/welldata.h"
#include "src/tasklist.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/io.h"
#include "src/waveletfilter.h"

Wavelet3D::Wavelet3D(const std::string                          & filterFile,
                     const std::vector<Surface *>               & estimInterval,
                     const NRLib::Grid2D<float>                 & refTimeGradX,
                     const NRLib::Grid2D<float>                 & refTimeGradY,
                     const std::vector<std::vector<double> >    & tGradX,
                     const std::vector<std::vector<double> >    & tGradY,
                     const SeismicStorage                       * seismic_data,
                     const ModelSettings                        * modelSettings,
                     std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                     const Simbox                               * simBox,
                     const float                                * reflCoef,
                     int                                          angle_index,
                     int                                        & errCode,
                     std::string                                & errText)
  : Wavelet(3),
    filter_(filterFile, errCode, errText)
{
  LogKit::LogFormatted(LogKit::Medium,"  Estimating 3D wavelet pulse from seismic data and (nonfiltered) blocked wells\n");

  //theta_      = seisCube->getTheta();;
  theta_      = seismic_data->GetAngle();
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

  //nzp_        = seisCube->getNzp();
  std::vector<float> tmp_trace = seismic_data->GetTraceData(0);
  nzp_ = tmp_trace.size(); ///H Correct? nzp_ is "size of padded FFT grid in depth (time)" All traces equal in length?

  cnzp_       = nzp_/2+1;
  rnzp_       = 2*cnzp_;

  unsigned int nWells = modelSettings->getNumberOfWells();
  float        v0     = modelSettings->getAverageVelocity();
  int nhalfWl         = static_cast<int> (0.5 * modelSettings->getWaveletTaperingL() / dz_);
  int nWl             = 2 * nhalfWl + 1;
  std::string angle   = NRLib::ToString((180.0/NRLib::Pi)*theta_, 1);

  std::vector<std::vector<fftw_real> > wellWavelets(nWells/*, std::vector<float>(rnzp_, 0.0)*/);
  std::vector<float>                   wellWeight(nWells, 0.0);
  std::vector<float>                   dzWell(nWells, 0.0);

  int w = 0;
  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs.begin(); it != mapped_blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if(blocked_log->GetUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::Medium, "  Well :  %s\n", blocked_log->GetWellName().c_str());

      const std::vector<int> iPos = blocked_log->GetIposVector();
      const std::vector<int> jPos = blocked_log->GetJposVector();

      std::vector<double> az(nz_);
      std::vector<double> bz(nz_);
      std::vector<double> at0(nz_);
      std::vector<double> bt0(nz_);
      unsigned int nBlocks = blocked_log->GetNumberOfBlocks();

      calculateGradients(blocked_log,
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
                         blocked_log,
                         seismic_data,
                         simBox,
                         az,
                         bz,
                         hasWellData);

      int start, length;
      blocked_log->FindContinuousPartOfData(hasWellData,
                                            nz_,
                                            start,
                                            length);

      dzWell[w]          = static_cast<float>(simBox->getRelThick(iPos[0],jPos[0]) * dz_);
      if (length > nWl) {
        std::string wellname(blocked_log->GetWellName());
        NRLib::Substitute(wellname,"/","_");
        NRLib::Substitute(wellname," ","_");
        std::vector<float> Halpha;
        std::vector<fftw_real> cppAdj = adjustCpp(blocked_log,
                                                  az,
                                                  bz,
                                                  Halpha,
                                                  start,
                                                  length,
                                                  wellname,
                                                  angle);

        if( ModelSettings::getDebugLevel() > 0 ) {
          std::string fileName = "xgrad_depth_" + wellname + "_" + angle;

          //std::vector<float> az_tmp(az.size());
          //for(size_t i = 0; i < az_tmp.size(); i++)
          //  az_tmp[i] = static_cast<float>(az[i]);
          //printVecToFile(fileName, &az_tmp[0], length);

          //printVecToFile(fileName, &az[0], length);
          printVecDoubleToFile(fileName, &az[0], length);


          fileName = "ygrad_depth_" + wellname + "_" + angle;
          //printVecToFile(fileName, &bz[0], length);
          printVecDoubleToFile(fileName, &bz[0], length);
          fileName = "cpp_adjust_" + wellname + "_" + angle;
          printVecToFile(fileName, &cppAdj[0], length);
        }

        std::vector<double> z_log(nBlocks);
        for (unsigned int b=0; b<nBlocks; b++) {
          double z_top     = static_cast<double> (simBox->getTop(iPos[b], jPos[b]));
//          zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(iPos[b], jPos[b]) * dz_);
          z_log[b]         = static_cast<double> (z_top + b * dzWell[w]);
        }
        std::vector<double> z_pos_well(nz_);
        blocked_log->GetVerticalTrend(z_log, z_pos_well);

        int nTracesX    = static_cast<int> (modelSettings->getEstRangeX(angle_index) / dx);
        int nTracesY    = static_cast<int> (modelSettings->getEstRangeY(angle_index) / dy);

        std::vector<std::vector<float> > gMat;
        std::vector<float> dVec;
        int nPoints = 0;
        for (int xTr = -nTracesX; xTr <= nTracesX; xTr++) {
          for (int yTr = -nTracesY; yTr <= nTracesY; yTr++) {
            std::vector<double> seis_log(nBlocks);
            blocked_log->GetBlockedGrid(seismic_data, simBox, seis_log, xTr, yTr);
            std::vector<double> seis_data(nz_);
            blocked_log->GetVerticalTrend(seis_log, seis_data);
            for (unsigned int b=0; b<nBlocks; b++) {
              //int xIndex      = iPos[b] + xTr; //NBNB Frode: Gir warning fordi den ikkje blir brukt. Slett om du ikkje skal bruka.
              //int yIndex      = jPos[b] + yTr; //NBNB Frode: Gir warning fordi den ikkje blir brukt. Slett om du ikkje skal bruka.
              float zTop     = static_cast<float> (simBox->getTop(iPos[b], jPos[b]));
//              zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(xIndex, yIndex) * dz_);
              z_log[b]         = static_cast<float> (zTop + b * dzWell[w]);
            }
            std::vector<double> z_pos_trace(nz_);
            blocked_log->GetVerticalTrend(z_log, z_pos_trace);
            for (int t=start; t < start+length; t++) {
              if (seis_data[t] != RMISSING) {
                dVec.push_back(static_cast<float>(seis_data[t]));
                std::vector<std::vector<float> > lambda(length, std::vector<float>(nWl,0.0)); //NBNB-Frode: Bør denne deklareres utenfor?
                for (int tau = start; tau < start+length; tau++) {
                  //Hva gjør vi hvis zData[t] er RMISSING. Kan det skje?
                  double at = at0[tau] + 2.0f*az[tau]/v0;
                  double bt = bt0[tau] + 2.0f*bz[tau]/v0;
                  double time_scale = cos(theta_)/sqrt(1+az[tau]*az[tau]+bz[tau]*bz[tau]);
//                  float u = static_cast<float> (zPosTrace[t] - zPosWell[tau] - at*xTr*dx - bt*yTr*dy);
                  float u = static_cast<float> ((z_pos_well[tau] + at*xTr*dx + bt*yTr*dy - z_pos_trace[t])*time_scale);

                  if ((filter_.hasHalpha() )&& (Halpha[tau-start]!=0)) {
                    float h = Halpha[tau-start];
                    for (int i=0; i<nWl; i++) {
                      int indexPlace = i;
                      if(i > nhalfWl)
                        indexPlace -=nWl;
                      float v = u - static_cast<float>(indexPlace*dzWell[w]*time_scale);
                      lambda[tau-start][i] =  std::min(1.0f,static_cast<float> (2*h*dzWell[w]*1e-3 / (NRLib::Pi*(h*h + 4*v*v*1e-6)))); // invers fouriertransform of exp(-pi*h*Omega*v)
                      // the minimum of 1 and the value isto avoid problem when halpha is very small i.e. 0.0005
                    }
                  }
                  else {
                    int tLow  = static_cast<int> (floor(u / (dzWell[w]*time_scale)));
                    int tHigh = tLow + 1;
                    float lambdaValue = static_cast<float>(u/(dzWell[w]*time_scale)) - static_cast<float> (tLow);
                    if ((tLow >= -nhalfWl) && (tHigh <= nhalfWl)) {
                      if (tLow >= 0)
                        lambda[tau-start][tLow] = 1-lambdaValue;
                      else
                        lambda[tau-start][tLow+nWl] = 1-lambdaValue;
                      if (tHigh < 0)
                        lambda[tau-start][tHigh+nWl] = lambdaValue;
                      else
                        lambda[tau-start][tHigh] = lambdaValue;
                    }
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

        wellWavelets[w] = calculateWellWavelet(gMat,
                                               dVec,
                                               modelSettings->getWavelet3DTuningFactor(),
                                               nWl,
                                               nhalfWl,
                                               nPoints);
        if( ModelSettings::getDebugLevel() > 0 ) {
          std::string fileName = "seismic_" + wellname + "_" + angle;
          printVecToFile(fileName, &dVec[0], nPoints);
          fileName = "gmat_" + wellname + "_" + angle;
          printMatToFile(fileName, gMat, nPoints, nWl);
        }
        wellWeight[w] = calculateWellWeight(nWl,
                                            nPoints,
                                            gMat,
                                            wellWavelets[w],
                                            dVec);
      } //if (length > nWl)
      else {
        LogKit::LogFormatted(LogKit::Medium,"     No enough data for 3D wavelet estimation in well %s\n", blocked_log->GetWellName().c_str());
      }
    } // if(wells->getUseForEstimation)
    w++;
  } // for (w=0...nWells)

  rAmp_ = averageWavelets(wellWavelets, nWells, nzp_, wellWeight, dzWell, dz_);
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);
  waveletLength_ = findWaveletLength(modelSettings->getMinRelWaveletAmp(),modelSettings->getWaveletTaperingL());
  LogKit::LogFormatted(LogKit::Low,"  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( ModelSettings::getDebugLevel() > 0 )
    writeWaveletToFile("estimated_wavelet_", 1.0f,false);

  norm_ = findNorm();

  fftw_real * trueAmp = rAmp_;
  rAmp_               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
  cAmp_               = reinterpret_cast<fftw_complex *>(rAmp_);

  w = 0;
  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs.begin(); it != mapped_blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if(blocked_log->GetUseForWaveletEstimation() &&
       ((modelSettings->getWaveletOutputFlag() & IO::WELL_WAVELETS)>0 || modelSettings->getEstimationMode()))
    {
      for(int i=0; i<nzp_; i++)
        rAmp_[i] = wellWavelets[w][i];
      std::string wellname(blocked_log->GetWellName());
      NRLib::Substitute(wellname,"/","_");
      NRLib::Substitute(wellname," ","_");
      std::string fileName = IO::PrefixWellWavelet() + wellname + "_";
      writeWaveletToFile(fileName, 1.0f,false);
    }
    w++;
  }
  fftw_free(rAmp_);
  rAmp_ = trueAmp;
  cAmp_ = reinterpret_cast<fftw_complex *>(rAmp_);
  averageWavelet_ = createAverageWavelet(simBox);
}


Wavelet3D::Wavelet3D(const std::string & fileName,
            int                 fileFormat,
            const ModelSettings * modelSettings,
            const float         * reflCoef,
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
  //The delete is needed to avoid memory leak, but it causes a segmentation fault when doing the
  // delete modelAVOdynamic in main.cpp
  //delete averageWavelet_;
}

float
Wavelet3D::getLocalStretch(int i,  int j) // Note: Not robust towards padding
{
  float gx=0.0f;
  float gy=0.0f;

  if(structureDepthGradX_.GetN()>0)
  {
    gx        = GetLocalDepthGradientX(i,j);
    gy        = GetLocalDepthGradientY(i,j);
  }

  float r     = sqrt(gx*gx + gy*gy + 1);
  double psi  = findPsi(r);

  float  sf   = float(cos(psi)*cos( getTheta() )); // First is local dip adjustment second is for anglestack
  return sf;
}

Wavelet1D*
Wavelet3D::createLocalWavelet1D(int i, int j)// note Not robust towards padding
{
  // returns the local wavelet exccept from the strech factor
  // The stretch factor is given by  getLocalStretch
  // this is done to avoid double resampling of wavelet.

  float gx=0.0f;
  float gy=0.0f;

  Wavelet1D* localWavelet;
  if(structureDepthGradX_.GetN()>0)
  {
    gx        = GetLocalDepthGradientX(i,j);
    gy        = GetLocalDepthGradientY(i,j);
  }

  // NBNB transform to local vs not local xy coordinates ??
  double phi       = findPhi(gx, gy);
  float r         = sqrt(gx*gx + gy*gy + 1);
  double psi       = findPsi(r);


  localWavelet = extractLocalWaveletByDip1D(phi,psi);// makes a new wavelet

  doLocalShiftAndScale1D(localWavelet,i,j);

  return localWavelet;
}

Wavelet1D*
Wavelet3D::createAverageWavelet(const Simbox * simBox)
{
  Wavelet1D*  w1;
  fftw_complex* average1;
  fftw_complex* average2;
  average1 = static_cast<fftw_complex*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)));
  average2 = static_cast<fftw_complex*>(fftw_malloc(2*(nzp_/2+1)*sizeof(fftw_real)));
  int k;
  for(k=0;k < (nzp_/2 +1);k++)
  {
    average1[k].re=0;
    average1[k].im=0;
  }
  int nx =simBox->getnx();
  int ny =simBox->getny();
  double divNx=static_cast<double>(1.0/nx);
  double divNy=static_cast<double>(1.0/ny);

  for(int i=0;i<nx;i++){
    for(k=0;k < (nzp_/2 +1);k++)
    {
      average2[k].re=0;
      average2[k].im=0;
    }
    for(int j=0;j<ny;j++)
    {
      w1= createLocalWavelet1D( i, j);
      w1->fft1DInPlace();
      double sfLoc =(simBox->getRelThick(i,j)*w1->getLocalStretch(i,j));// scale factor from thickness stretch + (local stretch when 3D wavelet)
      //double relT   = simBox->getRelThick(i,j);
      //double deltaF = static_cast<double>(nz_)*1000.0/(relT*simBox->getlz()*static_cast<double>(nzp_));
      for(int k=0;k < (nzp_/2 +1);k++)
      {
        fftw_complex amp = w1->getCAmp(k,static_cast<float>(sfLoc));
        average2[k].re+= static_cast<fftw_real>(amp.re*divNy);
        average2[k].im+= static_cast<fftw_real>(amp.im*divNy);
      }
      delete w1;
    }
    for(k=0;k < (nzp_/2 +1);k++)
    {
      average1[k].re+= static_cast<fftw_real>(average2[k].re*divNx);
      average1[k].im+= static_cast<fftw_real>(average2[k].im*divNx);
    }
  }
  w1= createLocalWavelet1D( 0, 0);
  w1->fft1DInPlace();
  for(k=0;k < (nzp_/2 +1);k++)
  {
    w1->setCAmp(average1[k],k);
  }

  fftw_free(average1);
  fftw_free(average2);
  return w1;
}


Wavelet1D *
Wavelet3D::extractLocalWaveletByDip1D(double phi,double psi)
{
  Wavelet1D* wavelet = createSourceWavelet();
  dipAdjustWavelet(wavelet,phi,psi);

  return wavelet;
}

Wavelet1D *
Wavelet3D::createSourceWavelet()
{
    Wavelet1D * sourceWavelet = new Wavelet1D(this);
    return sourceWavelet;
}

void
Wavelet3D::dipAdjustWavelet(Wavelet1D* wavelet, double phi, double psi)
{
  // NBNB open cos(psi)*cos(getTheta()) is  included
  double multiplyer = filter_.getAlpha1(phi, psi)*(cos(psi)*cos(getTheta()));

  double Halpha = filter_.getHalpha(phi, psi);
  wavelet->adjustForAmplitudeEffect(multiplyer,Halpha);
}



Wavelet1D*
Wavelet3D::createWavelet1DForErrorNorm(void)
{
  Wavelet1D* errorWavelet;
  errorWavelet = createLocalWavelet1D( 0, 0);
  errorWavelet->fft1DInPlace();
  for(int k=0;k < (nzp_/2 +1);k++)
    errorWavelet->setCAmp(averageWavelet_->getCAmp(k),k);
  errorWavelet->invFFT1DInPlace();
  errorWavelet->findNorm();
  // NBNB OK gjør om denne skaleringen til RMS gain.
  return errorWavelet;
}

float
Wavelet3D::calculateSNRatio(const Simbox                                     * simbox,
                            const SeismicStorage                             * seismic_data,
                            const std::map<std::string, BlockedLogsCommon *> & mapped_blocked_logs,
                            const ModelSettings                              * modelSettings,
                            std::string                                      & errText,
                            int                                              & error,
                            const NRLib::Grid2D<float>                       & refTimeGradX,
                            const NRLib::Grid2D<float>                       & refTimeGradY,
                            const std::vector<std::vector<double> >          & tGradX,
                            const std::vector<std::vector<double> >          & tGradY,
                            int                                                number,
                            float                                              SNRatio,
                            bool                                               estimateSNRatio,
                            bool                                               estimateWavelet)
{
  std::string angle    = NRLib::ToString((180.0/NRLib::Pi)*theta_, 1);
  unsigned int nWells  = modelSettings->getNumberOfWells();
  int          nhalfWl = static_cast<int> (0.5 * modelSettings->getWaveletTaperingL() / dz_);
  int          nWl     = 2 * nhalfWl + 1;
  float        v0      = modelSettings->getAverageVelocity();
  float        dataVar = 0.0f;
  float        errVar  = 0.0f;
  int          nData   = 0;
  std::vector<float> shiftWell  (nWells, 0.0f);
  std::vector<float> dzWell(nWells, 0.0);
  std::vector<float> dataVarWell(nWells, 0.0f);
  std::vector<float> errVarWell (nWells, 0.0f);
  std::vector<int>   nActiveData(nWells, 0);

  unsigned int w = 0;
  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs.begin(); it != mapped_blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if(blocked_log->GetUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::Medium, "  Well :  %s\n", blocked_log->GetWellName().c_str());

      const std::vector<int> iPos = blocked_log->GetIposVector();
      const std::vector<int> jPos = blocked_log->GetJposVector();

      std::vector<double> az(nz_);
      std::vector<double> bz(nz_);
      std::vector<double> at0(nz_);
      std::vector<double> bt0(nz_);
      unsigned int nBlocks = blocked_log->GetNumberOfBlocks();
      calculateGradients(blocked_log,
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
      std::vector<double> seis_log(nBlocks);
      blocked_log->GetBlockedGrid(seismic_data, simbox, seis_log);
      std::vector<double> seis_data(nz_);
      blocked_log->GetVerticalTrend(seis_log, seis_data);

      std::vector<double> vp(nz_);
      std::vector<double> vs(nz_);
      std::vector<double> rho(nz_);
      blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), vp);
      blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), vs);
      blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), rho);

      for (int k=0; k<nz_; k++)
        hasWellData[k] = (vp[k] != RMISSING && vs[k] != RMISSING && rho[k] != RMISSING && az[k] != RMISSING && bz[k] != RMISSING && seis_data[k] != RMISSING);

      int start, length;

      blocked_log->FindContinuousPartOfData(hasWellData,
                                            nz_,
                                            start,
                                            length);


      dzWell[w]          = static_cast<float>(simbox->getRelThick(iPos[0],jPos[0]) * dz_);
      if (length > nWl) {
        std::string wellname(blocked_log->GetWellName());
        NRLib::Substitute(wellname,"/","_");
        NRLib::Substitute(wellname," ","_");
        std::vector<float> Halpha;
        std::vector<fftw_real> cppAdj = adjustCpp(blocked_log,
                                                  az,
                                                  bz,
                                                  Halpha,
                                                  start,
                                                  length,
                                                  wellname,
                                                  angle);

        std::vector<double> zLog(nBlocks);
        for (unsigned int b=0; b<nBlocks; b++) {
          float zTop     = static_cast<float> (simbox->getTop(iPos[b], jPos[b]));
          //          zLog[b]         = static_cast<float> (zTop + b * simBox->getRelThick(iPos[b], jPos[b]) * dz_);
          zLog[b]         = static_cast<double> (zTop + b * dzWell[w]);
        }
        std::vector<double> zPosWell(nz_);
        blocked_log->GetVerticalTrend(zLog, zPosWell);

        std::vector<std::vector<float> > gMat;
        std::vector<float> dVec;

        for (int t=start; t < start+length; t++) {
          if (seis_data[t] != RMISSING) {
            dVec.push_back(static_cast<float>(seis_data[t]));
            std::vector<std::vector<float> > lambda(length, std::vector<float>(nWl,0.0)); //NBNB-Frode: Bør denne deklareres utenfor?
            for (int tau = start; tau < start+length; tau++) {
              double time_scale = cos(theta_)/sqrt(1+az[tau]*az[tau]+bz[tau]*bz[tau]);
              float u = static_cast<float> ((zPosWell[tau] - zPosWell[t])*time_scale);
              if ((filter_.hasHalpha() )&& (Halpha[tau-start]!=0)) {
                float h = Halpha[tau-start];
                for (int i=0; i<nWl; i++) {
                  int indexPlace = i;
                  if(i > nhalfWl)
                    indexPlace -=nWl;
                  float v = u - static_cast<float>(indexPlace*dzWell[w]*time_scale);
                  lambda[tau-start][i] = std::min(1.0f,static_cast<float> (2*h*dzWell[w]*1e-3 / (NRLib::Pi *(h*h + 4*v*v*1e-6)))); // ?
                }
              }
              else {
                int tLow  = static_cast<int> (floor(u / (dzWell[w]*time_scale)));
                int tHigh = tLow + 1;
                float lambdaValue = static_cast<float>(u/(dzWell[w]*time_scale)) - static_cast<float> (tLow);

                if ((tLow >= -nhalfWl) && (tHigh <= nhalfWl)) {
                  if (tLow >= 0)
                    lambda[tau-start][tLow] = 1-lambdaValue;
                  else
                    lambda[tau-start][tLow+nWl] = 1-lambdaValue;
                  if (tHigh < 0)
                    lambda[tau-start][tHigh+nWl] = lambdaValue;
                  else
                    lambda[tau-start][tHigh] = lambdaValue;
                }
              } // else
            } // for (tau=start...start+length)
            std::vector<float> gVec(nWl);
            for (int i=0; i<nWl; i++) {
              gVec[i] = 0.0f;
              for (int j=0; j<length; j++)
                gVec[i] += cppAdj[j] * lambda[j][i];
            }
            gMat.push_back(gVec);
          } //if (seisData[t] != RMISSING)
        }

        nActiveData[w] = length;
        fftw_real* syntSeisExt               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
        fftw_complex * syntSeisExt_c               = reinterpret_cast<fftw_complex *>(syntSeisExt);
        fftw_real* dataExt               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
        fftw_complex * dataExt_c               = reinterpret_cast<fftw_complex *>(dataExt);
        fftw_real* crossCorr               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
        fftw_complex * crossCorr_c               = reinterpret_cast<fftw_complex *>(crossCorr);


        std::vector<fftw_real> full_wavelet(rnzp_);
        fillInnWavelet(&full_wavelet[0], nzp_, dzWell[w]);


        std::vector<float> wavelet(nWl);
        wavelet[0] = static_cast<float> (full_wavelet[0]);
        for (int i=1; i<nhalfWl; i++) {
          wavelet[i]     = static_cast<float> (full_wavelet[i]);
          wavelet[nWl-i] = static_cast<float> (full_wavelet[nzp_-i]);
        }

        for (int i=0; i<length; i++) {
          syntSeisExt[i] = 0.0f;
          dataExt[i]=dVec[i];
          for (int j=0; j<nWl; j++)
            syntSeisExt[i] += gMat[i][j] * wavelet[j];
        }


        for (int i=length; i<nzp_; i++) {
          syntSeisExt[i] = 0.0f;
          dataExt[i]= 0.0f;
        }

        Utils::fft(syntSeisExt, syntSeisExt_c, nzp_);
        Utils::fft(dataExt, dataExt_c, nzp_);


        convolve(syntSeisExt_c, dataExt_c, crossCorr_c, cnzp_);
        Utils::fftInv(crossCorr_c, crossCorr, nzp_); //

        float shift = findBulkShift(crossCorr, dzWell[w], nzp_,  modelSettings->getMaxWaveletShift());


        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w] = shift;
        shiftReal(-shift/dzWell[w], crossCorr, nzp_);// NBNB just checking

        Utils::fftInv(syntSeisExt_c, syntSeisExt, nzp_);
        shiftReal(-shift/dzWell[w], syntSeisExt, nzp_);

        Utils::fftInv(dataExt_c, dataExt, nzp_);

        for (int i=0; i<length; i++) {
          float residual = dataExt[i] - syntSeisExt[i];
          errVarWell[w]  += residual * residual;
          dataVarWell[w] += dVec[i] * dVec[i];
        }
        errVar  += errVarWell[w];
        dataVar += dataVarWell[w];
        nData   += nActiveData[w];
        dataVarWell[w] /= static_cast<float>(nActiveData[w]);
        errVarWell[w]  /= static_cast<float>(nActiveData[w]);

        if(ModelSettings::getDebugLevel() > 0) {
          std::string fileName;
          //fileName = "seismic_" + wellname + "_" + angle;
          fileName = "seismic_Well_" + NRLib::ToString(w+1) + "_" + angle;
          printVecToFile(fileName, &dVec[0], length);
          //fileName = "synthetic_seismic_" + wellname + "_" + angle;
          fileName = "synthetic_seismic_Well_" + NRLib::ToString(w+1) + "_" + angle;
          printVecToFile(fileName, syntSeisExt, length);
        }
      }
      else {
        LogKit::LogFormatted(LogKit::Low, "\n  Not using vertical well %s for error estimation (length=%.1fms  required length=%.1fms).",
          blocked_log->GetWellName().c_str(), length*dz_, waveletLength_);
      }
    }
    w++;
  }

  dataVar /= nData;
  errVar  /= nData;
  float empSNRatio = dataVar/errVar;

  LogKit::LogFormatted(LogKit::Medium,"\n  Reporting errors (as standard deviations) estimated:\n\n");

  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"                                     SeisData            \n");
  LogKit::LogFormatted(LogKit::Low,"  Well                  shift[ms]     StdDev         S/N \n");
  LogKit::LogFormatted(LogKit::Low,"  -------------------------------------------------------\n");

  w = 0;
  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = mapped_blocked_logs.begin(); it != mapped_blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = mapped_blocked_logs.find(it->first);
    const BlockedLogsCommon * blocked_log = iter->second;

    if(nActiveData[w]>0) {
      float SNWell = dataVarWell[w]/errVarWell[w];
      LogKit::LogFormatted(LogKit::Low,"  %-20s   %6.2f     %9.2e      %6.2f \n",
            blocked_log->GetWellName().c_str(),shiftWell[w],sqrt(dataVarWell[w]),SNWell);
    }
    else
      LogKit::LogFormatted(LogKit::Low,"  %-20s      -            -             - \n",blocked_log->GetWellName().c_str());

    w++;
  }

  if(estimateSNRatio)
    LogKit::LogFormatted(LogKit::Low,"\n  The signal to noise ratio used for this angle stack is: %6.2f\n", empSNRatio);
  else {
    LogKit::LogFormatted(LogKit::Low,"\n  The signal to noise ratio given in the model file and used for this angle stack is : %6.2f\n", SNRatio);
    LogKit::LogFormatted(LogKit::Low,"  For comparison, the signal-to-noise ratio calculated from the available wells is   : %6.2f\n", empSNRatio);
    float minSN = 1.0f + (empSNRatio - 1.0f)/2.0f;
    float maxSN = 1.0f + (empSNRatio - 1.0f)*2.0f;
    if ((SNRatio<minSN || SNRatio>maxSN) && estimateWavelet) {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The difference between the SN ratio given in the model file and the calculated SN ratio is too large.\n");
      if (SNRatio < minSN)
        TaskList::addTask("Consider increasing the SN ratio for angle stack "+NRLib::ToString(number)+" to minimum "+NRLib::ToString(minSN,1)+"\n");
      else
        TaskList::addTask("Consider decreasing the SN ratio for angle stack "+NRLib::ToString(number)+" to maximum "+NRLib::ToString(minSN,1)+"\n");
    }
  }

  if (empSNRatio < 1.1f) {
    if (estimateSNRatio) {
      LogKit::LogFormatted(LogKit::Warning,"\nERROR: The empirical signal-to-noise ratio Var(data)/Var(noise) is %.2f. Ratios smaller",empSNRatio);
      LogKit::LogFormatted(LogKit::Warning,"\n       than 1.1 are not acceptable. The signal-to-noise ratio was not reliably estimated");
      LogKit::LogFormatted(LogKit::Warning,"\n       and you must give it as input in the model file.\n");
      LogKit::LogFormatted(LogKit::Warning,"\n       If the wavelet was estimated by CRAVA the solution may be to remove one or more wells");
      LogKit::LogFormatted(LogKit::Warning,"\n       from the wavelet estimation (compare shifts and SN-ratios for different wells).\n");

      errText += "Invalid signal-to-noise ratio obtained for the angle-gather of "+angle+" degrees.\n";
      error += 1;
    }
    else {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The empirical signal-to-noise ratio Var(data)/Var(noise) is %.2f. Ratios smaller",empSNRatio);
      LogKit::LogFormatted(LogKit::Warning,"\n         than 1.1 are not acceptable. If the low ratio is caused by one or more specific");
      LogKit::LogFormatted(LogKit::Warning,"\n         wells, these wells should not have been included in the wavelet estimation. If");
      LogKit::LogFormatted(LogKit::Warning,"\n         they indeed were omitted, you can avoid this warning by specifying");
      LogKit::LogFormatted(LogKit::Warning,"\n         <use-for-wavelet-estimation> no </use-...> for these wells in the model file.\n");
      TaskList::addTask("Check the signal-to-noise ratio given for angle stack "+NRLib::ToString(number)+" against that calculated by CRAVA.");
    }
  }

  return empSNRatio;
}

void
Wavelet3D::findLayersWithData(const std::vector<Surface *> & estimInterval,
                              BlockedLogsCommon            * blocked_log,
                              const SeismicStorage         * seismic_data,
                              const Simbox                 * simBox,
                              const std::vector<double>    & az,
                              const std::vector<double>    & bz,
                              std::vector<bool>            & hasWellData) const
{
  std::vector<double> seisLog(blocked_log->GetNumberOfBlocks());
  blocked_log->GetBlockedGrid(seismic_data, simBox, seisLog);
  std::vector<double> seis_data(nz_);
  blocked_log->GetVerticalTrend(seisLog, seis_data);

  std::vector<double> alpha(nz_);
  std::vector<double> beta(nz_);
  std::vector<double> rho(nz_);
  blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), alpha);
  blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), beta);
  blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), rho);

  for (int k=0; k<nz_; k++)
    hasWellData[k] = (alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING && az[k] != RMISSING && bz[k] != RMISSING && seis_data[k] != RMISSING);

  //Check that data are within wavelet estimation interval
  if (estimInterval.size() > 0) {
    const std::vector<double> xPos = blocked_log->GetXpos();
    const std::vector<double> yPos = blocked_log->GetYpos();
    const std::vector<double> zPos = blocked_log->GetZpos();
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
  double epsilon = 0.001;
  if (a > epsilon && b >= 0.0) //1. quadrant
    phi = atan(b/a);
  else if (a > epsilon && b < 0.0) //4. quadrant
    phi = 2*NRLib::Pi + atan(b/a);
  else if (a < - epsilon && b >= 0.0) //2. quadrant
    phi = NRLib::Pi + atan(b/a);
  else if (a < - epsilon && b < 0.0) //3. quadrant
    phi = NRLib::Pi + atan(b/a);
  else if (b  >= 0.0) //a very small
    phi = 0.5 * NRLib::Pi;
  else //a very small
    phi = 1.5 * NRLib::Pi;
  phi = phi;
  return(phi);
}

double
Wavelet3D::findPsi(float r) const
//Return value should be between 0 and PI/2
{
   double psi;
  if(r > 1)
  {
    psi  = acos(1.0 / r);
    psi = psi  ;
  }else
    psi=0;

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
Wavelet3D::adjustCpp(BlockedLogsCommon         * blocked_log,
                     const std::vector<double> & az,
                     const std::vector<double> & bz,
                     std::vector<float>        & Halpha,
                     int                         start,
                     int                         length,
                     const std::string         & wellname,
                     const std::string         & angle) const
{
  std::vector<fftw_real> cpp(rnzp_);
  blocked_log->FillInCpp(coeff_, start, length, &cpp[0], nzp_);

  std::vector<fftw_real> cppAdj(length, 0.0);
  std::vector<float>     alpha1_vec(length, 0.0);
  std::vector<float>     phi_vec(length, 0.0);
  std::vector<float>     psi_vec(length, 0.0);

  if (filter_.hasHalpha())
    Halpha.resize(length,0.0);

  for (int i=start; i < start+length-1; i++) {
    double phi       = findPhi(static_cast<float>(az[i]),static_cast<float>(bz[i]));
    phi_vec[i-start] = static_cast<float> (phi);
    float r          = sqrt(static_cast<float>(az[i]*az[i] + bz[i]*bz[i] + 1));
    double psi       = findPsi(r);
    psi_vec[i-start] = static_cast<float> (psi);
    float alpha1     = filter_.getAlpha1(phi, psi);
    alpha1_vec[i-start] = alpha1;
    float stretch = std::cos(theta_);
    cppAdj[i-start]     = static_cast<fftw_real> (cpp[i] * alpha1 * stretch / r);
    if (filter_.hasHalpha())
      Halpha[i-start]   = filter_.getHalpha(phi, psi);
  }

  if(ModelSettings::getDebugLevel() > 0) {
    std::string fileName = "cpp_" + wellname + "_" + angle;
    printVecToFile(fileName, &cpp[0], length);
    fileName = "phi_" + wellname + "_" + angle;
    printVecToFile(fileName, &phi_vec[0], length);
    fileName = "psi_" + wellname + "_" + angle;
    printVecToFile(fileName, &psi_vec[0], length);
    fileName = "alpha_" + wellname + "_" + angle;
    printVecToFile(fileName, &alpha1_vec[0], length);
  }

  return cppAdj;
}


void
Wavelet3D::calculateGradients(BlockedLogsCommon          * blocked_log,
                              const std::vector<int>     & iPos,
                              const std::vector<int>     & jPos,
                              const NRLib::Grid2D<float> & refTimeGradX,
                              const NRLib::Grid2D<float> & refTimeGradY,
                              const std::vector<double>  & tGradX,
                              const std::vector<double>  & tGradY,
                              float                        v0,
                              std::vector<double>        & az,
                              std::vector<double>        & bz,
                              std::vector<double>        & at0,
                              std::vector<double>        & bt0) const
{
  unsigned int nBlocks = blocked_log->GetNumberOfBlocks();
  std::vector<double> zGradX(nBlocks);
  std::vector<double> zGradY(nBlocks);
  std::vector<double> t0GradX(nBlocks);
  std::vector<double> t0GradY(nBlocks);
  for (unsigned int b = 0; b<nBlocks; b++) {
    t0GradX[b]    = refTimeGradX(iPos[b], jPos[b]);
    zGradX[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradX[b]) - t0GradX[b]); //0.001f is due to ms/s conversion
    t0GradY[b]    = refTimeGradY(iPos[b], jPos[b]);
    zGradY[b]     = 0.5f * v0 * 0.001f * (static_cast<float> (tGradY[b]) - t0GradY[b]);
  }

  blocked_log->GetVerticalTrend(zGradX, az);
  blocked_log->GetVerticalTrend(zGradY, bz);
  blocked_log->GetVerticalTrend(t0GradX, at0);
  blocked_log->GetVerticalTrend(t0GradY, bt0);
}


std::vector<fftw_real>
Wavelet3D::calculateWellWavelet(const std::vector<std::vector<float> > & gMat,
                                const std::vector<float>               & dVec,
                                double                                   SNR,
                                int                                      nWl,
                                int                                      nhalfWl,
                                int                                      nPoints) const
{
  std::vector<fftw_real> wellWavelet(rnzp_, 0.0);

  double maxG = 0.0;
  for (int i=0; i<nPoints; i++) {
    for (int j=0; j<nWl; j++)
      maxG = std::max(maxG, static_cast<double> (fabs(gMat[i][j])));
  }

  //double wScale = maxD/maxG;

  NRLib::SymmetricMatrix gTrg = NRLib::SymmetricZeroMatrix(nWl);

  for (int i=0; i<nWl; i++) {
    for (int j=0; j<=i; j++) {
      for (int k=0; k<nPoints; k++)
        gTrg(j,i) += gMat[k][i]*gMat[k][j];
    }
  }

  NRLib::Vector gTrd = NRLib::ZeroVector(nWl);
  NRLib::Vector x(nWl);

  double alpha = 4.0;
  double beta  = 0.5;

  for (int i=0; i<nWl; i++) {
    double a;
    if (i <= nhalfWl)
      a = static_cast<double> (i)/ static_cast<double> (nhalfWl);
    else
      a = static_cast<double> (nWl-i)/ static_cast<double> (nhalfWl);
    gTrg(i,i) += maxG * maxG * exp(alpha*a*a) / (SNR * beta * beta); //( GTG+ (maxG/maxD)^2*(Noise)^2*priorCov^-1)=( GTG+ maxG^2/SN*(priorCov^-1)

    for (int j=0; j<nPoints; j++)
      gTrd(i) += gMat[j][i] * dVec[j];
  }

  NRLib::CholeskySolve(gTrg, gTrd, x);

  wellWavelet[0] = static_cast<fftw_real> (x(0));
  for (int i=1; i<nhalfWl; i++) {
    wellWavelet[i]      = static_cast<fftw_real> (x(i));
    wellWavelet[nzp_-i] = static_cast<fftw_real> (x(nWl-i));
  }
  for (int i=nhalfWl+1; i<nzp_-nhalfWl; i++)
    wellWavelet[i] = 0.0f;
  for (int i=nzp_; i<rnzp_; i++)
    wellWavelet[i] = RMISSING;

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


NRLib::Grid2D<float> Wavelet3D::structureDepthGradX_    = NRLib::Grid2D<float>(0,0);
NRLib::Grid2D<float> Wavelet3D::structureDepthGradY_    = NRLib::Grid2D<float>(0,0);
