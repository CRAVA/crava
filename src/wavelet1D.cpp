#include <iostream>
#include <fstream>

#include <string.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

#include "lib/lib_matr.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/definitions.h"
#include "src/wavelet1D.h"
#include "src/welldata.h"
#include "src/tasklist.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/model.h"
#include "src/covgrid2d.h"
#include "src/kriging2d.h"
#include "src/krigingdata2d.h"
#include "src/io.h"

Wavelet1D::Wavelet1D(Simbox                       * simbox,
                     FFTGrid                      * seisCube,
                     WellData                    ** wells,
                     const std::vector<Surface *> & estimInterval,
                     ModelSettings                * modelSettings,
                     float                        * reflCoef,
                     int                            iAngle)
  : Wavelet(1)
{
  LogKit::LogFormatted(LogKit::Medium,"  Estimating 1D wavelet from seismic data and (nonfiltered) blocked wells\n");

  coeff_[0] = reflCoef[0];
  coeff_[1] = reflCoef[1];
  coeff_[2] = reflCoef[2];

  dz_                   = static_cast<float>(simbox->getdz());
  nz_                   = simbox->getnz();
  theta_                = seisCube->getTheta();
  nzp_                  = seisCube->getNzp();
  cnzp_                 = nzp_/2+1;
  rnzp_                 = 2*cnzp_;
  scale_                = 1.0f; 
  cz_                   = 0;
  inFFTorder_           = true;
  isReal_               = true; 
  formats_              = modelSettings->getWaveletFormatFlag();

  std::string fileName;
  int     nWells              = modelSettings->getNumberOfWells();
  float   waveletTaperLength  = modelSettings->getWaveletTaperingL();
 
  //Wavelet estimation
  fftw_real    ** cpp_r = new fftw_real*[nWells];
  fftw_complex ** cpp_c = reinterpret_cast<fftw_complex**>(cpp_r);
  
  fftw_real    ** seis_r = new fftw_real*[nWells];
  fftw_complex ** seis_c = reinterpret_cast<fftw_complex**>(seis_r);

  fftw_real    ** synt_seis_r = new fftw_real*[nWells];
  fftw_complex ** synt_seis_c = reinterpret_cast<fftw_complex**>(synt_seis_r);

  fftw_real    ** cor_cpp_r = new fftw_real*[nWells];
  fftw_complex ** cor_cpp_c = reinterpret_cast<fftw_complex**>(cor_cpp_r);
  
  fftw_real    ** ccor_seis_cpp_r = new fftw_real*[nWells];
  fftw_complex ** ccor_seis_cpp_c = reinterpret_cast<fftw_complex**>(ccor_seis_cpp_r);

  fftw_real    ** wavelet_r = new fftw_real*[nWells];
  fftw_complex ** wavelet_c = reinterpret_cast<fftw_complex**>(wavelet_r); 

  std::vector<float>  dzWell(nWells);
  for(int i=0;i<nWells;i++) {
    cpp_r[i]           = new fftw_real[rnzp_];
    synt_seis_r[i]     = new fftw_real[rnzp_];
    for(int j=0;j<rnzp_;j++) {
      cpp_r[i][j] = 0;
      synt_seis_r[i][j] = 0;
    }
    seis_r[i]          = new fftw_real[rnzp_];
    cor_cpp_r[i]       = new fftw_real[rnzp_];
    ccor_seis_cpp_r[i] = new fftw_real[rnzp_];
    wavelet_r[i]       = new fftw_real[rnzp_];

    const std::vector<int> ipos = wells[i]->getBlockedLogsOrigThick()->getIposVector();
    const std::vector<int> jpos = wells[i]->getBlockedLogsOrigThick()->getJposVector();
    dzWell[i]          = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])) * dz_;
  }

  std::vector<float>  z0(nWells);               // Needed to block syntSeis
  std::vector<int>    sampleStart(nWells,0);    // Needed to block syntSeis
  std::vector<int>    sampleStop(nWells,0);     // Needed to block syntSeis
  std::vector<float>  wellWeight(nWells,0.0f);  
  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::Medium,"  Well :  %s\n",wells[w]->getWellname().c_str());

      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
      std::vector<float> seisLog(bl->getNumberOfBlocks());
      bl->getBlockedGrid(seisCube,&seisLog[0]);
      //
      // Check seismic data outside estimation interval missing
      //
      if (estimInterval.size() > 0) {
        const std::vector<double> xPos = bl->getXposVector();
        const std::vector<double> yPos = bl->getYposVector();
        const std::vector<double> zPos = bl->getZposVector();
        for (int k = 0 ; k < bl->getNumberOfBlocks() ; k++) {
          const double zTop  = estimInterval[0]->GetZ(xPos[k],yPos[k]);
          const double zBase = estimInterval[1]->GetZ(xPos[k],yPos[k]);
          if ( (zPos[k] - 0.5*dz_) < zTop || (zPos[k] + 0.5*dz_) > zBase)
            seisLog[k] = RMISSING;
        }
      }
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      std::vector<float> alpha(nz_);
      bl->getVerticalTrend(bl->getAlpha(), &alpha[0]);
      std::vector<float> beta(nz_);
      bl->getVerticalTrend(bl->getBeta(), &beta[0]);
      std::vector<float> rho(nz_);
      bl->getVerticalTrend(bl->getRho(), &rho[0]);
      std::vector<float> seisData(nz_);
      bl->getVerticalTrend(&seisLog[0], &seisData[0]);
      std::vector<bool> hasData(nz_);
      for (int k = 0 ; k < nz_ ; k++)
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      //
      // Find continuous part of data
      //
      int start,length;
      bl->findContiniousPartOfData(hasData, nz_, start, length);
      if(length*dz_ > waveletTaperLength ) { // must have enough data
        bl->fillInCpp(coeff_, start, length, cpp_r[w], nzp_); 
        fileName = "cpp_1";
        printVecToFile(fileName, cpp_r[w], nzp_);  // Debug
        Utils::fft(cpp_r[w], cpp_c[w], nzp_);
        bl->fillInSeismic(&seisData[0], start, length, seis_r[w], nzp_);
        fileName = "seis_1";
        printVecToFile(fileName, seis_r[w], nzp_); // Debug
        Utils::fft(seis_r[w], seis_c[w], nzp_);
        bl->estimateCor(cpp_c[w], cpp_c[w], cor_cpp_c[w], cnzp_);
        Utils::fftInv(cor_cpp_c[w], cor_cpp_r[w], nzp_);
        bl->estimateCor(cpp_c[w], seis_c[w], ccor_seis_cpp_c[w], cnzp_);
        Utils::fftInv(ccor_seis_cpp_c[w], ccor_seis_cpp_r[w], nzp_);
        Utils::fftInv(cpp_c[w], cpp_r[w], nzp_);
        Utils::fftInv(seis_c[w], seis_r[w], nzp_);
        wellWeight[w] = length*dzWell[w]*(cor_cpp_r[w][0]+cor_cpp_r[w][1]);// Gives most weight to long datasets with  
                                                                       // large reflection coefficients
        z0[w] = static_cast<float> (bl->getZpos()[0]);
        sampleStart[w] = start;
        sampleStop[w]  = start + length;
      }
    }
  }

  std::vector<float> shiftWell(nWells);
  float shiftAvg = shiftOptimal(ccor_seis_cpp_r, wellWeight, dzWell, nWells, nzp_, shiftWell, modelSettings->getMaxWaveletShift());
  multiplyPapolouis(ccor_seis_cpp_r, dzWell, nWells, nzp_, waveletTaperLength, wellWeight);
  multiplyPapolouis(cor_cpp_r, dzWell, nWells, nzp_, waveletTaperLength, wellWeight);
  getWavelet(ccor_seis_cpp_r, cor_cpp_r, wavelet_r, wellWeight, nWells, nzp_);

  // Save estimated wavelet for each well, write to file later because waveletlength calculted for average wavelet is used for all well wavelets
  // to simplify comparison 
  std::vector<std::vector<fftw_real> > wellWavelets(nWells);
  for(int w=0; w<nWells; w++) {
    wellWavelets[w].resize(nzp_);
    for(int i=0;i<nzp_;i++)
      wellWavelets[w][i] = wavelet_r[w][i];
  }

  rAmp_ = averageWavelets(wellWavelets, nWells, nzp_, wellWeight, dzWell, dz_); // wavelet centered
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);

  for(int w=0;w<nWells;w++) { // gets syntetic seismic with estimated wavelet
    fillInnWavelet(wavelet_r[w], nzp_, dzWell[w]);
    shiftReal(shiftWell[w]/dzWell[w], wavelet_r[w], nzp_);
    fileName = "waveletShift";
    printVecToFile(fileName, wavelet_r[w], nzp_);
    Utils::fft(wavelet_r[w], wavelet_c[w], nzp_);
    fileName = "cpp";
    printVecToFile(fileName, cpp_r[w], nzp_);
    Utils::fft(cpp_r[w], cpp_c[w], nzp_);
    convolve(cpp_c[w], wavelet_c[w], synt_seis_c[w], cnzp_);
    Utils::fftInv(synt_seis_c[w], synt_seis_r[w], nzp_); // 
    fileName = "syntSeis";
    printVecToFile(fileName, synt_seis_r[w], nzp_);
    fileName = "seis";
    printVecToFile(fileName, seis_r[w], nzp_);

    std::vector<float> syntSeis(nz_, 0.0f); // Do not use RMISSING (fails in setLogFromVerticalTrend())
    if (wellWeight[w] > 0) {
      for (int i = sampleStart[w] ; i < sampleStop[w] ; i++)
        syntSeis[i] = synt_seis_r[w][i];
      wells[w]->getBlockedLogsOrigThick()->setLogFromVerticalTrend(&syntSeis[0], z0[w], dzWell[w], nz_,
                                                                   "WELL_SYNTHETIC_SEISMIC", iAngle);
    }
  }

  float scaleOpt = findOptimalWaveletScale(synt_seis_r, seis_r, nWells, nzp_, wellWeight);

  shiftAndScale(shiftAvg, scaleOpt);//shifts wavelet average from wells
  invFFT1DInPlace();
  waveletLength_ = findWaveletLength(modelSettings->getMinRelWaveletAmp());
  LogKit::LogFormatted(LogKit::Low,"  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( ModelSettings::getDebugLevel() > 0 ) 
    writeWaveletToFile("estimated_wavelet_", 1.0f);

  norm_ = findNorm();

  //Writing well wavelets to file. Using writeWaveletToFile, so manipulating rAmp_
  fftw_real * trueAmp = rAmp_;
  rAmp_               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
  cAmp_               = reinterpret_cast<fftw_complex *>(rAmp_);
  for(int w=0; w<nWells; w++) {
    if (wells[w]->getUseForWaveletEstimation() && 
       ((modelSettings->getWaveletOutputFlag() & IO::WELL_WAVELETS)>0 || modelSettings->getEstimationMode()))
    {
      for(int i=0;i<nzp_;i++)
        rAmp_[i] = wellWavelets[w][i];
      std::string wellname(wells[w]->getWellname());
      NRLib::Substitute(wellname,"/","_");
      NRLib::Substitute(wellname," ","_");
      fileName = IO::PrefixWellWavelet() + wellname + "_"; 
      writeWaveletToFile(fileName, 1.0f);
    }
  }
  fftw_free(rAmp_);
  rAmp_ = trueAmp;
  cAmp_ = reinterpret_cast<fftw_complex *>(rAmp_);

  if(ModelSettings::getDebugLevel() > 1)
    writeDebugInfo(seis_r, cor_cpp_r, ccor_seis_cpp_r, cpp_r, nWells);

  for(int i=0;i<nWells;i++) {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    delete [] synt_seis_r[i] ;
    delete [] cor_cpp_r[i] ;
    delete [] ccor_seis_cpp_r[i] ;
    delete [] wavelet_r[i];
  }
  delete [] cpp_r; 
  delete [] seis_r;
  delete [] synt_seis_r;
  delete [] cor_cpp_r;
  delete [] ccor_seis_cpp_r;
  delete [] wavelet_r;
}


Wavelet1D::Wavelet1D(const std::string & fileName, 
            int                 fileFormat, 
            ModelSettings     * modelSettings, 
            float             * reflCoef,
            float               theta,
            int               & errCode, 
            std::string       & errText)
  : Wavelet(fileName, fileFormat, modelSettings, reflCoef, theta, 1, errCode, errText)
{
}


Wavelet1D::Wavelet1D(Wavelet * wavelet)
  : Wavelet(1, wavelet)
{
}

Wavelet1D::~Wavelet1D()
{
}

Wavelet1D::Wavelet1D(ModelSettings * modelSettings, 
                     float         * reflCoef,
                     float           theta,
                     float          peakFrequency,
                     int           & errCode)
  : Wavelet(modelSettings, reflCoef, theta, 1, peakFrequency, errCode)
{
}

float 
Wavelet1D::findGlobalScaleForGivenWavelet(ModelSettings * modelSettings, 
                                          Simbox        * simbox,
                                          FFTGrid       * seisCube, 
                                          WellData     ** wells)
{
  int nWells          = modelSettings->getNumberOfWells();
//The reason why the next three variables are not taken from the class variables is that this function is called
//before resample in Model::processWavelets. This means that at this point nz_, nzp_ and rnzp_ are not properly set.
  int nz              = simbox->getnz();                     
  int nzp             = seisCube->getNzp();
  int rnzp            = 2*(nzp/2+1);

  fftw_real ** cpp_r  = new fftw_real*[nWells];
  fftw_real ** seis_r = new fftw_real*[nWells];

  for(int i=0;i<nWells;i++) {
    cpp_r[i]      = new fftw_real[rnzp];
    seis_r[i]     = new fftw_real[rnzp];
    for(int j=0; j<rnzp; j++) {
      cpp_r[i][j] = 0;
      seis_r[i][j] = 0;
    }
  }
  float maxSeis = 0.0;
  float maxCpp = 0.0;
  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++)  {
    if (wells[w]->getUseForWaveletEstimation()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
      std::vector<float> seisLog(wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks());
      bl->getBlockedGrid(seisCube, &seisLog[0]);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      std::vector<bool> hasData(nz);
      std::vector<float> seisData(nz);
      bl->getVerticalTrend(&seisLog[0], &seisData[0]);
      for (int k = 0 ; k < nz; k++)
        hasData[k] = seisData[k] != RMISSING;
      int start,length;
      bl->findContiniousPartOfData(hasData,nz,start,length);
      bl->fillInCpp(coeff_,start,length,cpp_r[w],nzp); 
      bl->fillInSeismic(&seisData[0],start,length,seis_r[w],nzp);

      for(int i=0;i<nzp;i++) {
        if(cpp_r[w][i] > maxCpp)
          maxCpp = cpp_r[w][i];
        if(seis_r[w][i] > maxSeis)
          maxSeis = seis_r[w][i];
      }
    }
  }

  for(int i=0;i<nWells;i++) {
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


float         
Wavelet1D::calculateSNRatioAndLocalWavelet(Simbox        * simbox, 
                                           FFTGrid       * seisCube, 
                                           WellData     ** wells, 
                                           ModelSettings * modelSettings,
                                           std::string   & errText, 
                                           int           & error, 
                                           int             number,
                                           Grid2D      *& noiseScaled, 
                                           Grid2D      *& shiftGrid, 
                                           Grid2D      *& gainGrid)
{
  LogKit::LogFormatted(LogKit::Medium,"\n  Estimating noise from seismic data and (nonfiltered) blocked wells");
  float errStd  = 0.0f;
  float dataVar = 0.0f;

  std::string angle                 = NRLib::ToString((180.0/M_PI)*theta_, 1);
  int         nWells                = modelSettings->getNumberOfWells();
  bool        doEstimateLocalShift  = modelSettings->getEstimateLocalShift(number);
  bool        doEstimateLocalScale  = modelSettings->getEstimateLocalScale(number);
  bool        doEstimateLocalNoise  = modelSettings->getEstimateLocalNoise(number);
  bool        doEstimateGlobalScale = modelSettings->getEstimateGlobalWaveletScale(number);
  bool        doEstimateSNRatio     = modelSettings->getEstimateSNRatio(number);
  bool        doEstimateWavelet     = modelSettings->getEstimateWavelet(number);

  bool        estimateSomething     = doEstimateLocalShift || doEstimateLocalScale || doEstimateLocalNoise 
                                      || doEstimateGlobalScale || doEstimateSNRatio || doEstimateWavelet;

  //Noise estimation
  fftw_real    ** cpp_r           = new fftw_real*[nWells];
  fftw_complex ** cpp_c           = reinterpret_cast<fftw_complex**>(cpp_r);
  
  fftw_real    ** seis_r          = new fftw_real*[nWells];
  fftw_complex ** seis_c          = reinterpret_cast<fftw_complex**>(seis_r); 

  fftw_real    ** synt_r          = new fftw_real*[nWells];
  fftw_complex ** synt_c          = reinterpret_cast<fftw_complex**>(synt_r);

  fftw_real    ** cor_seis_synt_r = new fftw_real*[nWells];
  fftw_complex ** cor_seis_synt_c = reinterpret_cast<fftw_complex**>(cor_seis_synt_r);

  fftw_real    ** wavelet_r       = new fftw_real*[nWells];
  fftw_complex ** wavelet_c       = reinterpret_cast<fftw_complex**>(wavelet_r);

  std::vector<float>  dzWell(nWells);

  for(int w=0; w<nWells; w++) {
    cpp_r[w]           = new fftw_real[rnzp_];
    seis_r[w]          = new fftw_real[rnzp_];
    cor_seis_synt_r[w] = new fftw_real[rnzp_];
    wavelet_r[w]       = new fftw_real[rnzp_];
    synt_r[w]          = new fftw_real[rnzp_];
    const std::vector<int> ipos   = wells[w]->getBlockedLogsOrigThick()->getIposVector();
    const std::vector<int> jpos   = wells[w]->getBlockedLogsOrigThick()->getJposVector();
    dzWell[w]          = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])) * dz_;
  }
  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  std::vector<float> dataVarWell(nWells, 0.0f);
  std::vector<float> errVarWell (nWells, 0.0f);
  std::vector<float> shiftWell  (nWells, 0.0f);
  std::vector<int>   nActiveData(nWells, 0);
  for (int w = 0 ; w < nWells ; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
      std::vector<float> seisLog(bl->getNumberOfBlocks());
      bl->getBlockedGrid(seisCube, &seisLog[0]);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      std::vector<float> alpha(nz_);
      bl->getVerticalTrend(bl->getAlpha(), &alpha[0]);
      std::vector<float> beta(nz_);
      bl->getVerticalTrend(bl->getBeta(), &beta[0]);
      std::vector<float> rho(nz_);
      bl->getVerticalTrend(bl->getRho(), &rho[0]);
      std::vector<float> seisData(nz_);
      bl->getVerticalTrend(&seisLog[0], &seisData[0]);
      std::vector<bool> hasData(nz_);
      for (int k = 0 ; k < nz_ ; k++)
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;

      int start,length;
      bl->findContiniousPartOfData(hasData, nz_, start, length);

      if(length * dz_ > waveletLength_) { // must have enough data
        bl->fillInCpp(coeff_, start, length, cpp_r[w], nzp_);  // fills in reflection coefficients
        Utils::fft(cpp_r[w], cpp_c[w], nzp_);
        fillInnWavelet(wavelet_r[w], nzp_, dzWell[w]); // fills inn wavelet
        Utils::fft(wavelet_r[w], wavelet_c[w], nzp_);
        convolve(cpp_c[w], wavelet_c[w], synt_c[w], cnzp_);
        bl->fillInSeismic(&seisData[0], start, length, seis_r[w], nzp_);
        Utils::fft(seis_r[w], seis_c[w], nzp_);
        bl->estimateCor(synt_c[w], seis_c[w], cor_seis_synt_c[w], cnzp_);
        Utils::fftInv(cor_seis_synt_c[w], cor_seis_synt_r[w], nzp_);
        //Estimate shift. Do not run if shift given, use given shift.        
        float shift = findBulkShift(cor_seis_synt_r[w], 
                                    dzWell[w], 
                                    nzp_, 
                                    modelSettings->getMaxWaveletShift());
        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w] = shift;
        Utils::fftInv(synt_c[w], synt_r[w], nzp_);
        shiftReal(-shift/dzWell[w], synt_r[w], nzp_);
        bl->fillInSeismic(&seisData[0], start, length, seis_r[w], nzp_);
        if(ModelSettings::getDebugLevel() > 0) {
          std::string fileName;
          fileName = "seismic_Well_" + NRLib::ToString(w) + "_" + angle;
          printVecToFile(fileName,seis_r[w], nzp_);
          fileName = "synthetic_seismic_Well_" + NRLib::ToString(w) + "_" + angle;
          printVecToFile(fileName,synt_r[w], nzp_);
        }
        for(int i=start;i<start+length;i++) { 
          float err=(seis_r[w][i] - synt_r[w][i]);
          errVarWell[w]+=err*err;
          dataVarWell[w]+=seis_r[w][i] *seis_r[w][i] ;
        }
        nActiveData[w]=length;
      }
      else {
        LogKit::LogFormatted(LogKit::Low, "\n  Not using vertical well %s for error estimation (length=%.1fms  required length=%.1fms).",
                             wells[w]->getWellname().c_str(), length*dz_, waveletLength_);
      }
    }
  }
  
  float globalScale = modelSettings->getWaveletScale(number);
  std::vector<float> scaleOptWell(nWells, -1.0f);
  std::vector<float> errWellOptScale(nWells);
  std::vector<float> errWell(nWells);
  float errOptScale = 1.0;
  //Estimate global scale, local scale and error.
  //If global scale given, do not use return value. Do kriging with global scale as mean.
  //If local scale given, run separate routine to find local noise if wanted.
  float optScale;
  if(doEstimateLocalScale==true || doEstimateGlobalScale==true) {
    optScale = findOptimalWaveletScale(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, scaleOptWell, errWellOptScale);
    if(doEstimateGlobalScale==false)
      optScale = globalScale;
    else {
      scale(optScale);
      for(int i=0;i<nWells;i++)
        scaleOptWell[i]/=optScale;
    }
  }
  // local scale given means gain !=NULL
  else if(doEstimateLocalNoise==true && gainGrid!=NULL) {
    optScale = globalScale; // Global scale must be given if local scale is given
    findLocalNoiseWithGainGiven(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, scaleOptWell, errWellOptScale, gainGrid, wells, simbox);
  }
  else {
    optScale = globalScale;
    // only for loging
    findOptimalWaveletScale(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, scaleOptWell, errWellOptScale);
  }
    
  for(int w=0; w<nWells; w++) {
    delete [] cpp_r[w]; 
    delete [] seis_r[w] ;
    delete [] synt_r[w] ;
    delete [] wavelet_r[w];
    delete [] cor_seis_synt_r[w];
  }
  delete [] cpp_r;
  delete [] seis_r;
  delete [] synt_r;
  delete [] wavelet_r;
  delete [] cor_seis_synt_r;

  int nData=0;
  for(int w=0; w<nWells; w++) {
    nData   += nActiveData[w];
    errStd  += errVarWell[w];
    dataVar += dataVarWell[w];
    if(nActiveData[w]>0) {    
      errVarWell[w]  /= nActiveData[w];
      dataVarWell[w] /= nActiveData[w];
    }
  }

  if (nData == 0) {
    errText += "Cannot estimate signal-to-noise ratio. No legal well data available.\n";
    error += 1;
  }

  dataVar /= static_cast<float> (nData);
  errStd  /= static_cast<float> (nData);
  errStd   = sqrt(errStd);

  LogKit::LogFormatted(LogKit::Medium,"\n  Reporting errors (as standard deviations) estimated in different ways:\n\n");

  LogKit::LogFormatted(LogKit::Low,"\n");
  LogKit::LogFormatted(LogKit::Low,"                                     SeisData       OptimalGlobal      OptimalLocal\n");
  LogKit::LogFormatted(LogKit::Low,"  Well                  shift[ms]     StdDev         Gain   S/N         Gain   S/N \n");
  LogKit::LogFormatted(LogKit::Low,"  ----------------------------------------------------------------------------------\n");
  for(int w=0; w<nWells; w++) {
    if(nActiveData[w]>0) {    
      float SNOptimalGlobal = dataVarWell[w]/(errWell[w]*errWell[w]);
      float SNOptimalLocal  = dataVarWell[w]/(errWellOptScale[w]*errWellOptScale[w]);   
      LogKit::LogFormatted(LogKit::Low,"  %-20s   %6.2f     %9.2e      %6.2f %6.2f      %6.2f %6.2f\n", 
            wells[w]->getWellname().c_str(),shiftWell[w],sqrt(dataVarWell[w]),
            optScale,SNOptimalGlobal,scaleOptWell[w],SNOptimalLocal);
    }
    else
      LogKit::LogFormatted(LogKit::Low,"  %-20s      -            -             -      -           -      -\n",wells[w]->getWellname().c_str());
  }

  if (estimateSomething) {
    float minLocalGain = 0.3334f;
    float maxLocalGain = 3.0f;
    for(int w=0; w<nWells; w++) {
      if((scaleOptWell[w] >= maxLocalGain || scaleOptWell[w] <= minLocalGain) && nActiveData[w]>0) {
        std::string text;
        text  = "\nWARNING: An optimal local gain cannot be established for well "+wells[w]->getWellname()+". The gain-value found ("+NRLib::ToString(scaleOptWell[w],3)+")\n";
        text += "         is outside the interval accepted by CRAVA which is <"+NRLib::ToString(minLocalGain,2)+", "+NRLib::ToString(maxLocalGain,2)+">.\n";
        LogKit::LogFormatted(LogKit::Warning, text);
        if (modelSettings->getEstimateWavelet(number)) {
          TaskList::addTask("Well "+wells[w]->getWellname()+" should not be used in the wavelet estimation for angle stack "+angle+".");
        }
      }
    }
  }

  if(doEstimateLocalScale==true) {
    // Estimate global noise with local waveletscale
    dataVar = 0.0;
    errStd = 0.0;
    for(int w=0; w<nWells; w++) {
      dataVar +=(dataVarWell[w]*nActiveData[w]);
      errStd  +=(errWellOptScale[w]*errWellOptScale[w]*nActiveData[w]);
    }
    dataVar/=nData;
    errStd/=nData;
    errStd = sqrt(errStd);
  }
  else if(doEstimateGlobalScale==true)
   errStd = errOptScale;

  if(doEstimateLocalShift || doEstimateLocalScale || doEstimateLocalNoise) {
    //
    // Pretabulate correlations
    //
    Vario  * localWaveletVario     = modelSettings->getLocalWaveletVario();
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
      estimateLocalShift(cov, shiftGrid, shiftWell, nActiveData, simbox,wells, nWells);
    
    if(doEstimateLocalScale)
      estimateLocalGain(cov, gainGrid, scaleOptWell, 1.0, nActiveData, simbox,wells, nWells);

    if(doEstimateLocalNoise) {
      float errStdLN;
      if(doEstimateSNRatio==true)
        errStdLN = errStd;
      else //SNRatio given in model file
        errStdLN = sqrt(dataVar/modelSettings->getSNRatio(number));
      if(gainGrid == NULL && doEstimateLocalScale==false && doEstimateGlobalScale==false) { // No local wavelet scale 
        for(int w=0; w<nWells; w++) 
          errVarWell[w] = sqrt(errVarWell[w]);
        estimateLocalNoise(cov, noiseScaled, errStdLN, errVarWell, nActiveData, simbox,wells, nWells); 
      }
      else if(doEstimateGlobalScale==true && doEstimateLocalScale==false) // global wavelet scale
        estimateLocalNoise(cov, noiseScaled, errStdLN,errWell, nActiveData, simbox,wells, nWells); 
      else
        estimateLocalNoise(cov, noiseScaled, errStdLN, errWellOptScale, nActiveData, simbox,wells, nWells); 
    }
  }

  float empSNRatio = dataVar/(errStd*errStd);
  if(doEstimateSNRatio==true)
    LogKit::LogFormatted(LogKit::Low,"\n  The signal to noise ratio used for this angle stack is: %6.2f\n", empSNRatio);
  else
  {
    float SNRatio = modelSettings->getSNRatio(number);
    LogKit::LogFormatted(LogKit::Low,"\n  The signal to noise ratio given in the model file and used for this angle stack is : %6.2f\n", SNRatio);
    LogKit::LogFormatted(LogKit::Low,"  For comparison, the signal-to-noise ratio calculated from the available wells is   : %6.2f\n", empSNRatio);
    float minSN = 1.0f + (empSNRatio - 1.0f)/2.0f;
    float maxSN = 1.0f + (empSNRatio - 1.0f)*2.0f;
    if ((SNRatio<minSN || SNRatio>maxSN) && modelSettings->getEstimateWavelet(number))
    {
      LogKit::LogFormatted(LogKit::Warning,"\nWARNING: The difference between the SN ratio given in the model file and the calculated SN ratio is too large.\n");
      if (SNRatio < minSN)
        TaskList::addTask("Consider increasing the SN ratio for angle stack "+NRLib::ToString(number)+" to minimum "+NRLib::ToString(minSN,1));
      else
        TaskList::addTask("Consider decreasing the SN ratio for angle stack "+NRLib::ToString(number)+" to maximum "+NRLib::ToString(minSN,1));
    }        
  }

  if (empSNRatio < 1.1f) {
    if (doEstimateSNRatio) {
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

float          
Wavelet1D::findOptimalWaveletScale(fftw_real               ** synt_seis_r,
                                   fftw_real               ** seis_r,
                                   int                        nWells,
                                   int                        nzp,
                                   const std::vector<float> & wellWeight,
                                   float                    & err, // NBNB-PAL: Det er uheldig å returnere err slik
                                   std::vector<float>       & errWell,
                                   std::vector<float>       & scaleOptWell,
                                   std::vector<float>       & errWellOptScale) const
{
  float   optScale   = 1.0f;
  float   scaleLimit = 3.0f;
  int     nScales    = 51; // should be odd to include 1.00
  float * scales     = new float[nScales];
  float * error      = new float[nScales];

  for(int i=0;i<nScales;i++) {
    scales[i] = exp(-log(scaleLimit)+i*2*(log(scaleLimit))/(nScales-1));
    error[i]  = 0.0f;
  }

  int    * counter  = new int[nWells];
  float  * seisNorm = new float[nWells];
  float ** resNorm  = new float*[nWells];

  for(int i=0;i<nWells;i++) {
    resNorm[i]  = new float[nScales];
    seisNorm[i] = 0.0f;
  }

  float minSeisAmp = static_cast<float> (1e-7);
  int totCount=0;

  for(int i=0;i<nWells;i++) {
    counter[i]=0;
    if(wellWeight[i]>0) {
      // Count number of layers with seismic data
      for(int k=0;k<nzp;k++)
        if(fabs(seis_r[i][k]) > minSeisAmp)
          counter[i]++;
      totCount+=counter[i];

      for(int j=0;j<nScales;j++) {
        resNorm[i][j]=0.0;
        for(int k=0;k<nzp;k++) {
          if(fabs(seis_r[i][k]) > minSeisAmp) {
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
    if(error[i]<optValue) {
      optValue=error[i];
      optInd=i;
    }

    delete [] error;

    err = sqrt(optValue/static_cast<float>(totCount));
    optScale = scales[optInd];

    for(int i=0;i<nWells;i++) {
      if(counter[i]>0)
        errWell[i] = sqrt(resNorm[i][optInd]/counter[i]);
      else
        errWell[i] = 0.0f;
    }

    for(int i=0;i<nWells;i++) {
      if(wellWeight[i]>0) {
        optValue = resNorm[i][0];
        optInd=0;
        for(int j=1;j<nScales;j++)
          if(resNorm[i][j]<optValue) {
            optValue=resNorm[i][j];
            optInd=j;
          }
          scaleOptWell[i]    = scales[optInd];
          errWellOptScale[i] = sqrt(optValue/float(counter[i]));
      }
      else {
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
Wavelet1D::findOptimalWaveletScale(fftw_real                ** synt_seis_r,
                                   fftw_real                ** seis_r,
                                   int                         nWells,
                                   int                         nzp,
                                   const std::vector<float>  & wellWeight) const
{
  //Wrapper used when we don't need err, errWell, scaleOptWell and errWellOptScale later. Used in first constructor
  float err;
  std::vector<float> scaleOptWell(nWells);
  std::vector<float> errWellOptScale(nWells);
  std::vector<float> errWell(nWells);
  float scaleOpt = findOptimalWaveletScale(synt_seis_r, seis_r, nWells, nzp, wellWeight, err, errWell, scaleOptWell, errWellOptScale);
  return scaleOpt;
}

void 
Wavelet1D::findLocalNoiseWithGainGiven(fftw_real              ** synt_seis_r,
                                     fftw_real                ** seis_r,
                                     int                         nWells,
                                     int                         nzp,
                                     const std::vector<float>  & wellWeight,
                                     float                     & err,
                                     std::vector<float>        & errWell,
                                     std::vector<float>        & scaleOptWell,
                                     std::vector<float>        & errWellOptScale, 
                                     Grid2D                    * gain, 
                                     WellData                 ** wells, 
                                     Simbox                    * simbox) const
{
  double *scale = new double[nWells];   
  float error = 0.0; 

  int    * counter  = new int[nWells];
  float  * seisNorm = new float[nWells];
  float * resNorm  = new float[nWells];

  for(int i=0;i<nWells;i++) {
    resNorm[i]  = 0.0f;
    seisNorm[i] = 0.0f;
  }

  float minSeisAmp = static_cast<float> (1e-7);
  int totCount=0;
  const double *x, *y;
  int nData;
  for(int i=0;i<nWells;i++) {
    x = wells[i]->getXpos(nData);
    y = wells[i]->getYpos(nData);
    int ix, iy;
    simbox->getIndexes(x[0],y[0],ix,iy);
    //scale[i] = gain->GetZ(x[0],y[0]);
    scale[i] = (*gain)(ix,iy);
    counter[i]=0;
    if(wellWeight[i]>0) {
      // Count number of layers with seismic data
      for(int k=0;k<nzp;k++)
        if(fabs(seis_r[i][k]) > minSeisAmp)
          counter[i]++;
      totCount+=counter[i];
      for(int k=0;k<nzp;k++) {
        if(fabs(seis_r[i][k]) > minSeisAmp) {
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
  for(int i=0;i<nWells;i++) {
    if(counter[i]>0)
      errWell[i] = sqrt(resNorm[i]/counter[i]);
    else
      errWell[i] = 0.0f;
  }

  for(int i=0;i<nWells;i++) {
    if(wellWeight[i]>0) {
      optValue = resNorm[i];
      scaleOptWell[i]    = float(scale[i]);
      errWellOptScale[i] = sqrt(optValue/float(counter[i]));
    }
    else {
      scaleOptWell[i]    = 0.0f;
      errWellOptScale[i] = 0.0f;
    }
  }

  delete [] resNorm;
  delete [] seisNorm;
  delete [] counter;
  delete [] scale;
}

void
Wavelet1D::estimateLocalShift(const CovGrid2D          & cov,
                              Grid2D                  *& shift,
                              const std::vector<float> & shiftWell,
                              const std::vector<int>   & nActiveData,
                              Simbox                   * simbox,
                              WellData                ** wells,
                              int                        nWells)
{
  //
  // NBNB-PAL: Since slightly deviated wells are accepted, we should
  // eventually make gain- and shift-cubes rather than single maps.
  //

  //
  // Collect data for kriging
  //
  KrigingData2D shiftData;
  
  for(int i=0;i<nWells;i++) {
    if(nActiveData[i]>0)  {
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
Wavelet1D::estimateLocalGain(const CovGrid2D           & cov,
                             Grid2D                   *& gain,
                             const std::vector<float>  & scaleOptWell,
                             float                       globalScale, 
                             const std::vector<int>    & nActiveData,
                             Simbox                    * simbox,
                             WellData                 ** wells,
                             int                         nWells)
{
  //
  // Collect data for kriging
  //
  KrigingData2D gainData;
  
  for(int i=0;i<nWells;i++) {
    if(nActiveData[i]>0)  {
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
Wavelet1D::estimateLocalNoise(const CovGrid2D          & cov,
                              Grid2D                  *& noiseScaled,
                              float                      globalNoise,
                              const std::vector<float> & errWellOptScale,
                              const std::vector<int>   & nActiveData,
                              Simbox                   * simbox,
                              WellData                ** wells,
                              int                        nWells)
{
  //
  // Collect data for kriging
  //
  KrigingData2D noiseData;
  
  for(int i=0;i<nWells;i++) {
    if(nActiveData[i]>0)  {
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


float
Wavelet1D::shiftOptimal(fftw_real               ** ccor_seis_cpp_r,
                        const std::vector<float> & wellWeight,
                        const std::vector<float> & dz,
                        int                        nWells,
                        int                        nzp,
                        std::vector<float>       & shiftWell,
                        float                      maxShift)
{
  float shift=0.0f;
  float totalWeight=0;
  float sum=0;
  int w,i,polarity;
  // if the sum from -maxShift_ to maxShift_ ms is 
  // positive then polarity is positive
  for(w=0;w<nWells;w++)
  {
    shiftWell[w]=0;
    if(wellWeight[w] > 0) // only include those estimated
    {
      for(i=0;i<ceil(maxShift/dz[w]);i++)//zero included
        sum+=ccor_seis_cpp_r[w][i];
      for(i=0;i<floor(maxShift/dz[w]);i++)
        sum+=ccor_seis_cpp_r[w][nzp-i-1];
    }
  }

  polarity=-1;
  if(sum > 0)
    polarity=1;

  // gets optimal shift
  float maxValue;
  float shiftF;
  int shiftI;
  float f1,f2,f3;

  for(w=0;w<nWells;w++) {
    if(wellWeight[w]>0) {
      maxValue = 0.0f;
      shiftI=0;
      for(i=0;i<ceil(maxShift/dz[w]);i++) {
        if(ccor_seis_cpp_r[w][i]*polarity > maxValue) {
          maxValue = ccor_seis_cpp_r[w][i]*polarity;
          shiftI = i;
        }
      }
      for(i=0;i<floor(maxShift/dz[w]);i++) {
        if(ccor_seis_cpp_r[w][nzp-1-i]*polarity > maxValue) {
          maxValue = ccor_seis_cpp_r[w][nzp-1-i]*polarity;
          shiftI = -1-i;
        }
      }
      if(shiftI < 0) {
        if(ccor_seis_cpp_r[w][nzp+shiftI-1]*polarity < maxValue) { //then local max 
          f1 = ccor_seis_cpp_r[w][nzp+shiftI-1];
          f2 = ccor_seis_cpp_r[w][nzp+shiftI];
          int ind3;
          if(shiftI==-1)
            ind3 = 0;
          else
            ind3=nzp+shiftI+1;
          f3 = ccor_seis_cpp_r[w][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      else { //positive or zero shift
        if(ccor_seis_cpp_r[w][shiftI+1]*polarity < maxValue) { //then local max
          f3 = ccor_seis_cpp_r[w][shiftI+1];
          f2 = ccor_seis_cpp_r[w][shiftI];
          int ind1;
          if(shiftI==0)
            ind1 = nzp-1;
          else
            ind1=shiftI-1;
          f1 = ccor_seis_cpp_r[w][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI);
      }
      shiftWell[w] = shiftF*dz[w];
      shiftReal(-shiftF, ccor_seis_cpp_r[w],nzp);// 
      shift += wellWeight[w]*shiftF*dz[w];//weigthing shift according to wellWeight
      totalWeight += wellWeight[w];
    }
  }
  shift/=totalWeight;
  return shift;
}

void
Wavelet1D::multiplyPapolouis(fftw_real                ** vec, 
                             const std::vector<float>  & dz,
                             int                         nWells,
                             int                         nzp, 
                             float                       waveletLength,
                             const std::vector<float>  & wellWeight) const
{
  float wHL=float( waveletLength/2.0);
  float weight,dist;
  for(int w=0;w<nWells;w++) {
    if(wellWeight[w] > 0) {
      for(int i=1;i<nzp;i++) {
        dist = std::min(i,nzp-i)*dz[w];
        if(dist < wHL) {
          weight  = float(1.0/M_PI*fabs(sin(M_PI*(dist)/wHL))); 
          weight += float((1-fabs(dist)/wHL)*cos(M_PI*dist/wHL));
        }
        else
          weight=0;
        vec[w][i]*=weight;
      }
    }
  }
}

void
Wavelet1D::getWavelet(fftw_real               ** ccor_seis_cpp_r,
                      fftw_real               ** cor_cpp_r,
                      fftw_real               ** wavelet_r,
                      const std::vector<float> & wellWeight,
                      int                        nWells,
                      int                        nt)
{
  fftw_complex* c_sc,*c_cc,*wav;

  int cnzp = nt/2+1;
  for(int w=0;w<nWells;w++)
  {
    if(wellWeight[w] > 0)
    {
      c_sc   = reinterpret_cast<fftw_complex*>(ccor_seis_cpp_r[w]);
      Utils::fft(ccor_seis_cpp_r[w],c_sc,nt);
      c_cc   = reinterpret_cast<fftw_complex*>(cor_cpp_r[w]);
      Utils::fft(cor_cpp_r[w],c_cc,nt);
      wav    = reinterpret_cast<fftw_complex*>(wavelet_r[w]);

      for(int i=0;i<cnzp;i++)
      {
        wav[i].re=c_sc[i].re/c_cc[i].re;//note c_cc[i].im =0
        wav[i].im=c_sc[i].im/c_cc[i].re;
      }
      Utils::fftInv(c_sc,ccor_seis_cpp_r[w],nt);
      Utils::fftInv(c_cc,cor_cpp_r[w],nt);
      Utils::fftInv(wav,wavelet_r[w],nt);
    }
  }
}


void
Wavelet1D::convolve(fftw_complex  * var1_c ,
                    fftw_complex  * var2_c, 
                    fftw_complex  * out_c,
                    int             cnzp) const
{
  for(int i=0;i<cnzp;i++) {
    out_c[i].re = var1_c[i].re*var2_c[i].re+var1_c[i].im*var2_c[i].im; 
    out_c[i].im = var1_c[i].im*var2_c[i].re - var1_c[i].re*var2_c[i].im;
  }
}

void
Wavelet1D::writeDebugInfo(fftw_real ** seis_r,
                          fftw_real ** cor_cpp_r,
                          fftw_real ** ccor_seis_cpp_r,
                          fftw_real ** cpp_r,
                          int          nWells) const
{
  std::string angle    = NRLib::ToString(theta_*180/M_PI, 1);
  std::string fileName = "estimated_wavelet_full_"+angle+".txt";
        
  std::ofstream fid;
  NRLib::OpenWrite(fid,fileName);
  fid << std::setprecision(6);
  for(int i=0;i<nzp_;i++)
    fid <<  rAmp_[i]  << "\n";
  fid.close();
  fid.clear();
    
  for(int j=0;j<nWells;j++) {
    std::string jp1 = NRLib::ToString(j+1);
    fileName = "seis_"+angle+"_well_"+jp1+".txt";
    NRLib::OpenWrite(fid,fileName);
    fid << std::setprecision(6);
    for(int i=0;i<nzp_;i++)
      fid <<  seis_r[0][i]  << "\n";
    fid.close();
    fid.clear();
        
    fileName = "cor_cpp_"+angle+"well_"+jp1+".txt";
    NRLib::OpenWrite(fid,fileName);
    fid << std::setprecision(6);
    for(int i=0;i<nzp_;i++)
      fid <<  cor_cpp_r[0][i]  << "\n";
    fid.close();
    fid.clear();
      
    fileName = "ccor_seis_cpp_"+angle+"_well_"+jp1;
    NRLib::OpenWrite(fid,fileName);
    fid << std::setprecision(6);
    for(int i=0;i<nzp_;i++)
      fid <<  ccor_seis_cpp_r[0][i]  << "\n";
    fid.close();
    fid.clear();

    fileName = "cpp_"+angle+"_well_"+jp1;
    NRLib::OpenWrite(fid,fileName);
    fid << std::setprecision(6);
    for(int i=0;i<nzp_;i++)
      fid << cpp_r[0][i] << "\n";
    fid.close();
    fid.clear();
  }
}
