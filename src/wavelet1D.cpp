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

#include "lib/global_def.h"
#include "lib/lib_matr.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/modelsettings.h"
#include "src/blockedlogs.h"
#include "src/definitions.h"
#include "src/wavelet1D.h"
#include "src/welldata.h"
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
  LogKit::LogFormatted(LogKit::MEDIUM,"  Estimating 1D wavelet from seismic data and (nonfiltered) blocked wells\n");

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

  std::string fileName;
  int     nWells        = modelSettings->getNumberOfWells();
  float   waveletLength = modelSettings->getWaveletTaperingL();
  float * dzWell        = new float[nWells];
  float * wellWeight    = new float[nWells];
  float * z0            = new float[nWells]; // Needed to block syntSeis
  int   * sampleStart   = new int[nWells];   // Needed to block syntSeis
  int   * sampleStop    = new int[nWells];   // Needed to block syntSeis

  float * alpha         = new float[nz_];
  float * beta          = new float[nz_];
  float * rho           = new float[nz_];
  float * seisData      = new float[nz_];
  bool  * hasData       = new bool[nz_];
 
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

  int maxBlocks = 0;
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
    wellWeight[i]      = 0;
    sampleStart[i]     = 0;
    sampleStop[i]      = 0;
    const int * ipos   = wells[i]->getBlockedLogsOrigThick()->getIpos();
    const int * jpos   = wells[i]->getBlockedLogsOrigThick()->getJpos();
    dzWell[i]          = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])) * dz_;
    int nBlocks        = wells[i]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
      LogKit::LogFormatted(LogKit::MEDIUM,"  Well :  %s\n",wells[w]->getWellname().c_str());

      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Check seismic data outside estimation interval missing
      //
      if (estimInterval.size() > 0) {
        const double * xPos  = bl->getXpos();
        const double * yPos  = bl->getYpos();
        const double * zPos  = bl->getZpos();
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
      bl->getVerticalTrend(bl->getAlpha(), alpha);
      bl->getVerticalTrend(bl->getBeta(), beta);
      bl->getVerticalTrend(bl->getRho(), rho);
      bl->getVerticalTrend(seisLog, seisData);
      
      for (int k = 0 ; k < nz_ ; k++)
        hasData[k] = seisData[k] != RMISSING && alpha[k] != RMISSING && beta[k] != RMISSING && rho[k] != RMISSING;
      //
      // Find continuous part of data
      //
      int start,length;
      bl->findContiniousPartOfData(hasData, nz_, start, length);
      if(length*dz_ > waveletLength ) { // must have enough data
        bl->fillInCpp(coeff_, start, length, cpp_r[w], nzp_); 
        fileName = "cpp_1";
        printVecToFile(fileName, cpp_r[w], nzp_);  // Debug
        Utils::fft(cpp_r[w], cpp_c[w], nzp_);
        bl->fillInSeismic(seisData, start, length, seis_r[w], nzp_);
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
  delete [] seisLog;
  float * shiftWell = new float[nWells];
  float shiftAvg = shiftOptimal(ccor_seis_cpp_r, wellWeight, dzWell, nWells, nzp_, shiftWell, modelSettings->getMaxWaveletShift());

  multiplyPapolouis(ccor_seis_cpp_r, dzWell, nWells, nzp_, waveletLength, wellWeight);
  multiplyPapolouis(cor_cpp_r, dzWell, nWells, nzp_, waveletLength, wellWeight);
  getWavelet(ccor_seis_cpp_r, cor_cpp_r, wavelet_r, wellWeight, nWells, nzp_);

  // Save estimated wavelet for each well, write to file later
  std::vector<std::vector<fftw_real> > wellWavelets(nWells);
  for(int w=0; w<nWells; w++) {
    wellWavelets[w].resize(nzp_);
    for(int i=0;i<nzp_;i++)
      wellWavelets[w][i] = wavelet_r[w][i];
  }

  rAmp_ = averageWavelets(wavelet_r, nWells, nzp_, wellWeight, dzWell, dz_); // wavelet centered
  cAmp_ = reinterpret_cast<fftw_complex*>(rAmp_);

  float * syntSeis = new float[nz_];
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

    if (wellWeight[w] > 0) {
      for (int i = 0 ; i < nz_ ; i++)
        syntSeis[i] = 0.0; // Do not use RMISSING (fails in setLogFromVerticalTrend())
      for (int i = sampleStart[w] ; i < sampleStop[w] ; i++)
        syntSeis[i] = synt_seis_r[w][i];
      wells[w]->getBlockedLogsOrigThick()->setLogFromVerticalTrend(syntSeis, z0[w], dzWell[w], nz_,
                                                                   "WELL_SYNTHETIC_SEISMIC", iAngle);
    }
  }
  delete [] syntSeis;
  delete [] shiftWell;

  float * scaleOptWell    = new float[nWells];
  float * errWellOptScale = new float[nWells];
  float * errWell         = new float[nWells];

  float err;
  float scaleOpt = findOptimalWaveletScale(synt_seis_r, seis_r, nWells, nzp_, wellWeight, err, errWell, scaleOptWell, errWellOptScale);

  delete [] scaleOptWell ;
  delete [] errWellOptScale;
  delete [] errWell ;
  
  //  shiftReal(shiftAvg/dz0,rAmp_,nzp);
  
  shiftAndScale(shiftAvg, scaleOpt);//shifts wavelet average from wells
  invFFT1DInPlace();
  waveletLength_ = findWaveletLength(modelSettings->getMinRelWaveletAmp());
  LogKit::LogFormatted(LogKit::LOW,"  Estimated wavelet length:  %.1fms\n",waveletLength_);

  if( ModelSettings::getDebugLevel() > 0 ) {
    fileName = "estimated_wavelet";
    float dzOut = 1.0; // sample at least as dense as this
    writeWaveletToFile(fileName, dzOut);
  }
  double norm2=0.0;
  for(int i=0; i < nzp_; i++ )
      norm2 += rAmp_[i]*rAmp_[i];
  norm_= float( sqrt(norm2));

  //Writing well wavelets to file. Using writeWaveletToFile, so manipulating rAmp_
  fftw_real * trueAmp = rAmp_;
  rAmp_               = static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));
  cAmp_               = reinterpret_cast<fftw_complex *>(rAmp_);

  for(int w=0; w<nWells; w++) {
    for(int i=0;i<nzp_;i++)
      rAmp_[i] = wellWavelets[w][i];
    fileName = "Wavelet"; 
    std::string wellname(wells[w]->getWellname());
    NRLib::Substitute(wellname,"/","_");
    NRLib::Substitute(wellname," ","_");

    fileName += "_"+wellname; 
    float dzOut = 1.0; 
    writeWaveletToFile(fileName, dzOut);
  }
  fftw_free(rAmp_);
  rAmp_ = trueAmp;
  cAmp_ = reinterpret_cast<fftw_complex *>(rAmp_);

  if(ModelSettings::getDebugLevel() > 1) {
    fileName = "estimated_wavelet_full_"+NRLib::ToString(int(ceil((theta_*180/PI)-0.5)))+".txt";
        
    std::ofstream fid;
    NRLib::OpenWrite(fid,fileName);
    fid << std::setprecision(6);
    for(int i=0;i<nzp_;i++)
      fid <<  rAmp_[i]  << "\n";
    fid.close();
    fid.clear();
    
    for(int j=0;j<nWells;j++) {
      fileName = "seis_"+NRLib::ToString(int(theta_*180/PI+0.5))+"_well_"+NRLib::ToString(j+1)+".txt";
      NRLib::OpenWrite(fid,fileName);
      fid << std::setprecision(6);
      for(int i=0;i<nzp_;i++)
        fid <<  seis_r[0][i]  << "\n";
      fid.close();
      fid.clear();
        
      fileName = "cor_cpp_"+NRLib::ToString(int(theta_*180/PI+0.5))+"well_"+NRLib::ToString(j+1)+".txt";
      NRLib::OpenWrite(fid,fileName);
      fid << std::setprecision(6);
      for(int i=0;i<nzp_;i++)
        fid <<  cor_cpp_r[0][i]  << "\n";
      fid.close();
      fid.clear();
      
      fileName = "ccor_seis_cpp_"+NRLib::ToString(int(theta_*180/PI+0.5))+"_well_"+NRLib::ToString(j+1);
      NRLib::OpenWrite(fid,fileName);
      fid << std::setprecision(6);
      for(int i=0;i<nzp_;i++)
        fid <<  ccor_seis_cpp_r[0][i]  << "\n";
      fid.close();
      fid.clear();

      fileName = "cpp_"+NRLib::ToString(int(theta_*180/PI+0.5))+"_well_"+NRLib::ToString(j+1);
      NRLib::OpenWrite(fid,fileName);
      fid << std::setprecision(6);
      for(int i=0;i<nzp_;i++)
        fid << cpp_r[0][i] << "\n";
      fid.close();
      fid.clear();
    }
  }

  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  delete [] hasData;
  delete [] dzWell;
  delete [] z0;
  delete [] wellWeight;
  delete [] sampleStart;
  delete [] sampleStop;
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

float 
Wavelet1D::findGlobalScaleForGivenWavelet(ModelSettings * modelSettings, 
                                          Simbox        * simbox,
                                          FFTGrid       * seisCube, 
                                          WellData     ** wells)
{
  int nWells          = modelSettings->getNumberOfWells();
  int nz              = simbox->getnz();                     
  int nzp             = seisCube->getNzp();
  int rnzp            = 2*(nzp/2+1);
  fftw_real ** cpp_r  = new fftw_real*[nWells];
  fftw_real ** seis_r = new fftw_real*[nWells];
  float *seisData     = new float[nz];
  bool *hasData       = new bool[nz];
  int maxBlocks       = 0;
  for(int i=0;i<nWells;i++) {
    cpp_r[i]      = new fftw_real[rnzp];
    seis_r[i]     = new fftw_real[rnzp];
    for(int j=0;j<rnzp;j++) {
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
  for (int w = 0 ; w < nWells ; w++)  {
    if (wells[w]->getUseForWaveletEstimation()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsOrigThick();
      //
      // Block seismic data for this well
      //
      bl->getBlockedGrid(seisCube,seisLog);
      //
      // Extract a one-value-for-each-layer array of blocked logs
      //
      bl->getVerticalTrend(seisLog, seisData);
      for (int k = 0 ; k < nz ; k++)
        hasData[k] = seisData[k] != RMISSING;
      int start,length;
      bl->findContiniousPartOfData(hasData,nz,start,length);
      
      bl->fillInCpp(coeff_,start,length,cpp_r[w],nzp); 
      bl->fillInSeismic(seisData,start,length,seis_r[w],nzp);

      for(int i=0;i<nzp;i++) {
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
                                           Grid2D       *& shift, 
                                           Grid2D       *& gain, 
                                           ModelSettings * modelSettings,
                                           std::string   & errText, 
                                           int           & error, 
                                           Grid2D       *& noiseScaled, 
                                           int             number, 
                                           float           globalScale)
{
  LogKit::LogFormatted(LogKit::MEDIUM,"\n  Estimating noise from seismic data and (nonfiltered) blocked wells");
  float errStd  = 0.0f;
  float dataVar = 0.0f;
 
  shiftGrid_ = NULL;  
  gainGrid_  = NULL; 
  
  Vario  * localWaveletVario     = modelSettings->getLocalWaveletVario();
  int      nWells                = modelSettings->getNumberOfWells();
  bool     doEstimateLocalShift  = modelSettings->getEstimateLocalShift(number);
  bool     doEstimateLocalScale  = modelSettings->getEstimateLocalScale(number);
  bool     doEstimateLocalNoise  = modelSettings->getEstimateLocalNoise(number);
  bool     doEstimateGlobalScale = modelSettings->getEstimateGlobalWaveletScale(number);
  bool     doEstimateSNRatio     = modelSettings->getEstimateSNRatio(number);

  float  * dzWell                = new float[nWells];

  float  * alpha                 = new float[nz_];
  float  * beta                  = new float[nz_];
  float  * rho                   = new float[nz_];
  float  * seisData              = new float[nz_];
  bool   * hasData               = new bool[nz_];

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

  float * shiftWell               = new float[nWells];
  float * errVarWell              = new float[nWells];
  float * dataVarWell             = new float[nWells];
  int   * nActiveData             = new int[nWells];

  int maxBlocks = 0;
  for(int i=0;i<nWells;i++) {
    cpp_r[i]           = new fftw_real[rnzp_];
    seis_r[i]          = new fftw_real[rnzp_];
    cor_seis_synt_r[i] = new fftw_real[rnzp_];
    wavelet_r[i]       = new fftw_real[rnzp_];
    synt_r[i]          = new fftw_real[rnzp_];
    nActiveData[i]     = 0;
    errVarWell[i]      = 0.0f;
    dataVarWell[i]     = 0.0f;
    const int * ipos   = wells[i]->getBlockedLogsOrigThick()->getIpos();
    const int * jpos   = wells[i]->getBlockedLogsOrigThick()->getJpos();
    dzWell[i]          = static_cast<float>(simbox->getRelThick(ipos[0],jpos[0])) * dz_;
    int nBlocks        = wells[i]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * seisLog = new float[maxBlocks];

  //
  // Loop over wells and create a blocked well and blocked seismic
  //
  for (int w = 0 ; w < nWells ; w++) {
    if (wells[w]->getUseForWaveletEstimation()) {
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
        bl->fillInSeismic(seisData, start, length, seis_r[w], nzp_);
        Utils::fft(seis_r[w], seis_c[w], nzp_);
        bl->estimateCor(synt_c[w], seis_c[w], cor_seis_synt_c[w], cnzp_);
        Utils::fftInv(cor_seis_synt_c[w], cor_seis_synt_r[w], nzp_);
        //Estimate shift. Do not run if shift given, use given shift.        
        float shift=findBulkShift(cor_seis_synt_r[w], dzWell[w], nzp_, modelSettings->getMaxWaveletShift());
        shift = floor(shift*10.0f+0.5f)/10.0f;//rounds to nearest 0.1 ms (don't have more accuracy)
        shiftWell[w]=shift;
        Utils::fftInv(synt_c[w], synt_r[w], nzp_);
        shiftReal(-shift/dzWell[w], synt_r[w], nzp_);
        bl->fillInSeismic(seisData, start, length, seis_r[w], nzp_);
        if(ModelSettings::getDebugLevel() > 0) {
          std::string angle = NRLib::ToString(theta_/(M_PI*180.0),1);
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
        LogKit::LogFormatted(LogKit::LOW, "\n  Not using vertical well %s for error estimation (length=%.1fms  required length=%.1fms).",
                             wells[w]->getWellname().c_str(), length*dz_, waveletLength_);
      }
    }
  }
  delete [] seisLog;
  delete [] dzWell;
  delete [] hasData;

  float * scaleOptWell    = new float[nWells];
  for(int i=0;i<nWells;i++)
    scaleOptWell[i] = -1.0;
  float * errWellOptScale = new float[nWells];
  float * errWell         = new float[nWells];
  bool writelog = false;
  float errOptScale = 1.0;
  //Estimate global scale, local scale and error.
  //If global scale given, do not use return value. Do kriging with global scale as mean.
  //If local scale given, run separate routine to find local noise if wanted.
  float optScale;
  if(doEstimateLocalScale==true || doEstimateGlobalScale==true) {
    optScale = findOptimalWaveletScale(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, scaleOptWell, errWellOptScale);
    writelog = true;
    if(doEstimateGlobalScale==false)
      optScale = globalScale;
    else {
      scale(optScale);
      for(int i=0;i<nWells;i++)
        scaleOptWell[i]/=optScale;
    }
  }
  // local scale given means gain !=NULL
  else if(doEstimateLocalNoise==true && gain!=NULL) {
    optScale = globalScale; // Global scale must be given if local scale is given
    findLocalNoiseWithGainGiven(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, errWellOptScale, scaleOptWell, gain, wells, simbox);
    writelog = true;
  }
  else {
    optScale = globalScale;
    // only for loging
    findOptimalWaveletScale(synt_r, seis_r, nWells, nzp_, dataVarWell, errOptScale, errWell, scaleOptWell, errWellOptScale);
  }
  int nData=0;
  for(int i=0;i<nWells;i++) {
    nData   += nActiveData[i];
    errStd  += errVarWell[i];
    dataVar += dataVarWell[i];
    if(nActiveData[i]>0) {    
      errVarWell[i]  /= nActiveData[i];
      dataVarWell[i] /= nActiveData[i];
    }
  }

  if (nData == 0) {
    errText += "Cannot estimate signal-to-noise ratio. No legal well data available.\n";
    error += 1;
  }

  dataVar /= float(nData);
  errStd  /= float(nData);
  errStd   = sqrt(errStd);

  LogKit::LogFormatted(LogKit::MEDIUM,"\n  Reporting errors (as standard deviations) estimated in different ways:\n\n");

  LogKit::LogFormatted(LogKit::LOW,"\n");
  LogKit::LogFormatted(LogKit::LOW,"                                     SeisData       OptimalGlobal      OptimalLocal\n");
  LogKit::LogFormatted(LogKit::LOW,"  Well                  shift[ms]     StdDev         Gain   S/N         Gain   S/N \n");
  LogKit::LogFormatted(LogKit::LOW,"  ----------------------------------------------------------------------------------\n");
  for(int i=0;i<nWells;i++) {
    if(nActiveData[i]>0) {
      float SNOptimalGlobal, SNOptimalLocal;    
      SNOptimalGlobal = dataVarWell[i]/(errWell[i]*errWell[i]);
      SNOptimalLocal  = dataVarWell[i]/(errWellOptScale[i]*errWellOptScale[i]);   
      LogKit::LogFormatted(LogKit::LOW,"  %-20s   %6.2f     %9.2e      %6.2f %6.2f      %6.2f %6.2f\n", 
            wells[i]->getWellname().c_str(),shiftWell[i],sqrt(dataVarWell[i]),
            optScale,SNOptimalGlobal,scaleOptWell[i],SNOptimalLocal);
    }
    else
      LogKit::LogFormatted(LogKit::LOW,"  %-20s      -            -             -      -           -      -\n",
      wells[i]->getWellname().c_str()); 
  }
  for(int i=0;i<nWells;i++) {
    if((scaleOptWell[i]>=3.0 || scaleOptWell[i]<=0.3334) && nActiveData[i]>0) {
      LogKit::LogFormatted(LogKit::WARNING,"\nWARNING: The well %s has a optimal local gain value indicating that this well should not be used for wavelet estimation\n",
          wells[i]->getWellname().c_str());
    }
  }

  if(doEstimateLocalScale==true) {
// Estimate global noise with local waveletscale
    dataVar = 0.0;
    errStd = 0.0;
    for(int i=0;i<nWells;i++) {
      dataVar+=(dataVarWell[i]*nActiveData[i]);
      errStd+=(errWellOptScale[i]*errWellOptScale[i]*nActiveData[i]);
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

    if(doEstimateLocalNoise) {
      float errStdLN;
      if(doEstimateSNRatio==true)
        errStdLN = errStd;
      else //SNRatio given in model file
        errStdLN = sqrt(dataVar/modelSettings->getSNRatio(number));
      if(gain==NULL && doEstimateLocalScale==false && doEstimateGlobalScale==false) { // No local wavelet scale 
        for(int i=0;i<nWells;i++) 
          errVarWell[i] = sqrt(errVarWell[i]);
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

  if (empSNRatio < 1.1f) {
    LogKit::LogFormatted(LogKit::WARNING,"\nERROR: The empirical signal-to-noise ratio Var(data)/Var(noise) is %.2f. Ratios smaller",empSNRatio);
    LogKit::LogFormatted(LogKit::WARNING,"\n       than 1.1 are not acceptable. The signal-to-noise ratio was not reliably estimated");
    LogKit::LogFormatted(LogKit::WARNING,"\n       and you must give it as input in the model file.\n");
    LogKit::LogFormatted(LogKit::WARNING,"\n       If the wavelet was estimated by CRAVA the solution may be to remove one or more wells");
    LogKit::LogFormatted(LogKit::WARNING,"\n       from the wavelet estimation (compare shifts and SN-ratios for different wells).\n");

    errText += "Invalid signal-to-noise ratio obtained for the angle-gather of "+NRLib::ToString(static_cast<float>(180.0/M_PI)*seisCube->getTheta())+" degrees.\n";
    error += 1;
  }
 
  delete [] alpha;
  delete [] beta;
  delete [] rho;
  delete [] seisData;
  for(int i=0;i<nWells;i++) {
    delete [] cpp_r[i]; 
    delete [] seis_r[i] ;
    delete [] synt_r[i] ;
    delete [] wavelet_r[i];
    delete [] cor_seis_synt_r[i];
  }
  delete [] cpp_r;
  delete [] seis_r;
  delete [] synt_r;
  delete [] wavelet_r;
  delete [] cor_seis_synt_r;
   
  return empSNRatio;
}



float          
Wavelet1D::findOptimalWaveletScale(fftw_real ** synt_seis_r,
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

  for(int i=0;i<nScales;i++) {
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

void 
Wavelet1D::findLocalNoiseWithGainGiven(fftw_real ** synt_seis_r,
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
Wavelet1D::estimateLocalShift(const CovGrid2D  & cov,
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
Wavelet1D::estimateLocalGain(const CovGrid2D  & cov,
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
Wavelet1D::estimateLocalNoise(const CovGrid2D  & cov,
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
Wavelet1D::findBulkShift(fftw_real* vec_r,
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
Wavelet1D::fillInnWavelet(fftw_real* wavelet_r,
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
  //wavelet_r[nzp/2]    = 0;
  //wavelet_r[nzp/2+1]    = 0;
  int counterForWavelet   = 1;

  float w;
  int i;
  for(i=1;i<= nzp/2;i++) {
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
Wavelet1D::shiftOptimal(fftw_real** ccor_seis_cpp_r,
                        float* wellWeight,
                        float* dz,
                        int nWells,
                        int nzp,
                        float* shiftWell,
                        float maxShift)
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
Wavelet1D::multiplyPapolouis(fftw_real** vec, 
                             float* dz,
                             int nWells,
                             int nzp, 
                             float waveletLength,
                             float* wellWeight) const
{
  int i,w;
  float wHL=float( waveletLength/2.0);
  float weight,dist;
  for(w=0;w<nWells;w++)
  {
  if(wellWeight[w] > 0)
    {
      for(i=1;i<nzp;i++)
      {
        dist =MINIM(i,nzp-i)*dz[w];
        if(dist < wHL)
        {
          weight  = float(1.0/PI*fabs(sin(PI*(dist)/wHL))); 
          weight += float((1-fabs(dist)/wHL)*cos(PI*dist/wHL));
        }
        else
          weight=0;

        vec[w][i]*=weight;
      }
    }
  }
}

void
Wavelet1D::getWavelet(fftw_real** ccor_seis_cpp_r,
                      fftw_real** cor_cpp_r,
                      fftw_real** wavelet_r,
                      float* wellWeight,
                      int nWells,
                      int nt)
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

fftw_real* 
Wavelet1D::averageWavelets(fftw_real** wavelet_r,
                           int nWells,
                           int nzp,
                           float* wellWeight,
                           float* dz,
                           float dzOut) const
{
  // assumes dz[w] < dzOut for all w
  fftw_real* wave= static_cast<fftw_real*>(fftw_malloc(rnzp_*sizeof(fftw_real)));  

  int w,i;

  float* weight= new float[nWells];// weight is length of data interval in well
  float sum = 0;
  for(i=0;i<nWells;i++)
    sum+=wellWeight[i];

  for(i=0;i<nWells;i++)
    weight[i]=wellWeight[i]/sum;

  
  for(i=0;i<nzp;i++)
    wave[i] = 0.0; // initialize
  float ww;

  for( w=0;w<nWells;w++)
  {
    if(wellWeight[w]>0)
    {
      wave[0] += weight[w]* wavelet_r[w][0]; // central

      for(i=1;i <nzp/2-1;i++) // positive time
      {
        int ind     = int(floor(i*dzOut/dz[w]));
        float t     = (i*dzOut-ind*dz[w])/dz[w];  // fraction of distance to ind*dz[w]
        if(ind<nzp/2)
            ww= (1-t)* wavelet_r[w][ind]+t*(wavelet_r[w][ind+1]);
        else
          ww=0;
        wave[i]    += weight[w]*ww;

      }
      for(i=1;i < nzp/2-1;i++)// negative time Note dz[w] <= dzOut for all w
      {
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
  fileName = "wavelet_"+NRLib::ToString(int(floor(theta_/PI*180+0.5)))+"_fftOrder_noshift";
  Wavelet::printVecToFile(fileName,wave,nzp_);

  delete [] weight;
  return wave;
}

void 
Wavelet1D::shiftReal(float shift, 
                     fftw_real* rAmp,
                     int nt)
{
  fftw_complex* cAmp = reinterpret_cast<fftw_complex*>(rAmp);
  Utils::fft(rAmp,cAmp, nt);
  int cnzp= nt/2+1;
  float expo;
  fftw_complex tmp,mult;
  for(int i=0;i<cnzp;i++) {
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
Wavelet1D::convolve(fftw_complex* var1_c ,
                    fftw_complex* var2_c, 
                    fftw_complex* out_c,
                    int cnzp) const
{
  for(int i=0;i<cnzp;i++) {
    out_c[i].re = var1_c[i].re*var2_c[i].re+var1_c[i].im*var2_c[i].im; 
    out_c[i].im = var1_c[i].im*var2_c[i].re - var1_c[i].re*var2_c[i].im;
  }
}
