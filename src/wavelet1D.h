#ifndef WAVELET1D_H
#define WAVELET1D_H

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/wavelet.h"

class Wavelet1D : public Wavelet {
public:
  void           fft1DInPlace();
  void           invFFT1DInPlace();

  //Constructors and destructor
  Wavelet1D(Simbox         * simbox, 
            FFTGrid        * seisCube, 
            WellData      ** wells, 
            Surface       ** estimInterval,
            ModelSettings  * modelSettings, 
            float          * reflCoef);
  Wavelet1D(char          * fileName, 
            int             fileFormat, 
            ModelSettings * modelSettings, 
            float         * reflCoef,
            int           & errCode, 
            char          * errText);
  Wavelet1D(Wavelet * wavelet);
  Wavelet1D(Wavelet * wavelet,int difftype);
  Wavelet1D(int difftype, int nz, int nzp);
  virtual ~Wavelet1D();

  void           resample(float dz, int nz, float pz,float theta);
  bool           consistentSize(int nzp, int, int) const;
  fftw_complex   getCAmp(int k, int j=0, int i=0) const;
  fftw_real      getRAmp(int k, int j=0, int i=0);
  fftw_complex   getCAmp(int k, float scale, int j=0, int i=0) const;
  void           setRAmp(float value, int k, int j=0, int i=0);
  void           scale(float gain);
  void           printToFile(char* fileName, bool overrideDebug = false);
  void           writeWaveletToFile(char* fileName, float approxDz, Simbox *simbox = NULL);
  void           write1DWLas3DWL();
  void           write3DWLfrom1DWL();

private:
  void           flipUpDown();
  float          getWaveletValue(float z, float * Wavelet, int center,int nx, float dz);
  void           shiftAndScale(float shift,float gain);
  void           WaveletReadOld(char * fileName, int &errCode, char *errText);
  void           WaveletReadJason(char * fileName, int &errCode, char *errText);
  float          shiftOptimal(fftw_real** ccor_seis_cpp_r,float* wellWeight,float* dz,int nWells,int nzp,float* shiftWell);
  void           multiplyPapolouis(fftw_real** vec, float* dz,int nWells,int nzp, float waveletLength, float *wellWeight) const;
  void           getWavelet(fftw_real** ccor_seis_cpp_r,fftw_real** cor_cpp_r,fftw_real** wavelet_r,float* wellWeight,int nWells,int nt);
  fftw_real*     averageWavelets(fftw_real** wavelet_r,int nWells,int nzp,float* wellWeight,float* dz,float dzOut) const;
  float          getArrayValueOrZero(int i ,float * Wavelet, int nz) const;
//  void           fillInnWavelet(fftw_real* wavelet_r,int nzp,float dz);
  
  int	         cnzp_;	                 // size in z direction for storage inplace algorithm (complex grid) nzp_/2+1
  int	         rnzp_;                  // expansion in z direction for storage inplace algorithm (real grid) 2*(nzp_/2+1)
  fftw_real*     rAmp_;                  // The amplitude of the wavelet  
  fftw_complex*  cAmp_;                  // The amplitude of the wavelet complex (if fourier transformed )
};

#endif
