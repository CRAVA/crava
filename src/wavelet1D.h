#ifndef WAVELET1D_H
#define WAVELET1D_H

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/wavelet.h"

class Wavelet1D : public Wavelet {
public:
//  void fft1DInPlace();

  //Constructors and destructor
  Wavelet1D(Simbox * simbox, FFTGrid * seisCube, WellData ** wells, ModelSettings * modelSettings, float * coef, int dim = 1);
  Wavelet1D(char * fileName, ModelSettings * modelSettings, int fileFormat, int dim = 1);
  Wavelet1D(Wavelet * wavelet, int dim = 1);
  Wavelet1D(Wavelet * wavelet,int difftype, int dim = 1);
  Wavelet1D(int difftype, int nz, int nzp, int dim = 1);
  virtual ~Wavelet1D() {}

  void           resample(float dz, int nz, float pz,float theta);
  bool           consistentSize(int nzp, int, int) const; 

private:
  void           flipUpDown();
  float          getWaveletValue(float z, float * Wavelet, int center,int nx, float dz);
  int            getWaveletLengthI();
  float          getWaveletLengthF();
  void           WaveletReadOld(char * fileName);
  void           WaveletReadJason(char * fileName);
  float          shiftOptimal(fftw_real** ccor_seis_cpp_r,float* wellWeight,float* dz,int nWells,int nzp,float* shiftWell);
  void           multiplyPapolouis(fftw_real** vec, float* dz,int nWells,int nzp, float waveletLength) const;
  void           getWavelet(fftw_real** ccor_seis_cpp_r,fftw_real** cor_cpp_r,fftw_real** wavelet_r,float* wellWeight,int nWells,int nt);
  fftw_real*     averageWavelets(fftw_real** wavelet_r,int nWells,int nzp,float* wellWeight,float* dz,float dzOut) const;
  float          getArrayValueOrZero(int i ,float * Wavelet, int nz) const;
};

#endif
