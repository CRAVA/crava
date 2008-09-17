#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/wavelet.h"
#include "src/fftgrid.h"
#include "lib/sgri.h"

class Wavelet3D : public Wavelet {
public:
  void           fft1DInPlace();
  void           invFFT1DInPlace();
  //Constructors and destructor
  Wavelet3D(char * fileName, ModelSettings * modelSettings, Simbox *simBox, float theta, int &errCode, char *errText, int dim = 3);
  Wavelet3D(Wavelet * wavelet, int dim = 3);
  virtual ~Wavelet3D() {}

  bool           consistentSize(int nzp, int nyp, int nxp) const; 
  fftw_complex   getCAmp(int k, int j, int i) const;
  fftw_real      getRAmp(int k, int j, int i);
  fftw_complex   getCAmp(int k, float scale, int j, int i) const;
  void           setRAmp(float value, int k, int j, int i);
  void           scale(float gain);

  void           writeWaveletToFile(char* fileName, float, Simbox *simbox);
  void           printToFile(char* fileName, bool overrideDebug = false);
  
private:
  void           shiftFFTGrid(FFTGrid *shiftAmp);

  FFTGrid        ampCube_;
  float          dy_, dx_;
  int            ny_, nx_;
  int            nyp_, nxp_;
};

#endif
