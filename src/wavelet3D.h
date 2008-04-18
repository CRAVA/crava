#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/wavelet.h"

class Wavelet3D : public Wavelet {
public:
 // void fft1DInPlace();

  //Constructors and destructor
  Wavelet3D(char * fileName, ModelSettings * modelSettings, int dim = 3);
  Wavelet3D(Wavelet * wavelet, int dim = 3);
//  Wavelet(Wavelet * wavelet,int difftype);
//  Wavelet(int difftype, int nz, int nzp);
  virtual ~Wavelet3D() {}

  void           resample(float dz, int nz, float pz,float theta);
  bool           consistentSize(int nzp, int nyp, int nxp) const; 
  
private:
  void			     WaveletReadSgri(char *fileName);
  float          getWaveletValue(float, float *, int, int, float) {return(0.0);}

  float          dy_, dx_;
  int            ny_, nx_;
  int            nyp_, nxp_;
};

#endif
