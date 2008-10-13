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
  Wavelet3D(char          * fileName, 
            ModelSettings * modelSettings, 
            Simbox        * simBox, 
            float           theta, 
            float         * reflCoef,  
            int           & errCode, 
            char          * errText);
  Wavelet3D(Wavelet * wavelet, int difftype);
  Wavelet3D(Wavelet * wavelet);
  virtual ~Wavelet3D() {}

  bool           consistentSize(int nzp, int nyp, int nxp) const; 
  fftw_complex   getCAmp(int k, int j, int i) const;
  fftw_real      getRAmp(int k, int j, int i);
  fftw_complex   getCAmp(int k, float scale, int j, int i) const;
  void           setRAmp(float value, int k, int j, int i);
  void           setCAmp(fftw_complex value, int k, int j, int i);
  void           scale(float gain);
  int            getNx() const {return nx_;}
  int            getNy() const {return ny_;}
  int            getNxp() const {return nxp_;}
  int            getNyp() const {return nyp_;}
  float          getDx() const {return dx_;}
  float          getDy() const {return dy_;}
  FFTGrid      * getAmpCube() {return &ampCube_;}
  void           multiplyByR(float p);

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
