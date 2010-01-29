#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fft/include/fftw.h"

#include "lib/global_def.h"

#include "nrlib/surface/regularsurfacerotated.hpp"

#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/waveletfilter.h"
#include "src/fftgrid.h"

class BlockedLogs;

class Wavelet3D : public Wavelet {
public:
  void           fft1DInPlace();
  void           invFFT1DInPlace();
  //Constructors and destructor
  Wavelet3D(const std::string            & filterFile,
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
            std::string                  & errText);
  
  Wavelet3D(Wavelet1D                    * wavelet1d,
            const std::string            & filterfile,
            ModelSettings                * modelSettings,
            int                            angle_index,
            Simbox                       * simBox,
            float                          theta,
            int                          & errCode,
            std::string                  & errText);

  Wavelet3D(Wavelet                      * wavelet, 
            int                            difftype);

  Wavelet3D(Wavelet *                      wavelet);

  virtual ~Wavelet3D(); 

  bool           consistentSize(int nzp, int nyp, int nxp)          const; 
  fftw_complex   getCAmp(int k, int j, int i)                       const;
  fftw_real      getRAmp(int k, int j, int i);
  fftw_complex   getCAmp(int k, float scale, int j, int i)          const;
  void           setRAmp(float value, int k, int j, int i);
  void           setCAmp(fftw_complex value, int k, int j, int i);
  void           scale(float gain);
  int            getNx()                                            const {return nx_;}
  int            getNy()                                            const {return ny_;}
  int            getNxp()                                           const {return nxp_;}
  int            getNyp()                                           const {return nyp_;}
  float          getDx()                                            const {return dx_;}
  float          getDy()                                            const {return dy_;}
  Wavelet1D    * getWavelet1D()                                     const {return wavelet1D_;}
  WaveletFilter  getFilter()                                        const {return filter_;}
  FFTGrid      * getAmpCube()                                             {return &ampCube_;}
  void           multiplyByR(float p);

  void           writeWaveletToFile(const std::string & fileName, float approxDzIn);
  void           printToFile(const std::string & filename, bool overrideDebug = false);
  
private:
  void           shiftFFTGrid(FFTGrid *shiftAmp);
  void           findLayersWithData(const std::vector<Surface *> & estimInterval,
                                    BlockedLogs                  * bl,
                                    FFTGrid                      * seisCube,
                                    float                        * az,
                                    float                        * bz,
                                    bool                         * hasWellData) const;
  double         findPhi(float a, float b)                                    const;
  double         findPsi(float radius, float kz)                                const;
  fftw_complex   findWLvalue(Wavelet1D       * wavelet1d,
                             float             omega)                           const;

  FFTGrid        ampCube_;
  float          dy_, dx_;
  int            ny_, nx_;
  int            nyp_, nxp_;
  Wavelet1D    * wavelet1D_;
  WaveletFilter  filter_;

};

#endif
