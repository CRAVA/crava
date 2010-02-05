#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fft/include/fftw.h"

#include "lib/global_def.h"

#include "nrlib/surface/regularsurfacerotated.hpp"

#include "src/wavelet.h"
#include "src/waveletfilter.h"

class BlockedLogs;

class Wavelet3D : public Wavelet {
public:
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
            int                          & errCode,
            std::string                  & errText);
  
  Wavelet3D(const std::string            & fileName, 
            int                            fileFormat, 
            ModelSettings                * modelSettings, 
            float                        * reflCoef,
            float                          theta,
            int                          & errCode, 
            std::string                  & errText,
            const std::string            & filterFile);

  Wavelet3D(Wavelet *                      wavelet);

  virtual ~Wavelet3D(); 

  WaveletFilter  getFilter() const {return filter_;}
  
private:
  void          findLayersWithData(const std::vector<Surface *> & estimInterval,
                                   BlockedLogs                  * bl,
                                   FFTGrid                      * seisCube,
                                   float                        * az,
                                   float                        * bz,
                                   bool                         * hasWellData) const;

  double         findPhi(float                                    a, 
                         float                                    b)           const;

  double         findPsi(float                                    r)           const;

  fftw_complex   findWLvalue(float                                omega)       const;

  WaveletFilter  filter_;
};

#endif
