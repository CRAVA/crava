#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fftw.h"

#include "nrlib/surface/regularsurfacerotated.hpp"

#include "src/wavelet.h"
#include "src/waveletfilter.h"

class BlockedLogs;
class Wavelet1D;

class Wavelet3D : public Wavelet {
public:
  //Constructors and destructor
  Wavelet3D(const std::string                          & filterFile,
            const std::vector<Surface *>               & estimInterval,
            const NRLib::Grid2D<float>                 & refTimeGradX,
            const NRLib::Grid2D<float>                 & refTimeGradY,
            const std::vector<std::vector<double> >    & tGradX,
            const std::vector<std::vector<double> >    & tGradY,
            FFTGrid                                    * seisCube,
            ModelSettings                              * modelSettings,
            WellData                                  ** wells,
            Simbox                                     * simBox,
            float                                      * reflCoef,
            int                                          angle_index,
            int                                        & errCode,
            std::string                                & errText);

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
  static void setGradientMaps( NRLib::Grid2D<float> gradX , NRLib::Grid2D<float> gradY ){gradX_=gradX; gradY_=gradY;}

  // Methods that are virtual in Wavelet
//  float         calculateSNRatioAndLocalWavelet(Simbox        * /*simbox*/,
//                                                FFTGrid       * /*seisCube*/,
//                                                WellData     ** /*wells*/,
//                                                ModelSettings * /*modelSettings*/,
//                                                std::string   & errText,
//                                                int           & error,
//                                                int             /*number*/,
//                                                Grid2D       *& /*noiseScaled*/,
//                                                Grid2D       *& /*shift*/,
//                                                Grid2D       *& /*gain*/);
Wavelet1D*  createWavelet1DForErrorNorm(void);
Wavelet1D * createLocalWavelet1D( int i, int j);
float       getLocalStretch(int i,int j);

Wavelet1D*  getSourceWavelet();
Wavelet1D*  extractLocalWaveletByDip1D(double phi, double psi);
void        dipAdjustWavelet(Wavelet1D* Wavelet, double phi, double psi);
float      GetLocalDepthGradientX(int i, int j){ return gradX_(i,j);}
float      GetLocalDepthGradientY(int i, int j){ return gradY_(i,j);}

  float         calculateSNRatio(Simbox                                   * simbox,
                                 FFTGrid                                  * seisCube,
                                 WellData                                ** wells,
                                 ModelSettings                            * modelSettings,
                                 std::string                              & errText,
                                 int                                      & error,
                                 const NRLib::Grid2D<float>               & refTimeGradX,
                                 const NRLib::Grid2D<float>               & refTimeGradY,
                                 const std::vector<std::vector<double> >  & tGradX,
                                 const std::vector<std::vector<double> >  & tGradY,
                                 int                                        number);

private:
  void          findLayersWithData(const std::vector<Surface *> & estimInterval,
                                   BlockedLogs                  * bl,
                                   FFTGrid                      * seisCube,
                                   const std::vector<float>     & az,
                                   const std::vector<float>     & bz,
                                   std::vector<bool>            & hasWellData) const;

  double         findPhi(float                                    a,
                         float                                    b)           const;

  double         findPsi(float                                    r)           const;

  fftw_complex   findWLvalue(float                                omega)       const;

  float          calculateWellWeight(int                                      nWl,
                                     int                                      nPoints,
                                     const std::vector<std::vector<float> > & gMat,
                                     const std::vector<float>               & wellWavelet,
                                     const std::vector<float>               & dVec) const;

  std::vector<fftw_real> adjustCpp(BlockedLogs              * bl,
                                   const std::vector<float> & az,
                                   const std::vector<float> & bz,
                                   std::vector<float>       & Halpha,
                                   int                        start,
                                   int                        length,
                                   const std::string        & wellname,
                                   const std::string        & angle) const;

  void           calculateGradients(BlockedLogs                * bl,
                                    const std::vector<int>     & iPos,
                                    const std::vector<int>     & jPos,
                                    const NRLib::Grid2D<float> & refTimeGradX,
                                    const NRLib::Grid2D<float> & refTimeGradY,
                                    const std::vector<double>  & tGradX,
                                    const std::vector<double>  & tGradY,
                                    float                        v0,
                                    std::vector<float>         & az,
                                    std::vector<float>         & bz,
                                    std::vector<float>         & at0,
                                    std::vector<float>         & bt0) const;



  std::vector<fftw_real> calculateWellWavelet(const std::vector<std::vector<float> > & gMat,
                                              const std::vector<float>               & dVec,
                                              int                                      nWl,
                                              int                                      nhalfWl,
                                              int                                      nPoints) const;

  void           printMatToFile(const std::string                       & fileName,
                                const std::vector<std::vector<float> >  & mat,
                                int                                       n,
                                int                                       m) const;

  WaveletFilter  filter_;
  static NRLib::Grid2D<float>  gradX_;
  static NRLib::Grid2D<float>  gradY_;
};

#endif
