/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef WAVELET3D_H
#define WAVELET3D_H

#include "fftw.h"

#include "nrlib/surface/regularsurfacerotated.hpp"

#include "src/wavelet.h"
#include "src/waveletfilter.h"
#include "src/wavelet1D.h"

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
            const FFTGrid                              * seisCube,
            const ModelSettings                        * modelSettings,
            const std::vector<WellData *>              & wells,
            const Simbox                               * simBox,
            const float                                * reflCoef,
            int                                          angle_index,
            int                                        & errCode,
            std::string                                & errText);

  Wavelet3D(const std::string            & fileName,
            int                            fileFormat,
            const ModelSettings          * modelSettings,
            const float                  * reflCoef,
            float                          theta,
            int                          & errCode,
            std::string                  & errText,
            const std::string            & filterFile);

  Wavelet3D(Wavelet *                      wavelet);

  virtual ~Wavelet3D();

  WaveletFilter          getFilter(void) const {return filter_;}
  static void            setGradientMaps( NRLib::Grid2D<float> gradX , NRLib::Grid2D<float> gradY ){structureDepthGradX_=gradX; structureDepthGradY_=gradY;}
  Wavelet1D            * createWavelet1DForErrorNorm(void);
  Wavelet1D            * createLocalWavelet1D( int i, int j);
  Wavelet1D            * getGlobalWavelet(){ return averageWavelet_;}
  float                  getLocalStretch(int i,int j);

  Wavelet1D            * createSourceWavelet();
  Wavelet1D            * createAverageWavelet(const Simbox * simBox);
  Wavelet1D            * extractLocalWaveletByDip1D(double phi, double psi);
  void                   dipAdjustWavelet(Wavelet1D* Wavelet, double phi, double psi);
  float                  GetLocalDepthGradientX(int i, int j){ return structureDepthGradX_(i,j);}
  float                  GetLocalDepthGradientY(int i, int j){ return structureDepthGradY_(i,j);}

  float                  calculateSNRatio(const Simbox                             * simbox,
                                          const FFTGrid                            * seisCube,
                                          const std::vector<WellData *>            & wells,
                                          const ModelSettings                      * modelSettings,
                                          std::string                              & errText,
                                          int                                      & error,
                                          const NRLib::Grid2D<float>               & refTimeGradX,
                                          const NRLib::Grid2D<float>               & refTimeGradY,
                                          const std::vector<std::vector<double> >  & tGradX,
                                          const std::vector<std::vector<double> >  & tGradY,
                                          int                                        number,
                                          float                                      SNRatio,
                                          bool                                       estimateSNRatio,
                                          bool                                       estimateWavelet);

private:
  void                   findLayersWithData(const std::vector<Surface *> & estimInterval,
                                            BlockedLogs                  * bl,
                                            const FFTGrid                * seisCube,
                                            const std::vector<float>     & az,
                                            const std::vector<float>     & bz,
                                            std::vector<bool>            & hasWellData) const;

  double                 findPhi(float a,
                                 float b) const;

  double                 findPsi(float r) const;

  fftw_complex           findWLvalue(float omega) const;

  float                  calculateWellWeight(int                                      nWl,
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

  void                   calculateGradients(BlockedLogs                * bl,
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
                                              double                                   SNR,
                                              int                                      nWl,
                                              int                                      nhalfWl,
                                              int                                      nPoints) const;

  void                   printMatToFile(const std::string                       & fileName,
                                        const std::vector<std::vector<float> >  & mat,
                                        int                                       n,
                                        int                                       m) const;

  WaveletFilter                  filter_;
  Wavelet1D                    * averageWavelet_;
  static NRLib::Grid2D<float>    structureDepthGradX_;// gradX_
  static NRLib::Grid2D<float>    structureDepthGradY_;// gradY_
};

#endif
