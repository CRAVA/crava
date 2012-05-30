#ifndef CORR_H
#define CORR_H

#include "src/definitions.h"
#include "nrlib/surface/regularsurface.hpp"
#include "src/fftgrid.h"

class ModelSettings;
class Simbox;
class SeismicParametersHolder;
class SpatialWellFilter;
class WellData;

class Corr
{
public:
  Corr(float   ** pointVar0,
       float   ** priorVar0,
       float    * priorCorrT,
       int        n,
       float      dt,
       Surface  * priorCorrXY);

  Corr(int                       n,
       float                     dt,
       SeismicParametersHolder & seismicParameters);

  ~Corr(void);

  float   ** getPriorVar0(void)                const { return priorVar0_                              ;}
  float    * getPriorCorrTFiltered(void)       const { return priorCorrTFiltered_                     ;}
  Surface  * getPriorCorrXY(void)              const { return priorCorrXY_                            ;}

  int        getnx(void)                       const { return static_cast<int>(priorCorrXY_->GetNI()) ;}
  int        getny(void)                       const { return static_cast<int>(priorCorrXY_->GetNJ()) ;}
  int        getn(void)                        const { return n_                                      ;}
  float      getdt(void)                       const { return dt_                                     ;}

  float   ** getPostVar0(void)                 const { return postVar0_                               ;}
  float      getPostCovAlpha00(int k)          const { return postCovAlpha00_[k]                      ;}
  float      getPostCovBeta00(int k)           const { return postCovBeta00_[k]                       ;}
  float      getPostCovRho00(int k)            const { return postCovRho00_[k]                        ;}
  float      getPostCrCovAlphaBeta00(int k)    const { return postCrCovAlphaBeta00_[k]                ;}
  float      getPostCrCovAlphaRho00(int k)     const { return postCrCovAlphaRho00_[k]                 ;}
  float      getPostCrCovBetaRho00(int k)      const { return postCrCovBetaRho00_[k]                  ;}
  FFTGrid *& getPostCovAlpha(void)                   { return postCovAlpha_                           ;}
  FFTGrid *& getPostCovBeta(void)                    { return postCovBeta_                            ;}
  FFTGrid *& getPostCovRho(void)                     { return postCovRho_                             ;}
  FFTGrid *& getPostCrCovAlphaBeta(void)             { return postCrCovAlphaBeta_                     ;}
  FFTGrid *& getPostCrCovAlphaRho(void)              { return postCrCovAlphaRho_                      ;}
  FFTGrid *& getPostCrCovBetaRho(void)               { return postCrCovBetaRho_                       ;}

  void       getNextParameterCovariance(fftw_complex **& parVar);
  void       getNextErrorVariance(fftw_complex **& errVar,
                                  fftw_complex * errMult1,
                                  fftw_complex * errMult2,
                                  fftw_complex * errMult3,
                                  int              ntheta,
                                  float            wnc,
                                  double ** errThetaCov,
                                  bool invert_frequency);

  fftw_real* initializeCorrelations(SpatialWellFilter     * spatwellfilter,
                                    std::vector<WellData *> wells,
                                    float                   corrGradI,
                                    float                   corrGradJ,
                                    int                     lowIntCut,
                                    int                     nWells,
                                    int                     nz,
                                    int                     nzp);
  void       terminateAccess(void);
  void       initializeAccess(void);

  void       createPostVariances(void);

  void       computeCircCorrT(fftw_real* CircCorrT, int nzp);
  void       makeCircCorrTPosDef(fftw_real* CircCorrT, int minIntFq, int nzp);

  void       setPriorVar0(float ** priorVar0);
  void       setPriorCorrTFiltered(float * corrT, int nz, int nzp);

  void       createPostGrids(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid);

  void       FFT(void);              // Transform all posterior grids to fourier domain
  void       invFFT(void);           // Transform all posterior grids to time domain

  void       printPriorVariances(void) const;
  void       printPostVariances(void) const;

  void       writeFilePriorCorrT(float* corrT, int nzp) const;
  void       writeFilePriorVariances(ModelSettings * modelSettings) const;
  void       writeFilePostVariances(void) const;
  void       writeFilePostCovGrids(Simbox const * simbox) const;

private:
  FFTGrid  * createFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid);
  void       writeFilePostCorrT(float * postCov, int nz, const std::string & subDir, const std::string & fileName) const;
  float      getOrigin(FFTGrid * grid) const;
  float    * createPostCov00(FFTGrid * postCov);
  void       createPriorVar0(void);

  fftw_real* extractParamCorr(FFTGrid * covAlpha, int nzp);

  float   ** pointVar0_;             // Point variance calculated from using well log resolution

  float   ** priorVar0_;             // Blocked variance this is used in the computations
  Surface  * priorCorrXY_;
  float    * priorCorrT_;
  float    * priorCorrTFiltered_;

  int        n_;                     // priorCorrT - length
  float      dt_;                    // priorCorrT - time step

  bool       common_correlation_;    // True:  Correlations made by multiplying pointwise variances with common correlation grids
                                     // False: Get correlations directly from covariance grids


  float   ** postVar0_;

  float    * postCovAlpha00_;        // Posterior covariance in (i,j) = (0,0)
  float    * postCovBeta00_;
  float    * postCovRho00_;
  float    * postCrCovAlphaBeta00_;
  float    * postCrCovAlphaRho00_;
  float    * postCrCovBetaRho00_;

  FFTGrid  * postCovAlpha_;
  FFTGrid  * postCovBeta_;
  FFTGrid  * postCovRho_;
  FFTGrid  * postCrCovAlphaBeta_;
  FFTGrid  * postCrCovAlphaRho_;
  FFTGrid  * postCrCovBetaRho_;
  FFTGrid  * errCorr_;

};
#endif
