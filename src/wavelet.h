/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef WAVELET_H
#define WAVELET_H

#include "nrlib/surface/regularsurface.hpp"

#include "fftw.h"
#include "src/definitions.h"
#include "lib/utils.h"

#include "src/seismicstorage.h"
#include "src/blockedlogscommon.h"

class Simbox;
class FFTGrid;
class ModelSettings;
class Wavelet1D;

class Wavelet {
public:
  enum           difftypes{FIRSTORDERFORWARDDIFF, FIRSTORDERCENTRALDIFF, FIRSTORDERBACKWARDDIFF, FOURIER};
  enum           waveletDims{ONE_D = 0, THREE_D = 1};
  enum           waveletreadtypes{JASON, NORSAR};

//Constructors and destructor
  Wavelet(int                 dim);

  Wavelet(int                 dim,
          Wavelet           * wavelet);

  Wavelet(Wavelet              * wavelet ,
          int              /*   difftype */);

  Wavelet(int                /*  difftype*/,
          int                /*  nz */ ,
          int                 /* nzp */);

  Wavelet(const std::string & fileName,
          int                 fileFormat,
          const ModelSettings     * modelSettings,
          const float             * reflCoef,
          float               theta,
          int                 dim,
          int               & errCode,
          std::string       & errText);


  Wavelet(const ModelSettings     * modelSettings,
          const float             * reflCoef,
          float               theta,
          int                 dim,
          float               peakFrequency,
          int               & errCode);

  virtual ~Wavelet();

// Fast Fourier Transform for inplace storage
  void          fft1DInPlace();
  void          invFFT1DInPlace();

// Access methods for wavelet values
  fftw_real     getRAmp(int k);
  fftw_real*    getRAmp(){ return rAmp_;}

  fftw_complex  getCAmp(int k) const;

  fftw_complex  getCAmp(int   k,
                        float scale) const;

  void          setRAmp(float value,
                        int   k);


  void          setCAmp(fftw_complex value,
                        int k);

  void          resample(float dz,
                         int   nz,
                         int   nzp);

  void          shiftFromFFTOrder();

  void          multiplyRAmpByConstant(float c);
  void          setupAsVector(int nz, int nzp);


  void          scale(float gain);

  //Note: Function below is mainly controlled by debugflag. Set overrideDebug = true to force.
  void          printToFile(const std::string                 & fileName,
                            bool                                overrideDebug = false);

  void          writeWaveletToFile(const std::string          & fileName,
                                   float                        approxDz,
                                   bool     makePrintedWaveletIntegralZero);

  void          setShiftGrid(Grid2D                           * grid);

  void          setGainGrid(Grid2D                            * grid);

  float         getNorm()     const {return norm_;}
  void          setNorm(float norm) {norm_ = norm;}
  bool          getIsReal()   const {return(isReal_);}
  int           getDim()      const {return dim_;}
  int           getNz()       const {return nz_;}
  int           getNzp()      const {return nzp_;}
  float         getDz()       const {return dz_;}
  float         getScale()    const {return scale_;}
  virtual float getLocalStretch(int /*i*/,
                                int /*j*/) {return 1.0f;} // note Not robust towards padding


  virtual Wavelet1D * createLocalWavelet1D(int /*i*/,
                                           int /*j*/) {return 0;} // note Not robust towards padding
  virtual Wavelet1D * createWavelet1DForErrorNorm() {return 0;}
  virtual Wavelet1D * getGlobalWavelet(){return 0;}



  virtual float findGlobalScaleForGivenWavelet(const ModelSettings         * /*modelSettings*/,
                                               const Simbox                * /*simbox*/,
                                               SeismicStorage              * /*seismic_data*/,
                                               const std::map<std::string, BlockedLogsCommon *> & /*mapped_blocked_logs*/) {return 1.0f;}

  // for noise estimation
  virtual float calculateSNRatioAndLocalWavelet(const Simbox          * /*estimation_simbox*/,
                                                const Simbox          * /*inversion_simbox*/,
                                                const std::vector<std::vector<double> >          & /*seis_logs*/,
                                                const std::map<std::string, BlockedLogsCommon *> & /*mapped_blocked_logs*/,
                                                const ModelSettings   * /*modelSettings*/,
                                                std::string           & /*errText*/,
                                                int                   & /*error*/,
                                                int                     /*number*/,
                                                Grid2D               *& /*noiseScaled*/,
                                                Grid2D               *& /*shift*/,
                                                Grid2D               *& /*gain*/,
                                                float                   /*SNRatio*/,
                                                float                   /*waveletScale*/,
                                                bool                    /*estimateSNRatio*/,
                                                bool                    /*estimateGlobalScale*/,
                                                bool                    /*estimateLocalNoise*/,
                                                bool                    /*estimateLocalShift*/,
                                                bool                    /*estimateLocalScale*/,
                                                bool                    /*estimateWavelet*/) {return 1.0f;}

  virtual float calculateSNRatio(const Simbox                             * /*simbox*/,
                                 SeismicStorage                           * /*seismic_data*/,
                                 const std::map<std::string, BlockedLogsCommon *> & /*mapped_blocked_logs*/,
                                 const ModelSettings                      * /*modelSettings*/,
                                 std::string                              & /*errText*/,
                                 int                                      & /*error*/,
                                 const NRLib::Grid2D<float>               & /*refTimeGradX*/,
                                 const NRLib::Grid2D<float>               & /*refTimeGradY*/,
                                 const std::vector<std::vector<double> >  & /*tGradX*/,
                                 const std::vector<std::vector<double> >  & /*tGradY*/,
                                 int                                        /*i*/,
                                 float                                      /*SNRatio*/,
                                 bool                                       /*estimateSNRatio*/,
                                 bool                                       /*estimateWavelet*/) {return 1.0f;}

  float          findNormWithinFrequencyBand(float loCut ,float hiCut ) const;
  void           nullOutsideFrequencyBand(float loCut ,float hiCut );
  float          findNorm() const;
  void           SetReflectionCoeffs(const NRLib::Matrix & reflCoef, int i);
  bool           getInFFTOrder()     const { return inFFTorder_ ;}
  //Grid2D       * getGainGrid()       const { return gainGrid_   ;}
  //Grid2D       * getShiftGrid()      const { return shiftGrid_  ;}

protected:
  float          getTheta()          const {return theta_;}
  int            getCz()             const {return cz_;}
  float          getWaveletLength()  const {return waveletLength_;}

  void           doLocalShiftAndScale1D(Wavelet1D* localWavelet,// wavelet to shift and scale
                                        int                   i,
                                        int                   j);

  float          getLocalTimeshift(int                          i,
                                   int                          j) const;

  float          getLocalGainFactor(int                         i,
                                    int                         j) const;

  float          findWaveletLength(float                        minRelativeAmp,float minimumLength);

  void           convolve(fftw_complex                       * var1_c,
                          fftw_complex                       * var2_c,
                          fftw_complex                       * out_c,
                          int                                  cnzp)           const;

  void           printVecToFile(const std::string                       & fileName,
                                fftw_real                               * vec ,
                                int                                       nzp) const;

  void           printVecDoubleToFile(const std::string                       & fileName,
                                      double                                  * vec,
                                      int                                       nzp) const;


  fftw_real*     averageWavelets(const std::vector<std::vector<float> > & wavelet_r,
                                 int                                      nWells,
                                 int                                      nzp,
                                 const std::vector<float>               & wellWeight,
                                 const std::vector<float>               & dz,
                                 float                                    dzOut)         const;

  void          fillInnWavelet(fftw_real                     * wavelet_r,
                               int                             nzp,
                               float                           dz);


  float         findBulkShift(fftw_real                      * vec_r,
                              float                            dz,
                              int                              nzp,
                              float                            maxShift);

  void           shiftReal(float                               shift,
                           fftw_real                         * rAmp,
                           int                                 nt);



  double         Ricker(double t, float peakF);

  int            cnzp_;                  // size in z direction for storage inplace algorithm (complex grid) nzp_/2+1
  int            rnzp_;                  // expansion in z direction for storage inplace algorithm (real grid) 2*(nzp_/2+1)
  fftw_real*     rAmp_;                  // The amplitude of the wavelet
  fftw_complex*  cAmp_;                  // The amplitude of the wavelet complex (if fourier transformed )

  float          theta_;                 // the reflection angle that the wavelet correspond to

  float          dz_;                    // Sampling interval of wavelet, unit [ ms ]
  int            nz_;                    // length of wavelet
  int            nzp_;                   // length of padded wavelet
  int            cz_;                    // position of central point
  int            formats_;                // formats for output of wavelet
  bool           inFFTorder_;            // is true if the wavelet is ordred with the central point at the start
                                             // false if the central point is in the middle
  bool           isReal_;                // is true if the wavlet is real, false if it is fourier transformed
  float          norm_;                  // The (vector) norm of the wavelet (not function norm that divides by dz)
  float          waveletLength_;         // Length of wavelet estimated as is amplitudes larger than 1/1000 * max amplitude

  float          coeff_[3];              //Reflection coefficients.

  const int      dim_;

//NBNB The following parameters are NOT copied in copy constructor.
  float          scale_;
  Grid2D       * shiftGrid_;             // 2D grid of shift
  Grid2D       * gainGrid_;              // 2D grid of gain factors.

private:
  float          getWaveletValue(float                           z,
                                 float                         * Wavelet,
                                 int                             center,
                                 int                             nx,
                                 float                           dz);

  float          getArrayValueOrZero(int                        i,
                                     float                    * Wavelet,
                                     int                        nz) const;



  void           WaveletReadJason(const std::string           & fileName,
                                  int                         & errCode,
                                  fftw_real                  *& rAmp,
                                  fftw_complex               *& cAmp,
                                  float                       & dz,
                                  int                         & nz,
                                  int                         & cz,
                                  int                         & nzp,
                                  int                         & cnzp,
                                  int                         & rnzp,
                                  float                       & norm,
                                  std::string                 & errText) const;

  void           WaveletReadNorsar(const std::string          & fileName,
                                   int                        & errCode,
                                   std::string                & errText);

};

#endif
