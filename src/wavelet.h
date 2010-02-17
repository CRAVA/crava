#ifndef WAVELET_H
#define WAVELET_H

#include "nrlib/surface/regularsurface.hpp"

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/definitions.h"
#include "lib/utils.h"

class Simbox;
class FFTGrid;
class WellData;
class ModelSettings;

class Wavelet {
public:
  enum           difftypes{FIRSTORDERFORWARDDIFF, FIRSTORDERCENTRALDIFF, FIRSTORDERBACKWARDDIFF, FOURIER};
  enum           waveletDims{ONE_D = 0, THREE_D = 1};
  enum           waveletreadtypes{OLD, JASON, NORSAR};

//Constructors and destructor
  Wavelet(int                 dim);

  Wavelet(int                 dim, 
          Wavelet           * wavelet);
  
  Wavelet(const std::string & fileName, 
          int                 fileFormat, 
          ModelSettings     * modelSettings, 
          float             * reflCoef,
          float               theta,
          int                 dim,
          int               & errCode, 
          std::string       & errText);
  
  Wavelet(Wavelet           * wavelet,
          int                 difftype);

  Wavelet(int                 difftype, 
          int                 nz, 
          int                 nzp);

  virtual ~Wavelet();
  
// Fast Fourier Transform for inplace storage     
  void          fft1DInPlace();
  void          invFFT1DInPlace();

// Access methods for wavelet values
  fftw_real     getRAmp(int                                     k);

  fftw_complex  getCAmp(int                                     k) const;
  
  fftw_complex  getCAmp(int                                     k, 
                        float                                   scale) const;
  
  void          setRAmp(float                                   value, 
                        int                                     k);

  Wavelet *     getLocalWavelet(int                             i,
                                int                             j);

  void          resample(float                                  dz, 
                         int                                    nz, 
                         float                                  pz,
                         bool                                   flip);

  bool          consistentSize(int                              nzp) const;

  void          multiplyRAmpByConstant(float                    c);

  void          scale(float                                     gain);

  //Note: Function below is mainly controlled by debugflag. Set overrideDebug = true to force.
  void          printToFile(const std::string                 & fileName, 
                            bool                                overrideDebug = false);

  void          writeWaveletToFile(const std::string          & fileName, 
                                   float                        approxDz);

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


  virtual float findGlobalScaleForGivenWavelet(ModelSettings * /*modelSettings*/, 
                                               Simbox        * /*simbox*/, 
                                               FFTGrid       * /*seisCube*/, 
                                               WellData     ** /*wells*/) {return 0.0f;}

  // for noise estimation
  virtual float calculateSNRatioAndLocalWavelet(Simbox        * /*simbox*/, 
                                                FFTGrid       * /*seisCube*/, 
                                                WellData     ** /*wells*/, 
                                                Grid2D       *& /*shift*/, 
                                                Grid2D       *& /*gain*/, 
                                                ModelSettings * /*modelSettings*/,
                                                std::string   & /*errText*/, 
                                                int           & /*error*/,
                                                Grid2D       *& /*noiseScaled*/, 
                                                int             /*number*/, 
                                                float           /*globalScale*/) {return 0.0f;} 

protected:
  float          getTheta()          const {return theta_;}
  int            getCz()             const {return cz_;}
  bool           getInFFTOrder()     const {return inFFTorder_;}
  float          getWaveletLength()  const {return waveletLength_;}

  void           shiftAndScale(float                            shift,
                               float                            gain);

  float          findWaveletLength(float                        minRelativeAmp);

  
  fftw_real*     averageWavelets(const std::vector<std::vector<float> > & wavelet_r,
                                 int                                      nWells,
                                 int                                      nzp,
                                 const std::vector<float>               & wellWeight,
                                 const std::vector<float>               & dz,
                                 float                                    dzOut)         const;

  void           printVecToFile(const std::string                       & fileName, 
                                fftw_real                               * vec ,
                                int                                       nzp) const;

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
  
  void           flipUpDown();
  
  float          getArrayValueOrZero(int                        i,
                                     float                    * Wavelet, 
                                     int                        nz) const;
  
  float          getLocalTimeshift(int                          i, 
                                   int                          j) const;

  float          getLocalGainFactor(int                         i, 
                                    int                         j) const;

  void           WaveletReadOld(const std::string             & fileName,
                                int                           & errCode, 
                                std::string                   & errText);

  void           WaveletReadJason(const std::string           & fileName,
                                  int                         & errCode, 
                                  std::string                 & errText);


};

#endif
