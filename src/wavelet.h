#ifndef WAVELET_H
#define WAVELET_H

#include "nrlib/surface/regularsurface.hpp"

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/definitions.h"
#include "lib/utils.h"

class Vario;
class Simbox;
class FFTGrid;
class CovGrid2D;
class WellData;
class KrigingData2D;
class ModelSettings;

class Wavelet {
public:
  enum           difftypes{FIRSTORDERFORWARDDIFF, FIRSTORDERCENTRALDIFF, FIRSTORDERBACKWARDDIFF, FOURIER};
  enum           waveletDims{ONE_D = 0, THREE_D = 1};
  enum           waveletreadtypes{OLD, JASON, ESTIMATE, SGRI};
  
  virtual void   fft1DInPlace() = 0;
  virtual void   invFFT1DInPlace() = 0;

  //Constructors and destructor
  Wavelet(int dim);
  Wavelet(int dim, float * reflCoef);
  Wavelet(ModelSettings * modelSettings, int dim, float * reflCoef);
  Wavelet(Wavelet * wavelet, int dim);
  virtual ~Wavelet();

  Wavelet             * getLocalWavelet(int i,int j);
  virtual void          resample(float, int, float, float) {};
  virtual bool          consistentSize(int nzp, int nyp, int nxp) const = 0;
  virtual fftw_complex  getCAmp(int k, int j=0, int i=0) const = 0;
  virtual fftw_real     getRAmp(int k, int j=0, int i=0) = 0;
  virtual fftw_complex  getCAmp(int k, float scale, int j=0, int i=0) const = 0;
  virtual void          setRAmp(float value, int k, int j=0, int i=0) = 0;
  float                 getNorm() const {return norm_;}
  void                  setNorm(float norm) {norm_ = norm;}
  bool                  getIsReal() const {return(isReal_);} 
  virtual void          scale(float gain);
  int                   getDim() const {return dim_;}
  int                   getNz() const {return nz_;}
  float                 getDz() const {return dz_;}
  virtual int           getNx() const {return 0;}
  virtual int           getNy() const {return 0;}
  virtual int           getNxp() const {return 0;}
  virtual int           getNyp() const {return 0;}
  virtual float         getDx() const {return 0.0;}
  virtual float         getDy() const {return 0.0;}
  virtual FFTGrid     * getAmpCube() {return NULL;}
  virtual void          multiplyByR(float) {};
  
  //Note: Function below is mainly controlled by debugflag. Set overrideDebug = true to force.
  virtual void          printToFile(const std::string & fileName, bool overrideDebug = false) = 0;
  virtual void          writeWaveletToFile(const std::string & fileName, float, Simbox * simbox = NULL) = 0;
  void                  setShiftGrid(Grid2D * grid);
  void                  setGainGrid(Grid2D * grid);
  float                 getScale() const {return scale_;}
  
  // for noise estimation
  float                 calculateSNRatioAndLocalWavelet(Simbox        * simbox, 
                                                        FFTGrid       * seisCube, 
                                                        WellData     ** wells, 
                                                        Grid2D       *& shift, 
                                                        Grid2D       *& gain, 
                                                        ModelSettings * modelSettings,
                                                        char          * errText, 
                                                        int           & error,
                                                        Grid2D       *& noiseScaled, 
                                                        int             number, 
                                                        float           globalScale); 
  void                  printVecToFile(const std::string & fileName, fftw_real* vec ,int nzp) const;

  virtual void          write1DWLas3DWL() {};
  virtual void          write3DWLfrom1DWL() {};

protected:
//virtual float         getWaveletValue(float z, float * Wavelet, int center,int nx, float dz) = 0;
  virtual void          shiftAndScale(float, float) {};

  // for wavelet estimation
  void                  shiftReal(float shift, fftw_real* rAmp,int nt);
  void                  convolve(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* out_c,int cnzp) const;
  float                 findOptimalWaveletScale(fftw_real** synt_seis_r,fftw_real** seis_r,int nWells,int nzp,
                                                float* wellWeight,float& err,float* errWell,float* scaleOptWell,
                                                float* errWellOptScale) const;
  void                  findLocalNoiseWithGainGiven(fftw_real ** synt_r,
                                                    fftw_real ** seis_r,
                                                    int nWells,
                                                    int nzp,
                                                    float * wellWeight,
                                                    float & err,
                                                    float * errWell, 
                                                    float * errWellOptScale,
                                                    float * scaleOptWell,
                                                    Grid2D * gain, 
                                                    WellData **wells, Simbox *simbox) const;
  void                  estimateLocalGain(const CovGrid2D & cov,
                                          Grid2D         *& gain,
                                          float           * scaleOptWell,
                                          float             globalScale,
                                          int             * nActiveData,
                                          Simbox          * simbox,
                                          WellData       ** wells,
                                          int               nWells);
  
  void                  estimateLocalShift(const CovGrid2D & cov,
                                           Grid2D         *& shift,
                                           float           * shiftWell,
                                           int             * nActiveData,
                                           Simbox          * simbox,
                                           WellData       ** wells,
                                           int               nWells);
  
  void                  estimateLocalNoise(const CovGrid2D & cov,
                                           Grid2D         *& noiseScaled,
                                           float             globalNoise,
                                           float           * errWellOptScale,
                                           int             * nActiveData,
                                           Simbox          * simbox,
                                           WellData       ** wells,
                                           int               nWells);

  //void                flipVec(fftw_real* vec, int n);

  void                  fillInnWavelet(fftw_real* wavelet_r,int nzp,float dz);
  float                 findBulkShift(fftw_real* vec_r,float dz,int nzp);
  float                 getLocalTimeshift(int i, int j) const;
  float                 getLocalGainFactor(int i, int j) const;
  int                   getWaveletLengthI();
  float                 getWaveletLengthF();
  
  float                 theta_;                 // the reflection angle that the wavelet correspond to
  int                   readtype_;              // how is wavelet obtained? read from file[OLD JASON SGRI] or ESTIMATE

  float                 dz_;                    // Sampling interval of wavelet, unit [ ms ]
  int                   nz_;                    // length of wavelet
  int                   nzp_;                   // length of padded wavelet
  
  int                   cz_;                    // position of central point   
  bool                  inFFTorder_;            // is true if the wavelet is ordred with the central point at the start
                                         // false if the central point is in the middle
  bool                  isReal_;                // is true if the wavlet is real, false if it is fourier transformed 
  float                 norm_;                  // The (vector) norm of the wavelet (not function norm that divides by dz)
  float                 waveletLength_;         // Length of wavelet estimated as is amplitudes larger than 1/1000 * max amplitude

  float                 coeff_[3];              //Reflection coefficients.

  float                 maxShift_;//maximum shift of wavelet in ms
  float                 minRelativeAmp_;

  //NBNB The following parameters are NOT copied in copy constructor.
  
  float                 scale_;
  Grid2D               * shiftGrid_;             // 2D grid of shift
  Grid2D               * gainGrid_;              // 2D grid of gain factors.
 
  const int             dim_;
};

#endif
