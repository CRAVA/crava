#ifndef WAVELET_H
#define WAVELET_H

#include "nrlib/surface/regularsurface.hpp"

#include "fft/include/fftw.h"
#include "lib/global_def.h"
#include "src/definitions.h"

class Simbox;
class FFTGrid;
class WellData;
class ModelSettings;

class Wavelet {
public:
  enum           difftypes{FIRSTORDERFORWARDDIFF, FIRSTORDERCENTRALDIFF, FIRSTORDERBACKWARDDIFF, FOURIER};
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
  bool                  getIsReal() const {return(isReal_);} 
  virtual void          scale(float gain);
  int                   getDim() const {return dim_;}
  
  //Note: Function below is mainly controlled by debugflag. Set overrideDebug = true to force.
  virtual void          printToFile(char* fileName, bool overrideDebug = false) = 0;
  virtual void          writeWaveletToFile(char*, float, Simbox * simbox = NULL) = 0;
  void                  setShiftGrid(Surface * grid, Simbox * simbox);
  void                  setGainGrid(Surface * grid, Simbox * simbox);
  float                 getScale() const {return scale_;}
  
  // for noise estimation
  float                 getNoiseStandardDeviation(Simbox * simbox,
                                                  FFTGrid * seisCube, 
                                                  WellData ** wells, 
                                                  int nWells, 
                                                  char *errText, 
                                                  int &error); 
  void                  setReflCoeff(float * coeff) {for(int i=0;i<3;i++) coeff_[i] = coeff[i];}
  void                  printVecToFile(char* fileName, fftw_real* vec ,int nzp) const;

  virtual void          write1DWLas3DWL() {};
  virtual void          write3DWLfrom1DWL() {};

protected:
//virtual float         getWaveletValue(float z, float * Wavelet, int center,int nx, float dz) = 0;
  virtual void          shiftAndScale(float, float) {};

  // for wavelet estimation
  void                  fft(fftw_real* rAmp,fftw_complex* cAmp,int nt);     
  void                  fftInv(fftw_complex* cAmp,fftw_real* rAmp,int nt);    
  void                  shiftReal(float shift, fftw_real* rAmp,int nt);
  void                  fillInCpp(float* alpha,float* beta,float* rho,int start,int length,fftw_real* cpp_r,int nzp);
  void                  fillInSeismic(float* seismicData,int start,int length,fftw_real* seis_r,int nzp) const;
  void                  estimateCor(fftw_complex* var1_c,fftw_complex* var2_c,fftw_complex* ccor_1_2_c,int cnzp) const;
  void                  convolve(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* out_c,int cnzp) const;
  float                 computeElasticImpedance(float vp, float vs, float rho) const;
  void                  findContiniousPartOfData(bool* hasData,int nz,int &start,int &length) const;
  float                 findOptimalWaveletScale(fftw_real** synt_seis_r,fftw_real** seis_r,int nWells,int nzp,
                                                float* wellWeight,float& err,float* errWell,float* scaleOptWell,
                                                float* errWellOptScale) const;
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
  int                   gridNI_;                // number of elements for  i in shiftGrid_ and gainGrid_;
  int                   gridNJ_;
  float                 scale_;
  float               * shiftGrid_;             // 2D grid of shift
  float               * gainGrid_;              // 2D grid of gain factors.
  const int             dim_;
};

#endif
