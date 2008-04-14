#ifndef WAVELET_H
#define WAVELET_H

#include "fft/include/fftw.h"
#include "lib/global_def.h"

class Simbox;
class FFTGrid;
class WellData;
struct irapgrid;
class ModelSettings;

class Wavelet {
public:
  enum          difftypes{FIRSTORDERFORWARDDIFF,FIRSTORDERCENTRALDIFF,FIRSTORDERBACKWARDDIFF,FOURIER};
  
  void fft1DInPlace();

  //Constructors and destructor
  Wavelet(Simbox * simbox, FFTGrid * seisCube, WellData ** wells, ModelSettings * modelSettings, float * coef);
  Wavelet(char * fileName, ModelSettings * modelSettings);
  Wavelet(Wavelet * wavelet);
  Wavelet(Wavelet * wavelet,int difftype);
  Wavelet(int difftype, int nz, int nzp);
  ~Wavelet();

  Wavelet*       getLocalWavelet(int i,int j);
  void           resample(float dz, int nz, float pz,float theta);
  bool           consistentSize(int nzp) const { if (nzp!=nzp_) printf("nzp=%d  nzp_wavelet=%d\n",nzp,nzp_); return (nzp==nzp_);}
  fftw_complex   getCAmp(int k) const;
  fftw_real      getRAmp(int k);
  fftw_complex   getCAmp(int k, float scale) const;
  float          getNorm() const {return norm_;}
  bool           getIsReal() const {return(isReal_);} 
  void           scale(float gain);
   
  //Note: Function below is mainly controlled by debugflag. Set overrideDebug = true to force.
  void           printToFile(char* fileName, bool overrideDebug = false) const;
  void           writeWaveletToFile(char* fileName,float approxDz);
  void           setShiftGrid(irapgrid * grid, Simbox * simbox);
  void           setGainGrid(irapgrid * grid, Simbox * simbox);
  float          getScale() const {return scale_;}
  
  // for noise estimation
  float          getNoiseStandardDeviation(Simbox * simbox, FFTGrid * seisCube, WellData ** wells, int nWells); 
  void           setReflCoeff(float * coeff) {for(int i=0;i<3;i++) coeff_[i] = coeff[i];}

private:
  enum           waveletreadtypes{OLD,JASON,ESTIMATE,SGRI};
  
  void           invFFT1DInPlace();
  void           flipUpDown();
  float          getWaveletValue(float z, float * Wavelet, int center,int nx, float dz);
//  void           printRealWaveletToScreen(){for(int i=0;i < nzp_;i++) printf("value %i: %f \n", i,rAmp_[i]) ;}
//  void           printComplexWaveletToScreen(){for(int i=0;i < cnzp_;i++) printf("re: %f  im: %f \n", cAmp_[i].re,cAmp_[i].im) ;}
//  void           printCompNormWaveletToScreen(){for(int i=0;i < cnzp_;i++) printf("q.norm %i: %f \n", i,(cAmp_[i].re*cAmp_[i].re+cAmp_[i].im*cAmp_[i].im) );}
  int            getWaveletLengthI();
  float          getWaveletLengthF();
  void           shiftAndScale(float shift,float gain);
  void           printToFile(char* fileName, fftw_real* vec ,int nzp) const;

  // Functions for reading wavelet from file
  int            getWaveletFileFormat(char * fileName);
  void           WaveletReadOld(char * fileName);
  void           WaveletReadJason(char * fileName);
  void			 WaveletReadSgri(char *fileName);
  
  // for wavelet estimation
  void           fft(fftw_real* rAmp,fftw_complex* cAmp,int nt);     
  void           fftInv(fftw_complex* cAmp,fftw_real* rAmp,int nt);    
  void           shiftReal(float shift, fftw_real* rAmp,int nt);
//  void           shiftReal(int shift, fftw_real* rAmp,int nt);
  void           fillInCpp(float* alpha,float* beta,float* rho,int start,int length,fftw_real* cpp_r,int nzp);
  void           fillInSeismic(float* seismicData,int start,int length,fftw_real* seis_r,int nzp) const;
  void           estimateCor(fftw_complex* var1_c,fftw_complex* var2_c,fftw_complex* ccor_1_2_c,int cnzp) const;
  void           convolve(fftw_complex* var1_c ,fftw_complex* var2_c, fftw_complex* out_c,int cnzp) const;
  float          computeElasticImpedance(float vp, float vs, float rho) const;
  void           findContiniousPartOfData(bool* hasData,int nz,int &start,int &length) const;
//  int*           getIndexPrior(int start,int nInd,int nzp);
//  int*           getIndexPost(int start,int nInd,int nzp);
  float          shiftOptimal(fftw_real** ccor_seis_cpp_r,float* wellWeight,float* dz,int nWells,int nzp,float* shiftWell);
  void           multiplyPapolouis(fftw_real** vec, float* dz,int nWells,int nzp, float waveletLength) const;
  void           getWavelet(fftw_real** ccor_seis_cpp_r,fftw_real** cor_cpp_r,fftw_real** wavelet_r,float* wellWeight,int nWells,int nt);
  fftw_real*     averageWavelets(fftw_real** wavelet_r,int nWells,int nzp,float* wellWeight,float* dz,float dzOut) const;
  float          findOptimalWaveletScale(fftw_real** synt_seis_r,fftw_real** seis_r,int nWells,int nzp,
                    float* wellWeight,float& err,float* errWell,float* scaleOptWell,float* errWellOptScale) const;
 // void           sampleWaveletBackToWell(fftw_real** wavelet_r,int nWells,int nzp,float* dz,float dz0);

  void           fillInnWavelet(fftw_real* wavelet_r,int nzp,float dz);
  float          findBulkShift(fftw_real* vec_r,float dz,int nzp);
  float          getArrayValueOrZero(int i ,float * Wavelet, int nz) const;
  float          getLocalTimeshift(int i, int j) const;
  float          getLocalGainFactor(int i, int j) const;

  float          theta_;                 // the reflection angle that the wavelet correspond to
  int            readtype_;              // how is wavelet obtained? read from file[OLD JASON SGRI] or ESTIMATE

  float          dz_;                    // Sampling interval of wavelet, unit [ ms ]
  int            nz_;                    // length of wavelet
  int            nzp_;                   // length of padded wavelet
  int	         cnzp_;	                 // size in z direction for storage inplace algorithm (complex grid) nzp_/2+1
  int	         rnzp_	;                // expansion in z direction for storage inplace algorithm (real grid) 2*(nzp_/2+1)

  int            cz_;                    // position of central point   
  bool           inFFTorder_;            // is true if the wavelet is ordred with the central point at the start
                                         // false if the central point is in the middle
  bool           isReal_;                // is true if the wavlet is real, false if it is fourier transformed 
  fftw_real*     rAmp_;                  // The amplitude of the wavelet  
  fftw_complex*  cAmp_;                  // The amplitude of the wavelet complex (if fourier transformed )
  float          norm_;                  // The (vector) norm of the wavelet (not function norm that divides by dz) 
  float          waveletLength_;         // Length of wavelet estimated as is amplitudes larger than 1/1000 * max amplitude

  float          coeff_[3];              //Reflection coefficients.

  int            errCode_;               // Code of error message
  char           errText_[6*MAX_STRING]; // Error message in case of mistakes in the file format

  float          maxShift_;//maximum shift of wavelet in ms
  float          minRelativeAmp_;

  //NBNB The following parameters are NOT copied in copy constructor.
  int            gridNI_;                // number of elements for  i in shiftGrid_ and gainGrid_;
  int            gridNJ_;
  float          scale_;
  float *        shiftGrid_;             // 2D grid of shift
  float *        gainGrid_;              // 2D grid of gain factors.

  // for wavelet estimation
 // float *        waveletshift_; //shift of wavelet in well
 // float *        errorLocalWaveStd_;
  //float *        errorGlobalWaveStd_;
  //float *        weight_; // importance of data in well for estimation
};

#endif
