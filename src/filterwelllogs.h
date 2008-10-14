#ifndef FILTERWELLLOGS_H
#define FILTERWELLLOGS_H
#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

class FFTGrid;
class Simbox;
class RandomGen;
class ModelSettings;
class CKrigingAdmin;
class KrigingData;

class FilterWellLogs
{
public:
  FilterWellLogs(FFTGrid *postCovAlpha, FFTGrid *postCovBeta, FFTGrid *postCovRho, FFTGrid *postCrCovAlphaBeta,
                 FFTGrid *postCrCovAlphaRho, FFTGrid *postCrCovBetaRho, int nzp, int nz, WellData ** wells, int nWells, 
                 float lowCut, float highCut, int relative,
                 const Simbox * timeSimboxConstThick, 
                 const Simbox * timeSimboxPropThick, 
                 RandomGen *random, float *corrprior, float ** sigma0, bool faciesprob);
  ~FilterWellLogs();

  float * getAlphaFiltered() { return alphafiltered_ ;}
  float * getBetaFiltered()  { return betafiltered_  ;}
  float * getRhoFiltered()   { return rhofiltered_   ;}
  float * getAlphaBlock()    { return alphablock_    ;}
  float * getBetaBlock()     { return betablock_     ;}
  float * getRhoBlock()      { return rhoblock_      ;}
  int   * getFaciesLog()     { return facieslog_     ;}
  int     getNdata()         { return ndata_         ;} 

private:
  void doFiltering(WellData **wells, int nWells,
                   fftw_real *postcova, fftw_real *postcovb, fftw_real *postcovr,
                   fftw_real *postcrab, fftw_real *postcrar, fftw_real *postcrbr, 
                   float lowCut, float highCut, int relative, int nz, int nzp, 
                   const Simbox * timeSimboxConstThick, 
                   const Simbox * timeSimboxPropThick, 
                   RandomGen *random, float ** sigma0);
  void extrapolate(float * log,
                   int     nz) ;
  void calcFilter(fftw_complex **sigmaK, fftw_complex **sigmaE, double **F);


  //float         * vtAlphaFiltered_;
  //float         * vtBetaFiltered_; 
  //float         * vtRhoFiltered_;

  float         * alphafiltered_;
  float         * betafiltered_; 
  float         * rhofiltered_;

  float         * alphablock_;
  float         * betablock_;
  float         * rhoblock_;

  int           * facieslog_;
  fftw_real     * corrprior_;
  bool            faciesprob_;

  int             ndata_;
};
#endif

