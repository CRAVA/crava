#ifndef FACIESPROB_H
#define FACIESPROB_H

#include "fft/include/fftw.h"
#include "fft/include/rfftw.h"
#include "fft/include/fftw-int.h"
#include "fft/include/f77_func.h"

class FFTGrid;
class Crava;
class Simbox;
class RandomGen;
class CKrigingAdmin;
class KrigingData;
class ModelSettings;
class FilterWellLogs;

class FaciesProb
{
public:
  FaciesProb(ModelSettings * modelSettings,
             int filegrid, float ** sigma0, float *corrprior, 
             int nzp, int nz, 
             FFTGrid* bgAlpha, FFTGrid* bgBeta, FFTGrid* bgRho, 
             RandomGen *random, float p_undef, float *priorFacies);
  ~FaciesProb();
  float        ** makeFaciesHistAndSetPriorProb(float* alphafiltered_,float* betafiltered_,float* rhofiltered_,int* facieslog_);
  void            makeFaciesDens(int nfac);
  void            makeFaciesProb(int nfac, FFTGrid *postAlpha, FFTGrid *postBeta, FFTGrid *postRho);
  FFTGrid       * getFaciesProb(int i){return faciesProb_[i];};
 
  // FOR COMPARISION TO PCUBE
  float        ** makeFaciesHistAndSetPriorProb2(float* alpha, float* beta, float* rho,int* faciesL); // PCUBE
  void            makeFaciesDens2(int nfac);// PCUBE

  void            calculateConditionalFaciesProb(WellData **wells, int nwells);
  void            setFilteredLogs(FilterWellLogs *filterlogs);
  void            setSigmaPost(FFTGrid *postCovAlpha, FFTGrid *postCovBeta, FFTGrid *postCovRho, FFTGrid *postCrCovAlphaBeta,
                               FFTGrid *postCrCovAlphaRho, FFTGrid *postCrCovBetaRho);
private:
  
  void            getMinMax(float* alpha,float* beta,float* rho,int* facies);
  void            calculateVariances(float* alpha,float* beta,float* rho,int* facies);
  
  float           findDensity(float alpha, float beta, float rho, int facies);
  void            calculateFaciesProb( FFTGrid *alphagrid, FFTGrid *betagrid, FFTGrid *rhogrid);

  ModelSettings * modelSettings_;
  const Simbox  * simbox_;
  fftw_real     * corrprior_;
  float        ** sigma0_;
  int             nFacies_, nbins_, nbinsr_, nobs_;
  int           * nData_;
  FFTGrid      ** density_;  // Note: normalizes to the number of observations of facies type. 
  FFTGrid      ** faciesProb_;
  float           dalpha_, dbeta_, drho_;
  float           alphamin_, betamin_, rhomin_;
  float           alphamax_, betamax_, rhomax_;
  int             nzp_, nz_, ndata_;
  int             fileGrid_;

  FFTGrid       * bgAlpha_;
  FFTGrid       * bgBeta_;
  FFTGrid       * bgRho_;

  float         * alphafiltered_;
  float         * betafiltered_; 
  float         * rhofiltered_;

  float         * alphablock_;
  float         * betablock_;
  float         * rhoblock_;

  int           * facieslog_;

  float           varAlpha_;
  float           varBeta_;
  float           varRho_;  // variances

  float           meanA_;
  float           meanB_; 
  float           meanR_;   // means

  RandomGen     * random_;

  float           p_undefined_;
  float         * priorFacies_;
  float        ** sigmaPost_; // for Pcube
};
#endif


