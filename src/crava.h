#ifndef CRAVA_H
#define CRAVA_H

#include "fft/include/fftw.h"
#include "definitions.h"

class Model;
class FFTGrid;
class FFTFileGrid;
class Wavelet;
class Simbox;
class RandomGen;
class WellData;
class KrigingData3D;
class CKrigingAdmin;
class CovGridSeparated;
class FaciesProb;
class GridMapping;
class FilterWellLogs;

class Crava {
public:
  Crava(Model * model);
  ~Crava();
  int                computePostMeanResidAndFFTCov();
  int                simulate( RandomGen * randomGen );
  int                computePostCov();
  int                computeSyntSeismic(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho);

  FFTGrid          * getpostAlpha()           {return postAlpha_;}
  FFTGrid          * getpostBeta()            {return postBeta_;}
  FFTGrid          * getpostRho()             {return postRho_;}

  const FFTGrid    * getpostAlpha()           const {return postAlpha_;}
  const FFTGrid    * getpostBeta()            const {return postBeta_;}
  const FFTGrid    * getpostRho()             const {return postRho_;}

  int                getNTheta()              const {return ntheta_;}
  int                getWarning(char* wText)  const {if(scaleWarning_>0) sprintf(wText,"%s",scaleWarningText_); return scaleWarning_;}

  void               printEnergyToScreen();
  void               computeFaciesProb(FilterWellLogs *filteredlogs);
  void               filterLogs(Simbox          * timeSimboxConstThick,
                                FilterWellLogs *& filterlogs);

private: 
  int                computeAcousticImpedance(FFTGrid * Alpha, FFTGrid * Rho, char * fileName);
  int                computeShearImpedance(FFTGrid * Beta, FFTGrid * Rho, char * fileName);
  int                computeVpVsRatio(FFTGrid * Alpha, FFTGrid * Beta, char * fileName);
  int                computePoissonRatio(FFTGrid * Alpha, FFTGrid * Beta, char * fileName);
  int                computeLameMu(FFTGrid * Beta, FFTGrid * Rho , char * FileName);
  int                computeLameLambda(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, char * fileName);
  int                computeMuRho(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, char * fileName);
  int                computeLambdaRho(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho, char * fileName);
  void               computeDataVariance(void);
  void               setupErrorCorrelation(Model * model);
  void               computeVariances(fftw_real* corrT, Model * model);
  const FFTGrid    * getpostCovAlpha()        const {return postCovAlpha_;}
  const FFTGrid    * getpostCovBeta()         const {return postCovBeta_;}       
  const FFTGrid    * getpostCovRho()          const {return postCovRho_;}        
  const FFTGrid    * getpostCrCovAlphaBeta()  const {return postCrCovAlphaBeta_;}
  const FFTGrid    * getpostCrCovAlphaRho()   const {return postCrCovAlphaRho_;}
  const FFTGrid    * getpostCrCovBetaRho()    const {return postCrCovBetaRho_;}
  float              getTheta(int l)          const {return theta_[l];}
  float              getEmpSNRatio(int l)     const {return empSNRatio_[l];}
  float              getTheoSNRatio(int l)    const {return theoSNRatio_[l];}
  float              getSignalVariance(int l) const {return signalVariance_[l];}
  float              getErrorVariance(int l)  const {return errorVariance_[l];}
  float              getDataVariance(int l)   const {return dataVariance_[l];}
  const Simbox     * getSimbox()              const {return(simbox_);}
  int                checkScale(void);
    //Conventions for writePars:
  // simNum = -1 indicates prediction, otherwise filename ends with n+1.
  // All grids are in normal domain, and on log scale.
  void               writePars(FFTGrid * alpha, FFTGrid * beta, FFTGrid * rho, int simNum); 

  void               fillkW(int k, fftw_complex* kW );
  void               fillInverseAbskWRobust(int k, fftw_complex* invkW );
  void               fillkWNorm(int k, fftw_complex* kWNorm, Wavelet** wavelet);
  FFTGrid          * createFFTGrid();
  FFTGrid          * copyFFTGrid(FFTGrid * fftGridOld);
  FFTFileGrid      * copyFFTGrid(FFTFileGrid * fftGridOld);

  float              computeWDCorrMVar (Wavelet* WD, fftw_real* corrT);
  float              computeWDCorrMVar (Wavelet* WD);

  void               divideDataByScaleWavelet();
  void               multiplyDataByScaleWaveletAndWriteToFile(const char* typeName);
  void               dumpCorrT(float* corrT,float dt);
  void               initPostKriging();          
  void               writeToFile(char * fileName, FFTGrid * grid, std::string sgriLabel = "NO_LABEL");

  int                fileGrid_;        // is true if is storage is on file 
  const Simbox     * simbox_;          // the simbox
  //const Simbox     * depthSimbox_;     // simbox with depth surfaces
  int                nx_;              // dimensions of the problem
  int                ny_;
  int                nz_;
  int                nxp_;             // padded dimensions
  int                nyp_;
  int                nzp_; 

  int                ntheta_;          // number of seismic cubes and number of wavelets
  float              lowCut_;          // lowest frequency that is inverted
  float              highCut_;         // highest frequency that is inverted

  int                nSim_;            // number of simulations
  float            * theta_;           // in radians

  FFTGrid          * meanAlpha_;       // mean values
  FFTGrid          * meanBeta_;
  FFTGrid          * meanRho_;
  FFTGrid          * meanAlpha2_;       // copy of mean values, to be used for facies prob, new method 
  FFTGrid          * meanBeta2_;
  FFTGrid          * meanRho2_;
  FFTGrid          * parSpatialCorr_;   // parameter correlation
  float           ** parPointCov_; 

  Wavelet         ** seisWavelet_;      // wavelet operator that define the forward map.
  FFTGrid         ** seisData_;         // Data
  FFTGrid          * errCorrUnsmooth_;  // Error correlation
  float           ** errThetaCov_;      //
  float              wnc_ ;          // if wnc=0.01 1% of the error wariance is white this has largest effect on
                                     // high frequency components. It makes everything run smoother we
                                     // avoid ill posed problems.
  float           ** A_;             // 

  float            * empSNRatio_;    // signal noise ratio empirical
  float            * theoSNRatio_;   // signal noise ratio from model
  float            * modelVariance_;
  float            * signalVariance_;
  float            * errorVariance_;
  float            * dataVariance_;
  FFTGrid          * postAlpha_;     // posterior values 
  FFTGrid          * postBeta_;
  FFTGrid          * postRho_;

  FFTGrid          * postCovAlpha_;
  FFTGrid          * postCovBeta_;
  FFTGrid          * postCovRho_;
  FFTGrid          * postCrCovAlphaBeta_;
  FFTGrid          * postCrCovAlphaRho_;
  FFTGrid          * postCrCovBetaRho_;

  float            * krigingParams_;
  WellData        ** wells_;
  int                nWells_;
  
  int                scaleWarning_;
  char             * scaleWarningText_;

  int                outputFlag_; //See model.h for bit interpretation.
  
  float              energyTreshold_; //If energy in reflection trace divided by mean energy
                                      //in reflection trace is lower than this, the reflections
                                      //will be interpolated. Default 0, set from model.
  RandomGen        * random_;
  CKrigingAdmin    * pKriging_;
  CovGridSeparated * covAlpha_, *covBeta_, *covRho_, *covCrAlphaBeta_, *covCrAlphaRho_, *covCrBetaRho_;
  KrigingData3D    * kd_;        
  FaciesProb       * fprob_;
  Model            * model_;
  float            * corrprior_;
  
};
#endif
