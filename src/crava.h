/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef CRAVA_H
#define CRAVA_H

#include "fftw.h"
#include "definitions.h"
#include "libs/nrlib/flens/nrlib_flens.hpp"

class ModelGeneral;
class ModelAVOStatic;
class ModelAVODynamic;
class FFTGrid;
class FFTFileGrid;
class Wavelet;
class Wavelet1D;
class Simbox;
class RandomGen;
class WellData;
class CKrigingAdmin;
class CovGridSeparated;
class KrigingData3D;
class FaciesProb;
class GridMapping;
class ModelSettings;
class Corr;
class SpatialWellFilter;


class Crava
{
public:
  Crava(ModelSettings     * modelSettings,
        ModelGeneral      * modelGeneral,
        ModelAVOStatic    * modelAVOstatic,
        ModelAVODynamic   * modelAVOdynamic,
        SpatialWellFilter * spatwellfilter);
  ~Crava();

  int                    computePostMeanResidAndFFTCov();
  int                    simulate( RandomGen * randomGen );
  void                   computeSyntSeismic(FFTGrid * alpha, FFTGrid * beta, FFTGrid * rho);
  int                    computeSyntSeismicOld(FFTGrid * Alpha, FFTGrid * Beta, FFTGrid * Rho);
  FFTGrid              * getPostAlpha() { return postAlpha_ ;}
  FFTGrid              * getPostBeta()  { return postBeta_  ;}
  FFTGrid              * getPostRho()   { return postRho_   ;}

  int                    getWarning(std::string & wText)  const {if(scaleWarning_>0) wText=scaleWarningText_; return scaleWarning_;}

  void                   printEnergyToScreen();
  void                   computeFaciesProb(SpatialWellFilter *filteredlogs, bool useFilter);
  int                    getRelative();


  void                   computeG(NRLib::Matrix & G) const;
  NRLib::Matrix          getPriorVar0() const;
  NRLib::Matrix          getPostVar0() const;
  NRLib::SymmetricMatrix getSymmetricPriorVar0() const;
  NRLib::SymmetricMatrix getSymmetricPostVar0() const;
  void                   newPosteriorCovPointwise(NRLib::Matrix & sigmanew,
                                                  NRLib::Matrix & G,
                                                  NRLib::Vector & scales,
                                                  NRLib::Matrix & sigmamdnew) const;
  NRLib::Matrix          computeFilter(NRLib::SymmetricMatrix & priorCov,
                                   NRLib::SymmetricMatrix & posteriorCov) const;

  void                   doPredictionKriging();

private:
  void                   computeDataVariance(void);
  void                   setupErrorCorrelation(ModelSettings * modelSttings, const std::vector<Grid2D *> & noiseScale);
  void                   computeVariances(fftw_real* corrT, ModelSettings * modelSettings);
  void                   writeBWPredicted(void);
  float                  getEmpSNRatio(int l)     const { return empSNRatio_[l]     ;}
  float                  getTheoSNRatio(int l)    const { return theoSNRatio_[l]    ;}
  float                  getSignalVariance(int l) const { return signalVariance_[l] ;}
  float                  getErrorVariance(int l)  const { return errorVariance_[l]  ;}
  float                  getDataVariance(int l)   const { return dataVariance_[l]   ;}
  void                   computeElasticImpedanceTimeCovariance(fftw_real* eiCovT,const float* corrT,float** Var0,float * A) const;
  void                   computeReflectionCoefficientTimeCovariance(fftw_real* refCovT,const float* corrT,float** Var0,float * A ) const ;

  int                    checkScale(void);
  void                   fillkW(int k, fftw_complex* kW, Wavelet** seisWavelet);
  void                   fillInverseAbskWRobust(int k, fftw_complex* invkW ,Wavelet1D** seisWaveletForNorm);
  void                   fillkWNorm(int k, fftw_complex* kWNorm, Wavelet1D** wavelet);

  void                   computeAdjustmentFactor(fftw_complex* relativeWeights, Wavelet1D* wLocal, double scaleF, Wavelet * wGlobal, const Corr* corr,float * A, float errorVar);

  FFTGrid              * createFFTGrid();
  FFTGrid              * copyFFTGrid(FFTGrid * fftGridOld);
  FFTFileGrid          * copyFFTGrid(FFTFileGrid * fftGridOld);

  float                  computeWDCorrMVar (Wavelet1D* WD, fftw_real* corrT);
  //float                computeWDCorrMVar (Wavelet* WD);

  void                   divideDataByScaleWavelet();
  void                   multiplyDataByScaleWaveletAndWriteToFile(const std::string & typeName);
  void                   doPostKriging(FFTGrid & postAlpha, FFTGrid & postBeta, FFTGrid & postRho);

  void                   correctAlphaBetaRho(ModelSettings *modelSettings);

  FFTGrid *              computeSeismicImpedance(FFTGrid * alpha,
                                                 FFTGrid * beta,
                                                 FFTGrid * rho,
                                                 int angle);

  bool               fileGrid_;         // is true if is storage is on file
  const Simbox     * simbox_;           // the simbox
  int                nx_;               // dimensions of the problem
  int                ny_;
  int                nz_;
  int                nxp_;              // padded dimensions
  int                nyp_;
  int                nzp_;

  Corr             * correlations_;     //

  int                ntheta_;           // number of seismic cubes and number of wavelets
  float              lowCut_;           // lowest frequency that is inverted
  float              highCut_;          // highest frequency that is inverted

  int                nSim_;             // number of simulations
  float            * thetaDeg_;         // in degrees

  FFTGrid          * meanAlpha_;        // mean values
  FFTGrid          * meanBeta_;
  FFTGrid          * meanRho_;
  FFTGrid          * meanAlpha2_;       // copy of mean values, to be used for facies prob, new method
  FFTGrid          * meanBeta2_;
  FFTGrid          * meanRho2_;
  float           ** parPointCov_;

  Wavelet         ** seisWavelet_;      // wavelet operator that define the forward map.
  FFTGrid         ** seisData_;         // Data
  double          ** errThetaCov_;      //
  float              wnc_ ;             // if wnc=0.01 1% of the error wariance is white this has largest effect on
                                        // high frequency components. It makes everything run smoother we
                                        // avoid ill posed problems.
  float           ** A_;                // coefficients in Aki-Richards 3 term reflection coefficients

  float            * empSNRatio_;       // signal noise ratio empirical
  float            * theoSNRatio_;      // signal noise ratio from model
  float            * modelVariance_;
  float            * signalVariance_;
  float            * errorVariance_;
  float            * dataVariance_;

  FFTGrid          * postAlpha_;        // posterior values
  FFTGrid          * postBeta_;
  FFTGrid          * postRho_;

  int                krigingParameter_;
  WellData        ** wells_;
  int                nWells_;

  int                scaleWarning_;
  std::string        scaleWarningText_;

  int                outputGridsSeismic_; // See modelsettings.h for bit interpretation.
  int                outputGridsElastic_;
  bool               writePrediction_;  // Write prediction grids?

  float              energyTreshold_;   // If energy in reflection trace divided by mean energy
                                        // in reflection trace is lower than this, the reflections
                                        // will be interpolated. Default 0, set from model.
  RandomGen        * random_;
  FaciesProb       * fprob_;

  ModelSettings    * modelSettings_;
  ModelGeneral     * modelGeneral_;
  ModelAVOStatic   * modelAVOstatic_;
  ModelAVODynamic  * modelAVOdynamic_;


  NRLib::Grid2D<double **> * sigmamdnew_;
};

#endif
