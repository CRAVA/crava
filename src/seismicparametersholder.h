/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SEISMIC_PARAMETERS_HOLDER
#define SEISMIC_PARAMETERS_HOLDER

#include "nrlib/flens/nrlib_flens.hpp"
#include <src/fftgrid.h>

class ModelSettings;

// A class holding the pointers for the seismic parameters
// for easy parameter transmission of the pointers to the class TimeEvolution.

class SeismicParametersHolder
{
public:
  SeismicParametersHolder(void);

  ~SeismicParametersHolder(void);

  //void getMeanReferenceVector(std::vector<FFTGrid *> mu);
  //void getCovReferenceVector(std::vector<FFTGrid *> sigma);

  FFTGrid                     * GetMuAlpha()                            { return muAlpha_        ;}
  FFTGrid                     * GetMuBeta()                             { return muBeta_         ;}
  FFTGrid                     * GetMuRho()                              { return muRho_          ;}
  FFTGrid                     * GetCovAlpha()                           { return covAlpha_       ;}
  FFTGrid                     * GetCovBeta()                            { return covBeta_        ;}
  FFTGrid                     * GetCovRho()                             { return covRho_         ;}
  FFTGrid                     * GetCrCovAlphaBeta()                     { return crCovAlphaBeta_ ;}
  FFTGrid                     * GetCrCovAlphaRho()                      { return crCovAlphaRho_  ;}
  FFTGrid                     * GetCrCovBetaRho()                       { return crCovBetaRho_   ;}

  void                          invFFTAllGrids();
  void                          invFFTCovGrids();
  void                          FFTCovGrids();
  void                          FFTAllGrids();
  void                          updatePriorVar();

  void                          setPriorVar0(NRLib::Matrix priorVar0)   { priorVar0_ = priorVar0 ;}

  void                          setBackgroundParameters(FFTGrid  * muAlpha,
                                                        FFTGrid  * muBeta,
                                                        FFTGrid  * muRho);

  void                          setBackgroundParametersIntervals(const std::vector<std::vector<NRLib::Grid<double> > > & mu_parameters);

  void                          setBackgroundParametersInterval(const std::vector<NRLib::Grid<double> > & mu_parameters);

  void                          setCovParameters(const std::vector<NRLib::Grid<double> > & cov_parameters);

  void                          setCrCovParameters(const NRLib::Grid<double> & cr_cov_vp_vs,
                                                   const NRLib::Grid<double> & cr_cov_vp_rho,
                                                   const NRLib::Grid<double> & cr_cov_vs_rho);

  void                          copyBackgroundParameters(FFTGrid  * muAlpha,
                                                         FFTGrid  * muBeta,
                                                         FFTGrid  * muRho);

  void                          setCorrelationParameters(float                    ** priorVar0,
                                                         const std::vector<float>  & priorCorrT,
                                                         Surface                   * priorCorrXY,
                                                         const int                 & minIntFq,
                                                         const float               & corrGradI,
                                                         const float               & corrGradJ,
                                                         const int                 & nx,
                                                         const int                 & ny,
                                                         const int                 & nz,
                                                         const int                 & nxPad,
                                                         const int                 & nyPad,
                                                         const int                 & nzPad);

  void                          allocateGrids(const int nx,
                                              const int ny,
                                              const int nz,
                                              const int nxPad,
                                              const int nyPad,
                                              const int nzPad);


  NRLib::Matrix                 getPriorVar0(void) const;

  float                       * getPriorCorrTFiltered(int nz, int nzp) const;

  fftw_real                   * computeCircCorrT(const std::vector<float> & priorCorrT,
                                                 const int                & minIntFq,
                                                 const int                & nzp) const;

  fftw_real                   * extractParamCorrFromCovAlpha(int nzp) const;

  void                          printPriorVariances(void) const;

  void                          printPostVariances(const NRLib::Matrix & postVar0) const;

  void                          getNextParameterCovariance(fftw_complex **& parVar) const;

  void                          writeFilePriorCorrT(fftw_real   * priorCorrT,
                                                    const int   & nzp,
                                                    const float & dt) const;

  void                          writeFilePriorVariances(const ModelSettings      * modelSettings,
                                                        const std::vector<float> & priorCorrT,
                                                        const Surface            * priorCorrXY,
                                                        const float              & dt) const;

  void                          writeFilePostVariances(const NRLib::Matrix      & postVar0,
                                                       const std::vector<float> & postCovAlpha00,
                                                       const std::vector<float> & postCovBeta00,
                                                       const std::vector<float> & postCovRho00) const;

  void                          writeFilePostCovGrids(Simbox const * simbox) const;

  std::vector<float>            createPostCov00(FFTGrid * postCov) const;

private:
  void                          createCorrGrids(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid);

  void                          initializeCorrelations(const Surface            * priorCorrXY,
                                                       const std::vector<float> & priorCorrT,
                                                       const float              & corrGradI,
                                                       const float              & corrGradJ,
                                                       const int                & lowIntCut,
                                                       const int                & nzp);

  FFTGrid                     * createFFTGrid(int nx,  int ny,  int nz,
                                              int nxp, int nyp, int nzp,
                                              bool fileGrid);


  void                          makeCircCorrTPosDef(fftw_real * circCorrT,
                                                    const int & minIntFq,
                                                    const int & nzp) const;

  void                         makeCorrXYPosDef(Surface         & priorCorrXY);

  float                         getOrigin(FFTGrid * grid) const;

  void                          writeFilePostCorrT(const std::vector<float> & postCov,
                                                   const std::string        & subDir,
                                                   const std::string        & baseName) const;

  FFTGrid * muAlpha_;
  FFTGrid * muBeta_ ;
  FFTGrid * muRho_  ;
  FFTGrid * covAlpha_;
  FFTGrid * covBeta_ ;
  FFTGrid * covRho_  ;
  FFTGrid * crCovAlphaBeta_;
  FFTGrid * crCovAlphaRho_ ;
  FFTGrid * crCovBetaRho_  ;

  NRLib::Grid<double> mu_vp_;
  NRLib::Grid<double> mu_vs_;
  NRLib::Grid<double> mu_rho_;

  NRLib::Grid<double> cov_vp_;
  NRLib::Grid<double> cov_vs_;
  NRLib::Grid<double> cov_rho_;

  NRLib::Grid<double> cr_cov_vp_vs_;
  NRLib::Grid<double> cr_cov_vp_rho_;
  NRLib::Grid<double> cr_cov_vs_rho_;

  std::vector<NRLib::Grid<double> > mu_vp_intervals_; //Vector over intervals.
  std::vector<NRLib::Grid<double> > mu_vs_intervals_;
  std::vector<NRLib::Grid<double> > mu_rho_intervals_;

  NRLib::Matrix priorVar0_;

};
#endif
