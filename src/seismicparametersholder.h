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

  FFTGrid                     * GetMuVp()                          { return muVp_       ;}
  FFTGrid                     * GetMuVs()                          { return muVs_       ;}
  FFTGrid                     * GetMuRho()                         { return muRho_      ;}
  FFTGrid                     * GetCovVp()                         { return covVp_      ;}
  FFTGrid                     * GetCovVs()                         { return covVs_      ;}
  FFTGrid                     * GetCovRho()                        { return covRho_     ;}
  FFTGrid                     * GetCrCovVpVs()                     { return crCovVpVs_  ;}
  FFTGrid                     * GetCrCovVpRho()                    { return crCovVpRho_ ;}
  FFTGrid                     * GetCrCovVsRho()                    { return crCovVsRho_ ;}

  void                          invFFTAllGrids();
  void                          invFFTCovGrids();
  void                          FFTCovGrids();
  void                          FFTAllGrids();
  void                          updatePriorVar();

  void                          setPriorVar0(NRLib::Matrix priorVar0)   { priorVar0_ = priorVar0 ;}

  void                          setBackgroundParameters(FFTGrid  * muVp,
                                                        FFTGrid  * muVs,
                                                        FFTGrid  * muRho);

  void                          setBackgroundParametersInterval(const std::vector<NRLib::Grid<double> > & mu_parameters,
                                                                int                                       nx_pad,
                                                                int                                       ny_pad,
                                                                int                                       nz_pad);

  void                          setCovParameters(const std::vector<NRLib::Grid<double> > & cov_parameters,
                                                 int                                       nx_pad,
                                                 int                                       ny_pad,
                                                 int                                       nz_pad);

  void                          setCrCovParameters(const NRLib::Grid<double> & cr_cov_vp_vs,
                                                   const NRLib::Grid<double> & cr_cov_vp_rho,
                                                   const NRLib::Grid<double> & cr_cov_vs_rho,
                                                   int                         nx_pad,
                                                   int                         ny_pad,
                                                   int                         nz_pad);

  void                          copyBackgroundParameters(FFTGrid  * muVp,
                                                         FFTGrid  * muVs,
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

  fftw_real                   * extractParamCorrFromCovVp(int nzp) const;

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
                                                       const std::vector<float> & postCovVp00,
                                                       const std::vector<float> & postCovVs00,
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

  FFTGrid * muVp_;
  FFTGrid * muVs_ ;
  FFTGrid * muRho_  ;
  FFTGrid * covVp_;
  FFTGrid * covVs_ ;
  FFTGrid * covRho_  ;
  FFTGrid * crCovVpVs_;
  FFTGrid * crCovVpRho_ ;
  FFTGrid * crCovVsRho_  ;

  NRLib::Matrix priorVar0_;

};
#endif
