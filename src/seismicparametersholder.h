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

  FFTGrid                     * GetMeanVp()                        { return meanVp_     ;}
  FFTGrid                     * GetMeanVs()                        { return meanVs_     ;}
  FFTGrid                     * GetMeanRho()                       { return meanRho_    ;}
  FFTGrid                     * GetCovVp()                         { return covVp_      ;}
  FFTGrid                     * GetCovVs()                         { return covVs_      ;}
  FFTGrid                     * GetCovRho()                        { return covRho_     ;}
  FFTGrid                     * GetCrCovVpVs()                     { return crCovVpVs_  ;}
  FFTGrid                     * GetCrCovVpRho()                    { return crCovVpRho_ ;}
  FFTGrid                     * GetCrCovVsRho()                    { return crCovVsRho_ ;}

  FFTGrid                     * GetPostVp()                        { return postVp_     ;}
  FFTGrid                     * GetPostVs()                        { return postVs_     ;}
  FFTGrid                     * GetPostRho()                       { return postRho_    ;}

  FFTGrid                     * GetPostVpKriging()                 { return postVpKriging_  ;}
  FFTGrid                     * GetPostVsKriging()                 { return postVsKriging_  ;}
  FFTGrid                     * GetPostRhoKriging()                { return postRhoKriging_ ;}

  NRLib::Matrix               & GetPostVar0()                      { return postVar0_     ;}
  std::vector<float>          & GetPostCovVp00()                   { return postCovVp00_  ;}
  std::vector<float>          & GetPostCovVs00()                   { return postCovVs00_  ;}
  std::vector<float>          & GetPostCovRho00()                  { return postCovRho00_ ;}

  fftw_real                   * GetCorrT()                         { return corr_T_          ;}
  float                       * GetCorrTFiltered()                 { return corr_T_filtered_ ;}

  FFTGrid                     * GetSimulationSeed0(int i)          { return simulations_seed0_[i] ;}
  FFTGrid                     * GetSimulationSeed1(int i)          { return simulations_seed1_[i] ;}
  FFTGrid                     * GetSimulationSeed2(int i)          { return simulations_seed2_[i] ;}

  std::vector<FFTGrid *>        GetSimulationsSeed0()              { return simulations_seed0_ ;}
  std::vector<FFTGrid *>        GetSimulationsSeed1()              { return simulations_seed1_ ;}
  std::vector<FFTGrid *>        GetSimulationsSeed2()              { return simulations_seed2_ ;}

  std::vector<FFTGrid *>        GetSyntSeismicData()               { return synt_seismic_data_ ;}
  std::vector<FFTGrid *>        GetSyntResiduals()                 { return synt_residuals_    ;}

  FFTGrid                     * GetBlockGrid()                     { return block_grid_        ;}

  std::vector<FFTGrid *>        GetFaciesProb()                    { return facies_prob_       ;}
  FFTGrid                     * GetFaciesProbUndefined()           { return facies_prob_undef_ ;}

  std::vector<FFTGrid *>        GetFaciesProbGeomodel()            { return facies_prob_geo_   ;}

  std::vector<FFTGrid *>        GetLHCube()                        { return lh_cube_           ;}

  FFTGrid                     * GetQualityGrid()                   { return quality_grid_      ;}

  void                          invFFTAllGrids();
  void                          invFFTCovGrids();
  void                          FFTCovGrids();
  void                          FFTAllGrids();
  void                          updatePriorVar();

  void SetPostVp(FFTGrid * vp)   { postVp_ = new FFTGrid(vp)   ;}
  void SetPostVs(FFTGrid * vs)   { postVs_ = new FFTGrid(vs)   ;}
  void SetPostRho(FFTGrid * rho) { postRho_ = new FFTGrid(rho) ;}

  void SetPostVpKriging(FFTGrid * vp)   { postVpKriging_ = new FFTGrid(vp)   ;}
  void SetPostVsKriging(FFTGrid * vs)   { postVsKriging_ = new FFTGrid(vs)   ;}
  void SetPostRhoKriging(FFTGrid * rho) { postRhoKriging_ = new FFTGrid(rho) ;}

  void SetPostVar0(NRLib::Matrix & postVar0)              { postVar0_     = postVar0     ;}
  void SetPostCovVp00(std::vector<float> & postCovVp00)   { postCovVp00_  = postCovVp00  ;}
  void SetPostCovVs00(std::vector<float> & postCovVs00)   { postCovVs00_  = postCovVs00  ;}
  void SetPostCovRho00(std::vector<float> & postCovRho00) { postCovRho00_ = postCovRho00 ;}

  void SetCorrT(fftw_real * corr_T)              { corr_T_          = corr_T          ;}
  void SetCorrTFiltered(float * corr_T_filtered) { corr_T_filtered_ = corr_T_filtered ;}

  void AddSimulationSeed0(FFTGrid * seed0) { simulations_seed0_.push_back(new FFTGrid(*seed0)) ;}
  void AddSimulationSeed1(FFTGrid * seed1) { simulations_seed1_.push_back(new FFTGrid(*seed1)) ;}
  void AddSimulationSeed2(FFTGrid * seed2) { simulations_seed2_.push_back(new FFTGrid(*seed2)) ;}

  void AddSyntSeismicData(FFTGrid * grid) { synt_seismic_data_.push_back(new FFTGrid(*grid)) ;}
  void AddSyntResiduals(FFTGrid * grid)   { synt_residuals_.push_back(new FFTGrid(*grid))    ;}

  void SetBlockGrid(FFTGrid * grid) { block_grid_ = new FFTGrid(*grid) ;}

  void AddFaciesProb(FFTGrid * facies_prob) { facies_prob_.push_back(new FFTGrid(facies_prob)) ;}
  void AddFaciesProbUndef(FFTGrid * facies_prob_undef) { facies_prob_undef_ = new FFTGrid(facies_prob_undef) ;}

  void AddFaciesProbGeomodel(FFTGrid * facies_prob) { facies_prob_geo_.push_back(new FFTGrid(facies_prob)) ;}

  void AddLHCube(FFTGrid * lh_cube) { lh_cube_.push_back(new FFTGrid(lh_cube)) ;}

  void SetQualityGrid(FFTGrid * grid) { quality_grid_ = new FFTGrid(*grid) ;}


  void                          setBackgroundParameters(FFTGrid  * meanVp,
                                                        FFTGrid  * meanVs,
                                                        FFTGrid  * meanRho);

  void                          setBackgroundParametersInterval(const std::vector<NRLib::Grid<float> *>  & mean_parameters,
                                                                int                                        nx_pad,
                                                                int                                        ny_pad,
                                                                int                                        nz_pad);

  void                          copyBackgroundParameters(FFTGrid  * meanVp,
                                                         FFTGrid  * meanVs,
                                                         FFTGrid  * meanRho);

  void                          setCorrelationParameters(const NRLib::Matrix       & priorVar0,
                                                         const std::vector<double> & priorCorrT,
                                                         const Surface             * priorCorrXY,
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

  fftw_real                   * computeCircCorrT(const std::vector<double> & priorCorrT,
                                                 const int                 & minIntFq,
                                                 const int                 & nzp) const;

  fftw_real                   * extractParamCorrFromCovVp(int nzp) const;

  void                          printPriorVariances(void) const;

  void                          printPostVariances(const NRLib::Matrix & postVar0) const;

  void                          getNextParameterCovariance(fftw_complex **& parVar) const;

  //void                          writeFilePriorCorrT(fftw_real   * priorCorrT,
  //                                                  const int   & nzp,
  //                                                  const float & dt) const;

  void                          writeFilePriorVariances(const ModelSettings      * modelSettings,
                                                        const std::vector<float> & priorCorrT,
                                                        const Surface            * priorCorrXY,
                                                        const float              & dt) const;

  //void                          writeFilePostVariances(const NRLib::Matrix      & postVar0,
  //                                                     const std::vector<float> & postCovVp00,
  //                                                     const std::vector<float> & postCovVs00,
  //                                                     const std::vector<float> & postCovRho00) const;

  //void                          writeFilePostCovGrids(Simbox const * simbox) const;

  std::vector<float>            createPostCov00(FFTGrid * postCov) const;

private:
  void                          createCorrGrids(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid);

  void                          initializeCorrelations(const Surface             * priorCorrXY,
                                                       const std::vector<double> & priorCorrT,
                                                       const float               & corrGradI,
                                                       const float               & corrGradJ,
                                                       const int                 & lowIntCut,
                                                       const int                 & nzp);

  FFTGrid                     * createFFTGrid(int nx,  int ny,  int nz,
                                              int nxp, int nyp, int nzp,
                                              bool fileGrid);


  void                          makeCircCorrTPosDef(fftw_real * circCorrT,
                                                    const int & minIntFq,
                                                    const int & nzp) const;

  void                         makeCorrXYPosDef(Surface         & priorCorrXY);

  float                         getOrigin(FFTGrid * grid) const;

  //void                          writeFilePostCorrT(const std::vector<float> & postCov,
  //                                                 const std::string        & subDir,
  //                                                 const std::string        & baseName) const;

  FFTGrid * meanVp_;
  FFTGrid * meanVs_;
  FFTGrid * meanRho_;
  FFTGrid * covVp_;
  FFTGrid * covVs_;
  FFTGrid * covRho_;
  FFTGrid * crCovVpVs_;
  FFTGrid * crCovVpRho_;
  FFTGrid * crCovVsRho_;

  NRLib::Matrix priorVar0_;

  //Stored variables for writing:

  FFTGrid * postVp_; //From avoinversion computePostMeanResidAndFFTCov()
  FFTGrid * postVs_;
  FFTGrid * postRho_;

  FFTGrid * postVpKriging_; //From avoinversion doPredictionKriging()
  FFTGrid * postVsKriging_;
  FFTGrid * postRhoKriging_;

  NRLib::Matrix      postVar0_;
  std::vector<float> postCovVp00_;        // Posterior covariance in (i,j) = (0,0)
  std::vector<float> postCovVs00_;
  std::vector<float> postCovRho00_;

  fftw_real * corr_T_;
  float     * corr_T_filtered_;

  std::vector<FFTGrid *> simulations_seed0_; //Vector over number of simulations;
  std::vector<FFTGrid *> simulations_seed1_;
  std::vector<FFTGrid *> simulations_seed2_;

  std::vector<FFTGrid *> synt_seismic_data_; //Vector over angles
  std::vector<FFTGrid *> synt_residuals_;

  FFTGrid * block_grid_;

  std::vector<FFTGrid *> facies_prob_;
  FFTGrid * facies_prob_undef_;
  std::vector<FFTGrid *> facies_prob_geo_;

  std::vector<FFTGrid *> lh_cube_;

  FFTGrid * quality_grid_;


};
#endif
