/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RESULT_H
#define RESULT_H

#include <math.h>
#include <string>

//#include "nrlib/segy/segy.hpp"
#include "src/definitions.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/modelsettings.h"

#include "src/multiintervalgrid.h"
#include "src/blockedlogscommon.h"
#include "src/commondata.h"
#include "src/seismicparametersholder.h"
#include "src/krigingdata3d.h"
#include "src/parameteroutput.h"

//#include "src/wavelet.h"
#include "src/wavelet1D.h"

class FFTGrid;
class Simbox;
class CommonData;
class ParameterOutput;
class Wavelet1D;

class CravaResult
{
public:
  CravaResult();

  ~CravaResult();

  void CombineResults(ModelSettings                        * model_settings,
                      CommonData                           * common_data,
                      std::vector<SeismicParametersHolder> & seismic_paramters_intervals);

  void CombineResult(StormContGrid         *& final_grid,
                     std::vector<FFTGrid *> & post_vp_intervals,
                     MultiIntervalGrid      * multi_interval_grid,
                     const std::vector<int> & erosion_priorities,
                     double                   dz_min);

  void WriteResults(ModelSettings * model_settings,
                    CommonData    * common_data);

  void WriteFilePriorCorrT(fftw_real   * prior_corr_T,
                           const int   & nzp,
                           const float & dt) const;

  void WriteFilePostVariances(const NRLib::Matrix      & post_var0,
                              const std::vector<float> & post_cov_vp00,
                              const std::vector<float> & post_cov_vs00,
                              const std::vector<float> & post_cov_rho00) const;

  void WriteFilePostCorrT(const std::vector<float> & post_cov,
                          const std::string        & sub_dir,
                          const std::string        & base_name) const;

  void WriteFilePostCovGrids(const Simbox & simbox) const;

  void WriteBlockedWells(const std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                         const ModelSettings                        * model_settings,
                         std::vector<std::string>                     facies_name,
                         std::vector<int>                             facies_label);

  void FindDzMin(const Simbox      & full_inversion_simbox,
                 MultiIntervalGrid * multi_interval_grid,
                 double            & dz);

  void ResampleSimple(std::vector<float>       & new_trace,
                      const std::vector<float> & old_trace);

  StormContGrid * CreateStormGrid(const Simbox & simbox,
                                  FFTGrid      * fft_grid);

  void WriteBackgrounds(const ModelSettings     * model_settings,
                        const Simbox            * simbox,
                        GridMapping             * depth_mapping,
                        const TraceHeaderFormat & thf);

  void ExpTransf(StormContGrid * grid);

  void ComputeSyntSeismic(const ModelSettings * model_settings,
                          const Simbox * simbox,
                          StormContGrid * vp,
                          StormContGrid * vs,
                          StormContGrid * rho);

  StormContGrid * ComputeSeismicImpedance(StormContGrid * vp,
                                          StormContGrid * vs,
                                          StormContGrid * rho,
                                          float        ** reflection_matrix,
                                          int             angle);

  //GET FUNCTIONS

  //SET FUNCTIONS
  void AddBackgroundVp(FFTGrid * vp)   { background_vp_intervals_.push_back(new FFTGrid(vp))   ;}
  void AddBackgroundVs(FFTGrid * vs)   { background_vs_intervals_.push_back(new FFTGrid(vs))   ;}
  void AddBackgroundRho(FFTGrid * rho) { background_rho_intervals_.push_back(new FFTGrid(rho)) ;}

private:

  //Resuls per interval
  std::vector<FFTGrid *> background_vp_intervals_;
  std::vector<FFTGrid *> background_vs_intervals_;
  std::vector<FFTGrid *> background_rho_intervals_;

  //Results combined
  FFTGrid                        * cov_vp_;
  FFTGrid                        * cov_vs_;
  FFTGrid                        * cov_rho_;
  FFTGrid                        * cr_cov_vp_vs_;
  FFTGrid                        * cr_cov_vp_rho_;
  FFTGrid                        * cr_cov_vs_rho_;

  //FFTGrid                        * post_vp_; //From avoinversion computePostMeanResidAndFFTCov()
  //FFTGrid                        * post_vs_;
  //FFTGrid                        * post_rho_;

  StormContGrid                  * post_vp_; //From avoinversion computePostMeanResidAndFFTCov()
  StormContGrid                  * post_vs_;
  StormContGrid                  * post_rho_;

  StormContGrid                  * post_vp_kriged_; //From avoinversion doPredictionKriging()
  StormContGrid                  * post_vs_kriged_;
  StormContGrid                  * post_rho_kriged_;

  StormContGrid                  * background_vp_;
  StormContGrid                  * background_vs_;
  StormContGrid                  * background_rho_;

  std::vector<FFTGrid *>           simulations_seed0_; //Vector over number of simulations
  std::vector<FFTGrid *>           simulations_seed1_;
  std::vector<FFTGrid *>           simulations_seed2_;

  fftw_real                      * corr_T_;
  float                          * corr_T_filtered_;

  NRLib::Matrix                    post_var0_;
  std::vector<float>               post_cov_vp00_;        // Posterior covariance in (i,j) = (0,0)
  std::vector<float>               post_cov_vs00_;
  std::vector<float>               post_cov_rho00_;

  std::vector<StormContGrid *>     synt_seismic_data_; //Vector angles
  std::vector<StormContGrid *>     synt_residuals_;

  StormContGrid                  * block_grid_;

  std::vector<StormContGrid *>     facies_prob_;
  StormContGrid                  * facies_prob_undef_;

  std::vector<StormContGrid *>     facies_prob_geo_;

  std::vector<StormContGrid*>      lh_cubes_;

  StormContGrid                  * quality_grid_;

  std::vector<Wavelet *>           wavelets_; //Vector angles
  float                         ** reflection_matrix_;

  int                              n_intervals_;

};

#endif
