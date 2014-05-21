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

class FFTGrid;
class Simbox;
class CommonData;
class ParameterOutput;

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
                    CommonData    * common_data,
                    const Simbox  & simbox);

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

  //GET FUNCTIONS

  //SET FUNCTIONS
  void AddBackgroundVp(std::string interval_name, FFTGrid * vp)   { background_vp_intervals_[interval_name]  = new FFTGrid(vp)  ;}
  void AddBackgroundVs(std::string interval_name, FFTGrid * vs)   { background_vs_intervals_[interval_name]  = new FFTGrid(vs)  ;}
  void AddBackgroundRho(std::string interval_name, FFTGrid * rho) { background_rho_intervals_[interval_name] = new FFTGrid(rho) ;}

private:

  //Resuls per interval
  std::map<std::string, FFTGrid *> background_vp_intervals_;
  std::map<std::string, FFTGrid *> background_vs_intervals_;
  std::map<std::string, FFTGrid *> background_rho_intervals_;

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

  FFTGrid                        * post_vp_kriged_test_; //From avoinversion doPredictionKriging()
  //FFTGrid                        * post_vs_kriged_;
  //FFTGrid                        * post_rho_kriged_;

  StormContGrid                  * post_vp_kriged_; //From avoinversion doPredictionKriging()
  StormContGrid                  * post_vs_kriged_;
  StormContGrid                  * post_rho_kriged_;

  FFTGrid                        * background_vp_;
  FFTGrid                        * background_vs_;
  FFTGrid                        * background_rho_;

  std::vector<FFTGrid *>           simulations_seed0_; //Vector over number of simulations
  std::vector<FFTGrid *>           simulations_seed1_;
  std::vector<FFTGrid *>           simulations_seed2_;

  fftw_real                      * corr_T_;
  float                          * corr_T_filtered_;

  NRLib::Matrix                    post_var0_;
  std::vector<float>               post_cov_vp00_;        // Posterior covariance in (i,j) = (0,0)
  std::vector<float>               post_cov_vs00_;
  std::vector<float>               post_cov_rho00_;

  std::vector<FFTGrid *>           synt_seismic_data_; //Vector angles
  std::vector<FFTGrid *>           synt_residuals_;

  FFTGrid                        * block_grid_;

  std::vector<FFTGrid *>           facies_prob_;
  FFTGrid                        * facies_prob_undef_;

  std::vector<FFTGrid *>           facies_prob_geo_;

  std::vector<FFTGrid *>           lh_cubes_;

  FFTGrid                        * quality_grid_;

  int                              n_intervals_;

















    //std::map<std::string, FFTGrid *> post_vp_intervals_;
  //std::map<std::string, FFTGrid *> post_vs_intervals_;
  //std::map<std::string, FFTGrid *> post_rho_intervals_;

  //std::map<std::string, FFTGrid *> post_vp_kriging_intervals_;
  //std::map<std::string, FFTGrid *> post_vs_kriging_intervals_;
  //std::map<std::string, FFTGrid *> post_rho_kriging_intervals_;

  //std::map<std::string, fftw_real *> corr_T_intervals_;
  //std::map<std::string, fftw_real *> corr_T_filtered_intervals_;

  //std::map<std::string, NRLib::Matrix>       post_var0_intervals_;
  //std::map<std::string, std::vector<float> > post_cov_vp00_intervals_;        // Posterior covariance in (i,j) = (0,0)
  //std::map<std::string, std::vector<float> > post_cov_vs00_intervals_;
  //std::map<std::string, std::vector<float> > post_cov_rho00_intervals_;

  //std::map<std::string, FFTGrid *> cov_vp_intervals_;
  //std::map<std::string, FFTGrid *> cov_vs_intervals_;
  //std::map<std::string, FFTGrid *> cov_rho_intervals_;
  //std::map<std::string, FFTGrid *> cr_cov_vp_vs_intervals_;
  //std::map<std::string, FFTGrid *> cr_cov_vp_rho_intervals_;
  //std::map<std::string, FFTGrid *> cr_cov_vs_rho_intervals_;

  //std::map<std::string, std::vector<FFTGrid *> > simulations_seed0_intervals_; //Vector over number of simulations
  //std::map<std::string, std::vector<FFTGrid *> > simulations_seed1_intervals_;
  //std::map<std::string, std::vector<FFTGrid *> > simulations_seed2_intervals_;

  //std::map<std::string, std::vector<FFTGrid *> > synt_seismic_data_intervals_; //Vector over angles
  //std::map<std::string, std::vector<FFTGrid *> > synt_residuals_intervals_;

  //std::map<std::string, FFTGrid *> block_grid_intervals_;

  //std::map<std::string, std::vector<FFTGrid *> > facies_prob_intervals_; //Vector over facies
  //std::map<std::string, FFTGrid *> facies_prob_undef_intervals_;

  //std::map<std::string, std::vector<FFTGrid *> > facies_prob_geo_intervals_; //Vector over facies

  //std::map<std::string, std::vector<FFTGrid *> > lh_cubes_intervals_;

  //std::map<std::string, FFTGrid *> quality_grid_intervals_;


};

#endif
