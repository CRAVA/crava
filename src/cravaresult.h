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

class FFTGrid;
class Simbox;

class CravaResult
{
public:
  CravaResult();

  ~CravaResult();

  void CombineResults(MultiIntervalGrid * multi_interval_grid);

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

  //GET FUNCTIONS

  //SET FUNCTIONS

  void AddPostVp(std::string interval_name, FFTGrid * vp)     { post_vp_intervals_[interval_name]  = new FFTGrid(vp)  ;}
  void AddPostVs(std::string interval_name, FFTGrid * vs)     { post_vs_intervals_[interval_name]  = new FFTGrid(vs)  ;}
  void AddPostRho(std::string interval_name, FFTGrid * rho)   { post_rho_intervals_[interval_name] = new FFTGrid(rho) ;}

  void AddCorrT(std::string interval_name, fftw_real * corr_T)              { corr_T_intervals_[interval_name]          = corr_T          ;}
  void AddCorrTFiltered(std::string interval_name, float * corr_T_filtered) { corr_T_filtered_intervals_[interval_name] = corr_T_filtered ;}

  void AddPostVar0(std::string interval_name, NRLib::Matrix & post_var0)    { post_var0_intervals_[interval_name]       = post_var0       ;}
  void AddPostCovVp00(std::string interval_name, std::vector<float> & post_cov_vp00) { post_cov_vp00_intervals_[interval_name] = post_cov_vp00 ;}
  void AddPostCovVs00(std::string interval_name, std::vector<float> & post_cov_vs00) { post_cov_vs00_intervals_[interval_name] = post_cov_vs00 ;}
  void AddPostCovRho00(std::string interval_name, std::vector<float> & post_cov_rho00) { post_cov_rho00_intervals_[interval_name] = post_cov_rho00 ;}

  void AddCovVp(std::string interval_name, FFTGrid * cov_vp) { cov_vp_intervals_[interval_name] = new FFTGrid(cov_vp) ;}
  void AddCovVs(std::string interval_name, FFTGrid * cov_vs) { cov_vs_intervals_[interval_name] = new FFTGrid(cov_vs) ;}
  void AddCovRho(std::string interval_name, FFTGrid * cov_rho) { cov_rho_intervals_[interval_name] = new FFTGrid(cov_rho) ;}
  void AddCrCovVpVs(std::string interval_name, FFTGrid * cr_cov_vp_vs) { cr_cov_vp_vs_intervals_[interval_name] = new FFTGrid(cr_cov_vp_vs) ;}
  void AddCrCovVpRho(std::string interval_name, FFTGrid * cr_cov_vp_rho) { cr_cov_vp_rho_intervals_[interval_name] = new FFTGrid(cr_cov_vp_rho) ;}
  void AddCrCovVsRho(std::string interval_name, FFTGrid * cr_cov_vs_rho) { cr_cov_vs_rho_intervals_[interval_name] = new FFTGrid(cr_cov_vs_rho) ;}

private:

  std::map<std::string, FFTGrid *> post_vp_intervals_;
  std::map<std::string, FFTGrid *> post_vs_intervals_;
  std::map<std::string, FFTGrid *> post_rho_intervals_;

  std::map<std::string, fftw_real *> corr_T_intervals_;
  std::map<std::string, fftw_real *> corr_T_filtered_intervals_;

  std::map<std::string, NRLib::Matrix>       post_var0_intervals_;
  std::map<std::string, std::vector<float> > post_cov_vp00_intervals_;        // Posterior covariance in (i,j) = (0,0)
  std::map<std::string, std::vector<float> > post_cov_vs00_intervals_;
  std::map<std::string, std::vector<float> > post_cov_rho00_intervals_;

  std::map<std::string, FFTGrid *> cov_vp_intervals_;
  std::map<std::string, FFTGrid *> cov_vs_intervals_;
  std::map<std::string, FFTGrid *> cov_rho_intervals_;
  std::map<std::string, FFTGrid *> cr_cov_vp_vs_intervals_;
  std::map<std::string, FFTGrid *> cr_cov_vp_rho_intervals_;
  std::map<std::string, FFTGrid *> cr_cov_vs_rho_intervals_;


  FFTGrid * cov_vp_;
  FFTGrid * cov_vs_;
  FFTGrid * cov_rho_;
  FFTGrid * cr_cov_vp_vs_;
  FFTGrid * cr_cov_vp_rho_;
  FFTGrid * cr_cov_vs_rho_;

  FFTGrid * post_vp_;
  FFTGrid * post_vs_;
  FFTGrid * post_rho_;

  fftw_real * corr_T_;
  float     * corr_T_filtered_;

  NRLib::Matrix      post_var0_;
  std::vector<float> post_cov_vp00_;        // Posterior covariance in (i,j) = (0,0)
  std::vector<float> post_cov_vs00_;
  std::vector<float> post_cov_rho00_;

  int n_intervals_;

};

#endif
