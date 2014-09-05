/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RESULT_H
#define RESULT_H

#include <math.h>
#include <string>

#include "src/definitions.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/modelsettings.h"
#include "nrlib/well/well.hpp"

class FFTGrid;
class Simbox;
class CommonData;
class ParameterOutput;
class Wavelet1D;
class MultiIntervalGrid;
class BlockedLogsCommon;

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

  float GetResampledTraceValue(const std::vector<float> & resampled_trace,
                               const double             & dz_resampled,
                               const double             & top,
                               const double             & global_z,
                               const double             & dz_final);

  double GetResampledTraceValue(const std::vector<double> & resampled_trace,
                                const std::vector<double> & z_pos_resampled,
                                //const double              & dz_resampled,
                                //const double              & top,
                                const double              & global_z); //z-value for this cell in the final blocked log
                                //const double              & dz_final);

  void CombineBlockedLogs(std::map<std::string, BlockedLogsCommon *>                     & blocked_logs_output,
                          const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                          MultiIntervalGrid                                              * multi_interval_grid,
                          const Simbox                                                   * output_simbox);

  void GetWellLogContributed(std::vector<double>       & log_out,
                             const std::vector<double> & log_old,
                             int                         start,
                             int                         end);

  void InterpolateMissing(std::vector<double> & well_log);

  void ExtrapolateLog(std::vector<double> & well_log);

  void ResampleLog(std::vector<double>                                            & final_log,
                   std::vector<std::vector<double> >                              & old_log_interval, //vector interval
                   const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                   MultiIntervalGrid                                              * multi_interval_grid,
                   const BlockedLogsCommon                                        * blocked_log_final,
                   std::string                                                      well_name,
                   float                                                            res_fac);

  void ResampleTrace(std::vector<double> & old_trace,
                     std::vector<double> & new_trace,
                     const float           res_fac);

  void CombineTraces(std::vector<double>                                            & final_log,
                     const BlockedLogsCommon                                        * blocked_log_final,
                     //const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                     //std::string                                                    & well_name,
                     MultiIntervalGrid                                              * multiple_interval_grid,
                     const std::vector<std::vector<double> >                        & resampled_logs,
                     const std::vector<std::vector<double> >                        & z_pos_resampled);

  void WriteResults(ModelSettings           * model_settings,
                    CommonData              * common_data,
                    SeismicParametersHolder & seismic_parameters);

  void WriteFilePriorCorrT(fftw_real   * prior_corr_T,
                           const int   & nzp,
                           const float & dt,
                           std::string   interval_name = "") const;

  void WriteFilePostVariances(const NRLib::Matrix      & post_var0,
                              const std::vector<float> & post_cov_vp00,
                              const std::vector<float> & post_cov_vs00,
                              const std::vector<float> & post_cov_rho00,
                              std::string                interval_name = "") const;

  void WriteFilePostCorrT(const std::vector<float> & post_cov,
                          const std::string        & sub_dir,
                          const std::string        & base_name) const;

  void WriteFilePostCovGrids(const ModelSettings * model_settings,
                             const Simbox        & simbox,
                             std::string           interval_name = "") const;

  void WriteBlockedWells(const std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                         const ModelSettings                              * model_settings,
                         std::vector<std::string>                           facies_name,
                         std::vector<int>                                   facies_label);

  void WriteWells(const std::vector<NRLib::Well> & wells,
                  const ModelSettings      * model_settings);

  StormContGrid * CreateStormGrid(const Simbox & simbox,
                                  FFTGrid      * fft_grid);

  //StormContGrid * CreateStormGrid(const Simbox       & simbox,
  //                                NRLib::Grid<float> * grid);

  void CreateStormGrid(StormContGrid & grid_new,
                       FFTGrid       * fft_grid);

  void WriteBackgrounds(const ModelSettings     * model_settings,
                        const Simbox            * simbox,
                        StormContGrid           * background_vp,
                        StormContGrid           * background_vs,
                        StormContGrid           * background_rho,
                        GridMapping             * depth_mapping,
                        const TraceHeaderFormat & thf);

  void ExpTransf(StormContGrid * grid);

  void ExpTransf(FFTGrid * grid);

  void LogTransf(FFTGrid * grid);

  void ComputeSyntSeismic(const ModelSettings          * model_settings,
                          const Simbox                 * simbox,
                          std::vector<Wavelet *>       & wavelets,
                          StormContGrid                * vp,
                          StormContGrid                * vs,
                          StormContGrid                * rho,
                          std::vector<StormContGrid *> & synt_seis_data,
                          std::vector<StormContGrid *> & synt_residuals);

  StormContGrid * ComputeSeismicImpedance(StormContGrid       * vp,
                                          StormContGrid       * vs,
                                          StormContGrid       * rho,
                                          const NRLib::Matrix & reflection_matrix,
                                          int                   angle) const;

  void GenerateSyntheticSeismicLogs(std::vector<Wavelet *>                     & wavelet,
                                    std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                    const NRLib::Matrix                        & reflection_matrix,
                                    const Simbox                               & simbox);

  void GenerateWellOptSyntSeis(ModelSettings                              * model_settings,
                               CommonData                                 * common_data,
                               std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                               const Simbox                               & simbox,
                               const NRLib::Matrix                        & reflection_matrix);

  //void WriteCravaGrids(ModelSettings           * model_settings,
  //                     CommonData              * common_data,
  //                     const Simbox            & output_simbox,
  //                     SeismicParametersHolder & seismic_parameters_intervals);

  //GET FUNCTIONS

  //SET FUNCTIONS
  void AddBackgroundVp(FFTGrid * vp)   { background_vp_intervals_.push_back(new FFTGrid(vp))   ;}
  void AddBackgroundVs(FFTGrid * vs)   { background_vs_intervals_.push_back(new FFTGrid(vs))   ;}
  void AddBackgroundRho(FFTGrid * rho) { background_rho_intervals_.push_back(new FFTGrid(rho)) ;}

  void AddBlockedLogs(const std::map<std::string, BlockedLogsCommon *> & blocked_logs);

  void SetBgBlockedLogs(const std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs) { bg_blocked_logs_ = bg_blocked_logs ;}

private:

  //Resuls per interval
  std::vector<FFTGrid *>                                   background_vp_intervals_;
  std::vector<FFTGrid *>                                   background_vs_intervals_;
  std::vector<FFTGrid *>                                   background_rho_intervals_;

  std::vector<std::map<std::string, BlockedLogsCommon *> > blocked_logs_intervals_;

  std::vector<fftw_real *>                                 corr_T_;
  std::vector<float *>                                     corr_T_filtered_;

  std::vector<NRLib::Matrix>                               post_var0_;
  std::vector<std::vector<float> >                         post_cov_vp00_;        // Posterior covariance in (i,j) = (0,0)
  std::vector<std::vector<float> >                         post_cov_vs00_;
  std::vector<std::vector<float> >                         post_cov_rho00_;

  //Results combined
  StormContGrid                                          * cov_vp_;
  StormContGrid                                          * cov_vs_;
  StormContGrid                                          * cov_rho_;
  StormContGrid                                          * cr_cov_vp_vs_;
  StormContGrid                                          * cr_cov_vp_rho_;
  StormContGrid                                          * cr_cov_vs_rho_;

  StormContGrid                                          * post_vp_; //From avoinversion computePostMeanResidAndFFTCov()
  StormContGrid                                          * post_vs_;
  StormContGrid                                          * post_rho_;

  StormContGrid                                          * post_vp_kriged_; //From avoinversion doPredictionKriging()
  StormContGrid                                          * post_vs_kriged_;
  StormContGrid                                          * post_rho_kriged_;

  StormContGrid                                          * background_vp_;
  StormContGrid                                          * background_vs_;
  StormContGrid                                          * background_rho_;

  std::map<std::string, BlockedLogsCommon *>               blocked_logs_;
  std::map<std::string, BlockedLogsCommon *>               bg_blocked_logs_;

  std::vector<StormContGrid *>                             simulations_seed0_; //Vector over number of simulations
  std::vector<StormContGrid *>                             simulations_seed1_;
  std::vector<StormContGrid *>                             simulations_seed2_;

  std::vector<StormContGrid *>                             synt_seismic_data_; //Vector angles
  std::vector<StormContGrid *>                             synt_residuals_;

  StormContGrid                                          * block_grid_;

  std::vector<StormContGrid *>                             facies_prob_;
  StormContGrid                                          * facies_prob_undef_;

  std::vector<StormContGrid *>                             facies_prob_geo_;

  std::vector<StormContGrid *>                              lh_cubes_;

  StormContGrid                                          * quality_grid_;

  std::vector<StormContGrid *>                             trend_cubes_; //vector trend_parameters

  std::vector<Wavelet *>                                   wavelets_; //Vector angles //Wavelet from common_data based on estimation simbox
  NRLib::Matrix                                            reflection_matrix_;

  bool                                                     write_crava_;
  int                                                      n_intervals_;
};

#endif
