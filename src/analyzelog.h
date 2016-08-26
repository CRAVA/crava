/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef ANALYZELOG_H
#define ANALYZELOG_H

#include "src/simbox.h"
#include "src/modelsettings.h"
#include "nrlib/grid/grid.hpp"
#include "nrlib/well/well.hpp"
#include "src/multiintervalgrid.h"
#include "src/commondata.h"

class BlockedLogsCommon;
class Background;

class Analyzelog {

public:

  Analyzelog(const std::vector<NRLib::Well *>                        & wells,
             const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs,
             const std::vector<std::vector<NRLib::Grid<float> *> >   & background,
             const std::vector<Simbox *>                             & interval_simboxes,
             double                                                    dz_min,
             const ModelSettings                                     * model_settings,
             bool                                                      multi_zone_available,
             std::string                                             & err_txt);


  ~Analyzelog(void);

  const NRLib::Matrix                    & GetVar0(void)                     const { return var_0_                           ;}
  int                                      GetNumberOfLags(void)             const { return n_lags_                          ;}
  int                                      GetMaxLagWithNonZeroAutoCovData() const { return max_lag_with_data_               ;}
  bool                                     GetEnoughData()                   const { return enough_data_for_corr_estimation_ ;}
  const std::vector<NRLib::Matrix>       & GetAutoCovariance()               const { return auto_cov_                        ;}

private:

  void            EstimateLnData(std::map<std::string, std::vector<double> >             & log_data,
                                 const std::vector<std::vector<NRLib::Grid<float> *> >  & background,
                                 const std::vector<std::string>                         & well_names,
                                 const std::map<std::string, BlockedLogsCommon *>       & mapped_blocked_logs,
                                 const std::vector<Simbox *>                            & interval_simboxes,
                                 const std::string                                      & log_name,
                                 std::string                                            & err_txt);

    void          EstimateCorrelation(const ModelSettings                                      * model_settings,
                                      const std::vector<std::string>                           & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>         & mapped_blocked_logs,
                                      const std::vector<Simbox *>                              & interval_simboxes,
                                      bool                                                     & enough_data_for_corr_estimation,
                                      NRLib::Vector                                            & regression_coef,
                                      std::vector<double>                                      & residual_variance_vs,
                                      double                                                     dz_min,
                                      std::vector<double>                                        dz_rel,
                                      const std::vector<std::vector<NRLib::Grid<float> *> >    & background,
                                      int                                                      & n_lags,
                                      bool                                                       multi_zone_available,
                                      std::string                                              & errTxt);

  void            EstimateAutoCovarianceFunction(std::vector<NRLib::Matrix >                        & auto_cov,
                                                 const std::vector<std::string>                     & well_names,
                                                 const std::map<std::string, BlockedLogsCommon *>   & mapped_blocked_logs,
                                                 const std::vector<Simbox *>                        & interval_simboxes,
                                                 const std::map<std::string, std::vector<double> >   & log_data_vp,
                                                 const std::map<std::string, std::vector<double> >   & log_data_vs,
                                                 const std::map<std::string, std::vector<double> >   & log_data_rho,
                                                 bool                                                 all_Vs_logs_synthetic,
                                                 bool                                                 all_Vs_logs_non_synthetic,
                                                 NRLib::Vector                                      & regression_coef,
                                                 std::vector<double>                                & residual_variance_vs,
                                                 float                                                dt,
                                                 int                                                  max_nd,
                                                 int                                                  min_blocks_with_data_for_corr_estim,
                                                 int                                                & max_lag_with_data,
                                                 std::string                                        & err_text);

  void            SetParameterCov(const NRLib::Matrix                           & auto_cov,
                                  NRLib::Matrix                                 & var_0,
                                  int                                             n_params);

  void            estimate(const ModelSettings  * modelSettings,
                           Background           * background,
                           std::string          & errTxt);

  void            estimatePointVar0(float       ** Var0,
                                    float       ** lnDataAlpha,
                                    float       ** lnDataBeta,
                                    float       ** lnDataRho,
                                    std::string & errTxt);

  void            estimateCorrTAndVar0(float       * CorrT,
                                       float      ** Var0,
                                       float      ** lnDataAlpha,
                                       float      ** lnDataBeta,
                                       float      ** lnDataRho,
                                       bool          allVsLogsAreSynthetic,
                                       float         dt,
                                       int           n,
                                       int           maxnd,
                                       std::string & errTxt);

  void            readMeanData(FFTGrid *cube, int nd, const double *xpos, const double *ypos,
                               const double *zpos, float *meanValue);

  void            CalculateNumberOfLags(int                                                   & max_nd,
                                        double                                                  dz_min,
                                        const std::vector<double>                             & dz_rel,
                                        const std::vector<Simbox *>                           & simboxes);

  void            findConstructedVsLogs(void);

  void            CheckVariances(const ModelSettings    * model_settings,
                                 const NRLib::Matrix    & var_0,
                                 double                   dz,
                                 std::string            & err_txt);

  bool            CheckConsistencyBackground(const std::vector<double>                              & ln_data_blocked,
                                              const std::vector<double>                              & background,
                                              const std::vector<double>                              & low_freq_log,
                                              int                                                     nd_tot);

  NRLib::Vector   Regress(const NRLib::Matrix                             & A,
                          const NRLib::Vector                             & b);

  //--------------------------------------------------------------------------------------------

  int                                                   min_blocks_with_data_for_corr_estim_;
  bool                                                  enough_data_for_corr_estimation_;       ///< Vector with size == n_intervals
  std::vector<std::string>                              well_names_;
  int                                                   n_wells_;             // Number of wells
  NRLib::Vector                                         regression_coef_vs_;  // Regression coefficients for Vs ~ Vp + Rho in case of mixture of synthetic and non-synt wells (CRA-257)
  std::vector<double>                                   residual_variance_vs_;

  std::string                                           interval_name_;
  std::vector<NRLib::Matrix >                           auto_cov_;
  NRLib::Matrix                                         var_0_;
  int                                                   n_lags_;
  int                                                   max_lag_with_data_;

};
#endif
