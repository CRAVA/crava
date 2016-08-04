/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "nrlib/iotools/logkit.hpp"
#include "src/definitions.h"
#include "src/analyzelog.h"
#include "src/modelsettings.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/blockedlogscommon.h"
#include "src/simbox.h"
#include "src/io.h"

Analyzelog::Analyzelog(const std::vector<NRLib::Well *>                        & wells,
                       const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs,
                       const std::vector<std::vector<NRLib::Grid<float> *> >   & background,
                       const std::vector<Simbox *>                             & interval_simboxes,
                       double                                                    dz_min,
                       const ModelSettings                                     * model_settings,
                       bool                                                      multi_zone_available,
                       std::string                                             & err_txt):
min_blocks_with_data_for_corr_estim_(model_settings->getMinBlocksForCorrEstimation()){

  enough_data_for_corr_estimation_  = false;
  n_lags_                           = 0;
  interval_name_                    = "";
  var_0_.resize(3,3);

  n_wells_    = static_cast<int>(mapped_blocked_logs.size());
  well_names_ = std::vector<std::string>();

  for (int i = 0; i < static_cast<int>(wells.size()); i++) {
    std::string well_name_tmp = wells[i]->GetWellName();
    for (int j = 0; j < n_wells_; j++) {

      if (mapped_blocked_logs.find(well_name_tmp) != mapped_blocked_logs.end()) {
        well_names_.push_back(well_name_tmp);
        j = n_wells_ - 1;
      }
    }
  }

  // Relative dz
  std::vector<double> dz_rel(interval_simboxes.size());
  for (size_t i = 0; i < interval_simboxes.size(); i++){
    dz_rel[i] = interval_simboxes[i]->getdz() / dz_min;
  }

  EstimateCorrelation(model_settings,
                      well_names_,
                      mapped_blocked_logs,
                      interval_simboxes,
                      enough_data_for_corr_estimation_,
                      regression_coef_vs_,
                      residual_variance_vs_,
                      dz_min,
                      dz_rel,
                      background,
                      n_lags_,
                      multi_zone_available,
                      err_txt);
}


Analyzelog::~Analyzelog(void)
{
}

void  Analyzelog::EstimateCorrelation(const ModelSettings                                       * model_settings,
                                      const std::vector<std::string>                            & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>          & mapped_blocked_logs,
                                      const std::vector<Simbox *>                               & interval_simboxes,
                                      bool                                                      & enough_data_for_corr_estimation,
                                      NRLib::Vector                                             & regression_coef,
                                      std::vector<double>                                       & residual_variance_vs,
                                      double                                                      dz_min,
                                      std::vector<double>                                         dz_rel,
                                      const std::vector<std::vector<NRLib::Grid<float> *> >     & background,
                                      int                                                       & n_lags,
                                      bool                                                        multi_zone_available,
                                      std::string                                               & err_txt)
{
  //
  // Covariance and correlation estimation with blocked logs
  //
  std::map<std::string, std::vector<double> > log_data_vp;
  std::map<std::string, std::vector<double> > log_data_vs;
  std::map<std::string, std::vector<double> > log_data_rho;

  //
  // Find out if there is enough data for correlation estimation
  //

  int n_blocks_tot = 0;
  for (int i = 0; i < n_wells_; i++){
    std::string well_name_tmp = well_names[i];
    for (size_t j = 0; j < interval_simboxes.size(); j++){
      std::string interval_name = interval_simboxes[j]->GetIntervalName();
      n_blocks_tot += mapped_blocked_logs.find(well_name_tmp)->second->GetNBlocksWithData(interval_name);
    }
  }
  enough_data_for_corr_estimation = (min_blocks_with_data_for_corr_estim_ <= n_blocks_tot);

  if (!enough_data_for_corr_estimation) {
    if (interval_simboxes.size() > 1) {
      TaskList::addTask("For estimation of prior correlations in all intervals at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
        " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n");
      LogKit::LogFormatted(LogKit::Low,"\nThere is not enough well data for estimation of prior correlations in all intervals.\n");
      err_txt += "There is not enough well data for estimation of prior correlations in all intervals. For estimation of prior correlations in all intervals at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
        " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n";
    }
    else {
      if (multi_zone_available == true)
        return;
      else {
        TaskList::addTask("For estimation of prior correlations at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
          " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n");
        LogKit::LogFormatted(LogKit::Low,"\nThere is not enough well data for estimation of prior correlations.\n");
        err_txt += "There is not enough well data for estimation of prior correlations. For estimation of prior correlations at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
          " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n";
      }
    }
    return;
  }
  else{
    //
    // Check if all wells have synthetic Vs logs
    //

    bool all_Vs_logs_synthetic      = true;
    bool all_Vs_logs_non_synthetic  = true;

    for(int i = 0; i < n_wells_; i++) {
      std::string well_name_tmp = well_names[i];
      bool synt_vs_log          = mapped_blocked_logs.find(well_name_tmp)->second->HasSyntheticVsLog();
      all_Vs_logs_synthetic     = all_Vs_logs_synthetic && synt_vs_log;
      all_Vs_logs_non_synthetic = all_Vs_logs_non_synthetic && !synt_vs_log;
    }
    if (all_Vs_logs_synthetic)
    {
      LogKit::LogFormatted(LogKit::Low,"\nThere are no nonsynthetic Vs logs available. Corr(Vp,Vs) is set to 0.7 and Corr(Vs,Rho) is set to 0.\n");
    }

    //
    // Check that the background model is consistent with the low frequency log data (part of CRA-257)
    //

    if(background[0][0]->GetN() == 0){
      std::vector<std::vector<NRLib::Grid<float> *> > temp_null;
      temp_null.resize(n_wells_);
      for(int j = 0; j < n_wells_; j++)
        temp_null[j].resize(3);
      EstimateLnData(log_data_vp, temp_null, well_names_, mapped_blocked_logs, interval_simboxes, "Vp", err_txt);
      EstimateLnData(log_data_vs, temp_null, well_names_, mapped_blocked_logs, interval_simboxes, "Vs", err_txt);
      EstimateLnData(log_data_rho, temp_null, well_names_, mapped_blocked_logs, interval_simboxes, "Rho", err_txt);
    }
    else{
      EstimateLnData(log_data_vp, background, well_names_, mapped_blocked_logs, interval_simboxes, "Vp", err_txt);
      EstimateLnData(log_data_vs, background, well_names_, mapped_blocked_logs, interval_simboxes, "Vs", err_txt);
      EstimateLnData(log_data_rho, background, well_names_, mapped_blocked_logs, interval_simboxes, "Rho", err_txt);
    }

    //EstimatePointVar0(point_var_0_, log_data_vp, log_data_vs, log_data_rho, err_txt);

    if (err_txt != "")
      return;

    CalculateNumberOfLags(n_lags, dz_min, dz_rel, interval_simboxes);

    n_lags = n_lags + 1; //H Lags in EstimateAutoCovarianceFunction may exceed n_lags

    //corr_T_.resize(n_lags, 0);
    auto_cov_.resize(n_lags);

    EstimateAutoCovarianceFunction(auto_cov_, well_names, mapped_blocked_logs, interval_simboxes, log_data_vp, log_data_vs, log_data_rho,
      all_Vs_logs_synthetic, all_Vs_logs_non_synthetic, regression_coef, residual_variance_vs, static_cast<float>(dz_min), n_lags, min_blocks_with_data_for_corr_estim_, max_lag_with_data_, err_txt);

    SetParameterCov(auto_cov_[0], var_0_, 3);

    /*
    EstimateCorrTAndVar0(corr_T_, var_0_,
                         log_data_vp, log_data_vs, log_data_rho,
                         allVsLogsAreSynthetic, dt,
                         n_lags_, max_nd,
                         err_txt);
    */

    CheckVariances(model_settings, var_0_, dz_min, err_txt);

  }
}

void  Analyzelog::CalculateNumberOfLags(int                                                   & max_nd,
                                        double                                                  dz_min,
                                        const std::vector<double>                             & dz_rel,
                                        const std::vector<Simbox *>                           & simboxes){
  //
  // Find the largest n: divide largest simbox lz by the minimum simbox dz
  //
  (void) dz_rel;
  max_nd = 0;
  for (size_t i = 0; i < simboxes.size(); i++){
    int nd = static_cast<int>(simboxes[i]->getMaxLz()/(dz_min) + 0.5);
    if (nd > max_nd)
      max_nd = nd;
  }
}

void Analyzelog::EstimateLnData(std::map<std::string, std::vector<double> >           & ln_data_map, //map from well to log
                                const std::vector<std::vector<NRLib::Grid<float> *> > & background,
                                const std::vector<std::string>                        & well_names,
                                const std::map<std::string, BlockedLogsCommon *>      & mapped_blocked_logs,
                                const std::vector<Simbox *>                           & interval_simboxes,
                                const std::string                                     & log_name,
                                std::string                                           & err_txt){

  double global_mean    = 0.0;
  int   count           = 0;
  int   nd              = 0;
  int   nd_tot          = 0;
  int   log_nr          = -1;
  //
  // The three vectors below contain data from all wells and all intervals
  //
  std::vector<double>  mean_vector;
  std::vector<double>  ln_data_vector;
  std::vector<double>  low_freq_log;

  assert(log_name == "Vp" || log_name == "Vs" || log_name == "Rho");

  for (int i = 0; i < n_wells_; i++) {
    std::vector<double>              temp_log(0);
    ln_data_map.insert(std::pair<std::string, std::vector<double> >(well_names[i], temp_log));

    std::vector<double> temp_low_freq_log =  mapped_blocked_logs.find(well_names[i])->second->GetContLogHighCutBackground(log_name);

    std::vector<double>   mean;
    nd = mapped_blocked_logs.find(well_names[i])->second->GetNumberOfBlocks();

    const std::vector<int>    & i_pos = mapped_blocked_logs.find(well_names[i])->second->GetIposVector();
    const std::vector<int>    & j_pos = mapped_blocked_logs.find(well_names[i])->second->GetJposVector();
    const std::vector<int>    & k_pos = mapped_blocked_logs.find(well_names[i])->second->GetKposVector();

    const std::vector<double> & x_pos_blocked = mapped_blocked_logs.find(well_names[i])->second->GetXposBlocked();
    const std::vector<double> & y_pos_blocked = mapped_blocked_logs.find(well_names[i])->second->GetYposBlocked();
    const std::vector<double> & z_pos_blocked = mapped_blocked_logs.find(well_names[i])->second->GetZposBlocked();

    std::vector<double> blocked_well_log;

    if(log_name      == "Vp"){
      blocked_well_log     = mapped_blocked_logs.find(well_names[i])->second->GetVpBlocked();
      log_nr               = 0;
    }
    else if(log_name == "Vs"){
      blocked_well_log     = mapped_blocked_logs.find(well_names[i])->second->GetVsBlocked();
      log_nr               = 1;
    }
    else if(log_name == "Rho"){
      blocked_well_log     = mapped_blocked_logs.find(well_names[i])->second->GetRhoBlocked();
      log_nr               = 2;
    }

    if (interval_simboxes.size() == 1) {

      if (background[0][log_nr]->GetN() > 0) {
        for(int n = 0; n < nd; n++) {
          mean_vector.push_back( (*background[0][log_nr])(i_pos[n], j_pos[n], k_pos[n]));
          mean.push_back( (*background[0][log_nr])(i_pos[n], j_pos[n], k_pos[n]));
        }
      }
      else{
        for(int n = 0; n < nd; n++){
          mean_vector.push_back(0);
          mean.push_back(0);
        }
      }

      for (int n = 0; n < nd; n++) {
        low_freq_log.push_back(temp_low_freq_log[n]);
        if (blocked_well_log[n] != RMISSING && mean[n] != RMISSING) {
          ln_data_map.find(well_names[i])->second.push_back(blocked_well_log[n] - mean[n]);
          global_mean   += blocked_well_log[n] - mean[n];
          ln_data_vector.push_back(blocked_well_log[n] - mean[n]);
          count++;
        }
        else {
          ln_data_vector.push_back(RMISSING);
          ln_data_map.find(well_names[i])->second.push_back(RMISSING);
        }
      }
    }
    //For multiple intervals the blocked logs are common for all intervals with an adjusted resolution, but the background models are created in interval-simbox resolutions
    else {
      //Get correct background values
      for (int n = 0; n < nd; n++) {

        //Find simbox
        double x = x_pos_blocked[n];
        double y = y_pos_blocked[n];
        double z = z_pos_blocked[n];
        int simbox_num = -1;
        for (size_t i = 0; i < interval_simboxes.size(); i++) {
          if (interval_simboxes[i]->IsPointBetweenVisibleSurfaces(x,y,z)) {
            simbox_num = static_cast<int>(i);
            break;
          }
        }

        //Find correct kpos
        int i = 0;
        int j = 0;
        int k = 0;
        if (simbox_num > -1)
          interval_simboxes[simbox_num]->getIndexes(x, y, z, i, j, k);

        if (simbox_num > -1 && background[simbox_num][log_nr]->GetN() > 0) { //if simbox_num = -1 then well_log should be RMISSING
          mean_vector.push_back( (*background[simbox_num][log_nr])(i_pos[n], j_pos[n], k));
          mean.push_back( (*background[simbox_num][log_nr])(i_pos[n], j_pos[n], k));
        }
        else {
          mean_vector.push_back(0);
          mean.push_back(0);
        }
      }

      //Subtract background
      for (int n = 0; n < nd; n++) {
        low_freq_log.push_back(temp_low_freq_log[n]);
        if (blocked_well_log[n] != RMISSING && mean[n] != RMISSING) {
          ln_data_map.find(well_names[i])->second.push_back(blocked_well_log[n] - mean[n]);
          global_mean   += blocked_well_log[n] - mean[n];
          ln_data_vector.push_back(blocked_well_log[n] - mean[n]);
          count++;
        }
        else {
          ln_data_vector.push_back(RMISSING);
          ln_data_map.find(well_names[i])->second.push_back(RMISSING);
        }
      }
    }

    nd_tot += nd;
  }

  //
  // Check consistency with background model
  //

  bool background_ok = CheckConsistencyBackground(ln_data_vector, mean_vector, low_freq_log, nd_tot);
  if(background_ok == false){
    err_txt += "The background model does not correspond with well log '"+ log_name +"'; prior correlation estimation failed.\n";
    TaskList::addTask("Check that the background model does not deviate significantly from well log '" + log_name + "'.\n");
  }

  //
  // Subtract global mean from data.
  // Erik N: This is left out as part of the new covariance estimation; ref Odd Kolbjørnsen
  //

  /*
  if(count > 0){
    global_mean /= count;
    for(int i = 0 ; i < n_wells_ ; i++){
      for(size_t k = 0 ; k < ln_data_map.find(well_names[i])->second.size(); k++){
        if (ln_data_map.find(well_names[i])->second[k] != RMISSING){
          ln_data_map.find(well_names[i])->second[k] -= global_mean;
        }
      }
    }
  }
  else{
    err_txt += std::string("Could not estimate global mean for log '" + log_name + "' because there are no data within the inversion interval.\n");
  }
  */
}

bool   Analyzelog::CheckConsistencyBackground(const std::vector<double>                              & ln_data_blocked,
                                              const std::vector<double>                              & background,
                                              const std::vector<double>                              & low_freq_log,
                                              int                                                     nd_tot){

  bool    background_ok    = true;
  double  sum_top          = 0;
  double  sum_base         = 0;
  double  ratio            = 0;


  for (int i = 0; i < nd_tot; i++){
    if (ln_data_blocked[i] != RMISSING){
      sum_top   += std::pow(ln_data_blocked[i] - background[i], 2);
      sum_base  += std::pow(ln_data_blocked[i] - low_freq_log[i], 2);
    }
  }

  ratio = sum_top / sum_base;

  //
  // The choice of using log(ratio) < 1 as a test statistic was made by
  // Ragnar Hauge and Odd Kolbjørnsen
  //

  if (log(ratio) > 1) {
    background_ok = false;
  }

  return background_ok;
}

//
// CRA-257: new implementation of estimation of autocovariance function
//
void Analyzelog::EstimateAutoCovarianceFunction(std::vector<NRLib::Matrix >                        & auto_cov,
                                                const std::vector<std::string>                     & well_names,
                                                const std::map<std::string, BlockedLogsCommon *>   & mapped_blocked_logs,
                                                const std::vector<Simbox *>                        & interval_simboxes,
                                                const std::map<std::string, std::vector<double> >  & log_data_vp,
                                                const std::map<std::string, std::vector<double> >  & log_data_vs,
                                                const std::map<std::string, std::vector<double> >  & log_data_rho,
                                                bool                                                 all_Vs_logs_synthetic,
                                                bool                                                 all_Vs_logs_non_synthetic,
                                                NRLib::Vector                                      & regression_coef,
                                                std::vector<double>                                & residual_variance_vs,
                                                float                                                min_dz,
                                                int                                                  max_nd,
                                                int                                                  min_blocks_with_data_for_corr_estim,
                                                int                                                & max_lag_with_data,
                                                std::string                                        & err_text)
{
  time_t timestart_tot, timeend_tot;
  time(&timestart_tot);
  max_lag_with_data = 0;
  // 1: Estimation of mean: Assumption: The data has mean 0, so there is no need to estimate it.

  // 2: Estimation of covariance
  auto_cov.resize(max_nd);
  std::vector<NRLib::Matrix> temp_auto_cov(max_nd);
  std::vector<NRLib::Matrix> count(max_nd);
  for(int i = 0; i < max_nd; i++){
    auto_cov[i].resize(3, 3);
    temp_auto_cov[i].resize(3,3);
    count[i].resize(3,3);
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
        auto_cov[i](j,k)      = 0.0;
        temp_auto_cov[i](j,k) = 0.0;
        count[i](j,k)         = 0;
      }
    }
  }
  std::vector<double> cov_vp_vp(max_nd, 0);
  std::vector<double> cov_vs_vs(max_nd, 0);
  std::vector<double> cov_rho_rho(max_nd, 0);

  //
  // 2.1: If there is a mixture of wells with synthetic and non-synthetic Vs:
  // 2.1.1 Use the wells with non-synthetic vs to calculate the linear regression of Vs on Vp and Rho.
  // We use a regression model with a constant term = 0 since we have already subtracted
  // the background model. This SHOULD result in elastic parameter vectors with mean ~ 0
  //
  if(!all_Vs_logs_synthetic && !all_Vs_logs_non_synthetic){
    int                                 n_data = 0;
    NRLib::Matrix                       A;
    NRLib::Vector                       b;
    std::vector<double>                 vp;
    std::vector<double>                 rho;
    std::vector<double>                 vs;
    for (size_t i = 0; i < well_names.size(); i++){
      if(mapped_blocked_logs.find(well_names[i])->second->HasSyntheticVsLog() == false){
        //int n_lags;
        for (size_t j = 0; j<interval_simboxes.size(); j++){
          const std::vector<double> well_log_vp     = log_data_vp.find(well_names[i])->second;
          const std::vector<double> well_log_rho    = log_data_rho.find(well_names[i])->second;
          const std::vector<double> well_log_vs     = log_data_vs.find(well_names[i])->second;
          //const std::vector<double> z_pos         = mapped_blocked_logs.find(well_names[i])->second->GetZposBlocked();
          // Loop over all blocks
          for (int k = 0; k < static_cast<float>(well_log_vp.size()); k++){
            if(well_log_vp[k] != RMISSING && well_log_rho[k] != RMISSING && well_log_vs[k] != RMISSING){
              //n_lags = static_cast<int>(std::floor(std::abs(z_pos[l] - z_pos[k])/min_dz + 0.5));
              n_data++;
              vp.push_back(well_log_vp[k]);
              rho.push_back(well_log_rho[k]);
              vs.push_back(well_log_vs[k]);
            }
          }
        }
      }
    }

    A.resize(n_data, 2);
    b.resize(n_data);
    for(int j = 0; j < n_data; j++){
      A(j, 0) = vp[j];
      A(j, 1) = rho[j];
      b(j)    = vs[j];
    }
    if (n_data > 0){
      regression_coef = Regress(A, b);

      //
      // 2.1.2 calculate residuals for vs and the resulting residual variance for the non-synthetic vs wells
      //
      residual_variance_vs.resize(max_nd, 0);
      std::vector<int> count_obs(vp.size(), 0);
      //double sum_resid_square = 0.0;
      double residual_k, residual_l;
      int last_obs = 0;
      int lag = 0;
      for (size_t i = 0; i < well_names.size(); i++){
        if(mapped_blocked_logs.find(well_names[i])->second->HasSyntheticVsLog() == false){
        //int n_lags;
          const std::vector<double> well_log_vp     = log_data_vp.find(well_names[i])->second;
          const std::vector<double> well_log_rho    = log_data_rho.find(well_names[i])->second;
          const std::vector<double> well_log_vs     = log_data_vs.find(well_names[i])->second;
          const std::vector<double> & x_pos         = mapped_blocked_logs.find(well_names[i])->second->GetXposBlocked();
          const std::vector<double> & y_pos         = mapped_blocked_logs.find(well_names[i])->second->GetYposBlocked();
          const std::vector<double> & z_pos         = mapped_blocked_logs.find(well_names[i])->second->GetZposBlocked();

          for (size_t j = 0; j<interval_simboxes.size(); j++){

            for (size_t k = 0; k < well_log_vp.size(); k++){
              for (size_t l = k; l < well_log_vp.size(); l++){
                if (well_log_vp[k] != RMISSING && well_log_vp[l] != RMISSING
                  && well_log_rho[k] != RMISSING && well_log_rho[l] != RMISSING
                  && well_log_vs[k] != RMISSING && well_log_vs[l] != RMISSING){
                  double z_k_rel = (z_pos[k] - interval_simboxes[j]->getTop(x_pos[k], y_pos[k]))/interval_simboxes[j]->getRelThick(x_pos[k], y_pos[k]);
                  double z_l_rel = (z_pos[l] - interval_simboxes[j]->getTop(x_pos[l], y_pos[l]))/interval_simboxes[j]->getRelThick(x_pos[l], y_pos[l]);
                  lag = static_cast<int>(std::floor(std::abs(z_k_rel - z_l_rel)/min_dz + 0.5));
                  residual_k = regression_coef(0)*well_log_vp[k] + regression_coef(1)*well_log_rho[k] - well_log_vs[k];
                  residual_l = regression_coef(0)*well_log_vp[l] + regression_coef(1)*well_log_rho[l] - well_log_vs[l];
                  residual_variance_vs[lag] += residual_k*residual_l;
                  count_obs[lag] += 1;
                  if (lag > last_obs)
                    last_obs = lag;
                }
              }
            }
          }
        }
      }
      for (int i = 0; i <= last_obs; i++){
        if (count_obs[i] > 2)
          residual_variance_vs[i] /= (count_obs[i] - 2);
        else
          residual_variance_vs[i] = 0;
      }

      // 2.1.3 Linear downscaling of error covariance
      for (int i = 0; i <= last_obs; i++){
        residual_variance_vs[i] *= static_cast<float>(last_obs - i)/(last_obs);
      }
      for (int i = last_obs + 1; i < max_nd; i++)
        residual_variance_vs[i] = 0.0;
    }
  }

  //
  // 2.2: Estimate autocovariance using both synthetic and non-synthetic well logs
  // OBS The autocovariance function does not necessarily result in symmetric
  // matrices for each time lag, i.e. cov(h)(vp, vs) != cov(h)(vs, vp)
  // but cov(h)(vp,vs) = cov(-h)(vs,vp) and cov(h)(vs,vp) = cov(-h)(vp,vs)
  //
  for (size_t i = 0; i < well_names.size(); i++){
    time_t timestart, timeend;
    time(&timestart);
    for (size_t j = 0; j < interval_simboxes.size(); j++){
      std::string interval_name = interval_simboxes[j]->GetIntervalName();
      size_t nd = mapped_blocked_logs.find(well_names[i])->second->GetNBlocksWithData(interval_name);
      const std::vector<double> & x_pos = mapped_blocked_logs.find(well_names[i])->second->GetXposBlocked();
      const std::vector<double> & y_pos = mapped_blocked_logs.find(well_names[i])->second->GetYposBlocked();
      const std::vector<double> & z_pos = mapped_blocked_logs.find(well_names[i])->second->GetZposBlocked();
      std::vector<double> z(nd);
      const std::vector<double> log_vp   = log_data_vp.find(well_names[i])->second;
      const std::vector<double> log_vs   = log_data_vs.find(well_names[i])->second;
      const std::vector<double> log_rho  = log_data_rho.find(well_names[i])->second;
      //
      // 2.2.1 Add autocovariance data
      //
      int lag = 0;

      for (size_t k = 0; k < nd; k++){
        for (size_t l = k; l < nd; l++){

          if(log_vp[k] != RMISSING || log_vs[k] != RMISSING || log_rho[k] != RMISSING){
            if (log_vp[l] != RMISSING || log_vs[l] != RMISSING || log_rho[l] != RMISSING) {
              double xk = x_pos[k];
              double yk = y_pos[k];
              double xl = x_pos[l];
              double yl = y_pos[l];
              double z_k_rel = (z_pos[k] - interval_simboxes[j]->getTop(xk, yk))/interval_simboxes[j]->getRelThick(xk, yk);
              double z_l_rel = (z_pos[l] - interval_simboxes[j]->getTop(xl, yl))/interval_simboxes[j]->getRelThick(xl, yl);
              lag = static_cast<int>(std::floor(std::abs(z_k_rel - z_l_rel)/min_dz + 0.5));

              // cov(vp_k, vp_l)
              if(log_vp[k] != RMISSING && log_vp[l] != RMISSING){
                if (lag > max_lag_with_data)
                  max_lag_with_data = lag;
                temp_auto_cov[lag](0,0) += log_vp[k]*log_vp[l];
                count[lag](0,0) += 1;
              }
              // cov(rho_k, rho_l)
              if(log_rho[k] != RMISSING && log_rho[l] != RMISSING){
                if (lag > max_lag_with_data)
                  max_lag_with_data = lag;
                temp_auto_cov[lag](2,2) += log_rho[k]*log_rho[l];
                count[lag](2,2) += 1;
              }
              // cov(vp_k, rho_l)
              if(log_vp[k] != RMISSING && log_rho[l] != RMISSING){
                temp_auto_cov[lag](0,2) += log_vp[k]*log_rho[l];
                count[lag](0,2) += 1;
                if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                  temp_auto_cov[lag](2,0) += log_vp[k]*log_rho[l];
                  count[lag](2,0)         += 1;
                }
              }
              // cov(rho_k, vp_l)
              if(log_rho[k] != RMISSING && log_vp[l] != RMISSING){
                temp_auto_cov[lag](2,0) += log_rho[k]*log_vp[l];
                count[lag](2,0) += 1;
                if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                  temp_auto_cov[lag](0,2) += log_rho[k]*log_vp[l];
                  count[lag](0,2)         += 1;
                }
              }
              //
              // If this Vs log is synthetic and there exist real Vs logs: use regression coefficients
              //
              if(!all_Vs_logs_synthetic && mapped_blocked_logs.find(well_names[i])->second->HasSyntheticVsLog() == true){
                // Use the relation Vs = a*Vp + b*Rho + e, where e is iid
                // cov[t](vs_i, vs_j) = cov[t](a*vp_k + b*rho_k + e_k, a*vp_l + b*rho_l + e_l) = a*a*cov(vp_k,vp_l) + a*b*cov(vp_k, rho_l) + a*b*cov(vp_l, rho_k) + b*b*cov(rho_k, rho_l) + I(k = l) var(e)
                if(log_vp[k] != RMISSING && log_rho[k] != RMISSING && log_rho[l] != RMISSING && log_vp[l] != RMISSING){
                  double vs_k = regression_coef(0)*log_vp[k] + regression_coef(1)*log_rho[k];
                  double vs_l = regression_coef(0)*log_vp[l] + regression_coef(1)*log_rho[l];

                  temp_auto_cov[lag](1,1) += vs_k*vs_l + residual_variance_vs[lag];
                  //if (k == l)
                  //  temp_auto_cov[lag](1,1) += var_vs_resid;
                  count[lag](1,1) += 1;
                }
                // cov[l-k](vp, vs) = cov[l-k](vp_k, a*vp_l + b*rho_l + e_l) = a*autocov[l-k](vp_k, vp_l) + b*autocov[l-k](vp_k, rho_l)
                if (log_vp[k] != RMISSING && log_vp[l] != RMISSING && log_rho[l] != RMISSING){
                  temp_auto_cov[lag](0,1) += regression_coef(0)*log_vp[k]*log_vp[l] + regression_coef(1)*log_vp[k]*log_rho[l];
                  count[lag](0,1)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](1,0) += regression_coef(0)*log_vp[k]*log_vp[l] + regression_coef(1)*log_vp[k]*log_rho[l];
                    count[lag](1,0)         += 1;
                  }
                }
                // cov[l-k](vs, vp) = a*cov[l-k](vp_k, vp_l) + b*cov[l-k](rho_k, vp_l)
                if (log_vp[k] != RMISSING && log_rho[k] != RMISSING && log_vp[l] != RMISSING){
                  temp_auto_cov[lag](1,0) += regression_coef(0)*log_vp[k]*log_vp[l] + regression_coef(1)*log_rho[k]*log_vp[l];
                  count[lag](1,0)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](0,1) += regression_coef(0)*log_vp[k]*log_vp[l] + regression_coef(1)*log_rho[k]*log_vp[l];
                    count[lag](0,1)         += 1;
                  }
                }
                // cov[l-k](rho_k, vs_l) = cov[l-k](rho_k, a*vp_l + b*rho_l + e_l) = a*cov[l-k](rho_k,vp_l) + b*cov[l-k](rho_k, rho_l)
                if (log_rho[k] != RMISSING && log_vp[l] != RMISSING && log_rho[l] != RMISSING){
                  temp_auto_cov[lag](1,2) += regression_coef(0)*log_rho[k]*log_vp[l] + regression_coef(1)*log_rho[k]*log_rho[l];
                  count[lag](1,2)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](2,1) += regression_coef(0)*log_rho[k]*log_vp[l] + regression_coef(1)*log_rho[k]*log_rho[l];
                    count[lag](2,1)         += 1;
                  }
                }
                // cov[l-k](vs_k, rho_l) = a*cov[l-k](vp_k, rho_l) + b*cov[l-k](rho_k, rho_l)
                if (log_rho[k] != RMISSING && log_vp[k] != RMISSING && log_rho[l] != RMISSING){
                  temp_auto_cov[lag](2,1) += regression_coef(0)*log_vp[k]*log_rho[l] + regression_coef(1)*log_rho[k]*log_rho[l];
                  count[lag](2,1)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](1,2) += regression_coef(0)*log_vp[k]*log_rho[l] + regression_coef(1)*log_rho[k]*log_rho[l];
                    count[lag](1,2)         += 1;
                  }
                }
              }
              //
              // Non-synthetic Vs log
              //
              else if(mapped_blocked_logs.find(well_names[i])->second->HasSyntheticVsLog() == false){
                // cov[t](vs, vs)
                if(log_vs[k] != RMISSING && log_vs[l] != RMISSING){
                  if (lag > max_lag_with_data)
                    max_lag_with_data = lag;
                  temp_auto_cov[lag](1,1) += log_vs[k]*log_vs[l];
                  count[lag](1,1)         += 1;
                }
                // cov[t](vp, vs)
                if(log_vp[k] != RMISSING && log_vs[l] != RMISSING){
                  temp_auto_cov[lag](0,1) += log_vp[k]*log_vs[l];
                  count[lag](0,1)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](1,0) += log_vp[k]*log_vs[l];
                    count[lag](1,0)         += 1;
                  }
                }
                // cov[t](vs, vp)
                if(log_vs[k] != RMISSING && log_vp[l] != RMISSING){
                  temp_auto_cov[lag](1,0) += log_vs[k]*log_vp[l];
                  count[lag](1,0)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](0,1) += log_vs[k]*log_vp[l];
                    count[lag](0,1)         += 1;
                  }
                }
                // cov[t](vs, rho)
                if(log_vs[k] != RMISSING && log_rho[l] != RMISSING){
                  temp_auto_cov[lag](1,2) += log_vs[k]*log_rho[l];
                  count[lag](1,2)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](2,1) += log_vs[k]*log_rho[l];
                    count[lag](2,1)         += 1;
                  }
                }
                // cov[t](rho, vs)
                if(log_rho[k] != RMISSING && log_vs[l] != RMISSING){
                  temp_auto_cov[lag](2,1) += log_rho[k]*log_vs[l];
                  count[lag](2,1)         += 1;
                  if (lag == 0){ // In lag 0, the autocov matrix is symmetric
                    temp_auto_cov[lag](1,2) += log_rho[k]*log_vs[l];
                    count[lag](1,2)         += 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    time(&timeend);
    //long int time = static_cast<long int>(timeend - timestart);
    //printf("\nWell %s processed in %ld seconds.",well_names[i].c_str(),time);
  }

  //
  // 2.2.2 Calculate diagonal autocovariances
  //
  if(all_Vs_logs_synthetic){
    LogKit::LogFormatted(LogKit::Low,"\n Estimating Vs autocovariance as 2 * Vp autocovariance, corr(Vp, Vs) = 0.7 and corr(Vs, Rho) = 0.\n");
  }
  bool vp_fail  = false;
  bool vs_fail  = false;
  bool rho_fail = false;
  for (int i = 0; i < max_lag_with_data && err_text==""; i++){
    // If there is not enough data to estimate autocovariance within the first ~50 lags, return an error
    //double n_vp = count[i](0,0);
    //double n_vs = count[i](1,1);
    //double n_rho = count[i](2,2);

    //double c_vp = temp_auto_cov[i](0,0);
    //double c_vs = temp_auto_cov[i](1,1);
    //double c_rho = temp_auto_cov[i](2,2);

    if (i < static_cast<int>(min_blocks_with_data_for_corr_estim/2)){
      if(count[i](0,0)<2)
      {
        err_text += "Not enough well data within simulation area to estimate autocovariance of Vp for time lag = "+ CommonData::ConvertIntToString(i) +", dz = " + CommonData::ConvertFloatToString(min_dz) +" .\n";
      }
      if(count[i](2,2)<2)
      {
        err_text += "Not enough well data within simulation area to estimate autocovariance of Rho for time lag "+ CommonData::ConvertIntToString(i) +", dz = "+ CommonData::ConvertFloatToString(min_dz) + ".\n";
      }
      if (err_text == ""){
        temp_auto_cov[i](0,0) = temp_auto_cov[i](0,0) / (count[i](0,0) - 1);
        temp_auto_cov[i](2,2) = temp_auto_cov[i](2,2) / (count[i](2,2) - 1);
      }
      if(all_Vs_logs_synthetic){
        temp_auto_cov[i](1,1) = 2*temp_auto_cov[i](0,0);
      } else if (count[i](1,1) > 1){
        temp_auto_cov[i](1,1) = temp_auto_cov[i](1,1) / (count[i](1,1) - 1);
      } else{
        err_text += "Not enough well data within simulation area to estimate autocovariance of Vs for time lag "+ CommonData::ConvertIntToString(i) +", dz = "+ CommonData::ConvertFloatToString(min_dz) + ".\n";
      }
    }
    // If there is not enough data to estimate autocovariance for lags > 50, set to 0 and return a warning
    else{
      if(count[i](0,0)<2 || vp_fail == true){
        temp_auto_cov[i](0,0) = 0.0;
        if(vp_fail == false) {
          LogKit::LogFormatted(LogKit::Low,"\n WARNING: Vp autocovariance for time lag " + CommonData::ConvertIntToString(i) +", dz = " +
                                            CommonData::ConvertFloatToString(min_dz) +", could not be estimated and is set to 0.\n", i, min_dz);
          vp_fail = true;
        }
      }
      else{
        temp_auto_cov[i](0,0) = temp_auto_cov[i](0,0) / (count[i](0,0) - 1);
      }
      if(count[i](2,2)<2 || rho_fail == true){
        temp_auto_cov[i](2,2) = 0.0;
        if(rho_fail == false) {
          LogKit::LogFormatted(LogKit::Low,"\n WARNING: Rho autocovariance for time lag " + CommonData::ConvertIntToString(i) +", dz = " +
                                            CommonData::ConvertFloatToString(min_dz) +", could not be estimated and is set to 0.\n", i, min_dz);
          rho_fail = true;
        }
      }
      else{
        temp_auto_cov[i](2,2) = temp_auto_cov[i](2,2) / (count[i](2,2) - 1);
      }
      if(all_Vs_logs_synthetic){
        temp_auto_cov[i](1,1) = 2*temp_auto_cov[i](0,0);
      } else if (count[i](1,1) > 1 && vs_fail == false){
        temp_auto_cov[i](1,1) = temp_auto_cov[i](1,1) / (count[i](1,1) - 1);
      } else{
        temp_auto_cov[i](1,1) = 0.0;
        if(vs_fail == false) {
          LogKit::LogFormatted(LogKit::Low,"\n WARNING: Vs autocovariance for time lag " + CommonData::ConvertIntToString(i) +", dz = " +
                                            CommonData::ConvertFloatToString(min_dz) +", could not be estimated and is set to 0.\n", i, min_dz);
          vs_fail = true;
        }
      }
    }
  }

  //
  // 2.2.3 Calculate cross autocovariances
  //
  if(err_text == ""){

    // Default covariance between Vp and Vs. Value is explained in Jira issue CRA-220.
    for (int i = 0; i < max_lag_with_data; i++){

      // cov[t](vp, vs) = cov[-t](vs,vp)
      if(count[i](0,1) > 1){
        temp_auto_cov[i](0,1) /= (count[i](0,1) - 1);
      }
      else{
        temp_auto_cov[i](0,1) = std::sqrt(2.0)*temp_auto_cov[i](0,0)*0.7;
      }

      // cov[t](vs, vp) = cov[-t](vp,vs)
      if (count[i](1,0) > 1){
        temp_auto_cov[i](1,0) /= (count[i](1,0) - 1);
      }
      else{
        temp_auto_cov[i](1,0) = std::sqrt(2.0)*temp_auto_cov[i](0,0)*0.7;
      }

      // cov[t](vp, rho) = cov[-t](rho,vp)
      if(count[i](0,2) > 1){
        temp_auto_cov[i](0,2) = temp_auto_cov[i](0,2)/(count[i](0,2) - 1);
      }
      else{
        temp_auto_cov[i](0,2)  = 0.0;
      }

      // cov[t](rho, vp) = cov[-t](vp,rho)
      if (count[i](2,0) > 1){
        temp_auto_cov[i](2,0) /= (count[i](2,0) - 1);
      }
      else{
        temp_auto_cov[i](2,0)  = 0.0;
      }

      // cov[t](vs, rho) = cov[-t](rho, vs)
      if(count[i](1,2) > 1){
        temp_auto_cov[i](1,2) = temp_auto_cov[i](1,2)/(count[i](1,2) - 1);
      }
      else{
        temp_auto_cov[i](1,2)  = 0.0;
      }

      // cov[t](rho, vs) = cov[-t](vs, rho)
      if(count[i](2,1) > 1){
        temp_auto_cov[i](2,1) /= (count[i](2,1) - 1);
      }
      else{
        temp_auto_cov[i](2,1)  = 0.0;
      }
    }

  }



  // CRA-257, Erik N: Should this still be done?
  /************************************************************
  //
  // x) Scale corTT correlation to account for blocking effects
  //
  for(i=0;i<n+1;i++)
    corTT[i]/=corTT[0];
  //
  // Replace missing-values after last nonzero element with zero
  //
  int   nend;
  int   iprev,inext,nipol;
  float cprev,cnext,cipol;
  i=n;
  while(nTT[i]==0 && i>0)
  {
    CorrT[i]=0.0;
    i--;
  }
  **************************************************************/

  //int   nend;
  std::vector<std::vector<int> > i_prev(3), i_next(3);
  std::vector<std::vector<double> > c_prev(3), c_next(3);
  for (int i = 0; i < 3; i++){
    i_prev[i].resize(3,0);
    i_next[i].resize(3,0);
    c_prev[i].resize(3,0.0);
    c_next[i].resize(3,0.0);
  }

  double   nipol;
  double   cipol;

  //
  // 3) Smooth observations by making a weighted average of each observation
  //    and the linear interpolation of its two nearest neighbours.
  //

  for(int j=0 ; j<3; j++){
    for(int k=0; k<3; k++){
      auto_cov[0](j,k) = temp_auto_cov[0](j,k);
    }
  }

  for(int i=1; i<max_lag_with_data; i++){
    for(int j=0 ; j<3; j++){
      for(int k=0; k<3; k++){
        //
        // Vp-Vp, Vp-Rho, Rho-Rho or non-synthetic logs
        //
        if ((j == 0 && k == 0) || (j == 0 && k ==2) || (j ==2 && k == 2) || (j ==2 && k==0) || !all_Vs_logs_synthetic){
          if(count[i-1](j,k)>0)                 // Find previous estimated autocov with data
            i_prev[j][k]=i-1;
          if(i_next[j][k]<i+1){
            i_next[j][k]++;
            while(count[i_next[j][k]](j,k) == 0 && i_next[j][k]<max_lag_with_data)   // Find next estimated autocov with data
              i_next[j][k]++;
          }
          if(i_next[j][k] == max_lag_with_data)
            i_next[j][k]--;

          c_prev[j][k]      = temp_auto_cov[i_prev[j][k]](j,k);
          c_next[j][k]      = temp_auto_cov[i_next[j][k]](j,k);
          nipol             = 0.5*(count[i_prev[j][k]](j,k)+count[i_next[j][k]](j,k));
          cipol             = ((i_next[j][k]-i)*c_prev[j][k] + (i-i_prev[j][k])*c_next[j][k])/(i_next[j][k]-i_prev[j][k]);
          if((count[i](j,k) + nipol) > 0)
            auto_cov[i](j,k)  = (count[i](j,k)*temp_auto_cov[i](j,k) + nipol*cipol)/(count[i](j,k)+nipol);
          //auto_cov[i](k,j)  = auto_cov[i](j,k);

          //LogKit::LogFormatted(LogKit::Low," i nTT,corTT  iprev,cprev  inext,cnext  nipol,cipol  %4d %4d%8.3f %4d%8.3f %4d%8.3f    %4d%8.3f\n",
          //                   i,nTT[i],corTT[i],iprev,cprev,inext,cnext,nipol,cipol);
        }
        //
        // Vp-Vs, Vs-Vs and Vs-Rho if all Vs logs are synthetic
        //
        else {
          if ((j == 0 && k == 1) || (j == 1 && k == 1) || (j == 1 && k == 0)){
            // use same weights as for Vp correlations
            i_prev[j][k]      = i_prev[0][0];
            i_next[j][k]      = i_next[0][0];
            c_prev[j][k]      = temp_auto_cov[i_prev[j][k]](j,k);
            c_next[j][k]      = temp_auto_cov[i_next[j][k]](j,k);
            nipol             = 0.5*(count[i_prev[j][k]](0,0)+count[i_next[j][k]](0,0));
            cipol             = ((i_next[j][k]-i)*c_prev[j][k] + (i-i_prev[j][k])*c_next[j][k])/(i_next[j][k]-i_prev[j][k]);
            if((count[i](j,k) + nipol) > 0)
              auto_cov[i](j,k)  = (count[i](0,0)*temp_auto_cov[i](j,k) + nipol*cipol)/(count[i](0,0)+nipol);
            //auto_cov[i](k,j)  = auto_cov[i](j,k);
          }
          else if((j == 1 && k ==2) || (k == 1 && j == 1)){
            i_prev[j][k]      = i-1;
            i_next[j][k]      = i;
            auto_cov[i](j,k)  = 0.0;
            //auto_cov[i](k,j)  = 0.0;
          }
        }
      }
    }
  }

  // Last time step with data
  for(int j=0 ; j<3; j++){
    for(int k=0; k<3; k++){
      auto_cov[max_lag_with_data-1](j,k) = (count[max_lag_with_data-1](j,k)*temp_auto_cov[max_lag_with_data-1](j,k) + i_prev[j][k]*c_prev[j][k]*0.5)/(count[max_lag_with_data-1](j,k) + i_prev[j][k]);
      //auto_cov[max_lag_with_data-1](k,j) = auto_cov[max_lag_with_data-1](j,k);
    }
  }

  // Set remaining autocovariances to 0
  for(int i=max_lag_with_data; i<max_nd; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
        auto_cov[i](j,k) = 0.0;
      }
    }
  }


  //
  // 4) Downscale correlations linearly
  //

  int last_element = max_nd - 1;
  if (max_lag_with_data < max_nd){
    last_element = max_lag_with_data;
  }

  for(int j=0 ; j<3; j++){
    for(int k=0; k<3; k++){
      for(int i=0; i <= last_element; i++){
        auto_cov[i](j,k) *= 1 - static_cast<float>(i)/(last_element);
        //auto_cov[i](k,j) = auto_cov[i](j,k);
      }
    }
  }

  //
  // 5) Print to file autocov.dat
  //
  if(ModelSettings::getDebugLevel() > 0) {
    std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), std::string("Autocov.dat"));
    std::ofstream file;
    NRLib::OpenWrite(file, fileName);
    file << " Time increment is: " << min_dz << "\n"
         << "   i      n_VpVp    cov(Vp_0,Vp_i)      n_VsVs    cov(Vs_0,Vs_i)      nRR    cov(R_0,R_i)  \n"
         << "----------------------------------------------------------------------------------------\n";
    for(int i=0; i<max_nd; i++){
      file << std::fixed
           << std::right
           << std::setprecision(3)
           << std::setw(4) << i        << " "
           << std::setw(8) << count[i](0,0)   << " "
           << std::setw(8) << auto_cov[i](0,0) << " "
           << std::setw(8) << count[i](1,1)   << " "
           << std::setw(8) << auto_cov[i](1,1) << " "
           << std::setw(8) << count[i](2,2)   << " "
           << std::setw(8) << auto_cov[i](2,2) << "    "
           << std::endl;
    }
      file.close();
  }

  //
  // Return replace n by nonzero elements nend?
  //
  // n=nend;
  //
  time(&timeend_tot);
  //LogKit::LogFormatted(LogKit::Low,"\nEstimated parameter autocovariance in %d seconds.\n",
  //                 static_cast<int>(timeend_tot-timestart_tot));

}

void  Analyzelog::SetParameterCov(const NRLib::Matrix                           & auto_cov,
                                  NRLib::Matrix                                 & var_0,
                                  int                                             n_params){
  for (int i=0; i<n_params; i++){
    for (int j=0; j<n_params; j++){
      var_0(i,j) = auto_cov(i,j);
    }
  }
}

void  Analyzelog::CheckVariances(const ModelSettings      * model_settings,
                                 const NRLib::Matrix      & var_0,
                                 double                     dz,
                                 std::string              & err_txt)
{
  //| These min and max values below are used for consistency check. If a variance
  //| is outside these ranges there is probably a problem with the log.
  //|
  //| The limits are for point variances. The minimum allowed variance
  //| for parameters will be scaled with 1/dt*dt
  float minVarVp  = model_settings->getVarVpMin();
  float maxVarVp  = model_settings->getVarVpMax();
  float minVarVs  = model_settings->getVarVsMin();
  float maxVarVs  = model_settings->getVarVsMax();
  float minVarRho = model_settings->getVarRhoMin();
  float maxVarRho = model_settings->getVarRhoMax();

  //
  // We scale the minimum variances allowed with 1/dt since the variance decreases with increasing dt..
  //
  std::string err_txt_tmp = "";
  if (var_0(0,0) < minVarVp/dz || var_0(0,0) > maxVarVp)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vp variance "  << var_0(0,0)
      << " is outside allowed interval Min=" << minVarVp/dz << " Max=" << maxVarVp << "\n";
    err_txt_tmp += o.str();
  }
  if (var_0(1,1) < minVarVs/dz || var_0(1,1)  > maxVarVs)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vs variance "  << var_0(1,1)
      << " is outside allowed interval Min=" << minVarVs/dz << " Max=" << maxVarVs << "\n";
    err_txt_tmp += o.str();
  }
  if (var_0(2,2)  < minVarRho/dz || var_0(2,2)  > maxVarRho)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Rho variance "  << var_0(2,2)
      << " is outside allowed interval Min=" << minVarRho/dz << " Max=" << maxVarRho << "\n";
    err_txt_tmp += o.str();
  }

  if (err_txt_tmp != "") {
    err_txt += err_txt_tmp;
    err_txt += "Allowed minimum and maximum variance limits can adjusted under <well> and <allowed-parameter-values> in the model file.\n";
  }

  if (err_txt != "")
  {
    LogKit::LogFormatted(LogKit::Low,"\n --------------------------------------------------------------------");
    LogKit::LogFormatted(LogKit::Low,"\n                           ln Vp     ln Vs    ln Rho ");
    LogKit::LogFormatted(LogKit::Low,"\n Parameter variances:   %.2e  %.2e  %.2e (used by program)",var_0(0,0) ,var_0(1,1) ,var_0(2,2) );
    LogKit::LogFormatted(LogKit::Low,"\n --------------------------------------------------------------------\n");
  }
}

NRLib::Vector Analyzelog::Regress(const NRLib::Matrix & A,
                                  const NRLib::Vector & b){

  NRLib::Matrix AtA           = NRLib::transpose(A)*A;
  assert(AtA.numCols() == AtA.numRows());

  //Make symmetric matrix
  NRLib::SymmetricMatrix AtA_sym(A.numCols());
  for(int i = 0; i < AtA.numCols(); i++){
    for(int j = i; j < AtA.numCols(); j++)
      AtA_sym(i,j) = AtA(i,j);
  }

  NRLib::Vector Atb           = NRLib::transpose(A)*b;
  NRLib::Vector x(2);

  NRLib::CholeskySolve(AtA_sym, Atb, x);

  return x;
}
