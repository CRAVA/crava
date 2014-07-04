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

// CRA-257: New correlation estimation routine
/*
Analyzelog::Analyzelog(const std::vector<NRLib::Well>                          & wells,
                       const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs_for_correlation,
                       const std::vector<NRLib::Grid<float> *>                 & background,
                       const Simbox                                            * interval_simbox,
                       const ModelSettings                                     * model_settings,
                       std::string                                             & err_txt):
min_blocks_with_data_for_corr_estim_(model_settings->getMinBlocksForCorrEstimation()){

  enough_data_for_corr_estimation_  = false;
  n_lags_                           = 0;
  interval_name_                    = interval_simbox->GetIntervalName();
  point_var_0_.resize(3,3);
  var_0_.resize(3,3);

  //mapped_blocked_logs_for_correlation_ = mapped_blocked_logs_for_correlation;
  //simbox_       = interval_simbox;
  //wells_        = wells;
  n_wells_       = static_cast<int>(mapped_blocked_logs_for_correlation.size());
  well_names_.resize(n_wells_);
  for (int i = 0; i < n_wells_; i++){
    well_names_[i] = wells[i].GetWellName();
  }


  EstimateCorrelation(model_settings,
                      wells,
                      well_names_,
                      mapped_blocked_logs_for_correlation,
                      interval_name_,
                      interval_simbox,
                      enough_data_for_corr_estimation_,
                      background,
                      err_txt);

}
*/

// CRA-257: New correlation estimation routine
Analyzelog::Analyzelog(const std::vector<NRLib::Well>                          & wells,
                       const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs_for_correlation,
                       const std::vector<std::vector<NRLib::Grid<float> *> >   & background,
                       const std::vector<Simbox *>                             & interval_simboxes,
                       double                                                    dz_min,
                       const ModelSettings                                     * model_settings,
                       std::string                                             & err_txt):
min_blocks_with_data_for_corr_estim_(model_settings->getMinBlocksForCorrEstimation()){

  enough_data_for_corr_estimation_  = false;
  n_lags_                           = 0;
  interval_name_                    = "";
  //point_var_0_.resize(3,3);
  var_0_.resize(3,3);

  //mapped_blocked_logs_for_correlation_ = mapped_blocked_logs_for_correlation;
  //simbox_       = interval_simbox;
  //wells_        = wells;
  n_wells_       = static_cast<int>(mapped_blocked_logs_for_correlation.size());
  well_names_.resize(n_wells_);
  for (int i = 0; i < n_wells_; i++){
    well_names_[i] = wells[i].GetWellName();
  }

  // Relative dz
  std::vector<double> dz_rel(interval_simboxes.size());
  for (size_t i = 0; i < interval_simboxes.size(); i++){
    dz_rel[i] = interval_simboxes[i]->getdz() / dz_min;
  }

  EstimateCorrelation(model_settings,
                      wells,
                      well_names_,
                      mapped_blocked_logs_for_correlation,
                      interval_simboxes,
                      enough_data_for_corr_estimation_,
                      regression_coef_vs_,
                      var_vs_resid_,
                      dz_min,
                      dz_rel,
                      background,
                      n_lags_,
                      err_txt);
}

/*
Analyzelog::Analyzelog(std::vector<WellData *>  & wells,
                       Background               * background,
                       const Simbox             * simbox,
                       const ModelSettings      * modelSettings,
                       std::string              & errTxt)
{
  pointVar0_ = new float*[3];
  for(int i=0 ; i<3 ; i++)
    pointVar0_[i] = new float[3];

  Var0_ = new float*[3];
  for(int i=0 ; i<3 ; i++)
    Var0_[i] = new float[3];

  CorrT_        = NULL;
  numberOfLags_ = 0;

  simbox_       = simbox;
  wells_        = wells;
  nwells_       = modelSettings->getNumberOfWells();

  estimate(modelSettings,
           background,
           errTxt);
}
*/

Analyzelog::~Analyzelog(void)
{
}

// CRA-257: New implementation of correlation estimation for intervals
void  Analyzelog::EstimateCorrelation(const ModelSettings                                       * model_settings,
                                      const std::vector<NRLib::Well>                            & wells,
                                      const std::vector<std::string>                            & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>          & mapped_blocked_logs_for_correlation,
                                      const std::vector<Simbox *>                               & interval_simboxes,
                                      bool                                                      & enough_data_for_corr_estimation,
                                      NRLib::Vector                                             & regression_coef,
                                      double                                                    & var_vs_resid,
                                      double                                                      dz_min,
                                      std::vector<double>                                         dz_rel,
                                      const std::vector<std::vector<NRLib::Grid<float> *> >     & background,
                                      int                                                       & n_lags,
                                      std::string                                               & err_txt){
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
  for (size_t i = 0; i < wells.size(); i++){
    std::string well_name_tmp = wells[i].GetWellName();
    for (size_t j = 0; j < interval_simboxes.size(); j++){
      std::string interval_name = interval_simboxes[j]->GetIntervalName();
      n_blocks_tot += mapped_blocked_logs_for_correlation.find(well_name_tmp)->second->GetNBlocksWithData(interval_name);
    }
  }
  enough_data_for_corr_estimation = (min_blocks_with_data_for_corr_estim_ <= n_blocks_tot);

  if(!enough_data_for_corr_estimation){
    if (interval_simboxes.size() > 1){
      TaskList::addTask("For estimation of prior correlations in all intervals at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
        " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n");
      LogKit::LogFormatted(LogKit::Low,"\nThere is not enough well data for estimation of prior correlations in all intervals.\n");
    }
    else{
      std::string interval_name = interval_simboxes[0]->GetIntervalName();
      TaskList::addTask("For estimation of prior correlations in interval '" + interval_name +  "', at least "+ CommonData::ConvertIntToString(min_blocks_with_data_for_corr_estim_) +
        " blocks are needed; currently there are "+ CommonData::ConvertIntToString(n_blocks_tot) +". Either increase the number of layers or add wells.\n");
      LogKit::LogFormatted(LogKit::Low,"\nThere is not enough well data for estimation of prior correlations in interval '"+ interval_name + "'.\n");
    }
    return;
  }
  else{
    //
    // Check if all wells have synthetic Vs logs
    //

    bool all_Vs_logs_synthetic      = true;
    bool all_Vs_logs_non_synthetic  = true;

    for(int i=0 ; i<n_wells_ ; i++)
    {
      std::string well_name_tmp = wells[i].GetWellName();
      bool synt_vs_log          = mapped_blocked_logs_for_correlation.find(well_name_tmp)->second->HasSyntheticVsLog();
      all_Vs_logs_synthetic     = all_Vs_logs_synthetic && synt_vs_log;
      all_Vs_logs_non_synthetic = all_Vs_logs_non_synthetic && !synt_vs_log;
    }
    if (all_Vs_logs_synthetic)
    {
      LogKit::LogFormatted(LogKit::Low,"\nThere are no nonsynthetic Vs logs available. Corr(Vp,Vs) and Corr(Vs,Rho) are set to 0.7.\n");
    }

    //
    // Check that the background model is consistent with the low frequency log data (part of CRA-257)
    //

    if(background[0][0]->GetN() == 0){
      std::vector<std::vector<NRLib::Grid<float> *> > temp_null;
      temp_null.resize(n_wells_);
      for(int j = 0; j < n_wells_; j++)
        temp_null[j].resize(3);
      EstimateLnData(log_data_vp, temp_null, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Vp", err_txt);
      EstimateLnData(log_data_vs, temp_null, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Vs", err_txt);
      EstimateLnData(log_data_rho, temp_null, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Rho", err_txt);
    }
    else{
      EstimateLnData(log_data_vp, background, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Vp", err_txt);
      EstimateLnData(log_data_vs, background, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Vs", err_txt);
      EstimateLnData(log_data_rho, background, well_names_, mapped_blocked_logs_for_correlation, interval_simboxes, "Rho", err_txt);
    }

    //EstimatePointVar0(point_var_0_, log_data_vp, log_data_vs, log_data_rho, err_txt);

    if (err_txt != "")
      return;

    CalculateNumberOfLags(n_lags, dz_min, dz_rel, interval_simboxes);

    //corr_T_.resize(n_lags, 0);
    auto_cov_.resize(n_lags);

    EstimateAutoCovarianceFunction(auto_cov_, well_names, mapped_blocked_logs_for_correlation, interval_simboxes, log_data_vp, log_data_vs, log_data_rho,
      all_Vs_logs_synthetic, all_Vs_logs_non_synthetic, regression_coef, var_vs_resid, static_cast<float>(dz_min), n_lags, min_blocks_with_data_for_corr_estim_, max_lag_with_data_, err_txt);

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

//
// Do the covariance and correlation estimation
//
/*
void
Analyzelog::estimate(const ModelSettings * modelSettings,
                     Background          * background,
                     std::string         & errTxt)
{
  float ** lnDataAlpha = new float*[nwells_];
  float ** lnDataBeta  = new float*[nwells_];
  float ** lnDataRho   = new float*[nwells_];

  // We do not want syntethic Vs logs to part of the estimation of Var(Vs),
  // Cov(Vp,Vs) and Cov(Vs,Rho) if non-synthetic logs exists. If all Vs logs
  // are synthetic we use them for Var(Vs).
  bool allVsLogsAreSynthetic = true;
  for(int i=0 ; i<nwells_ ; i++)
  {
    allVsLogsAreSynthetic = allVsLogsAreSynthetic && wells_[i]->hasSyntheticVsLog();
  }
  if (allVsLogsAreSynthetic)
  {
    LogKit::LogFormatted(LogKit::Low,"\nThere are no nonsynthetic Vs logs available. Corr(Vp,Vs) and Corr(Vs,Rho) are set 0.7.\n");
  }

  if(background == NULL){
    estimateLnData(lnDataAlpha, NULL, 0, errTxt);
    estimateLnData(lnDataBeta,  NULL, 1, errTxt);
    estimateLnData(lnDataRho,   NULL, 2, errTxt);
  }
  else{
    estimateLnData(lnDataAlpha, background->getAlpha(), 0, errTxt);
    estimateLnData(lnDataBeta, background->getBeta(), 1, errTxt);
    estimateLnData(lnDataRho, background->getRho(), 2, errTxt);
  }
  if (errTxt != "")
    return;

  estimatePointVar0(pointVar0_, lnDataAlpha, lnDataBeta, lnDataRho, errTxt);

  int maxnd;
  calculateNumberOfLags(numberOfLags_, maxnd, errTxt);

  CorrT_ = new float[numberOfLags_+1];
  float dt = static_cast<float>(simbox_->getdz());

  estimateCorrTAndVar0(CorrT_, Var0_,
                       lnDataAlpha, lnDataBeta, lnDataRho,
                       allVsLogsAreSynthetic, dt,
                       numberOfLags_, maxnd,
                       errTxt);

  checkVariances(modelSettings, pointVar0_, Var0_, dt, errTxt);

  for(int i=0 ; i<nwells_ ; i++)
  {
    delete [] lnDataAlpha[i];
    delete [] lnDataBeta[i];
    delete [] lnDataRho[i];
  }
  delete [] lnDataAlpha;
  delete [] lnDataBeta;
  delete [] lnDataRho;
}
*/

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

/*
void
Analyzelog::CalculateNumberOfLags(int                                                   & n_lags,
                                  const std::map<std::string, BlockedLogsCommon *>      & mapped_blocked_logs_for_correlation,
                                  int                                                   & max_nd,
                                  const std::vector<Simbox *>                           & simboxes,
                                  std::string                                           & err_txt)
{
  //
  // Find n (number of lags)
  //
  // Note: The old approach that didn't involve simbox->getTop and
  // simbox->getRelThick produced out-of-bounds lags in estimateCorrTAndVar0.
  //
  // Maybe we can calculate 'n' in a more slick way using min(topsurface)
  // and max(bottomsurface)? This will give a faster and nicer code, but
  // also lags that are larger than actually needed.
  //
  std::map<std::string, BlockedLogsCommon *>::const_iterator it;
  double dt;//      = simbox->getdz();
  double maxdist = 0.0;
  max_nd = 0;
  for(it = mapped_blocked_logs_for_correlation.begin() ; it != mapped_blocked_logs_for_correlation.end() ; it++)
  {
    int   nd                        = it->second->GetNumberOfBlocks();
    const std::vector<int>    s_pos = it->second->GetSposVector();
    const std::vector<double> x_pos = it->second->GetXposBlocked();
    const std::vector<double> y_pos = it->second->GetYposBlocked();
    const std::vector<double> z_pos = it->second->GetZposBlocked();

    std::vector<float> z(nd);
    if(nd>max_nd)
      max_nd=nd;
    for(int j=0 ; j<nd ; j++)
    {
      float xj = static_cast<float>(x_pos[j]);
      float yj = static_cast<float>(y_pos[j]);
      int   s  = s_pos[j];
      z[j] = static_cast<float>((z_pos[j]-simboxes[s_pos[j]]->GetTopErodedSurface(xj,yj))/simboxes[s_pos[s]]->getRelThick(xj,yj));
      for (size_t t = s - 1; s >= 0; s--){
        z[j] += simboxes[t]->getbas
      }
    }
    std::string well_err_txt = "";
    for(int j=0 ; j<nd ; j++)
    {
      for(int k=j+1 ; k<nd ; k++)
      {
        double dist = z[k] - z[j];
        if(dist > maxdist)
          maxdist = dist;
        if(dist < 0) {             // Small negative lags were observed on Smorbukk Sor
          if (floor(dist/dt+0.5) < 0 && wellErrTxt == "") // Check that error is numerically significant
          {
            wellErrTxt += std::string("Negative lags in well \'") + wells[i].GetWellName() +  + "\'.z[k]=";
            wellErrTxt += NRLib::ToString(z[k],3) + std::string(" z[j]=") + NRLib::ToString(z[j],3) + "\n";
          }
        }
      }
    }
    err_txt += wellErrTxt;
    delete [] z;
  }
  n_lags = int(maxdist/dt)+1;
  }
}
*/

void Analyzelog::EstimateLnData(std::map<std::string, std::vector<double> >           & ln_data_map, //map from well to log
                                const std::vector<std::vector<NRLib::Grid<float> *> > & background,
                                const std::vector<std::string>                        & well_names,
                                const std::map<std::string, BlockedLogsCommon *>      & mapped_blocked_logs_for_correlation,
                                const std::vector<Simbox *>                           & interval_simboxes,
                                const std::string                                     & log_name,
                                std::string                                           & err_txt){

  double global_mean   = 0.0;
  int   count         = 0;
  int   nd            = 0;
  int   nd_tot        = 0;
  int   log_nr        = -1;
  //
  // The three vectors below contain data from all wells and all intervals
  //
  std::vector<double>  mean_vector;
  std::vector<double>  ln_data_vector;
  std::vector<double>  low_freq_log;

  assert(log_name == "Vp" || log_name == "Vs" || log_name == "Rho");

  for (int i = 0; i < n_wells_; i++){
    std::vector<double>              temp_log(0);
    ln_data_map.insert(std::pair<std::string, std::vector<double> >(well_names[i], temp_log));

    std::vector<double> temp_low_freq_log =  mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetContLogHighCutBackground(log_name);


    for (size_t j = 0; j < interval_simboxes.size(); j++){
      std::vector<double>   mean;
      std::string           interval_name = interval_simboxes[j]->GetIntervalName();
      nd                    = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNumberOfBlocks();

      const std::vector<int>    & i_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetIposVector();
      const std::vector<int>    & j_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetJposVector();
      const std::vector<int>    & k_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetKposVector();

      std::vector<double> blocked_well_log;

      if(log_name      == "Vp"){
        blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVpBlocked();
        log_nr               = 0;
      }
      else if(log_name == "Vs"){
        blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVsBlocked();
        log_nr               = 1;
      }
      else if(log_name == "Rho"){
        blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetRhoBlocked();
        log_nr               = 2;
      }

      if (background[j][log_nr]->GetN() > 0){
          for(int n = 0; n < nd; n++){
            mean_vector.push_back( (*background[j][log_nr])(i_pos[n], j_pos[n], k_pos[n]));
            mean.push_back( (*background[j][log_nr])(i_pos[n], j_pos[n], k_pos[n]));
          }
      }
      else{
        for(int n = 0; n < nd; n++){
          mean_vector.push_back(0);
          mean.push_back( (*background[j][log_nr])(i_pos[n], j_pos[n], k_pos[n]));
        }
      }


      for (int j = 0; j < nd; j++){
        low_freq_log.push_back(temp_low_freq_log[j]);
        if(blocked_well_log[j] != RMISSING && mean_vector[j] != RMISSING){
          ln_data_map.find(well_names[i])->second.push_back(blocked_well_log[j] - mean[j]);
          global_mean   += blocked_well_log[j] - mean[j];
          ln_data_vector.push_back(blocked_well_log[j] - mean[j]);
          count++;
        }
        else{
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
  //

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

/*
void
Analyzelog::estimateLnData(float      **& lnData,
                           FFTGrid      * background,
                           int            logNr,
                           std::string  & errTxt)
{
  float globalMean = 0.0f;
  int tell = 0;
  if(background !=NULL)
    background->setAccessMode(FFTGrid::RANDOMACCESS);

  for(int i=0 ; i<n_wells_ ; i++)
  {
    int nd;
    const double * xpos  = wells_[i]->getXpos(nd);
    const double * ypos  = wells_[i]->getYpos(nd);
    const double * zpos  = wells_[i]->getZpos(nd);

    const float * wLog = NULL;
    if (logNr == 0)
      wLog = wells_[i]->getAlpha(nd);
    else if (logNr == 1)
      wLog = wells_[i]->getBeta(nd);
    else if (logNr == 2)
      wLog = wells_[i]->getRho(nd);
    else
      errTxt += std::string("In Analyzelog::estimateLnData: Log number ") + NRLib::ToString(logNr)
              + std::string(" does not exist (Vp=1,Vs=2,Rho=3)\n");

    float * mean = new float[nd];
    int useBackgroundMean = 0;
    if(background!=NULL){
      readMeanData(background, nd, xpos, ypos, zpos, mean);
      useBackgroundMean = 1;
    }

    lnData[i] = new float[nd];
    for(int j=0 ; j<nd ; j++)
    {
      if(wLog[j]!= RMISSING && mean[j]!=RMISSING)
      {
        lnData[i][j] = static_cast<float>(log(wLog[j])-useBackgroundMean*mean[j]); //mean is ln(background)
        globalMean += lnData[i][j];
        tell++;
      }
      else
        lnData[i][j] = RMISSING;
    }
    delete [] mean;
  }
  if(background!=NULL)
    background->endAccess();

  //
  // Subtract global mean from data.
  //
  if(tell > 0)
  {
    globalMean /= tell;
    for(int i=0 ; i<n_wells_ ; i++)
    {
      for(int j=0 ; j<wells_[i]->getNd() ; j++)
      {
        if (lnData[i][j] != RMISSING)
        {
          lnData[i][j] -= globalMean;
        }
      }
    }
  }
  else
  {
    errTxt += std::string("Could not estimate globalMean for log") + NRLib::ToString(logNr) + " (Vp=0,Vs=1,Rho=2)\n";
  }
}
*/

//
// Read background model in well positions.
//
/*
void Analyzelog::readMeanData(FFTGrid *cube, int nd,
                              const double * xpos,
                              const double * ypos,
                              const double * zpos,
                              float * meanValue)
{
  int i, j, k, n;
  for(n=0;n<nd;n++)
  {
    simbox_->getIndexes(xpos[n], ypos[n], zpos[n], i,j,k);
    meanValue[n] = cube->getRealValue(i,j,k);     //changed
  }
}
*/

//
// Estimate covariance matrix for alpha, beta and rho.
//
/*
void
Analyzelog::estimatePointVar0(float      ** Var0,
                              float      ** lnDataAlpha,
                              float      ** lnDataBeta,
                              float      ** lnDataRho,
                              std::string & errTxt)
{
  int nd, i, j, tell1, tell2, tell3, tell4, tell5, tell6;
  double sum1, sum2, sum3, sum4, sum5, sum6;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Var0[i][j] = 0.0;

  tell1 = 0;
  tell2 = 0;
  tell3 = 0;
  tell4 = 0;
  tell5 = 0;
  tell6 = 0;
  sum1 = 0.0;
  sum2 = 0.0;
  sum3 = 0.0;
  sum4 = 0.0;
  sum5 = 0.0;
  sum6 = 0.0;

  for(i=0;i<nwells_;i++)
  {
    nd = wells_[i]->getNd();
    for(j=0;j<nd;j++)
    {
      if(lnDataAlpha[i][j]!=RMISSING)
      {
        sum1 += lnDataAlpha[i][j]*lnDataAlpha[i][j];
        tell1++;
      }
    }

    for(j=0;j<nd;j++)
    {
      if(lnDataBeta[i][j]!=RMISSING)
      {
        sum2 += lnDataBeta[i][j]*lnDataBeta[i][j];
        tell2++;
      }
    }

    for(j=0;j<nd;j++)
    {
      if(lnDataRho[i][j]!=RMISSING)
      {
        sum3 += lnDataRho[i][j]*lnDataRho[i][j];
        tell3++;
      }
    }

    for(j=0;j<nd;j++)
    {
      if(lnDataAlpha[i][j]!=RMISSING && lnDataBeta[i][j]!=RMISSING)
      {
        sum4 += lnDataAlpha[i][j]*lnDataBeta[i][j];
        tell4++;
      }
    }

    for(j=0;j<nd;j++)
    {
      if(lnDataAlpha[i][j]!=RMISSING && lnDataRho[i][j]!=RMISSING)
      {
        sum5 += lnDataAlpha[i][j]*lnDataRho[i][j];
        tell5++;
      }
    }

    for(j=0;j<nd;j++)
    {
      if(lnDataRho[i][j]!=RMISSING && lnDataBeta[i][j]!=RMISSING)
      {
        sum6 = lnDataBeta[i][j]*lnDataRho[i][j];
        tell6++;
      }
    }
  }

  if(tell1 < 2)
  {
    errTxt += "\nNot enough well data within simulation area to estimate variance of Vp.\n";
  }
  if(tell3 < 2)
  {
    errTxt += "\nNot enough well data within simulation area to estimate variance of Rho.\n";
  }

  if(errTxt == "")
  {
    Var0[0][0] = float (sum1/(tell1-1));
    if(tell2>1)
      Var0[1][1] = float (sum2/(tell2-1));
    else
    {
      Var0[1][1] = 2*Var0[0][0];
      LogKit::LogFormatted(LogKit::Low,"\nEstimating Vs variance as 2 * Vp variance.\n");
    }
    Var0[2][2] = float (sum3/(tell3-1));
    if(tell4>1)
      Var0[0][1] = Var0[1][0] = float (sum4/(tell4-1));
    if(tell5>1)
      Var0[0][2] = Var0[2][0] = float (sum5/(tell5-1));
    if(tell6>1)
      Var0[1][2] = Var0[2][1] = float (sum6/(tell6-1));
  }
}
*/

//
// CRA-257: new implementation of estimation of autocovariance function
//
void            Analyzelog::EstimateAutoCovarianceFunction(std::vector<NRLib::Matrix >                        & auto_cov,
                                                           const std::vector<std::string>                     & well_names,
                                                           const std::map<std::string, BlockedLogsCommon *>   & mapped_blocked_logs_for_correlation,
                                                           const std::vector<Simbox *>                        & interval_simboxes,
                                                           const std::map<std::string, std::vector<double> >  & log_data_vp,
                                                           const std::map<std::string, std::vector<double> >  & log_data_vs,
                                                           const std::map<std::string, std::vector<double> >  & log_data_rho,
                                                           bool                                                 all_Vs_logs_synthetic,
                                                           bool                                                 all_Vs_logs_non_synthetic,
                                                           NRLib::Vector                                      & regression_coef,
                                                           double                                             & var_vs_resid,
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
  // use the wells with non-synthetic vs to calculate the linear regression of Vs on Vp and Rho
  //
  if(!all_Vs_logs_synthetic && !all_Vs_logs_non_synthetic){
    int                                 n_data = 0;
    NRLib::Matrix                       A;
    NRLib::Vector                       b;
    std::vector<double>                 vp;
    std::vector<double>                 rho;
    std::vector<double>                 vs;
    for (size_t i = 0; i < well_names.size(); i++){
      if(mapped_blocked_logs_for_correlation.find(well_names[i])->second->HasSyntheticVsLog() == false){
        //int n_lags;
        for (size_t j = 0; j<interval_simboxes.size(); j++){
          const std::string interval_name         = interval_simboxes[j]->GetIntervalName();
          const std::vector<double> vp_blocked    = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVpBlocked();
          const std::vector<double> rho_blocked   = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetRhoBlocked();
          const std::vector<double> vs_blocked    = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVsBlocked();
          const std::vector<double> z_pos         = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetZposBlocked();
          // Loop over all blocks
          for (int k = 0; k < mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNBlocksWithData(interval_name); k++){
            if(vp_blocked[k] != RMISSING && rho_blocked[k] != RMISSING && vs_blocked[k] != RMISSING){
              //n_lags = static_cast<int>(std::floor(std::abs(z_pos[l] - z_pos[k])/min_dz + 0.5));
              n_data++;
              vp.push_back(vp_blocked[k]);
              rho.push_back(rho_blocked[k]);
              vs.push_back(vs_blocked[k]);
            }
          }
        }
      }
    }
    A.resize(n_data, 2);
    b.resize(n_data);
    for(size_t j = 0; j < vp.size(); j++){
      A(j, 0) = vp[j];
      A(j, 1) = rho[j];
      b(j)    = vs[j];
    }
    if (n_data > 0){
      regression_coef = Regress(A, b);

      //
      // calculate residuals for vs and the resulting residual variance for the non-synthetic vs wells
      //
      double sum_resid_square = 0.0;
      double residual;
      for (size_t i = 0; i < well_names.size(); i++){
        if(mapped_blocked_logs_for_correlation.find(well_names[i])->second->HasSyntheticVsLog() == false){
          //int n_lags;
          for (size_t j = 0; j < interval_simboxes.size(); j++){
            const std::string interval_name         = interval_simboxes[j]->GetIntervalName();
            const std::vector<double> vp_blocked    = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVpBlocked();
            const std::vector<double> rho_blocked   = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetRhoBlocked();
            const std::vector<double> vs_blocked    = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVsBlocked();
            const std::vector<double> z_pos         = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetZposBlocked();
            // Loop over all blocks
            for (int k = 0; k < mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNBlocksWithData(interval_name); k++){
              if(vp_blocked[k] != RMISSING && rho_blocked[k] != RMISSING && vs_blocked[k] != RMISSING){
                residual = (regression_coef(0)*vp_blocked[k] + regression_coef(1)*rho_blocked[k] - vs[k]);
                sum_resid_square += residual*residual;
              }
            }
          }
        }
      }
      var_vs_resid = sum_resid_square / (n_data - 2);
    }
  }

  //
  // 2.2: Estimate autocovariance using both synthetic and non-synthetic well logs
  //
  for (size_t i = 0; i < well_names.size(); i++){
    time_t timestart, timeend;
    time(&timestart);
    for (size_t j = 0; j < interval_simboxes.size(); j++){
      std::string interval_name = interval_simboxes[j]->GetIntervalName();
      size_t nd = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNBlocksWithData(interval_name);
      const std::vector<double> x_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetXposBlocked();
      const std::vector<double> y_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetYposBlocked();
      const std::vector<double> z_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetZposBlocked();
      std::vector<double> z(nd);
      const std::vector<double> log_vp   = log_data_vp.find(well_names[i])->second;
      const std::vector<double> log_vs   = log_data_vs.find(well_names[i])->second;
      const std::vector<double> log_rho  = log_data_rho.find(well_names[i])->second;
      //
      // Calculate autocovariance
      //
      int lag = 0;
      for (size_t k = 0; k < nd; k++){
        for (size_t l = k; l < nd; l++){

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
          }
          // cov(rho_k, vp_l)
          if(log_rho[k] != RMISSING && log_vp[l] != RMISSING){
            temp_auto_cov[lag](0,2) += log_rho[k]*log_vp[l];
            count[lag](0,2) += 1;
          }
          //
          // If this Vs log is synthetic and there exist real Vs logs: use regression coefficients
          //
          if(!all_Vs_logs_synthetic && mapped_blocked_logs_for_correlation.find(well_names[i])->second->HasSyntheticVsLog() == true){
            // Use the relation Vs = a*Vp + b*Rho + e, where e is iid
            // cov[t](vs_i, vs_j) = cov[t](a*vp_k + b*rho_k + e_k, a*vp_l + b*rho_l + e_l) = a*a*cov(vp_k,vp_l) + a*b*cov(vp_k, rho_l) + a*b*cov(vp_l, rho_k) + b*b*cov(rho_k, rho_l) + I(k = l) var(e)
            if(log_vp[k] != RMISSING && log_rho[l] != RMISSING && log_vp[l] != RMISSING && log_rho[k] != RMISSING){
              temp_auto_cov[lag](1,1) += regression_coef(0)*regression_coef(0)*log_vp[k]*log_vp[k] + regression_coef(1)*regression_coef(1)*log_rho[l]*log_rho[l]
                                    + regression_coef(0)*regression_coef(1)*(log_vp[k]*log_rho[l] + log_vp[l]*log_rho[k]);
              if (k == l)
                temp_auto_cov[lag](1,1) += var_vs_resid;
              count[lag](1,1) += 1;
            }
            // cov[t](vp, vs) = cov[t](vp, a*vp + b*rho + e) = a*autocov[t](vp)
            if (log_vp[k] != RMISSING && log_vp[l] != RMISSING){
              temp_auto_cov[lag](0,1) += regression_coef(0)*log_vp[k]*log_vp[l];
              count[lag](0,1) += 1;
            }
            // cov[t](rho, vs) = cov[t](rho, a*vp + b*rho + e) = b*autocov[t](rho)
            if (log_rho[k] != RMISSING && log_rho[l] != RMISSING){
              temp_auto_cov[lag](1,2) += regression_coef(1)*log_rho[k]*log_rho[l];
              count[lag](1,2) += 1;
            }
          }
          //
          // Non-synthetic Vs log
          //
          else if(mapped_blocked_logs_for_correlation.find(well_names[i])->second->HasSyntheticVsLog() == false){
            // cov[t](vs, vs)
            if(log_vs[k] != RMISSING && log_vs[l] != RMISSING){
              if (lag > max_lag_with_data)
                max_lag_with_data = lag;
              temp_auto_cov[lag](1,1) += log_vs[k]*log_vs[l];
              count[lag](1,1) += 1;
            }
            // cov[t](vp, vs)
            if(log_vp[k] != RMISSING && log_vs[l] != RMISSING){
              temp_auto_cov[lag](0,1) += log_vp[k]*log_vs[l];
              count[lag](0,1) += 1;
            }
            // cov[t](vs, vp)
            if(log_vs[k] != RMISSING && log_vp[l] != RMISSING){
              temp_auto_cov[lag](0,1) += log_vs[k]*log_vp[l];
              count[lag](0,1) += 1;
            }
            // cov[t](rho, vs)
            if(log_vs[k] != RMISSING && log_rho[l] != RMISSING){
              temp_auto_cov[lag](1,2) += log_vs[k]*log_rho[l];
              count[lag](1,2) += 1;
            }
            // cov[t](vs, rho)
            if(log_rho[k] != RMISSING && log_vs[l] != RMISSING){
              temp_auto_cov[lag](1,2) += log_rho[k]*log_vs[l];
              count[lag](1,2) += 1;
            }
          }
        }
      }

    }
    time(&timeend);
    long int time = static_cast<long int>(timeend - timestart);
    printf("\nWell %s processed in %ld seconds.",well_names[i].c_str(),time);
  }

  // Calculate autocovariances
  //
  if(all_Vs_logs_synthetic){
    LogKit::LogFormatted(LogKit::Low,"\nEstimating Vs autocovariance as 2 * Vp autocovariance and corr(Vp, Vs) and corr(Vs, Rho) are set to 0.7.\n");
  }
  for (int i = 0; i < max_lag_with_data; i++){
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
        err_text += "Not enough well data within simulation area to estimate autocovariance of Vp for time lag = "+ CommonData::ConvertIntToString(i) +", dz = " + CommonData::ConvertFloatToString(min_dz) +" +.\n";
      }
      if(count[i](2,2)<2)
      {
        err_text += "Not enough well data within simulation area to estimate autocovariance of Rho for time lag "+ CommonData::ConvertIntToString(i) +", dz = "+ CommonData::ConvertFloatToString(min_dz) + ".\n";
      }
      if (err_text == ""){
        temp_auto_cov[i](0,0) = temp_auto_cov[i](0,0) / (count[i](0,0) - 1);
        temp_auto_cov[i](2,2) = temp_auto_cov[i](2,2) / (count[i](2,2) - 1);
        if(all_Vs_logs_synthetic){
          temp_auto_cov[i](1,1) = 2*temp_auto_cov[i](0,0);
        } else if (count[i](1,1) > 1){
          temp_auto_cov[i](1,1) = temp_auto_cov[i](1,1) / (count[i](1,1) - 1);
        } else{
          err_text += "There is not enough data in the Vs logs.\n";
        }
      }
    }
    // If there is not enough data to estimate autocovariance for lags > 50, set to 0 and return a warning
    else{
      if(count[i](0,0)<2){
        temp_auto_cov[i](0,0) = 0.0;
        LogKit::LogFormatted(LogKit::Low,"\nWARNING: Vp autocovariance for time lag %i, dz = %d could not be estimated and is set to 0.\n", i, min_dz);
      }
      else{
        temp_auto_cov[i](0,0) = temp_auto_cov[i](0,0) / (count[i](0,0) - 1);
      }
      if(count[i](2,2)<2){
        temp_auto_cov[i](2,2) = 0.0;
        LogKit::LogFormatted(LogKit::Low,"\nWARNING: Rho autocovariance for time lag %i, dz = %d could not be estimated and is set to 0.\n", i, min_dz);
      }
      else{
        temp_auto_cov[i](2,2) = temp_auto_cov[i](2,2) / (count[i](2,2) - 1);
      }
      if(all_Vs_logs_synthetic){
        temp_auto_cov[i](1,1) = 2*temp_auto_cov[i](0,0);
      } else if (count[i](1,1) > 1){
        temp_auto_cov[i](1,1) = temp_auto_cov[i](1,1) / (count[i](1,1) - 1);
      } else{
        err_text += "There is not enough data in the Vs logs.\n";
      }
    }
  }

  if(err_text == ""){

    // Default covariance between Vp and Vs. Value is explained in Jira issue CRA-220.
    for (int i = 0; i < max_lag_with_data; i++){

      if(count[i](0,1) > 1){
        temp_auto_cov[i](0,1) = temp_auto_cov[i](0,1)/(count[i](0,1) - 1);
      }
      else{
        temp_auto_cov[i](0,1)  = std::sqrt(temp_auto_cov[i](0,0))*std::sqrt(temp_auto_cov[i](1,1))*0.70;
      }
      if(count[i](0,2) > 1){
        temp_auto_cov[i](0,2) = temp_auto_cov[i](0,2)/(count[i](0,2) - 1);
      }
      else{
        temp_auto_cov[i](0,2)  = 0.0;
      }
      if(count[i](1,2) > 1){
        temp_auto_cov[i](1,2) = temp_auto_cov[i](1,2)/(count[i](1,2) - 1);
      }
      else{
        temp_auto_cov[i](1,2)  = 0.0;
      }
      temp_auto_cov[i](1,0) = temp_auto_cov[i](0,1);
      temp_auto_cov[i](2,0) = temp_auto_cov[i](2,0);
      temp_auto_cov[i](2,1) = temp_auto_cov[i](1,2);
    }

  }


  // CRA-257, Erik N: Should this still be done?
  /************************************************************
  //
  // 3) Scale corTT correlation to account for blocking effects
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
    for(int k=j; k<3; k++){
      auto_cov[0](j,k) = temp_auto_cov[0](j,k);
      auto_cov[0](k,j) = auto_cov[0](j,k);
    }
  }

  //CorrT[0]=corTT[0];
  for(int i=1; i<max_lag_with_data; i++){
    for(int j=0 ; j<3; j++){
      for(int k=j; k<3; k++){
        if(count[i-1](j,k)>0)                 // Find previous estimated autocov with data
          i_prev[j][k]=i-1;
        if(i_next[j][k]<i+1){
          i_next[j][k]++;
        while(count[i_next[j][k]](j,k) == 0 && i_next[j][k]<max_lag_with_data)   // Find next estimated autocov with data
          i_next[j][k]++;
      }
      c_prev[j][k]    = temp_auto_cov[i_prev[j][k]](j,k);
      c_next[j][k]    = temp_auto_cov[i_next[j][k]](j,k);
      nipol    = std::min(count[i_prev[j][k]](j,k),count[i_next[j][k]](j,k));
      cipol    = ((i_next[j][k]-i)*c_prev[j][k] + (i-i_prev[j][k])*c_next[j][k])/(i_next[j][k]-i_prev[j][k]);
      auto_cov[i](j,k) = (count[i](j,k)*temp_auto_cov[i](j,k) + nipol*cipol)/(count[i](j,k)+nipol);
      auto_cov[i](k,j) = auto_cov[i](j,k);

      //LogKit::LogFormatted(LogKit::Low," i nTT,corTT  iprev,cprev  inext,cnext  nipol,cipol  %4d %4d%8.3f %4d%8.3f %4d%8.3f    %4d%8.3f\n",
      //                   i,nTT[i],corTT[i],iprev,cprev,inext,cnext,nipol,cipol);
      }
    }

  }

  // Last time step with data
  for(int j=0 ; j<3; j++){
    for(int k=j; k<3; k++){
      auto_cov[max_lag_with_data-1](j,k) = (count[max_lag_with_data-1](j,k)*temp_auto_cov[max_lag_with_data-1](j,k) + i_prev[j][k]*c_prev[j][k]*0.5)/(count[max_lag_with_data-1](j,k) + i_prev[j][k]);
      auto_cov[max_lag_with_data-1](k,j) = auto_cov[max_lag_with_data-1](j,k);
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

  // CRA-257 Erik N: Should this still be done?
  /*****************************************************************
    long int ntot=0;
    for(i=1;i<nend+1;i++)
      ntot+=nTT[i];
    float b = static_cast<float>((2*ntot)/static_cast<float>(nend+1));  // counts the zero-element as well
    float a = b/static_cast<float>(nend+1);
    for(i=1;i<nend+1;i++)
    {
      CorrT[i] *= (b - a*i)/b;
    }
  ********************************************************************/

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
    for(int i=0; i<max_nd+1; i++){
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
  LogKit::LogFormatted(LogKit::Low,"\nEstimated parameter variance and parameter temporal correlation in %d seconds.\n",
                   static_cast<int>(timeend_tot-timestart_tot));

}

//
// Estimate auto correlation and variances. Use lag 0 to estimate variance.
//
/*
void
Analyzelog::estimateCorrTAndVar0(float       * CorrT,
                                 float      ** Var0,
                                 float      ** lnDataAlpha,
                                 float      ** lnDataBeta,
                                 float      ** lnDataRho,
                                 bool          allVsLogsAreSynthetic,
                                 float         dt,
                                 int           n,
                                 int           maxnd,
                                 std::string & errTxt)
{
  time_t timestart, timeend;
  time(&timestart);

  std::string tmpErrTxt = "";

  float * corTT    = new float[n+1];
  float * varAj    = new float[n+1];
  float * varBj    = new float[n+1];
  float * varRj    = new float[n+1];
  float * varAk    = new float[n+1];
  float * varBk    = new float[n+1];
  float * varRk    = new float[n+1];
  float * covAA    = new float[n+1];
  float * covBB    = new float[n+1];
  float * covRR    = new float[n+1];
  int   * nAA      = new int[n+1];
  int   * nBB      = new int[n+1];
  int   * nRR      = new int[n+1];
  int   * nTT      = new int[n+1];
  int   * indA     = new int[maxnd];   // arrays pointing to non-missing log entries
  int   * indB     = new int[maxnd];
  int   * indR     = new int[maxnd];
  float * z        = new float[maxnd];

  int i,j,k;
  for(i=0;i<n+1;i++)
  {
    varAj[i] = 0.0;
    varBj[i] = 0.0;
    varRj[i] = 0.0;
    varAk[i] = 0.0;
    varBk[i] = 0.0;
    varRk[i] = 0.0;
    covAA[i] = 0.0;
    covBB[i] = 0.0;
    covRR[i] = 0.0;
    corTT[i] = 0.0;
    nAA[i] = 0;
    nBB[i] = 0;
    nRR[i] = 0;
  }
  float covAB = 0.0;
  float covAR = 0.0;
  float covBR = 0.0;
  int   nAB = 0;
  int   nAR = 0;
  int   nBR = 0;

  for(i=0;i<nwells_;i++)
  {
    time_t timestart, timeend;
    time(&timestart);
    //
    // Extract z-coordinates of each log entry
    //
    int nd;
    const double * xpos= wells_[i]->getXpos(nd);
    const double * ypos= wells_[i]->getYpos(nd);
    const double * zpos= wells_[i]->getZpos(nd);
    for(j=0;j<nd;j++)
    {
      float xj = static_cast<float>(xpos[j]);
      float yj = static_cast<float>(ypos[j]);
      z[j] = static_cast<float>((zpos[j]-simbox_->getTop(xj,yj))/simbox_->getRelThick(xj,yj));
      z[j] = z[j]/dt; // simplify lag calculation  floor((z2-z1)/dt+0.5)
    }
    //
    // 1) Create array of indices pointing to nonmissing data and
    //
    // 2) dLnAlpha = (lnAlpha - lnBGAlpha) - meanResLnAlpha   (NB! lnX-lnBGX is done in constructor)
    //    dLnBeta  = (lnBeta  - lnBGBeta)  - meanResLnBeta
    //    dLnRho   = (lnRho   - lnBGRho)   - meanResLnRho
    //
    //    where the residuals (denoted globalMean in code below) are given by
    //
    //    meanResLnX = sum(lnX - lnBGX)/n)
    //
    float * dLnAlpha = lnDataAlpha[i];
    float * dLnBeta  = lnDataBeta[i];
    float * dLnRho   = lnDataRho[i];
    int na=0;
    int nb=0;
    int nr=0;
    for(j=0;j<nd;j++)
    {
      if (dLnAlpha[j] != RMISSING)
      {
        indA[na] = j;
        na++;
      }
      if (dLnBeta[j] != RMISSING)
      {
        indB[nb] = j;
        nb++;
      }
      if (dLnRho[j] != RMISSING)
      {
        indR[nr] = j;
        nr++;
      }
    }
    //
    // Calculate Cov_t(A,A), Cov_t(B,B), and Cov_t(R,R)
    // ------------------------------------------------
    //
    int h;
    for(j=0;j<na;j++)
    {
      for(k=j;k<na;k++) // indA points to nonmissing entries of dLnAlpha[]
      {
        h = static_cast<int>(floor(z[indA[k]]-z[indA[j]]+0.5));
        if(h >= 0) { //Guard against small upwards movements in transformed domain.
          covAA[h] += dLnAlpha[indA[j]]*dLnAlpha[indA[k]];
          varAj[h] += dLnAlpha[indA[j]]*dLnAlpha[indA[j]];
          varAk[h] += dLnAlpha[indA[k]]*dLnAlpha[indA[k]];
          nAA[h]++;
        }
      }
    }
    if (allVsLogsAreSynthetic || !wells_[i]->hasSyntheticVsLog())
    {
      for(j=0;j<nb;j++)
      {
        for(k=j;k<nb;k++) // indB points to nonmissing entries of dLnBeta[]
        {
          h = static_cast<int>(floor(z[indB[k]]-z[indB[j]]+0.5));
          if(h >= 0) { //Guard against small upwards movements in transformed domain.
            covBB[h] += dLnBeta[indB[j]]*dLnBeta[indB[k]];
            varBj[h] += dLnBeta[indB[j]]*dLnBeta[indB[j]];
            varBk[h] += dLnBeta[indB[k]]*dLnBeta[indB[k]];
            nBB[h]++;
          }
        }
      }
    }
    for(j=0;j<nr;j++)
    {
      for(k=j;k<nr;k++) // indR points to nonmissing entries of dLnRho[]
      {
        h = static_cast<int>(floor(z[indR[k]]-z[indR[j]]+0.5));
        if(h >= 0) { //Guard against small upwards movements in transformed domain.
          covRR[h] += dLnRho[indR[j]]*dLnRho[indR[k]];
          varRj[h] += dLnRho[indR[j]]*dLnRho[indR[j]];
          varRk[h] += dLnRho[indR[k]]*dLnRho[indR[k]];
          nRR[h]++;
        }
      }
    }
    //
    //  Cov_0(A,B) Cov_0(A,R), and Cov_0(B,R)   (Cov_0(B,A),... follows by sym)
    //  -------------------------------------
    //
    //  We cannot use the reduced loop structure involving indA, indB, and indC
    //  when we mix A, B, or R. This is not a performance problem, however, since
    //  most of the calculation time will go into calculating Cov_t
    //
    for(j=0;j<nd;j++)
    {
      for(k=j;k<nd;k++)
      {
        h = static_cast<int>(floor(z[k]-z[j]+0.5));
        if (h==0)
        {
          if (!wells_[i]->hasSyntheticVsLog())
          {
            if(dLnAlpha[j]!=RMISSING && dLnBeta[k]!=RMISSING) //cov(Aj,Bk)
            {
              covAB += dLnAlpha[j]*dLnBeta[k];
              nAB++;
            }
            if(dLnAlpha[k]!=RMISSING && dLnBeta[j]!=RMISSING) //cov(Ak,Bj)
            {
              covAB += dLnAlpha[k]*dLnBeta[j];
              nAB++;
            }
            if(dLnBeta[j]!=RMISSING && dLnRho[k]!=RMISSING)   //cov(Bj,Rk)
            {
              covBR += dLnBeta[j]*dLnRho[k];
              nBR++;
            }
            if(dLnBeta[k]!=RMISSING && dLnRho[j]!=RMISSING)   //cov(Bk,Rj)
            {
              covBR += dLnBeta[k]*dLnRho[j];
              nBR++;
            }
          }
          if(dLnAlpha[j]!=RMISSING && dLnRho[k]!=RMISSING)  //cov(Aj,Rk)
          {
            covAR += dLnAlpha[j]*dLnRho[k];
            nAR++;
          }
          if(dLnAlpha[k]!=RMISSING && dLnRho[j]!=RMISSING)  //cov(Ak,Rj)
          {
            covAR += dLnAlpha[k]*dLnRho[j];
            nAR++;
          }
        }
      }
    }
    time(&timeend);
    printf("\nWell %s processed in %ld seconds.",
      wells_[i]->getWellname().c_str(),timeend-timestart);
  }
  printf("\n");

  if(nAA[0]<2)
  {
    tmpErrTxt += "Not enough well data within simulation area to estimate variance of Vp.\n";
  }
  if(nRR[0]<2)
  {
    tmpErrTxt += "Not enough well data within simulation area to estimate variance of Rho.\n";
  }

  if(tmpErrTxt == "")
  { //
    // Calculate Var0 (covariances at zero lag)
    // ==============
    //
    Var0[0][0] = covAA[0]/(nAA[0]-1);
    Var0[2][2] = covRR[0]/(nRR[0]-1);
    if(nBB[0]>1)
      Var0[1][1] = covBB[0]/(nBB[0]-1);
    else
    {
      LogKit::LogFormatted(LogKit::Low,"\nEstimating Vs variance as 2 * Vp variance.\n");
      Var0[1][1] = 2*Var0[0][0];
    }
    // Default covariance between Vp and Vs. Value is explained in Jira issue CRA-220.
    Var0[0][1] = Var0[1][0] = static_cast<float>(std::sqrt(Var0[0][0])*std::sqrt(Var0[1][1])*0.70);
    Var0[0][2] = Var0[2][0] = 0.0;
    Var0[1][2] = Var0[2][1] = 0.0;
    if (nAB>1)
      Var0[0][1] = Var0[1][0] = covAB/(nAB-1);
    if (nAR>1)
      Var0[0][2] = Var0[2][0] = covAR/(nAR-1);
    if (nBR>1)
      Var0[1][2] = Var0[2][1] = covBR/(nBR-1);

    //
    // Calculate CorrT
    // ===============
    // The joint correlation CorrT is calculated as a weighted
    // sum of correlations corAA, corBB, and corRR.
    //
    // 1) Individual correlations: covAA -> corAA
    // 2) weighted_sum(corAA,corBB,corRR) -> corTT
    // 3) Scale corTT[i] -> corTT[i]/corTT[0];
    // 4) Linear interpolate corTT -> CorrT
    // 5) Linear downscale CorrT with distance from 0
    //
    for(i=0;i<n+1;i++)
    {
      int naa = 0;       // These must be initialised to zero.
      int nbb = 0;
      int nrr = 0;
      float corAA = 0.0; // These must be initialised to zero.
      float corBB = 0.0;
      float corRR = 0.0;
      //
      // 1) + 2)
      //
      if(nAA[i]>1 && varAj[i]>0.0 && varAk[i]>0.0)
      {
        corAA = covAA[i]/sqrt(varAj[i]*varAk[i]);
        naa = nAA[i];
      }
      if(nBB[i]>1 && varBj[i]>0.0 && varBk[i]>0.0)
      {
        corBB = covBB[i]/sqrt(varBj[i]*varBk[i]);
        nbb = nBB[i];
      }
      if(nRR[i]>1 && varRj[i]>0.0 && varRk[i]>0.0)
      {
        corRR = covRR[i]/sqrt(varRj[i]*varRk[i]);
        nrr = nRR[i];
      }

      if (corAA<-1 || corAA>1)
      {
        tmpErrTxt += std::string("Correlation corAA=") + NRLib::ToString(corAA, 3)
                + std::string(" out of range for element ") + NRLib::ToString(i) + "\n";
      }
      if (corBB<-1 || corBB>1)
      {
        tmpErrTxt += std::string("Correlation corBB=") + NRLib::ToString(corBB, 3)
                + std::string(" out of range for element ") + NRLib::ToString(i) + "\n";
      }
      if (corRR<-1 || corRR>1)
      {
        tmpErrTxt += std::string("Correlation corRR=") + NRLib::ToString(corRR, 3)
                + std::string(" out of range for element ") + NRLib::ToString(i) + "\n";
      }
      if (naa>1 || nbb>1 || nrr>1)
        corTT[i] = (naa*corAA + nbb*corBB + nrr*corRR)/(static_cast<float>(naa+nbb+nrr));
      else
        corTT[i] = WELLMISSING;
      nTT[i] = naa+nbb+nrr;

      covAA[i] = corAA; // quickhack to get correlations printed to file
      covBB[i] = corBB;
      covRR[i] = corRR;
    }
    //
    // 3) Scale corTT correlation to account for blocking effects
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
    //
    // 4) Smooth observations by making a weighted average of each observation
    //    and the linear interpolation of its two nearest neighbours.
    //
    CorrT[0]=corTT[0];
    cprev = 0;
    iprev = 0;
    inext = 1;
    nend  = i;
    for(i=1;i<nend;i++)
    {
      if(nTT[i-1]>0)                         // Find previous correlation
        iprev=i-1;
      if(inext<i+1)
      {
        inext++;
        while(nTT[inext]==0 && inext<nend)   // Find next correlation
          inext++;
      }
      cprev    = corTT[iprev];
      cnext    = corTT[inext];
      nipol    = std::min(nTT[iprev],nTT[inext]);
      cipol    = ((inext-i)*cprev + (i-iprev)*cnext)/(inext-iprev);
      CorrT[i] = (nTT[i]*corTT[i] + nipol*cipol)/(nTT[i]+nipol);

      //LogKit::LogFormatted(LogKit::Low," i nTT,corTT  iprev,cprev  inext,cnext  nipol,cipol  %4d %4d%8.3f %4d%8.3f %4d%8.3f    %4d%8.3f\n",
      //                   i,nTT[i],corTT[i],iprev,cprev,inext,cnext,nipol,cipol);

    }
    CorrT[nend] = (nTT[nend]*corTT[nend] + iprev*cprev*0.5f)/(nTT[nend]+iprev);
    //
    // 5) Downscale correlations linearly
    //
    long int ntot=0;
    for(i=1;i<nend+1;i++)
      ntot+=nTT[i];
    float b = static_cast<float>((2*ntot)/static_cast<float>(nend+1));  // counts the zero-element as well
    float a = b/static_cast<float>(nend+1);
    for(i=1;i<nend+1;i++)
    {
      CorrT[i] *= (b - a*i)/b;
    }

    if(ModelSettings::getDebugLevel() > 0) {
      std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), std::string("Autocorr.dat"));
      std::ofstream file;
      NRLib::OpenWrite(file, fileName);
      file << "   i      nAA    corAA      nBB    corBB      nRR    corRR         nTT    corTT    CorrT\n"
           << "----------------------------------------------------------------------------------------\n";
      for(i=0;i<n+1;i++)
      {
        file << std::fixed
             << std::right
             << std::setprecision(3)
             << std::setw(4) << i        << " "
             << std::setw(8) << nAA[i]   << " "
             << std::setw(8) << covAA[i] << " "
             << std::setw(8) << nBB[i]   << " "
             << std::setw(8) << covBB[i] << " "
             << std::setw(8) << nRR[i]   << " "
             << std::setw(8) << covRR[i] << "    "
             << std::setw(8) << nTT[i]   << " "
             << std::setw(8) << corTT[i] << " "
             << std::setw(8) << CorrT[i]
             << std::endl;
      }
      file.close();
    }
  }

  errTxt += tmpErrTxt;
  //
  // Return replace n by nonzero elements nend?
  //
  // n=nend;
  //
  time(&timeend);
  LogKit::LogFormatted(LogKit::Low,"\nEstimate parameter variance and parameter temporal correlation in %d seconds.\n",
                   static_cast<int>(timeend-timestart));

  delete [] z;
  delete [] indA;
  delete [] indB;
  delete [] indR;
  delete [] varAj;
  delete [] varBj;
  delete [] varRj;
  delete [] varAk;
  delete [] varBk;
  delete [] varRk;
  delete [] covAA;
  delete [] covBB;
  delete [] covRR;
  delete [] corTT;
  delete [] nAA;
  delete [] nBB;
  delete [] nRR;
  delete [] nTT;
}
*/

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
  float minVarAlpha = model_settings->getVarVpMin();
  float maxVarAlpha = model_settings->getVarVpMax();
  float minVarBeta  = model_settings->getVarVsMin();
  float maxVarBeta  = model_settings->getVarVsMax();
  float minVarRho   = model_settings->getVarRhoMin();
  float maxVarRho   = model_settings->getVarRhoMax();

  /*
  if (pointVar0[0][0] < minVarAlpha || pointVar0[0][0] > maxVarAlpha)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vp point variance "  << pointVar0[0][0]
      << " is outside allowed interval Min=" << minVarAlpha << " Max=" << maxVarAlpha << "\n";
    errTxt += o.str();
  }
  if (pointVar0[1][1] < minVarBeta || pointVar0[1][1] > maxVarBeta)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vs point variance "  << pointVar0[1][1]
      << " is outside allowed interval Min=" << minVarBeta << " Max=" << maxVarBeta << "\n";
    errTxt += o.str();
  }
  if (pointVar0[2][2] < minVarRho || pointVar0[2][2] > maxVarRho)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Rho point variance "  << pointVar0[2][2]
      << " is outside allowed interval Min=" << minVarRho << " Max=" << maxVarRho << "\n";
    errTxt += o.str();
  }

  if (errTxt != "")
  {
    LogKit::LogFormatted(LogKit::Low,"\n\n---------------------------------------------------");
    LogKit::LogFormatted(LogKit::Low,"\n                         ln Vp     ln Vs    ln Rho ");
    LogKit::LogFormatted(LogKit::Low,"\nWell log variances:   %.2e  %.2e  %.2e ",pointVar0[0][0],pointVar0[1][1],pointVar0[2][2]);
    LogKit::LogFormatted(LogKit::Low,"\n---------------------------------------------------\n");
  }
  */

  //
  // We scale the minimum variances allowed with 1/dt since the variance decreases with increasing dt..
  //
  if (var_0(0,0) < minVarAlpha/dz || var_0(0,0) > maxVarAlpha)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vp variance "  << var_0(0,0)
      << " is outside allowed interval Min=" << minVarAlpha/dz << " Max=" << maxVarAlpha << "\n";
    err_txt += o.str();
  }
  if (var_0(1,1) < minVarBeta/dz || var_0(1,1)  > maxVarBeta)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vs variance "  << var_0(1,1)
      << " is outside allowed interval Min=" << minVarBeta/dz << " Max=" << maxVarBeta << "\n";
    err_txt += o.str();
  }
  if (var_0(2,2)  < minVarRho/dz || var_0(2,2)  > maxVarRho)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Rho variance "  << var_0(2,2)
      << " is outside allowed interval Min=" << minVarRho/dz << " Max=" << maxVarRho << "\n";
    err_txt += o.str();
  }
  if (err_txt != "")
  {
    LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------------------------------------");
    LogKit::LogFormatted(LogKit::Low,"\n                          ln Vp     ln Vs    ln Rho ");
    //LogKit::LogFormatted(LogKit::Low,"\nWell log  variances:   %.2e  %.2e  %.2e",pointVar0[0][0],pointVar0[1][1],pointVar0[2][2]);
    LogKit::LogFormatted(LogKit::Low,"\nParameter variances:   %.2e  %.2e  %.2e (used by program)",var_0(0,0) ,var_0(1,1) ,var_0(2,2) );
    LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------------------------------------\n");
  }
}

NRLib::Vector Analyzelog::Regress(const NRLib::Matrix                             & A,
                                  const NRLib::Vector                             & b){

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
