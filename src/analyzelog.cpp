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
#include "src/welldata.h"
#include "src/background.h"
#include "src/simbox.h"
#include "src/io.h"

// CRA-257: New correlation estimation routine with blocked logs within one interval
Analyzelog::Analyzelog(const std::vector<NRLib::Well>                          & wells,
                       const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs_for_correlation,
                       const std::vector<NRLib::Grid<float> *>                 & background,
                       //const std::vector<NRLib::Grid<double> >                 & background_max_hz,
                       const Simbox                                            * interval_simbox,
                       const ModelSettings                                     * model_settings,
                       std::string                                             & err_txt):
min_blocks_with_data_for_corr_estim_(model_settings->getMinBlocksForCorrEstimation()){

  enough_data_for_corr_estimation_  = false;
  n_lags_                   = 0;
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

// CRA-257: New correlation estimation routine with blocked logs within multiple intervals
Analyzelog::Analyzelog(const std::vector<NRLib::Well>                          & wells,
                       const std::map<std::string, BlockedLogsCommon *>        & mapped_blocked_logs_for_correlation,
                       const std::vector<std::vector<NRLib::Grid<float> *> >   & background,
                       //const std::vector<std::vector<NRLib::Grid<double> > >   & background_max_hz,
                       const std::vector<Simbox>                               & interval_simboxes,
                       const ModelSettings                                     * model_settings,
                       std::string                                             & err_txt):
min_blocks_with_data_for_corr_estim_(model_settings->getMinBlocksForCorrEstimation()){

  enough_data_for_corr_estimation_  = false;
  n_lags_                           = 0;
  interval_name_                    = "";
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

  /*
  EstimateCorrelation(model_settings,
                      wells,
                      well_names_,
                      mapped_blocked_logs_for_correlation,
                      interval_simboxes,
                      enough_data_for_corr_estimation_,
                      background,
                      err_txt);
                      */
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

// CRA-257: New implementation of correlation estimation for an interval
void  Analyzelog::EstimateCorrelation(const ModelSettings                                       * model_settings,
                                      const std::vector<NRLib::Well>                            & wells,
                                      const std::vector<std::string>                            & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>          & mapped_blocked_logs_for_correlation,
                                      std::string                                               & interval_name,
                                      const Simbox                                              * interval_simbox,
                                      bool                                                      & enough_data_for_corr_estimation,
                                      const std::vector<NRLib::Grid<float> *>                   & background,
                                      std::string                                               & err_txt){

  // Covariance and correlation estimation with blocked logs
  std::vector<std::vector<float> > log_data_vp(n_wells_);
  std::vector<std::vector<float> > log_data_vs(n_wells_);
  std::vector<std::vector<float> > log_data_rho(n_wells_);

  //
  // Find out if there is enough data for correlation estimation
  //

  int n_blocks_tot = 0;
  for (size_t i = 0; i < wells.size(); i++){
    std::string well_name_tmp = wells[i].GetWellName();
    n_blocks_tot += mapped_blocked_logs_for_correlation.find(well_name_tmp)->second->GetNBlocksWithData(interval_name);
  }
  enough_data_for_corr_estimation = (min_blocks_with_data_for_corr_estim_ <= n_blocks_tot);

  //
  // Check if all wells have synthetic Vs logs
  //

  bool allVsLogsAreSynthetic = true;
  for(int i=0 ; i<n_wells_ ; i++)
  {
    std::string well_name_tmp = wells[i].GetWellName();
    bool synt_vs_log =mapped_blocked_logs_for_correlation.find(well_name_tmp)->second->HasSyntheticVsLog();
    allVsLogsAreSynthetic = allVsLogsAreSynthetic && synt_vs_log;
  }
  if (allVsLogsAreSynthetic)
  {
    LogKit::LogFormatted(LogKit::Low,"\nThere are no nonsynthetic Vs logs available. Corr(Vp,Vs) and Corr(Vs,Rho) are set to 0.7.\n");
  }

  if(!enough_data_for_corr_estimation){
    LogKit::LogFormatted(LogKit::Low,"\nThere is not enough well data for estimation of prior correlations in interval:"+ interval_name + "\n");
    return;
  }

  size_t test = background.size();
  if(background[0]->GetN() == 0){
    //int ni = interval_simbox->getnx();
    //int nj = interval_simbox->getny();
    //int nk = interval_simbox->getnz();
    //background[0][0].Resize(ni, nj, nk);
    model_settings->getMaxHzBackground();

    EstimateLnData(log_data_vp, NULL, well_names_, mapped_blocked_logs_for_correlation, interval_name, "vp", err_txt);
    EstimateLnData(log_data_vs, NULL, well_names_, mapped_blocked_logs_for_correlation, interval_name, "vs", err_txt);
    EstimateLnData(log_data_rho, NULL, well_names_, mapped_blocked_logs_for_correlation, interval_name, "rho", err_txt);
  }
  else{
    EstimateLnData(log_data_vp, background[0], well_names_, mapped_blocked_logs_for_correlation, interval_name, "vp", err_txt);
    EstimateLnData(log_data_vs, background[1], well_names_, mapped_blocked_logs_for_correlation, interval_name, "vs", err_txt);
    EstimateLnData(log_data_rho, background[2], well_names_, mapped_blocked_logs_for_correlation, interval_name, "rho", err_txt);
  }

  /*
  EstimatePointVar0(point_var_0_, log_data_vp, log_data_vs, log_data_rho, err_txt);

  if (err_txt != "")
    return;

  int max_nd;
  CalculateNumberOfLags(n_lags_, max_nd, interval_simbox, err_txt);

  corr_T_.resize(n_lags_+1, 0);
  float dt = static_cast<float>(interval_simbox->getdz());

  EstimateCorrTAndVar0(corr_T_, var_0_,
                       log_data_vp, log_data_vs, log_data_rho,
                       allVsLogsAreSynthetic, dt,
                       n_lags_, max_nd,
                       err_txt);

  CheckVariances(model_settings, point_var_0_, var_0_, dt, err_txt);
  */
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

void
Analyzelog::CalculateNumberOfLags(int                                                   & n_lags,
                                  const std::vector<NRLib::Well>                        & wells,
                                  const std::map<std::string, BlockedLogsCommon *>      & mapped_blocked_logs_for_correlation,
                                  int                                                   & max_nd,
                                  const Simbox                                          * simbox,
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
  double dt      = simbox->getdz();
  double maxdist = 0.0;
  max_nd = 0;
  for(int i=0 ; i<n_wells_ ; i++)
  {
    int   nd                        = wells[i].GetNData();
    const std::vector<double> x_pos = mapped_blocked_logs_for_correlation.find(wells[i].GetWellName())->second->GetXposRawLogs();
    const std::vector<double> y_pos = mapped_blocked_logs_for_correlation.find(wells[i].GetWellName())->second->GetYposRawLogs();
    const std::vector<double> z_pos = mapped_blocked_logs_for_correlation.find(wells[i].GetWellName())->second->GetZposRawLogs();
    float * z = new float[nd];
    if(nd>max_nd)
      max_nd=nd;
    for(int j=0 ; j<nd ; j++)
    {
      float xj = static_cast<float>(x_pos[j]);
      float yj = static_cast<float>(y_pos[j]);
      z[j] = static_cast<float>((z_pos[j]-simbox->getTop(xj,yj))/simbox->getRelThick(xj,yj));
    }
    std::string wellErrTxt = "";
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

void Analyzelog::EstimateLnData(std::vector<std::vector<float> >                    & ln_data,
                                const NRLib::Grid<float>                            * background,
                                const std::vector<std::string>                      & well_names,
                                const std::map<std::string, BlockedLogsCommon *>    & mapped_blocked_logs_for_correlation,
                                const std::string                                   & interval_name,
                                const std::string                                   & log_name,
                                std::string                                         & err_txt){

  float global_mean = 0.0;
  int   count       = 0;
  int   nd          = 0;
  int   log_nr      = -1;

  assert(log_name == "vp" || log_name == "vs" || log_name == "rho");

  for (int i = 0; i < n_wells_; i++){

    nd                              = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNBlocksWithData(interval_name);
    //const std::vector<double> x_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetXposBlocked();
    //const std::vector<double> y_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetYposBlocked();
    //const std::vector<double> z_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetZposBlocked();
    const std::vector<int>    & i_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetIposVector();
    const std::vector<int>    & j_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetJposVector();
    const std::vector<int>    & k_pos = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetKposVector();

    std::vector<double> blocked_well_log;

    if(log_name == "vp"){
      blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVpBlocked();
      log_nr               = 0;
    }
    else if(log_name == "vs"){
      blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetVsBlocked();
      log_nr               = 1;
    }
    else if(log_name == "rho"){
      blocked_well_log     = mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetRhoBlocked();
      log_nr               = 2;
    }

    std::vector<double> mean(nd);
    std::vector<double> mean_max_hz(nd);

    if (background != NULL){
        for(int n = 0; n < nd; n++){
          mean[n]         = (*background)(i_pos[n], j_pos[n], k_pos[n]);
          mean_max_hz[n]  = 0;
        }
    }
    else{
      for(int n = 0; n < nd; n++){
        mean[n]         = 0;
        mean_max_hz[n]  = 0;
      }
    }

    ln_data[i].resize(nd, 0.0);

    for (int j = 0; j < nd; j++){
      if(blocked_well_log[j] != RMISSING && mean[j] != RMISSING){
        ln_data[i][j]  = static_cast<float>(log(blocked_well_log[j]) - mean[j]);
        global_mean   += ln_data[i][j];
        count++;
      }
      else{
        ln_data[i][j] = RMISSING;
      }
    }
  }

  //
  // Subtract global mean from data.
  //
  if(count > 0){
    global_mean /= count;
    for(int i = 0 ; i < n_wells_ ; i++){
      for(int j = 0 ; j < mapped_blocked_logs_for_correlation.find(well_names[i])->second->GetNBlocksWithData(interval_name) ; j++){
        if (ln_data[i][j] != RMISSING){
          ln_data[i][j] -= global_mean;
        }
      }
    }
  }
  else{
    err_txt += std::string("Could not estimate globalMean for log") + NRLib::ToString(log_nr) + " (Vp=0,Vs=1,Rho=2)\n";
  }
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

/*
void
Analyzelog::checkVariances(const ModelSettings  * modelSettings,
                           const float  * const * pointVar0,
                           const float  * const * Var0,
                           float            dt,
                           std::string    & errTxt)
{
  //| These min and max values below are used for consistency check. If a variance
  //| is outside these ranges there is probably a problem with the log.
  //|
  //| The limits are for point variances. The minimum allowed variance
  //| for parameters will be scaled with 1/dt*dt

  float minVarAlpha = modelSettings->getVarAlphaMin();
  float maxVarAlpha = modelSettings->getVarAlphaMax();
  float minVarBeta  = modelSettings->getVarBetaMin();
  float maxVarBeta  = modelSettings->getVarBetaMax();
  float minVarRho   = modelSettings->getVarRhoMin();
  float maxVarRho   = modelSettings->getVarRhoMax();

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

  //
  // We scale the minimum variances allowed with 1/dt since the variance decreases with increasing dt..
  //
  if (Var0[0][0] < minVarAlpha/dt || Var0[0][0] > maxVarAlpha)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vp variance "  << Var0[0][0]
      << " is outside allowed interval Min=" << minVarAlpha/dt << " Max=" << maxVarAlpha << "\n";
    errTxt += o.str();
  }
  if (Var0[1][1] < minVarBeta/dt || Var0[1][1] > maxVarBeta)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Vs variance "  << Var0[1][1]
      << " is outside allowed interval Min=" << minVarBeta/dt << " Max=" << maxVarBeta << "\n";
    errTxt += o.str();
  }
  if (Var0[2][2] < minVarRho/dt || Var0[2][2] > maxVarRho)
  {
    std::ostringstream o;
    o << std::scientific << std::setprecision(2) << "The Rho variance "  << Var0[2][2]
      << " is outside allowed interval Min=" << minVarRho/dt << " Max=" << maxVarRho << "\n";
    errTxt += o.str();
  }
  if (errTxt != "")
  {
    LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------------------------------------");
    LogKit::LogFormatted(LogKit::Low,"\n                          ln Vp     ln Vs    ln Rho ");
    LogKit::LogFormatted(LogKit::Low,"\nWell log  variances:   %.2e  %.2e  %.2e",pointVar0[0][0],pointVar0[1][1],pointVar0[2][2]);
    LogKit::LogFormatted(LogKit::Low,"\nParameter variances:   %.2e  %.2e  %.2e (used by program)",Var0[0][0],Var0[1][1],Var0[2][2]);
    LogKit::LogFormatted(LogKit::Low,"\n--------------------------------------------------------------------\n");
  }
}
*/
