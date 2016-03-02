/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#ifdef PARALLEL
#include <omp.h>
#endif

#include "lib/kriging1d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/random/beta.hpp"
#include "nrlib/random/distribution.hpp"

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/background.h"
#include "src/covgrid2d.h"
#include "src/krigingdata2d.h"
#include "src/kriging2d.h"
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/io.h"
#include "src/multiintervalgrid.h"

void
Background::SetupBackground(std::vector<NRLib::Grid<float> *>                & parameters,
                            std::vector<std::vector<double> >                & vertical_trends,
                            NRLib::Grid<float>                               * velocity,
                            const Simbox                                     * simbox,
                            const Simbox                                     * bg_simbox,
                            const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                            const std::map<std::string, BlockedLogsCommon *> & bg_blocked_logs,
                            const ModelSettings                              * model_settings,
                            const std::string                                & interval_name,
                            std::string                                      & err_text)
{

  int nx_pad, ny_pad, nz_pad;

  if (bg_simbox == NULL) {
    nx_pad = simbox->getnx();
    ny_pad = simbox->getny();
    nz_pad = simbox->getnz();
  }
  else {
    nx_pad = bg_simbox->getnx();
    ny_pad = bg_simbox->getny();
    nz_pad = bg_simbox->getnz();
  }

  for (int i=0 ; i<3 ; i++)
    parameters[i]->Resize(nx_pad, ny_pad, nz_pad);

  if (bg_simbox == NULL) {
    GenerateBackgroundModel(parameters[0], parameters[1], parameters[2],
                            vertical_trends,
                            velocity,
                            simbox,
                            blocked_logs,
                            model_settings,
                            interval_name,
                            err_text);
  }
  else {
    GenerateBackgroundModel(parameters[0], parameters[1], parameters[2],
                            vertical_trends,
                            velocity,
                            bg_simbox,
                            bg_blocked_logs,
                            model_settings,
                            interval_name,
                            err_text);

    ResampleBackgroundModel(parameters[0], parameters[1], parameters[2],
                            bg_simbox,
                            simbox);

    float avg = 0.0f;
    float min = 0.0f;
    float max = 0.0f;
    for (int i = 0; i < 3; i++) {
      parameters[i]->GetAvgMinMaxWithMissing(avg, min, max, RMISSING);
      CommonData::SetUndefinedCellsToGlobalAverageGrid(parameters[i], avg);
    }

  }
}

//-------------------------------------------------------------------------------
void
Background::GenerateBackgroundModel(NRLib::Grid<float>                               * bg_vp,
                                    NRLib::Grid<float>                               * bg_vs,
                                    NRLib::Grid<float>                               * bg_rho,
                                    std::vector<std::vector<double> >                & vertical_trends,
                                    NRLib::Grid<float>                               * velocity,
                                    const Simbox                                     * simbox,
                                    const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                    const ModelSettings                              * model_settings,
                                    const std::string                                & interval_name,
                                    std::string                                      & err_text)
{
  const int   nz      = simbox->getnz();
  const int   n_wells = static_cast<int>(blocked_logs.size());
  const float dz      = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  std::string name_vp  = "Vp";
  std::string name_vs  = "Vs";
  std::string name_rho = "Rho";

  std::vector<std::vector<double> > well_trend_vp(n_wells);
  std::vector<std::vector<double> > well_trend_vs(n_wells);
  std::vector<std::vector<double> > well_trend_rho(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_vp(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_vs(n_wells);
  std::vector<std::vector<double> > high_cut_well_trend_rho(n_wells);

  std::string err_text_tmp = "";

  std::vector<std::string> well_names(n_wells);
  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    well_names[w] = iter->second->GetWellName();
    w++;
  }

  GetWellTrends(well_trend_vp,  high_cut_well_trend_vp,  blocked_logs, nz, name_vp,  err_text_tmp);
  GetWellTrends(well_trend_vs,  high_cut_well_trend_vs,  blocked_logs, nz, name_vs,  err_text_tmp);
  GetWellTrends(well_trend_rho, high_cut_well_trend_rho, blocked_logs, nz, name_rho, err_text_tmp);

  if (err_text_tmp != "") {
    err_text += err_text_tmp;
  }
  else {
    std::vector<double> trend_vp(nz);
    std::vector<double> trend_vs(nz);
    std::vector<double> trend_rho(nz);
    std::vector<double> trend_vel(nz);

    std::vector<double> avg_dev_vp(n_wells);
    std::vector<double> avg_dev_vs(n_wells);
    std::vector<double> avg_dev_rho(n_wells);
    std::vector<double> avg_dev_vel(n_wells);

    CalculateBackgroundTrend(trend_vp,
                             avg_dev_vp,
                             nz,
                             dz,
                             model_settings->getVpMin(),
                             model_settings->getVpMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_vp,
                             high_cut_well_trend_vp,
                             name_vp);
    CalculateBackgroundTrend(trend_vs,
                             avg_dev_vs,
                             nz,
                             dz,
                             model_settings->getVsMin(),
                             model_settings->getVsMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_vs,
                             high_cut_well_trend_vs,
                             name_vs);
    CalculateBackgroundTrend(trend_rho,
                             avg_dev_rho,
                             nz,
                             dz,
                             model_settings->getRhoMin(),
                             model_settings->getRhoMax(),
                             model_settings->getMaxHzBackground(),
                             well_trend_rho,
                             high_cut_well_trend_rho,
                             name_rho);

    if (velocity->GetN() != 0) {
      //
      // We still want CalculateBackgroundTrend() for alpha above. By calculating
      // avgDevAlpha we can check that the bgAlpha calculated from velocity is as
      // good as or better than that calculated by crava.
      //
      CalculateVelocityDeviations(velocity,
                                  simbox, blocked_logs,
                                  trend_vel, avg_dev_vel, avg_dev_vp,
                                  //model_settings->getOutputGridsElastic(),
                                  n_wells);

      velocity->LogTransform(RMISSING);

      delete bg_vp;
      bg_vp = velocity;
      velocity->Resize(0, 0, 0, 0.0);
      WriteDeviationsFromVerticalTrend(avg_dev_vel, avg_dev_vs, avg_dev_rho,
                                       trend_vel, trend_vs, trend_rho,
                                       blocked_logs, n_wells, nz);
    }
    else {
      WriteDeviationsFromVerticalTrend(avg_dev_vp, avg_dev_vs, avg_dev_rho,
                                       trend_vp, trend_vs, trend_rho,
                                       blocked_logs, n_wells, nz);
    }

    std::vector<KrigingData2D> kriging_data_vp(nz);
    std::vector<KrigingData2D> kriging_data_vs(nz);
    std::vector<KrigingData2D> kriging_data_rho(nz);

    std::vector<std::vector<double > > bl_vp(n_wells); // bl = blocked logs
    std::vector<std::vector<double > > bl_vs(n_wells);
    std::vector<std::vector<double > > bl_rho(n_wells);
    std::vector<std::vector<double > > vt_vp(n_wells);   // vt = vertical trend
    std::vector<std::vector<double > > vt_vs(n_wells);
    std::vector<std::vector<double > > vt_rho(n_wells);

    std::vector<const std::vector<int> *> ipos(n_wells);
    std::vector<const std::vector<int> *> jpos(n_wells);
    std::vector<const std::vector<int> *> kpos(n_wells);

    for (int i=0; i < n_wells; i++) {
      vt_vp[i]  = std::vector<double>(nz);
      vt_vs[i]  = std::vector<double>(nz);
      vt_rho[i] = std::vector<double>(nz);
    }

    std::vector<int> n_blocks(n_wells);
    int              tot_blocks;

    GetKrigingWellTrends(bl_vp,bl_vs,bl_rho,
                         vt_vp,vt_vs,vt_rho,
                         ipos,jpos,kpos,
                         n_blocks,tot_blocks,
                         blocked_logs, n_wells);

    SetupKrigingData2D(kriging_data_vp,kriging_data_vs,kriging_data_rho,
                       trend_vp,trend_vs,trend_rho,
                       model_settings->getOutputGridsElastic(),
                       nz,dz,tot_blocks,n_blocks,
                       bl_vp,bl_vs,bl_rho,
                       vt_vp,vt_vs,vt_rho,
                       ipos,jpos,kpos,
                       interval_name, well_names);

    const CovGrid2D & covGrid2D = MakeCovGrid2D(simbox,
                                                model_settings->getBackgroundVario(),
                                                model_settings->getDebugFlag(),
                                                interval_name);

    MakeKrigedBackground(kriging_data_vp,  bg_vp,  trend_vp,  simbox, covGrid2D, "Vp",  model_settings->getNumberOfThreads());
    MakeKrigedBackground(kriging_data_vs,  bg_vs,  trend_vs,  simbox, covGrid2D, "Vs",  model_settings->getNumberOfThreads());
    MakeKrigedBackground(kriging_data_rho, bg_rho, trend_rho, simbox, covGrid2D, "Rho", model_settings->getNumberOfThreads());

    vertical_trends[0] = trend_vp;
    vertical_trends[1] = trend_vs;
    vertical_trends[2] = trend_rho;

    delete &covGrid2D;
  }
}

void
Background::CalculateVelocityDeviations(NRLib::Grid<float>                               * velocity,
                                        const Simbox                                     * simbox,
                                        const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                        std::vector<double>                              & trend_vel,
                                        std::vector<double>                              & avg_dev_vel,
                                        std::vector<double>                              & avg_dev_vp,
                                        //int                                                output_flag,
                                        int                                                n_wells)
{
  //H Writing of NRLib::Grid missing.
  //if ((output_flag & IO::BACKGROUND_TREND) > 0) {
  //  std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
  //  velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
  //}

  //
  // Calculate deviation between well data and trend
  //
  int max_blocks = 0;
  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    int n_blocks = iter->second->GetNumberOfBlocks();

    if (n_blocks > max_blocks)
      max_blocks = n_blocks;

    w++;
  }

  std::vector<double> velocity_log(max_blocks);

  const int nz        = simbox->getnz();
  std::vector<double> vt_vp(nz);
  std::vector<double> vt_velocity(nz);

  for (int k=0 ; k<nz ; k++)
    trend_vel[k]=0.0;

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    const std::vector<double> & vp_log = blocked_log->GetVpHighCutBackground();
    blocked_log->GetVerticalTrend(vp_log, vt_vp);
    blocked_log->GetBlockedGrid(velocity, velocity_log);
    blocked_log->GetVerticalTrend(velocity_log, vt_velocity);
    double sumDev = 0.0;
    int count = 0;
    for (int k = 0; k < nz; k++) {
      if (vt_vp[k] != RMISSING) {
        trend_vel[k] += vt_velocity[k];
        double diff = exp(vt_vp[k]) - vt_velocity[k]; // Velocity trend is in exp-domain
        sumDev += diff*diff;
        count++;
      }
    }
    if (count > 0)
      sumDev /= count;
    avg_dev_vel[w] = sqrt(sumDev);

    w++;
  }

  for (int k = 0; k < nz; k++)
    trend_vel[k] /= n_wells;

  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
                         iter->second->GetWellName().c_str(),avg_dev_vel[w],avg_dev_vp[w]);

    w++;
  }
}

//---------------------------------------------------------------------------
void
Background::CalculateBackgroundTrend(std::vector<double>               & trend,
                                     std::vector<double>               & avg_dev,
                                     const int                           nz,
                                     const float                         dz,
                                     float                               log_min,
                                     float                               log_max,
                                     float                               max_hz,
                                     std::vector<std::vector<double> > & well_trend,
                                     std::vector<std::vector<double> > & high_cut_well_trend,
                                     const std::string                 & name)
{

  CalculateVerticalTrend(well_trend,
                         trend,
                         log_min,
                         log_max,
                         max_hz,
                         nz,
                         dz,
                         name);


  CalculateDeviationFromVerticalTrend(high_cut_well_trend, trend, avg_dev, nz);


}

//---------------------------------------------------------------------------
void
Background::GetKrigingWellTrends(std::vector<std::vector<double> >                & bl_vp,
                                 std::vector<std::vector<double> >                & bl_vs,
                                 std::vector<std::vector<double> >                & bl_rho,
                                 std::vector<std::vector<double> >                & vt_vp,
                                 std::vector<std::vector<double> >                & vt_vs,
                                 std::vector<std::vector<double> >                & vt_rho,
                                 std::vector<const std::vector<int> *>            & ipos,
                                 std::vector<const std::vector<int> *>            & jpos,
                                 std::vector<const std::vector<int> *>            & kpos,
                                 std::vector<int>                                 & n_blocks,
                                 int                                              & tot_blocks,
                                 const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                 const int                                        & n_wells)
{
  int max_blocks = 0;
  tot_blocks     = 0;

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);

    n_blocks[w] = iter->second->GetNumberOfBlocks();
    tot_blocks += n_blocks[w];
    if (n_blocks[w] > max_blocks)
      max_blocks = n_blocks[w];

    w++;
  }

  for (int i = 0; i < n_wells; i++) {
    bl_vp[i]  = std::vector<double>(max_blocks);
    bl_vs[i]  = std::vector<double>(max_blocks);
    bl_rho[i] = std::vector<double>(max_blocks);
  }

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    bl_vp[w]  = blocked_log->GetVpHighCutBackground();
    bl_vs[w]  = blocked_log->GetVsHighCutBackground();
    bl_rho[w] = blocked_log->GetRhoHighCutBackground();

    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    blocked_log->GetVerticalTrend(bl_vp[w],  vt_vp[w]);
    blocked_log->GetVerticalTrend(bl_vs[w],  vt_vs[w]);
    blocked_log->GetVerticalTrend(bl_rho[w], vt_rho[w]);

    ipos[w] = &(blocked_log->GetIposVector());
    jpos[w] = &(blocked_log->GetJposVector());
    kpos[w] = &(blocked_log->GetKposVector());

    w++;
  }
}

//---------------------------------------------------------------------------
void
Background::GetWellTrends(std::vector<std::vector<double> >                & well_trend,
                          std::vector<std::vector<double> >                & high_cut_well_trend,
                          const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                          const int                                        & nz,
                          const std::string                                & name,
                          std::string                                      & err_text)
{
  int i_wells = 0;

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log->GetUseForBackgroundTrend()) {

      if (blocked_log != NULL) {
        well_trend[w].resize(nz);

        if (name == "Vp")
          blocked_log->GetVerticalTrend(blocked_log->GetVpBlocked(), well_trend[w]);
        else if (name == "Vs")
          blocked_log->GetVerticalTrend(blocked_log->GetVsBlocked(), well_trend[w]);
        else if (name == "Rho")
          blocked_log->GetVerticalTrend(blocked_log->GetRhoBlocked(), well_trend[w]);
        else {
          err_text += "ERROR in Background::GetWellTrends(): ";
          err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
        }
        i_wells++;
      }
      else
        well_trend[w].resize(0);
    }
    else
      well_trend[w].resize(0);

    w++;
  }
  if (i_wells == 0) {
    err_text += "\nERROR in Background::GetWellTrends(): There are no wells\n";
    err_text += "available for the estimation of background trend.\n";
  }

  w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log != NULL) {
      high_cut_well_trend[w].resize(nz);
      if (name == "Vp")
        blocked_log->GetVerticalTrend(blocked_log->GetVpHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Vs")
        blocked_log->GetVerticalTrend(blocked_log->GetVsHighCutBackground(), high_cut_well_trend[w]);
      else if (name == "Rho")
        blocked_log->GetVerticalTrend(blocked_log->GetRhoHighCutBackground(), high_cut_well_trend[w]);
      else {
        err_text += "ERROR in Background::GetWellTrends(): ";
        err_text += "Log \'"+name+"\' requested, but no such log exists.\n";
      }
    }
    else
      high_cut_well_trend[w].resize(0);

    w++;
  }
}

//---------------------------------------------------------------------------
void
Background::WriteTrendsToFile(std::vector<double> & trend,
                              const Simbox        * simbox,
                              bool                  write1D,
                              const std::string   & name,
                              const std::string   & interval_name)
{
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());
  const int   nz = simbox->getnz();

  if (write1D == true) {
    WriteVerticalTrend(trend, dz, nz, name, interval_name);
  }
}

//-------------------------------------------------------------------------------
void
Background::SetupKrigingData2D(std::vector<KrigingData2D>                & kriging_data_vp,
                               std::vector<KrigingData2D>                & kriging_data_vs,
                               std::vector<KrigingData2D>                & kriging_data_rho,
                               std::vector<double>                       & trend_vp,
                               std::vector<double>                       & trend_vs,
                               std::vector<double>                       & trend_rho,
                               const int                                   output_flag,
                               const int                                 & nz,
                               const float                               & dz,
                               const int                                 & tot_blocks,
                               const std::vector<int>                    & n_blocks,
                               const std::vector<std::vector<double> >   & bl_vp,
                               const std::vector<std::vector<double> >   & bl_vs,
                               const std::vector<std::vector<double> >   & bl_rho,
                               const std::vector<std::vector<double> >   & vt_vp,
                               const std::vector<std::vector<double> >   & vt_vs,
                               const std::vector<std::vector<double> >   & vt_rho,
                               const std::vector<const std::vector<int> *> ipos,
                               const std::vector<const std::vector<int> *> jpos,
                               const std::vector<const std::vector<int> *> kpos,
                               const std::string                         & interval_name,
                               const std::vector<std::string>            & well_names)
{
  //
  // Although unnecessary, we have chosen to set up kriging data from
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  const int n_wells = static_cast<int>(bl_vp.size());

  KrigingData3D forLogging(tot_blocks);

  for (int w = 0; w < n_wells; w++) {

    if (vt_vp[w].size() > 0) {
      std::vector<double> vt_vp_well  = vt_vp[w];
      std::vector<double> vt_vs_well  = vt_vs[w];
      std::vector<double> vt_rho_well = vt_rho[w];

      std::vector<double> bl_vp_well  = bl_vp[w];
      std::vector<double> bl_vs_well  = bl_vs[w];
      std::vector<double> bl_rho_well = bl_rho[w];

      //Check if vertical trends only contains missing values
      CheckLogForOnlyMissing(vt_vp_well,  "vp",  well_names[w]);
      CheckLogForOnlyMissing(vt_vs_well,  "vs",  well_names[w]);
      CheckLogForOnlyMissing(vt_rho_well, "rho", well_names[w]);

      //
      // Kriging vertical trend (vt....) against global vertical trend (trend...)
      //
      Kriging1D::krigVector(&vt_vp_well[0],  &trend_vp[0],  nz, dz);
      Kriging1D::krigVector(&vt_vs_well[0],  &trend_vs[0],  nz, dz);
      Kriging1D::krigVector(&vt_rho_well[0], &trend_rho[0], nz, dz);
      //
      // Use kriged vertical trend where original log is not defined.
      //
      const std::vector<int> ipos_well = *(ipos[w]);
      const std::vector<int> jpos_well = *(jpos[w]);
      const std::vector<int> kpos_well = *(kpos[w]);

      for (int m = 0; m < n_blocks[w]; m++) {
        int i = ipos_well[m];
        int j = jpos_well[m];
        int k = kpos_well[m];

        if (bl_vp_well[m] == RMISSING)
          bl_vp_well[m] = vt_vp_well[k];
        if (bl_vs_well[m] == RMISSING)
          bl_vs_well[m] = vt_vs_well[k];
        if (bl_rho_well[m] == RMISSING)
          bl_rho_well[m] = vt_rho_well[k];

        kriging_data_vp[k].addData(i, j, static_cast<float>(bl_vp_well[m]));
        kriging_data_vs[k].addData(i, j, static_cast<float>(bl_vs_well[m]));
        kriging_data_rho[k].addData(i, j, static_cast<float>(bl_rho_well[m]));
      }

      forLogging.addData(bl_vp_well, bl_vs_well, bl_rho_well,
                         ipos_well, jpos_well, kpos_well,
                         n_blocks[w]);
    }

    for (int k = 0; k < nz; k++) {
      kriging_data_vp[k].findMeanValues();
      kriging_data_vs[k].findMeanValues();
      kriging_data_rho[k].findMeanValues();
    }
  }

  if ((output_flag & IO::BACKGROUND) > 0) {
    forLogging.divide();

    std::string interval_text = "";
    if (interval_name != "")
      interval_text = "_" + interval_name;

    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + interval_text + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
}

//---------------------------------------------------------------------------
const CovGrid2D &
Background::MakeCovGrid2D(const Simbox      * simbox,
                          Vario             * vario,
                          int                 debug_flag,
                          const std::string & interval_name)
{
  //
  // Pretabulate all needed covariances
  //
  const int    nx = simbox->getnx();
  const int    ny = simbox->getny();

  const float  dx = static_cast<float>(simbox->getdx());
  const float  dy = static_cast<float>(simbox->getdy());

  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);

  if (debug_flag == 1) {

    std::string interval_text = "";
    if (interval_name != "")
      interval_text = "_" + interval_name;

    std::string base_name = IO::PrefixBackground() + "covGrid2D" + interval_text + IO::SuffixAsciiIrapClassic();
    std::string file_name = IO::makeFullFileName(IO::PathToBackground(), base_name);
    cov->writeToFile(file_name);
  }
  return (*cov);
}

void
Background::MakeKrigedBackground(const std::vector<KrigingData2D> & kriging_data,
                                 NRLib::Grid<float>               * bg_grid,
                                 std::vector<double>              & trend,
                                 const Simbox                     * simbox,
                                 const CovGrid2D                  & cov_grid_2D,
                                 const std::string                & type,
                                 int                                n_threads)
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

  //
  // Store a surface for each layer (needed for parallelization)
  //
  NRLib::RegularSurface<double> tmp(x0, y0, lx, ly, nx, ny, RMISSING);
  std::vector<NRLib::RegularSurface<double> > surfaces(0);
  surfaces.reserve(nz);
  for (int k=0 ; k<nz ; k++)
    surfaces.push_back(tmp);

  float monitor_size = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float next_monitor = monitor_size;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  bg_grid->Resize(nx, ny, nz);
#ifdef PARALLEL
  int  chunk_size = 1;
#pragma omp parallel for schedule(dynamic, chunk_size) num_threads(n_threads)
#endif

  for (int k=0 ; k<nz ; k++) {
    // Set trend for layer
    surfaces[k].Assign(trend[k]);

    // Kriging of layer
    Kriging2D::krigSurface(surfaces[k], kriging_data[k], cov_grid_2D);

    // Log progress
    if (k+1 >= static_cast<int>(next_monitor)) {
      next_monitor += monitor_size;
      std::cout << "^";
      fflush(stdout);
    }
  }

  // Set layers in background model from surface
  for (int k=0 ; k<nz ; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0 ; i < nx ; i++) {
        bg_grid->SetValue(i, j, k, static_cast<float>(surfaces[k](i,j)));
      }
    }
  }
}

//-------------------------------------------------------------------------------
void
Background::CalculateVerticalTrend(std::vector<std::vector<double> > & well_trend,
                                   std::vector<double>               & trend,
                                   float                               log_min,
                                   float                               log_max,
                                   float                               max_hz,
                                   int                                 nz,
                                   float                               dz,
                                   const std::string                 & name)
{
  int     n_wells      = static_cast<int>(well_trend.size());
  std::vector<double> filtered_log(nz);
  int   * count        = new int[nz];
  //
  // Calculate the average values of well log
  // ----------------------------------------
  // When calculating the vertical trend, we do not want each well
  // to contribute with more than one value to each layer 'k',
  // and therefore we calculate average values for each well
  // first. This way we avoid strange behaviour caused by
  // deviated/horizontal wells.
  //
  for (int k = 0; k < nz; k++) {
    trend[k] = 0.0f;
    count[k] = 0;
  }
  int i_wells = 0;
  for (int w = 0; w < n_wells; w++) {
    if (well_trend[w].size() > 0) {
      std::vector<double> & w_trend = well_trend[w];

      for (int k = 0; k < nz; k++) {
        if (w_trend[k] != RMISSING) {
          trend[k] += exp(w_trend[k]);
          count[k]++;
        }
      }
      i_wells++;
    }
  }
  if (i_wells > 0) {
    for (int k = 0; k < nz; k++) {
      if (count[k] > 0) {
        trend[k] = trend[k]/static_cast<float>(count[k]);
      }
    }
  }

  Utils::writeVectorToFile(std::string("trend_mean_values_") + name, trend);

  SmoothTrendWithLocalLinearRegression(trend, count,
                                       i_wells, nz, dz,
                                       log_min,
                                       log_max,
                                       name);

  Utils::writeVectorToFile(std::string("trend_after_linreg_") + name, trend);

  CommonData::ApplyFilter(filtered_log,
                          trend,
                          nz,
                          dz,
                          max_hz);

  for (int i = 0; i < nz; i++) {
    trend[i] = filtered_log[i];
  }

  Utils::writeVectorToFile(std::string("trend_after_filter_") + name, trend);

  delete [] count;
}

//-------------------------------------------------------------------------------
void
Background::SmoothTrendWithLocalLinearRegression(std::vector<double> & trend,
                                                 int                 * count,
                                                 int                   i_wells,
                                                 int                   nz,
                                                 float                 dz,
                                                 float                 min_value,
                                                 float                 max_value,
                                                 std::string           par_name)
{
  bool debug = false;
  //
  // 1. Center-parts of scatter plots
  //
  // In the center parts of the scatter plots the average value should be
  // accepted as a trend value if the number of data points behind each
  // average is fraction * iWells, where fraction is the acceptance
  // fraction, typically larger than one.
  //
  // Sometimes we have only one, two, or three wells available. In the
  // case of two and three wells, the logs for these may differ considerably,
  // in which case the trend need to be stabilised by requiring a minimum
  // number of data points behind each trend value. The denser the log is
  // sampled, the more prone it is to numerical instabilities. This
  // minimum value is therefore linked to the sampling density 'dz'.
  //
  // Finally, we must require a definite minimum number of data points
  // that should be behind each trend value.
  //
  // nDataMin = max(min_req_points, max(min_time_sample, fraction * iWells))
  //
  // 2. End-parts of scatter plot
  //
  // If the blocked wells do not contain any values for the lower layers of
  // the simbox, the end part of the scatter plot needs special attention.
  //
  // We must possibly require a larger min_req_points to avoid an
  // "arbitrary" end behaviour.
  //

  float fraction     = 5.0f;                      // Require minimum 5*iWells
  int   n_time_limit = static_cast<int>(50.0/dz); // The smaller sampling density, the more values are needed.
  int   n_low_limit  = 10;                        // Require minimum 10
  int   n_data_min   = std::max(n_low_limit, std::max(n_time_limit, int(fraction * i_wells)));

  bool  use_weights  = true;
  bool  error_mid    = false;
  bool  error_head   = false;
  bool  error_trail  = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  float * mean = new float[nz];
  for (int k = 0; k < nz; k++) {
    mean[k] = static_cast<float>(trend[k]);
  }

  //
  // Find first non-missing value
  //
  int first_nonmissing = 0;
  for (int k = 0; k < nz; k++) {
    if (trend[k] > 0.0f) {
      first_nonmissing = k;
      break;
    }
  }

  //
  // Find last non-missing value
  //
  int last_nonmissing = nz - 1;
  for (int k = nz - 1; k > 0; k--) {
    if (trend[k] > 0.0f) {
      last_nonmissing = k;
      break;
    }
  }

  float * x = new float[nz];  // Time indices
  float * y = new float[nz];  // Log values
  float * w = new float[nz];  // Weights (number of data behind each avg.)

  for (int k = 0; k < nz; k++) {
    int n_cur_data_min = n_data_min;
    if (k < first_nonmissing || k > last_nonmissing) {
      n_cur_data_min *= 2;
    }

    int n      = 0;
    int n_data = 0;
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"k=%d\n",k);
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<float>(count[k]);
      x[0]   = static_cast<float>(k);
      y[0]   = static_cast<float>(trend[k]);
      n_data += count[k];
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   A:t=%.2f   x[0] y[0]  %d   %.2f\n",dz*(x[0] + 0.5f),int(x[0]),y[0]);
      n++;
    }

    //
    // 2. Add local data points to get 'nCurDataMin' points behind each trend
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (n_data < n_cur_data_min) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]    = static_cast<float>(count[k - i]);
        x[n]    = static_cast<float>(k - i);
        y[n]    = mean [k - i];
        n_data += count[k - i];
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   B:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]    = static_cast<float>(count[k + i]);
        x[n]    = static_cast<float>(k + i);
        y[n]    = mean [k + i];
        n_data += count[k + i];
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   C:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k-i < 0 && k+i >= nz) { // We will never find enough data
        break;
      }
    }

    //
    // Calculate normalised weights
    //
    if (use_weights)
      for (i = 0; i < n; i++)
        w[i] /= n_data;
    else
      for (i = 0; i < n; i++)
        w[i] = 1.0f/n;

    //
    // We need at least two points to make a line.
    //
    if (n > 1) {
      //
      // Estimate local regression line: y = bx + a
      //
      float Sx  = x[0]*w[0];
      float Sy  = y[0]*w[0];
      float Sxx = x[0]*w[0]*x[0];
      float Sxy = x[0]*w[0]*y[0];
      for (i = 1 ; i < n ; i++) {
        Sx  += x[i]*w[i];
        Sy  += y[i]*w[i];
        Sxx += x[i]*w[i]*x[i];
        Sxy += x[i]*w[i]*y[i];
      }
      float b = (Sxy - Sx*Sy)/(Sxx - Sx*Sx);
      float a = (Sy - b*Sx);
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"Sx, Sy, Sxx, Sxy  : %.4f, %.4f, %.4f, %.4f     a, b : %.4f %.4f\n",Sx,Sy,Sxx,Sxy,a,b);

      //
      // Estimate value of regression line at requested point.
      //
      float value = a + b*static_cast<float>(k);
      if (value < min_value || value > max_value) {
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f\n",value);
        if (k < first_nonmissing)
          error_head = true;
        else if (k > last_nonmissing)
          error_trail = true;
        else {
          error_mid   = true;
          break;
        }
      }
      trend[k] = log(static_cast<double>(value));
    }
    else {
      trend[k] = log(static_cast<double>(y[0]));
    }
    if (debug)
      LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f        (minLog/maxLog = %.2f / %.2f)\n",exp(trend[k]),min_value,max_value);
  }

  if (error_mid) {
    // Big problem ...
    LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+par_name+" using local linear\n");
    LogKit::LogFormatted(LogKit::Low,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::Low,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::Low,"          2) There are too many layers in grid compared to well logs available.\n");
    float sum  = 0.0f;
    int n_data = 0;
    for (int k = 0; k < nz; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum    += mean[k]*count[k];
          n_data += count[k];
          if (debug)
            LogKit::LogFormatted(LogKit::Low,"k=%d  count[k], mean[k]  nData, sum  %d  %8.3f     %d  %8.3f\n",
                                 k,count[k],mean[k],n_data,sum);
        }
        else {
          sum    += mean[k];
          n_data += 1;
        }
      }
    }
    float global_mean = sum/n_data;
    for (int k = 0; k < nz; k++) {
      trend[k] = log(global_mean);
      if (debug)
        LogKit::LogFormatted(LogKit::Low,"   TREND: k = %d   trend[k] = %.2f\n",k,exp(trend[k]));
    }
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean for parameter %s = %.2f\n\n",par_name.c_str(),global_mean);
  }
  else {
    if (error_head) {
      // Fix first part of trend containing missing-values.
      double first_value = trend[first_nonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter \'"+par_name+"\' using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [0-%d] where the log is undefined. The first\n",first_nonmissing-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(first_value),first_nonmissing);
      for (int k = 0 ; k < first_nonmissing ; k++) {
        trend[k] = first_value;
      }
    }
    if (error_trail) {
      // Fix last part of trend containing missing-values.
      double last_value = trend[last_nonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+par_name+" using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",last_nonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d will be used throughout this region.\n",exp(last_value),last_nonmissing);
      for (int k = last_nonmissing + 1; k < nz; k++) {
        trend[k] = last_value;
      }
    }
  }
  delete [] x;
  delete [] y;
  delete [] w;
  delete [] mean;
}

//-------------------------------------------------------------------------------
void
Background::WriteVerticalTrend(std::vector<double> & trend,
                               float                 dz,
                               int                   nz,
                               std::string           name,
                               const std::string   & interval_name)
{
  float z0 = dz/2.0f;

  std::string interval_text = "";
  if (interval_name != "")
    interval_text = "_" + interval_name;

  std::string base_name = IO::PrefixBackground() + IO::PrefixTrend() + name + interval_text + IO::SuffixAsciiIrapClassic();
  std::string file_name = IO::makeFullFileName(IO::PathToBackground(), base_name);
  std::ofstream file;
  NRLib::OpenWrite(file, file_name);
  for (int i = 0; i < nz; i++) {
    file << std::fixed
         << std::setprecision(2)
         << std::setw(8) << (z0 + i*dz)     << " "
         << std::setprecision(3)
         << std::setw(8) << exp( trend[i] ) << " "
         << "0.00\n";
  }
  file.close();
}

//-------------------------------------------------------------------------------
void
Background::CalculateDeviationFromVerticalTrend(std::vector<std::vector<double> > & well_trend,
                                                const std::vector<double>         & global_trend,
                                                std::vector<double>               & avg_dev,
                                                const int                           nz)
{
  int n_wells = static_cast<int>(well_trend.size());

  for (int w = 0; w < n_wells; w++) {
    if (well_trend[w].size() > 0) {
      std::vector<double> & w_trend = well_trend[w];
      double sum_dev = 0.0f;
      int count = 0;
      for (int k = 0 ; k < nz ; k++) {
        if (w_trend[k] != RMISSING) {
          double diff = exp(w_trend[k]) - exp(global_trend[k]);
          sum_dev += diff*diff;
          count++;
        }
      }
      if (count > 0)
        sum_dev /= count;
      avg_dev[w] = sqrt(sum_dev);
    }
    else
      avg_dev[w] = RMISSING;
  }
}

//-------------------------------------------------------------------------------
void
Background::WriteDeviationsFromVerticalTrend(const std::vector<double>                        & avg_dev_vp,
                                             const std::vector<double>                        & avg_dev_vs,
                                             const std::vector<double>                        & avg_dev_rho,
                                             const std::vector<double>                        & trend_vp,
                                             const std::vector<double>                        & trend_vs,
                                             const std::vector<double>                        & trend_rho,
                                             const std::map<std::string, BlockedLogsCommon *> & blocked_logs,
                                             const int                                          n_wells,
                                             const int                                          nz)
{
  double global_mean_vp  = 0.0;
  double global_mean_vs  = 0.0;
  double global_mean_rho = 0.0;

  for (int k = 0; k < nz; k++) {
    global_mean_vp  += exp(trend_vp[k]);
    global_mean_vs  += exp(trend_vs[k]);
    global_mean_rho += exp(trend_rho[k]);
  }
  global_mean_vp  /= nz;
  global_mean_vs  /= nz;
  global_mean_rho /= nz;

  //
  // Find the relative average deviations (mean of Vp,Vs and Rho deviations).
  //
  double * rel_avg_dev = new double[n_wells];
  for (int i = 0; i < n_wells; i++) {
    double rel_dev_vp  = avg_dev_vp[i]/global_mean_vp;
    double rel_dev_vs  = avg_dev_vs[i]/global_mean_vs;
    double rel_dev_rho = avg_dev_rho[i]/global_mean_rho;
    rel_avg_dev[i] = (rel_dev_vp + rel_dev_vs + rel_dev_rho)/3;
  }
  //
  // Sort deviations to find worst well.
  //

  std::vector<int> index(n_wells);
  std::vector<std::string> well_names(n_wells);

  int w = 0;
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);

    index[w]      = w;
    well_names[w] = iter->second->GetWellName();
    w++;
  }

  for (int i = 0; i < n_wells; i++) {
    for (int j = i; j < n_wells; j++) {
      if (rel_avg_dev[index[j]] > rel_avg_dev[index[i]]) {
        int tmp = index[i];
        index[i] = index[j];
        index[j] = tmp;
      }
    }
  }
  //
  // Print results
  //
  if (n_wells > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
    LogKit::LogFormatted(LogKit::Low,"Well                        Vp       Vs      Rho\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
  }
  for (int i = 0; i < n_wells; i++) {

    int ii = index[i];
    if (avg_dev_vp[ii] != RMISSING) {
      LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", well_names[ii].c_str(),
                           avg_dev_vp[ii], avg_dev_vs[ii], avg_dev_rho[ii]);
    }
  }

  if (n_wells == 1) {
    LogKit::LogFormatted(LogKit::High,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::High,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
  }
  delete [] rel_avg_dev;
}

//-------------------------------------------------------------------------------
void
Background::FillInVerticalTrend(FFTGrid                   * grid,
                                const std::vector<double> & trend)
{
  const int nzp = grid->getNzp();  // equals nx
  const int nyp = grid->getNyp();  // equals ny
  const int nxp = grid->getNxp();  // equals nz

  grid->createRealGrid();
  grid->setType(FFTGrid::PARAMETER);
  grid->setAccessMode(FFTGrid::WRITE);

  int rnxp = 2*(nxp/2 + 1);
  for (int k = 0; k < nzp; k++)
    for (int j = 0; j < nyp; j++)
      for (int i = 0; i < rnxp; i++)
        grid->setNextReal(static_cast<float>(trend[k]));

  grid->endAccess();
}

void
Background::FillInVerticalTrend(NRLib::Grid<float>        * grid,
                                const std::vector<double> & trend)
{
  int nk = static_cast<int>(grid->GetNK());
  int nj = static_cast<int>(grid->GetNJ());
  int ni = static_cast<int>(grid->GetNI());
  grid->Resize(ni, nj, nk);

  for (int k = 0; k < nk; k++)
    for (int j = 0; j < nj; j++)
      for (int i = 0; i < ni; i++)
        grid->SetValue(i,j,k, static_cast<float>(trend[k]));
}

//-------------------------------------------------------------------------------
void
Background::ResampleBackgroundModel(NRLib::Grid<float>  * & bg_vp,
                                    NRLib::Grid<float>  * & bg_vs,
                                    NRLib::Grid<float>  * & bg_rho,
                                    const Simbox        *   bg_simbox,
                                    const Simbox        *   simbox)
{
  NRLib::Grid<float> * res_bg_vp  = new NRLib::Grid<float>();
  NRLib::Grid<float> * res_bg_vs  = new NRLib::Grid<float>();
  NRLib::Grid<float> * res_bg_rho = new NRLib::Grid<float>();

  LogKit::LogFormatted(LogKit::Low,"\nResampling background model...\n");
  ResampleParameter(res_bg_vp,  bg_vp,  simbox, bg_simbox);
  ResampleParameter(res_bg_vs,  bg_vs,  simbox, bg_simbox);
  ResampleParameter(res_bg_rho, bg_rho, simbox, bg_simbox);

  delete bg_vp;
  delete bg_vs;
  delete bg_rho;

  bg_vp  = res_bg_vp;
  bg_vs  = res_bg_vs;
  bg_rho = res_bg_rho;
}

void
Background::ResampleParameter(NRLib::Grid<float> *& p_new, // Resample to
                              NRLib::Grid<float> *  p_old, // Resample from
                              const Simbox       *  simbox_new,
                              const Simbox       *  simbox_old)
{
  int nx = simbox_new->getnx();
  int ny = simbox_new->getny();
  int nz = simbox_new->getnz();
  //
  // Use same padding as for nonresampled cubes
  //
  // NBNB-PAL: These grids are unpadded, so all nxp, nyp, ... would probably
  //           better be replaced by nx, ny, ... to avoid confusion...
  //
  //int nxp = nx + (pOld->getNxp() - pOld->getNxp());
  //int nyp = ny + (pOld->getNyp() - pOld->getNyp());
  //int nzp = nz + (pOld->getNzp() - pOld->getNzp());

  //
  // Set up relation between old layer index and new layer index using
  //
  // k2 = dz1/dz2 * k1 + (z01 - z02)/dz2    (from dz2*k2 + z02 = dz1*k1 + z01)
  //
  double * a       = new double[nx*ny];
  double * b       = new double[nx*ny];
  bool   * missing = new bool[nx*ny];

  int ij = 0;
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      double dz_new = simbox_new->getdz(i,j);
      double dz_old = simbox_old->getdz(i,j);
      double z0_new = simbox_new->getTop(i,j);
      double z0_old = simbox_old->getTop(i,j);
      double zn_new = simbox_new->getBot(i,j);
      double zn_old = simbox_old->getBot(i,j);

      if (dz_new == RMISSING || dz_old == RMISSING || z0_new == RMISSING || z0_old == RMISSING || zn_new == RMISSING || zn_old == RMISSING) {
        missing[ij] = true;
      }
      else {
        missing[ij] = false;
        a[ij] = dz_new/dz_old;
        b[ij] = (z0_new - z0_old)/dz_old;
      }
      ij++;
    }
  }

  //
  // Resample parameter
  //
  p_new->Resize(nx, ny, nz);

  double * layer = new double[nx*ny];

  for (int k = 0; k < nz; k++) {
    //
    // Map a layer
    //
    int ij=0;
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        if (missing[ij] == false) {
          int k_old = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
          layer[ij] = p_old->GetValue(i, j, k_old);
        }
        else
          layer[ij] = RMISSING;
        ij++;
      }
    }

    //
    // Smooth the layer (equal weighting of all neighbouring cells)
    //
    double value = 0.0;
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        if (layer[j*nx + i] != RMISSING) {
          int n = 1;
          double sum = layer[j*nx + i];
          if (i > 1) {
            if (layer[j*nx + i - 1] != RMISSING) {
              sum += layer[j*nx + i - 1];
              n++;
            }
          }
          if (j > 1) {
            if (layer[(j - 1)*nx + i] != RMISSING) {
              sum += layer[(j - 1)*nx + i];
              n++;
            }
          }
          if (i > 1 && j > 1) {
            if (layer[(j - 1)*nx + i - 1] != RMISSING) {
              sum += layer[(j - 1)*nx + i - 1];
              n++;
            }
          }
          if (i < nx-1) {
            if (layer[j*nx + i + 1] != RMISSING) {
              sum += layer[j*nx + i + 1];
              n++;
            }
          }
          if (j < ny-1) {
            if (layer[(j + 1)*nx + i] != RMISSING) {
              sum += layer[(j + 1)*nx + i];
              n++;
            }
          }
          if (i < nx-1 && j < ny-1) {
            if (layer[(j + 1)*nx + i + 1] != RMISSING) {
              sum += layer[(j + 1)*nx + i + 1];
              n++;
            }
          }
          value = sum/static_cast<double>(n);
        }
        else
          value = RMISSING;

        p_new->SetValue(i, j, k, static_cast<float>(value));
      }
    }
  }

  delete [] layer;
  delete [] a;
  delete [] b;
  delete [] missing;
}

FFTGrid *
Background::CopyFFTGrid(FFTGrid   * orig_grid,
                        const bool  exp_trans,
                        const bool  file_grid)
{
  FFTGrid * new_grid;

  if (file_grid)
    new_grid = new FFTFileGrid(static_cast<FFTFileGrid *>(orig_grid), exp_trans);
  else
    new_grid = new FFTGrid(orig_grid, exp_trans);

  return (new_grid);
}

void
Background::CheckLogForOnlyMissing(const std::vector<double> & log,
                                   const std::string         & log_name,
                                   const std::string         & well_name)
{
  int count = 0;
  for (int i = 0; i < static_cast<int>(log.size()); i++) {
    if (log[i] != RMISSING)
      count++;
  }
  if (count == 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWARNING in Kriging of background model: ");
    LogKit::LogFormatted(LogKit::Low,"Only missing values found in " + log_name + " log in well " + well_name + ".\n");
  }

}

