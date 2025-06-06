/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

#include "src/cravaresult.h"
#include "src/multiintervalgrid.h"
#include "src/blockedlogscommon.h"
#include "src/commondata.h"
#include "src/seismicparametersholder.h"
#include "src/krigingdata3d.h"
#include "src/parameteroutput.h"
#include "src/wavelet1D.h"
#include "src/modelavodynamic.h"
#include "src/modelgeneral.h"
#include "src/timings.h"
#include "lib/timekit.hpp"

CravaResult::CravaResult():
cov_vp_(NULL),
cov_vs_(NULL),
cov_rho_(NULL),
cr_cov_vp_vs_(NULL),
cr_cov_vp_rho_(NULL),
cr_cov_vs_rho_(NULL),
post_vp_(NULL),
post_vs_(NULL),
post_rho_(NULL),
post_vp_kriged_(NULL),
post_vs_kriged_(NULL),
post_rho_kriged_(NULL),
background_vp_(NULL),
background_vs_(NULL),
background_rho_(NULL),
background_trend_vp_(NULL),
background_trend_vs_(NULL),
background_trend_rho_(NULL),
block_grid_(NULL),
facies_prob_undef_(NULL),
quality_grid_(NULL),
write_crava_(false),
n_intervals_(1)
{
}

CravaResult::~CravaResult()
{
  for (size_t i = 0; i < corr_T_.size(); i++) {
    if (corr_T_[i] != NULL) {
      fftw_free(corr_T_[i]);
    }
  }
  for (size_t i = 0; i < corr_T_filtered_.size(); i++) {
    if (corr_T_filtered_[i] != NULL)
      delete [] corr_T_filtered_[i];
  }

  if (cov_vp_ != NULL) {
    delete cov_vp_;
    cov_vp_ = NULL;
  }
  if (cov_vs_ != NULL) {
    delete cov_vs_;
    cov_vs_ = NULL;
  }
  if (cov_rho_ != NULL) {
    delete cov_rho_;
    cov_rho_ = NULL;
  }
  if (post_vp_ != NULL) {
    delete post_vp_;
    post_vp_ = NULL;
  }
  if (post_vs_ != NULL) {
    delete post_vs_;
    post_vs_ = NULL;
  }
  if (post_rho_ != NULL) {
    delete post_rho_;
    post_rho_ = NULL;
  }

  if (cr_cov_vp_vs_ != NULL) {
    delete cr_cov_vp_vs_;
    cr_cov_vp_vs_ = NULL;
  }
  if (cr_cov_vp_rho_ != NULL) {
    delete cr_cov_vp_rho_;
    cr_cov_vp_rho_ = NULL;
  }
  if (cr_cov_vs_rho_ != NULL) {
    delete cr_cov_vs_rho_;
    cr_cov_vs_rho_ = NULL;
  }

  if (post_vp_kriged_ != NULL) {
    delete post_vp_kriged_;
    post_vp_kriged_ = NULL;
  }
  if (post_vs_kriged_ != NULL) {
    delete post_vs_kriged_;
    post_vs_kriged_ = NULL;
  }
  if (post_rho_kriged_ != NULL) {
    delete post_rho_kriged_;
    post_rho_kriged_ = NULL;
  }

  for (size_t i = 0; i < synt_seismic_data_.size(); i++) {
    delete synt_seismic_data_[i];
  }

  if (background_vp_ != NULL) {
    delete background_vp_;
    background_vp_ = NULL;
  }
  if (background_vs_ != NULL) {
    delete background_vs_;
    background_vs_ = NULL;
  }
  if (background_rho_ != NULL) {
    delete background_rho_;
    background_rho_ = NULL;
  }

  if (block_grid_ != NULL) {
    delete block_grid_;
    block_grid_ = NULL;
  }

  /*
  for (size_t i = 0; i < background_vp_intervals_.size(); i++) {
    if (background_vp_intervals_[i] != NULL) {
      delete background_vp_intervals_[i];
      background_vp_intervals_[i] = NULL;
    }
  }
  for (size_t i = 0; i < background_vs_intervals_.size(); i++) {
    if (background_vs_intervals_[i] != NULL) {
      delete background_vs_intervals_[i];
      background_vs_intervals_[i] = NULL;
    }
  }
  for (size_t i = 0; i < background_rho_intervals_.size(); i++) {
    if (background_rho_intervals_[i] != NULL) {
      delete background_rho_intervals_[i];
      background_rho_intervals_[i] = NULL;
    }
  }
  */

  for (size_t i = 0; i < simulations_seed0_.size(); i++) {
    delete simulations_seed0_[i];
  }
  for (size_t i = 0; i < simulations_seed1_.size(); i++) {
    delete simulations_seed1_[i];
  }
  for (size_t i = 0; i < simulations_seed2_.size(); i++) {
    delete simulations_seed2_[i];
  }

  for (size_t i = 0; i < facies_prob_.size(); i++) {
    delete facies_prob_[i];
  }

  if (facies_prob_undef_!= NULL) {
    delete facies_prob_undef_;
    facies_prob_undef_ = NULL;
  }
  for (size_t i = 0; i < facies_prob_geo_.size(); i++) {
    delete facies_prob_geo_[i];
  }
  for (size_t i = 0; i < lh_cubes_.size(); i++) {
    delete lh_cubes_[i];
  }
  if (quality_grid_!= NULL) {
    delete quality_grid_;
    quality_grid_ = NULL;
  }

  //for (size_t i = 0; i < wavelets_.size(); i++) {
  //  if (wavelets_[i] != NULL) {
  //    delete wavelets_[i];
  //    wavelets_[i] = NULL;
  //  }
  //}

  if (blocked_logs_intervals_.size() > 0) {
    for (int i = 0; i < n_intervals_; i++) {
      for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_intervals_[i].begin(); it != blocked_logs_intervals_[i].end(); it++) {
        if (blocked_logs_intervals_[i].find(it->first)->second != NULL) {
          delete blocked_logs_intervals_[i].find(it->first)->second;
        }
      }
    }
  }

}

void CravaResult::CombineResults(ModelSettings                        * model_settings,
                                 CommonData                           * common_data,
                                 std::vector<SeismicParametersHolder> & seismic_parameters_intervals)
{
  //Combine interval grids to one stormgrid
  //Compaction: Use the finest dz between intervals

  //Final grids are stored as StormContGrids
  // Results in seismic_parameters are stored as fftgrid
  //CreateStormGrid() (inside CombineResult) creates a storm grid from fft grid, and deletes the fft_grid

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  MultiIntervalGrid * multi_interval_grid     = common_data->GetMultipleIntervalGrid();
  Simbox & output_simbox                      = common_data->GetOutputSimbox();
  n_intervals_                                = multi_interval_grid->GetNIntervals();

  //Rapport
  if (n_intervals_ > 1 || output_simbox.getnz() != multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
    LogKit::LogFormatted(LogKit::Low,"\nThe results are combined and resampled to the final output (visualization) grid with " + NRLib::ToString(output_simbox.getnz())
                                      + " layers before they are written to file.\n");
  }
  else {
    LogKit::LogFormatted(LogKit::Low,"\nThe results are written to file based on the output (visualization) grid with " + NRLib::ToString(output_simbox.getnz()) + " layers.\n");
  }

  std::vector<NRLib::Grid<float> *> dummy_grids;

  int nx           = output_simbox.getnx();
  int ny           = output_simbox.getny();
  int nz_output    = output_simbox.getnz();

  blocked_logs_ = common_data->GetBlockedLogsOutput(); //Logs blocked to output_simbox

  CombineBlockedLogs(blocked_logs_, blocked_logs_intervals_, multi_interval_grid, &output_simbox); //Combine and resample logs that cannot be created directly to blocked_logs_output

  //Add seismic data to blocked logs

  if ((n_intervals_ > 1 || output_simbox.getnz() != multi_interval_grid->GetIntervalSimbox(0)->getnz()) && common_data->HasSeismicData()) {
    std::vector<SeismicStorage *> seismic_data = common_data->GetSeismicDataTimeLapse(0);
    int n_angles = static_cast<int>(seismic_data.size());
    ModelAVODynamic::AddSeismicLogsFromStorage(blocked_logs_,
                                               seismic_data,
                                               output_simbox,
                                               n_angles);
  }

  //If we want to write grids on CRAVA-format we do not delete the fft-grids in CombineResult so we can write the grids with padding in WriteResults
  if (n_intervals_ == 1 && ((model_settings->getOutputGridFormat() & IO::CRAVA) > 0))
    write_crava_ = true;

  //Create missing map
  NRLib::Grid2D<bool> * missing_map = CreateMissingGrid(output_simbox);

  LogKit::LogFormatted(LogKit::Low,"\nSetting up zone probabilities");
  std::vector<StormContGrid> zone_prob_grid(multi_interval_grid->GetNIntervals());
  for(size_t i=0; i<zone_prob_grid.size();i++)
    zone_prob_grid[i] = NRLib::StormContGrid(output_simbox, nx, ny, nz_output);
  multi_interval_grid->FindZoneProbGrid(zone_prob_grid);

  if (model_settings->getWritePrediction() && !model_settings->getEstimationMode()) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Prediction Grids");

    //Post vp, vs and rho from avoinversion computePostMeanResidAndFFTCov()
    post_vp_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    post_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    post_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    std::vector<FFTGrid *> post_vp_intervals(n_intervals_);
    std::vector<FFTGrid *> post_vs_intervals(n_intervals_);
    std::vector<FFTGrid *> post_rho_intervals(n_intervals_);

    // Erik N: Should these two cases be treated differently? I leave
    // the conditional statement
    if (model_settings->getForwardModeling()) {
      for (int i = 0; i < n_intervals_; i++) {
        post_vp_intervals[i]  = seismic_parameters_intervals[i].GetMeanVp();
        post_vs_intervals[i]  = seismic_parameters_intervals[i].GetMeanVs();
        post_rho_intervals[i] = seismic_parameters_intervals[i].GetMeanRho();
      }
    }
    else {
      for (int i = 0; i < n_intervals_; i++) {
        post_vp_intervals[i]  = seismic_parameters_intervals[i].GetPostVp();
        post_vs_intervals[i]  = seismic_parameters_intervals[i].GetPostVs();
        post_rho_intervals[i] = seismic_parameters_intervals[i].GetPostRho();
      }
    }

    LogKit::LogFormatted(LogKit::Low,"\n Vp ");
    CombineResult(post_vp_,  post_vp_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(post_vs_,  post_vs_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(post_rho_, post_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);

    //Add predicted logs from resampled grid
    if (n_intervals_ > 1 || output_simbox.getnz() != multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
      for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_.begin(); it != blocked_logs_.end(); it++) {
        std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_.find(it->first);

        BlockedLogsCommon * blocked_log = iter->second;
        blocked_log->SetLogFromGrid(post_vp_,0,1,"VP_PREDICTED");
        blocked_log->SetLogFromGrid(post_vs_,0,1,"VS_PREDICTED");
        blocked_log->SetLogFromGrid(post_rho_,0,1,"RHO_PREDICTED");
      }
    }

    //Post vp, vs and rho from avoinversion doPredictionKriging()
    if (model_settings->getKrigingParameter() > 0) {
      post_vp_kriged_  = new StormContGrid(output_simbox, nx, ny, nz_output);
      post_vs_kriged_  = new StormContGrid(output_simbox, nx, ny, nz_output);
      post_rho_kriged_ = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> post_vp_kriged_intervals(n_intervals_);
      std::vector<FFTGrid *> post_vs_kriged_intervals(n_intervals_);
      std::vector<FFTGrid *> post_rho_kriged_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        post_vp_kriged_intervals[i]  = seismic_parameters_intervals[i].GetPostVpKriged();
        post_vs_kriged_intervals[i]  = seismic_parameters_intervals[i].GetPostVsKriged();
        post_rho_kriged_intervals[i] = seismic_parameters_intervals[i].GetPostRhoKriged();
      }
      LogKit::LogFormatted(LogKit::Low,"\n Vp kriged ");
      CombineResult(post_vp_kriged_,  post_vp_kriged_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
      LogKit::LogFormatted(LogKit::Low,"\n Vs kriged ");
      CombineResult(post_vs_kriged_,  post_vs_kriged_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
      LogKit::LogFormatted(LogKit::Low,"\n Rho kriged ");
      CombineResult(post_rho_kriged_, post_rho_kriged_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    }
  }

  //Background models
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Background Model");
    background_vp_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    LogKit::LogFormatted(LogKit::Low,"\n Vp ");
    CombineResult(background_vp_,  background_vp_intervals_,  multi_interval_grid, zone_prob_grid, dummy_grids,
                  missing_map, model_settings->getFilterMultizoneModel(), model_settings->getMaxHzBackground());
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(background_vs_,  background_vs_intervals_,  multi_interval_grid, zone_prob_grid, dummy_grids,
                  missing_map, model_settings->getFilterMultizoneModel(), model_settings->getMaxHzBackground());
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(background_rho_, background_rho_intervals_, multi_interval_grid, zone_prob_grid, dummy_grids,
                  missing_map, model_settings->getFilterMultizoneModel(), model_settings->getMaxHzBackground());
  }
  else { //These background grids are not used later if we are not going to write them to file
    for (size_t i = 0; i < background_vp_intervals_.size(); i++) {
      if (background_vp_intervals_[i] != NULL)
        delete background_vp_intervals_[i];
      if (background_vs_intervals_[i] != NULL)
        delete background_vs_intervals_[i];
      if (background_rho_intervals_[i] != NULL)
        delete background_rho_intervals_[i];
    }
  }

  //Background trends
  if((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0) {
    const NRLib::Grid2D<std::vector<double> > & vertical_trends = common_data->GetBackgroundVerticalTrends();
    if (vertical_trends.GetNI() > 0) {

      background_trend_vp_  = new StormContGrid(output_simbox, nx, ny, nz_output);
      background_trend_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
      background_trend_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);

      CombineVerticalTrends(multi_interval_grid,
                            common_data,
                            vertical_trends,
                            zone_prob_grid,
                            background_trend_vp_,
                            background_trend_vs_,
                            background_trend_rho_,
                            missing_map);
    }
  }
  //Covariance grids
  if (!model_settings->getForwardModeling() && (model_settings->getOutputGridsOther() & IO::CORRELATION)) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Covariance grids");
    cov_vp_        = new StormContGrid(output_simbox, nx, ny, nz_output);
    cov_vs_        = new StormContGrid(output_simbox, nx, ny, nz_output);
    cov_rho_       = new StormContGrid(output_simbox, nx, ny, nz_output);
    cr_cov_vp_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    cr_cov_vp_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);
    cr_cov_vs_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);
    std::vector<FFTGrid *> cov_vp_intervals(n_intervals_);
    std::vector<FFTGrid *> cov_vs_intervals(n_intervals_);
    std::vector<FFTGrid *> cov_rho_intervals(n_intervals_);
    std::vector<FFTGrid *> cr_cov_vp_vs_intervals(n_intervals_);
    std::vector<FFTGrid *> cr_cov_vp_rho_intervals(n_intervals_);
    std::vector<FFTGrid *> cr_cov_vs_rho_intervals(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      //if (!model_settings->getForwardModeling())
      //  seismic_parameters_intervals[i].invFFTCovGrids();
      cov_vp_intervals[i]        = seismic_parameters_intervals[i].GetCovVp();
      cov_vs_intervals[i]        = seismic_parameters_intervals[i].GetCovVs();
      cov_rho_intervals[i]       = seismic_parameters_intervals[i].GetCovRho();
      cr_cov_vp_vs_intervals[i]  = seismic_parameters_intervals[i].GetCrCovVpVs();
      cr_cov_vp_rho_intervals[i] = seismic_parameters_intervals[i].GetCrCovVpRho();
      cr_cov_vs_rho_intervals[i] = seismic_parameters_intervals[i].GetCrCovVsRho();
      //if (!model_settings->getForwardModeling())
      //  seismic_parameters_intervals[i].FFTCovGrids();
    }
    LogKit::LogFormatted(LogKit::Low,"\n Vp ");
    CombineResult(cov_vp_,        cov_vp_intervals,        multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(cov_vs_,        cov_vs_intervals,        multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(cov_rho_,       cov_rho_intervals,       multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n VpVs ");
    CombineResult(cr_cov_vp_vs_,  cr_cov_vp_vs_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n VpRho ");
    CombineResult(cr_cov_vp_rho_, cr_cov_vp_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    LogKit::LogFormatted(LogKit::Low,"\n VsRho ");
    CombineResult(cr_cov_vs_rho_, cr_cov_vs_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
  }

  //Facies prob
  if (model_settings->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Facies prob:");
    int n_facies = static_cast<int>(seismic_parameters_intervals[0].GetFaciesProb().size());
    facies_prob_.resize(n_facies);

    for (int j = 0; j < n_facies; j++) {
      facies_prob_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> facies_prob_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProb()[j];
      }
      LogKit::LogFormatted(LogKit::Low,"\n " + NRLib::ToString(j) + " ");
      CombineResult(facies_prob_[j],  facies_prob_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    }

    //Set facies prob in blocked logs
    if (n_intervals_ > 1 || output_simbox.getnz() != multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
      for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_.begin(); it != blocked_logs_.end(); it++) {
        std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_.find(it->first);
        BlockedLogsCommon * blocked_log = iter->second;

        for (int j = 0; j < n_facies; j++) {
          blocked_log->SetLogFromGrid(facies_prob_[j], j , n_facies, "FACIES_PROB");
        }
      }
    }

    //Undef
    facies_prob_undef_ = new StormContGrid(output_simbox, nx, ny, nz_output);
    std::vector<FFTGrid *> facies_prob_intervals_undef(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      facies_prob_intervals_undef[i] = seismic_parameters_intervals[i].GetFaciesProbUndefined();
    }
    LogKit::LogFormatted(LogKit::Low,"\n Undef ");
    CombineResult(facies_prob_undef_,  facies_prob_intervals_undef,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);

  }
  if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Facies prob:");
    int n_facies = static_cast<int>(seismic_parameters_intervals[0].GetFaciesProbGeomodel().size());
    facies_prob_geo_.resize(n_facies);
    for (int j = 0; j < n_facies; j++) {
      facies_prob_geo_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> facies_prob_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProbGeomodel()[j];
      }
      LogKit::LogFormatted(LogKit::Low,"\n " + NRLib::ToString(j) + " ");
      CombineResult(facies_prob_geo_[j],  facies_prob_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    }

    //Set facies prob in blocked logs
    if (n_intervals_ > 1 || output_simbox.getnz() != multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
      for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_.begin(); it != blocked_logs_.end(); it++) {
        std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_.find(it->first);
        BlockedLogsCommon * blocked_log = iter->second;

        for (int j = 0; j < n_facies; j++) {
          blocked_log->SetLogFromGrid(facies_prob_geo_[j], j , n_facies, "FACIES_PROB");
        }
      }
    }

  }

  //LH Cube
  if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nFacies Likelihood:");
    int n_facies = static_cast<int>(seismic_parameters_intervals[0].GetLHCube().size());
    lh_cubes_.resize(n_facies);

    for (int j = 0; j < n_facies; j++) {
      lh_cubes_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> lh_cubes_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        lh_cubes_intervals[i] = seismic_parameters_intervals[i].GetLHCube()[j];
      }
      LogKit::LogFormatted(LogKit::Low,"\n " + NRLib::ToString(j) + " ");
      CombineResult(lh_cubes_[j], lh_cubes_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    }
  }

  //Quality grid
  if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID) {
    LogKit::LogFormatted(LogKit::Low,"\nSeismic Quality Grid: ");
    quality_grid_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    std::vector<FFTGrid *> quality_grid_intervals(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      quality_grid_intervals[i] = seismic_parameters_intervals[i].GetQualityGrid();
    }
    CombineResult(quality_grid_, quality_grid_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
  }

  //Simulation grids
  if (model_settings->getNumberOfSimulations() > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nSimulations Grids: ");
    int n_simulations = static_cast<int>(seismic_parameters_intervals[0].GetSimulationsSeed0().size());
    simulations_seed0_.resize(n_simulations);
    simulations_seed1_.resize(n_simulations);
    simulations_seed2_.resize(n_simulations);

    for (int j = 0; j < n_simulations; j++) {
      simulations_seed0_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);
      simulations_seed1_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);
      simulations_seed2_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> simulations_seed0_intervals(n_intervals_);
      std::vector<FFTGrid *> simulations_seed1_intervals(n_intervals_);
      std::vector<FFTGrid *> simulations_seed2_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        simulations_seed0_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed0(j);
        simulations_seed1_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed1(j);
        simulations_seed2_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed2(j);
      }
      LogKit::LogFormatted(LogKit::Low,"\n seed0 " + NRLib::ToString(j) + " ");
      CombineResult(simulations_seed0_[j], simulations_seed0_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
      LogKit::LogFormatted(LogKit::Low,"\n seed1 " + NRLib::ToString(j) + " ");
      CombineResult(simulations_seed1_[j], simulations_seed1_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
      LogKit::LogFormatted(LogKit::Low,"\n seed2 " + NRLib::ToString(j) + " ");
      CombineResult(simulations_seed2_[j], simulations_seed2_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);

    }
  }

  //Block grid
  if (model_settings->getWritePrediction() && model_settings->getDebugFlag()) {
    LogKit::LogFormatted(LogKit::Low,"\nBlock Grid ");
    block_grid_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    std::vector<FFTGrid *> block_grid_intervals(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      block_grid_intervals[i] = seismic_parameters_intervals[i].GetBlockGrid();
    }
    CombineResult(block_grid_, block_grid_intervals, multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
  }

  //Correlations and post variances
  //Do not combine, store and write per interval
  for (int i = 0; i < n_intervals_; i++) {
    corr_T_.push_back(seismic_parameters_intervals[i].GetCorrT());
    if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      corr_T_filtered_.push_back(seismic_parameters_intervals[i].GetCorrTFiltered());

    post_var0_.push_back(seismic_parameters_intervals[i].GetPostVar0());
    post_cov_vp00_.push_back(seismic_parameters_intervals[i].GetPostCovVp00());
    post_cov_vs00_.push_back(seismic_parameters_intervals[i].GetPostCovVs00());
    post_cov_rho00_.push_back(seismic_parameters_intervals[i].GetPostCovRho00());
  }

  //ComputeSeismicImpedance from avoinversion.cpp
  reflection_matrix_ = common_data->GetReflectionMatrixTimeLapse(0);
  wavelets_          = common_data->GetWavelet(0);

  //Resample wavelet to output simbox. These wavelets are from CommonData, and should not be on FFTFormat (and should not have added gainfactor).
  int n_wavelets = static_cast<int>(wavelets_.size());
  for (int i = 0; i < n_wavelets; i++) {
    wavelets_[i]->resample(static_cast<float>(output_simbox.getdz()),
                           output_simbox.getnz(),
                           output_simbox.GetNZpad());

    if (common_data->GetShiftGrid(0, i) != NULL)
      wavelets_[i]->setShiftGrid(new Grid2D(*common_data->GetShiftGrid(0, i)));
    if (common_data->GetGainGrid(0, i) != NULL) {
      wavelets_[i]->setGainGrid(new Grid2D(*common_data->GetGainGrid(0, i)));
      wavelets_[i]->invFFT1DInPlace();
    }
  }

  if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {
    LogKit::LogFormatted(LogKit::Low,"\nCompute Synthetic Seismic and Residuals...");
    if (model_settings->getKrigingParameter() > 0)
      ComputeSyntSeismic(model_settings, &output_simbox, wavelets_, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_, synt_seismic_data_);
    else
      ComputeSyntSeismic(model_settings, &output_simbox, wavelets_, post_vp_, post_vs_, post_rho_, synt_seismic_data_);
    LogKit::LogFormatted(LogKit::Low,"ok");
  }

  //From CommonData::WaveletHandling
  LogKit::LogFormatted(LogKit::Low,"\nGenerate Synthetic Seismic Logs...");
  CommonData::GenerateSyntheticSeismicLogs(wavelets_, blocked_logs_, reflection_matrix_, &output_simbox);

  //We need synt_seis from well wavelets, but it needs to be based on output simbox
  //Estimate a temp wavelet, which adds well_synt_seismic_data to blocked logs
  bool wavelet_estimated = false;
  for (int i = 0; i < model_settings->getNumberOfAngles(0); i++) {
    if (model_settings->getEstimateWavelet(0)[i] == true)
      wavelet_estimated = true;
  }
  if (wavelet_estimated == true)
    GenerateWellOptSyntSeis(model_settings, common_data, blocked_logs_, output_simbox, reflection_matrix_);

  LogKit::LogFormatted(LogKit::Low,"ok");

  //Trend Cubes
  if (model_settings->getOutputGridsOther() && IO::TREND_CUBES > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nTrend Cubes: ");

    int n_parameters = static_cast<int>(model_settings->getTrendCubeParameters().size());
    trend_cubes_.resize(n_parameters);

    for (int i = 0; i < n_parameters; i++) {
      trend_cubes_[i] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> trend_cubes_intervals(n_intervals_);

      for (int j = 0; j < n_intervals_; j++) {

        NRLib::Grid<float> * trend_cube = common_data->GetTrendCube(j).GetTrendCube(i);
        trend_cubes_intervals[j] = new FFTGrid(common_data->GetTrendCube(j).GetTrendCube(i), static_cast<int>(trend_cube->GetNI()), static_cast<int>(trend_cube->GetNJ()), static_cast<int>(trend_cube->GetNK()));

        if (trend_cube != NULL)
          delete trend_cube;
      }
      LogKit::LogFormatted(LogKit::Low,"\n " + NRLib::ToString(i) + ": ");
      CombineResult(trend_cubes_[i],  trend_cubes_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids, missing_map);
    }
    LogKit::LogFormatted(LogKit::Low,"ok");
  }

  //Delete grids from seismicparamtersholder
  if (model_settings->getEstimationMode() || !model_settings->getForwardModeling()) {
    for (int i = 0; i < n_intervals_; i++)
      seismic_parameters_intervals[i].releaseExpGrids();
  }

  Timings::setTimeCombineResults(wall,cpu);
}

void CravaResult::CombineResult(StormContGrid                    *& final_grid,
                                std::vector<FFTGrid *>            & interval_grids,
                                MultiIntervalGrid                 * multi_interval_grid,
                                const std::vector<StormContGrid>  & zone_probability,
                                std::vector<NRLib::Grid<float> *> & interval_grids_nrlib, //Optional, send in an empty vector if FFTGrids are used
                                NRLib::Grid2D<bool>               * missing_map, //NULL if no missing.
                                bool                                apply_filter, //filter combined grid to a maxHz
                                float                               max_hz)
{
  int nx = static_cast<int>(final_grid->GetNI());
  int ny = static_cast<int>(final_grid->GetNJ());
  int nz = static_cast<int>(final_grid->GetNK());

  //If output simbox has the same size as the result grid there is no need to resample
  if (n_intervals_ == 1 && static_cast<int>(final_grid->GetNK()) == multi_interval_grid->GetIntervalSimbox(0)->getnz()) {

    if (interval_grids_nrlib.size() > 0) {
      FFTGrid * tmp_grid = new FFTGrid(interval_grids_nrlib[0], nx, ny, nz);
      CreateStormGrid(*final_grid, tmp_grid);
    }
    else
      CreateStormGrid(*final_grid, interval_grids[0]);

    if (missing_map != NULL)
      SetMissingInGrid(*final_grid, missing_map);

    LogKit::LogFormatted(LogKit::Low,"ok");
    return;
  }

  bool use_nrlib_grids = false;
  if (interval_grids_nrlib.size() > 0)
    use_nrlib_grids = true;

  float monitor_size = std::max(1.0f, static_cast<float>(nx*ny)*0.02f);
  float next_monitor = monitor_size;
  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |");
  printf("\n  ^");

  std::vector<int> nz_old(n_intervals_);
  for(size_t zone=0;zone<nz_old.size();zone++) {
    if (use_nrlib_grids == false)
      nz_old[zone] = interval_grids[zone]->getNzp();
    else
      nz_old[zone] = CommonData::FindClosestFactorableNumber(multi_interval_grid->GetIntervalSimbox(static_cast<int>(zone))->getnz()+100);
  }

  std::vector<rfftwnd_plan> small_plans;
  std::vector<rfftwnd_plan> big_plans;
  int scale = 10; //How densely to sample "fine" values.
  CreateDownscalingPlans(nz_old,
                         scale,
                         small_plans,
                         big_plans);

  //Resample
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      //Combine vectors for each interval to one trace in stormgrid
      std::vector<float> combined_trace(nz);

      if (missing_map != NULL && (*missing_map)(i,j) == true) {
        for (int k = 0; k < nz; k++)
          combined_trace[k] = RMISSING;
      }
      else {
        std::vector<std::vector<float> > new_traces(n_intervals_); //Finely interpolated values

        //Resample each trace to new nz
        for (int zone = 0; zone < n_intervals_; zone++) {
          std::vector<float> old_trace;
          if (use_nrlib_grids == false) {
            FFTGrid * interval_grid  = interval_grids[zone];
            old_trace = interval_grid->getRealTrace(i, j); //old_trace is changed below.
          }
          else {
            old_trace = GetNRLibGridTrace(interval_grids_nrlib[zone], i, j);
          }
          int prepad_size = static_cast<int>(old_trace.size());
          AddPadding(old_trace, nz_old[zone]);

          new_traces[zone].resize(nz_old[zone]*scale);

          DownscaleTrace(old_trace,
                         new_traces[zone],
                         scale,
                         prepad_size,
                         small_plans[zone],
                         big_plans[zone]);

        } //n_intervals

        for (int k = 0; k < nz; k++) {
          double global_x = 0.0;
          double global_y = 0.0;
          double global_z = 0.0;

          double value = 0;
          final_grid->FindCenterOfCell(i, j, k, global_x, global_y, global_z);
          for (int zone = 0; zone < n_intervals_; zone++) {
            if(zone_probability[zone](i,j,k) > 0) {
              Simbox * z_simbox = multi_interval_grid->GetIntervalSimbox(zone);
              double dummy1, dummy2, rel_index;
              z_simbox->getInterpolationIndexes(global_x, global_y, global_z, dummy1, dummy2, rel_index);
              rel_index -= 0.5; //First half grid cell is outside interpolation vector.
              rel_index /= static_cast<double>(z_simbox->getnz()-1);
              if(rel_index < 0)
                rel_index = 0;
              else if(rel_index > 1)
                rel_index = 1;

              int index = static_cast<int>(floor(0.5+rel_index*(new_traces[zone].size()-1))); //0 to first item, 1 to last item.
              value += zone_probability[zone](i,j,k)*new_traces[zone][index];
            }
          }
          combined_trace[k] = static_cast<float>(value);
        }

        //Filter background grids
        if (apply_filter == true) {

          std::vector<float> filtered_trace(nz);
          double lz = final_grid->GetLZ();
          double dz = lz/nz;
          CommonData::ApplyFilter(filtered_trace,
                                  combined_trace,
                                  nz,
                                  dz,
                                  max_hz);

          for (int k = 0; k < nz; k++) {
            combined_trace[k] = filtered_trace[k];
          }
        }
      }

      for (int k = 0; k < nz; k++) {
        (*final_grid)(i,j,k) = combined_trace[k];
      }

      // Log progress
      if (ny*i + j + 1 >= static_cast<int>(next_monitor)) {
        next_monitor += monitor_size;
        printf("^");
        fflush(stdout);
      }

    } //ny
  } //nx

  if (use_nrlib_grids == false) {
    for (size_t i = 0; i < interval_grids.size(); i++) {
      delete interval_grids[i];
      interval_grids[i] = NULL;
    }
  }

  for(size_t zone=0;zone<small_plans.size();zone++) {
    fftwnd_destroy_plan(small_plans[zone]);
    fftwnd_destroy_plan(big_plans[zone]);
  }
}


void
CravaResult::CombineVerticalTrends(MultiIntervalGrid                         * multiple_interval_grid,
                                   CommonData                                * common_data,
                                   const NRLib::Grid2D<std::vector<double> > & vertical_trend,
                                   const std::vector<StormContGrid>          & zone_probability,
                                   StormContGrid                            *& background_trend_vp,
                                   StormContGrid                            *& background_trend_vs,
                                   StormContGrid                            *& background_trend_rho,
                                   NRLib::Grid2D<bool>                       * missing_map)
{
  int nx = multiple_interval_grid->GetIntervalSimbox(0)->getnx();
  int ny = multiple_interval_grid->GetIntervalSimbox(0)->getny();
  int n_intervals = multiple_interval_grid->GetNIntervals();
  std::vector<FFTGrid *> dummy_grids;
  LogKit::LogFormatted(LogKit::Low,"\nCombine Background Trend Grids");

  LogKit::LogFormatted(LogKit::Low,"\n Vp ");
  std::vector<NRLib::Grid<float> *> bg_trend_vp_intervals(n_intervals);
  for (int i = 0; i < n_intervals; i++) {

    Simbox * interval_simbox = multiple_interval_grid->GetIntervalSimbox(i);
    const Simbox * bg_simbox = common_data->GetBgSimbox(i);

    int nz = interval_simbox->getnz();
    int nz_extended = 0;
    if (bg_simbox != NULL)
      nz_extended = bg_simbox->getnz();

    NRLib::Grid<float> * trend_grid = new NRLib::Grid<float>(nx, ny, nz);
    if (nz_extended > 0 && nz_extended != nz) { //Vertical trend is estimated in an extended simbox, we resample it to interval simbox similar to the resampling of background grid.

      NRLib::Grid<float> * trend_grid_extended = new NRLib::Grid<float>(nx, ny, nz_extended);
      Background::FillInVerticalTrend(trend_grid_extended, vertical_trend(i,0));

      std::string err_text = "";
      Background::ResampleParameter(trend_grid, trend_grid_extended, interval_simbox, bg_simbox, err_text);

      if (err_text != "")
        LogKit::LogFormatted(LogKit::Error,"\n Problem combining vertical trends\n: " + err_text);


      delete trend_grid_extended;

    }
    else {
      Background::FillInVerticalTrend(trend_grid, vertical_trend(i,0));
    }
    bg_trend_vp_intervals[i] = trend_grid;
  }
  CombineResult(background_trend_vp, dummy_grids, multiple_interval_grid, zone_probability, bg_trend_vp_intervals, missing_map);

  LogKit::LogFormatted(LogKit::Low,"\n Vs ");
  std::vector<NRLib::Grid<float> *> bg_trend_vs_intervals(n_intervals);
  for (int i = 0; i < n_intervals; i++) {

    Simbox * interval_simbox = multiple_interval_grid->GetIntervalSimbox(i);
    const Simbox * bg_simbox = common_data->GetBgSimbox(i);

    int nz = interval_simbox->getnz();
    int nz_extended = 0;
    if (bg_simbox != NULL)
      nz_extended = bg_simbox->getnz();

    NRLib::Grid<float> * trend_grid = new NRLib::Grid<float>(nx, ny, nz);
    if (nz_extended > nz) { //Vertical trend is estimated in an extended simbox, we resample it to interval simbox similar to the resampling of background grid.

      NRLib::Grid<float> * trend_grid_extended = new NRLib::Grid<float>(nx, ny, nz_extended);
      Background::FillInVerticalTrend(trend_grid_extended, vertical_trend(i,1));

      std::string err_text = "";
      Background::ResampleParameter(trend_grid, trend_grid_extended, interval_simbox, bg_simbox, err_text);

      if (err_text != "")
        LogKit::LogFormatted(LogKit::Error,"\n Problem combining vertical trends\n: " + err_text);

      delete trend_grid_extended;

    }
    else {
      Background::FillInVerticalTrend(trend_grid, vertical_trend(i,1));
    }
    bg_trend_vs_intervals[i] = trend_grid;
  }
  CombineResult(background_trend_vs, dummy_grids, multiple_interval_grid, zone_probability, bg_trend_vs_intervals, missing_map);

  LogKit::LogFormatted(LogKit::Low,"\n Rho ");
  std::vector<NRLib::Grid<float> *> bg_trend_rho_intervals(n_intervals);
  for (int i = 0; i < n_intervals; i++) {

    Simbox * interval_simbox = multiple_interval_grid->GetIntervalSimbox(i);
    const Simbox * bg_simbox = common_data->GetBgSimbox(i);

    int nz = interval_simbox->getnz();
    int nz_extended = 0;
    if (bg_simbox != NULL)
      nz_extended = bg_simbox->getnz();

    NRLib::Grid<float> * trend_grid = new NRLib::Grid<float>(nx, ny, nz);
    if (nz_extended > nz) { //Vertical trend is estimated in an extended simbox, we resample it to interval simbox similar to the resampling of background grid.

      NRLib::Grid<float> * trend_grid_extended = new NRLib::Grid<float>(nx, ny, nz_extended);
      Background::FillInVerticalTrend(trend_grid_extended, vertical_trend(i,2));

      std::string err_text = "";
      Background::ResampleParameter(trend_grid, trend_grid_extended, interval_simbox, bg_simbox, err_text);

      if (err_text != "")
        LogKit::LogFormatted(LogKit::Error,"\n Problem combining vertical trends\n: " + err_text);
      delete trend_grid_extended;

    }
    else {
      Background::FillInVerticalTrend(trend_grid, vertical_trend(i,2));
    }
    bg_trend_rho_intervals[i] = trend_grid;
  }
  CombineResult(background_trend_rho, dummy_grids, multiple_interval_grid, zone_probability, bg_trend_rho_intervals, missing_map);

}


std::vector<float>
CravaResult::GetNRLibGridTrace(NRLib::Grid<float> * grid,
                               int i,
                               int j)
{
  std::vector<float> trace;

  for (int k = 0; k < static_cast<int>(grid->GetNK()); k++) {
    trace.push_back(grid->GetValue(i,j,k));
  }

  return trace;
}

float CravaResult::GetResampledTraceValue(const std::vector<float> & resampled_trace, //interval
                                          const double             & dz_resampled,
                                          const double             & top,
                                          const double             & global_z, //center of cell
                                          const double             & dz_final)
{
  int nz_resampled    = static_cast<int>(resampled_trace.size());
  double global_z_top = global_z - 0.5*dz_final; //Use top of cell

  //Start on top and pick the closest value in z
  double trace_z = top;
  int index      = 0;
  while (index < (nz_resampled-1)) {
    if (trace_z >= global_z_top && trace_z <= (global_z_top + dz_resampled))
      break;
    else {
      trace_z += dz_resampled;
      index++;
    }
  }

  float value = resampled_trace[index];
  return(value);
}

double CravaResult::GetResampledTraceValue(const std::vector<double> & resampled_trace,
                                           const std::vector<double> & z_pos_resampled,
                                           const double              & global_z) //z-value for this cell in the final blocked log
{
  int nz_resampled = static_cast<int>(resampled_trace.size());
  //Get trace value to global_z. We search from the top of this interval until we find the corresponding z-value

  int index = 0;
  if (global_z > z_pos_resampled[0]) { //z-value above this blocked log, we use the first element
    for (int i = 0; i < (nz_resampled-1); i++) {
      if (global_z >= z_pos_resampled[i] && global_z <= z_pos_resampled[i+1]) { //wanted z is inside this fine interval
        index = i;
        break;
      }
    }
  }

  double value = resampled_trace[index];
  return(value);
}

void CravaResult::CombineBlockedLogs(std::map<std::string, BlockedLogsCommon *>                     & blocked_logs_output, //blocked to ouput_simbox in CommonData
                                     const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                                     MultiIntervalGrid                                              * multi_interval_grid,
                                     const Simbox                                                   * output_simbox)
{
  //Resample blocked logs to output_simbox and combine if multiple intervals
  //blocked_logs_output are blocked to output_simbox in commondata
  //Need to resample blocked logs that are added after commondata, and cannot be added directly to blocked_logs_output

  float res_fac   = 10.0; //Degree of refinement, must be integer.

  //Do not resample if there is only one interval, and if output_simbox and has the same resolution as the interval_simbox
  if (n_intervals_ == 1 && output_simbox->getnz() == multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
    blocked_logs_output = blocked_logs_intervals[0];
    return;
  }


  //Loop over wells and combine
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_output.begin(); it != blocked_logs_output.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_output.find(it->first);

    BlockedLogsCommon * blocked_log_final = iter->second;

    std::string well_name = blocked_log_final->GetWellName();

    int n_blocks = blocked_log_final->GetNumberOfBlocks();

    bool got_vp_rho_fac_log = false;
    for (int i = 0; i < n_intervals_; i++) {
      if (blocked_logs_intervals[i].find(well_name) != blocked_logs_intervals[i].end()) {
        if (blocked_logs_intervals[i].find(well_name)->second->GetVpFaciesFiltered().size() > 0) {
          got_vp_rho_fac_log = true;
        }
      }
    }

    bool got_filtered_logs = false;
    for (int i = 0; i < n_intervals_; i++) {
      if (blocked_logs_intervals[i].find(well_name) != blocked_logs_intervals[i].end()) {
        if (blocked_logs_intervals[i].find(well_name)->second->GetContLogsSeismicRes().size() > 0)
          got_filtered_logs = true;
      }
    }

    if (got_vp_rho_fac_log || got_filtered_logs) {
      LogKit::LogFormatted(LogKit::Low,"\nCombine Blocked Logs for well " + well_name + "...");
    }

    //ForFacies logs (vp, rho)
    if (got_vp_rho_fac_log) {
      LogKit::LogFormatted(LogKit::Low,"\n    vp rho fac log... ");
      std::vector<double> vp_for_facies_final(n_blocks);
      std::vector<double> rho_for_facies_final(n_blocks);
      std::vector<std::vector<double> > vp_for_facies_intervals;
      std::vector<std::vector<double> > rho_for_facies_intervals;

      //Get well logs, missing values are interpolated
      std::vector<int> intervals_used; //Keep track of which intervals contain this well
      for (int i = 0; i < n_intervals_; i++) {

        if (blocked_logs_intervals[i].find(well_name) != blocked_logs_intervals[i].end()) {
          std::vector<double> tmp_log_vp;
          CopyWellLog(tmp_log_vp,  blocked_logs_intervals[i].find(well_name)->second->GetVpFaciesFiltered());
          vp_for_facies_intervals.push_back(tmp_log_vp);

          std::vector<double> tmp_log_rho;
          CopyWellLog(tmp_log_rho, blocked_logs_intervals[i].find(well_name)->second->GetRhoFaciesFiltered());
          rho_for_facies_intervals.push_back(tmp_log_rho);

          intervals_used.push_back(i);
        }
      }

      std::vector<int> interval_log(n_blocks);
      ResampleLog(vp_for_facies_final,  vp_for_facies_intervals,  blocked_logs_intervals, intervals_used, multi_interval_grid, blocked_log_final, well_name, res_fac, interval_log);
      ResampleLog(rho_for_facies_final, rho_for_facies_intervals, blocked_logs_intervals, intervals_used, multi_interval_grid, blocked_log_final, well_name, res_fac, interval_log);

      blocked_log_final->SetVpFaciesFiltered(vp_for_facies_final);
      blocked_log_final->SetRhoFaciesFiltered(rho_for_facies_final);
      if (n_intervals_ > 1)
        blocked_log_final->SetIntervalLog(interval_log);
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    //Filtered logs
    if (got_filtered_logs) {
      LogKit::LogFormatted(LogKit::Low,"\n    filtered logs... ");
      std::vector<double> vp_filtered_final(n_blocks);
      std::vector<double> vs_filtered_final(n_blocks);
      std::vector<double> rho_filtered_final(n_blocks);
      std::vector<std::vector<double> > vp_filtered_intervals;
      std::vector<std::vector<double> > vs_filtered_intervals;
      std::vector<std::vector<double> > rho_filtered_intervals;

      //Get well logs, missing values are interpolated
      std::vector<int> intervals_used;
      for (int i = 0; i < n_intervals_; i++) {

        if (blocked_logs_intervals[i].find(well_name) != blocked_logs_intervals[i].end()) {
          std::vector<double> tmp_log_vp;
          CopyWellLog(tmp_log_vp,  blocked_logs_intervals[i].find(well_name)->second->GetVpSeismicResolution());
          vp_filtered_intervals.push_back(tmp_log_vp);

          std::vector<double> tmp_log_vs;
          CopyWellLog(tmp_log_vs,  blocked_logs_intervals[i].find(well_name)->second->GetVsSeismicResolution());
          vs_filtered_intervals.push_back(tmp_log_vs);

          std::vector<double> tmp_log_rho;
          CopyWellLog(tmp_log_rho,  blocked_logs_intervals[i].find(well_name)->second->GetRhoSeismicResolution());
          rho_filtered_intervals.push_back(tmp_log_rho);

          intervals_used.push_back(i);
        }
      }

      std::vector<int> interval_log(n_blocks);
      ResampleLog(vp_filtered_final,  vp_filtered_intervals,  blocked_logs_intervals, intervals_used, multi_interval_grid, blocked_log_final, well_name, res_fac, interval_log);
      ResampleLog(vs_filtered_final,  vs_filtered_intervals,  blocked_logs_intervals, intervals_used, multi_interval_grid, blocked_log_final, well_name, res_fac, interval_log);
      ResampleLog(rho_filtered_final, rho_filtered_intervals, blocked_logs_intervals, intervals_used, multi_interval_grid, blocked_log_final, well_name, res_fac, interval_log);

      //If we vp and vp_bg are missing, seismic_resolution should also be missing. But missing values are interpolated over during the resampling, so we reset them as missing.
      for (int i = 0; i < n_blocks; i++) {
        if (blocked_log_final->GetVpBlocked()[i] == RMISSING && blocked_log_final->GetVpHighCutBackground()[i] == RMISSING)
          vp_filtered_final[i] = RMISSING;
        if (blocked_log_final->GetVsBlocked()[i] == RMISSING && blocked_log_final->GetVsHighCutBackground()[i] == RMISSING)
          vs_filtered_final[i] = RMISSING;
        if (blocked_log_final->GetRhoBlocked()[i] == RMISSING && blocked_log_final->GetRhoHighCutBackground()[i] == RMISSING)
          rho_filtered_final[i] = RMISSING;
      }

      blocked_log_final->SetVpSeismicResolution(vp_filtered_final);
      blocked_log_final->SetVsSeismicResolution(vs_filtered_final);
      blocked_log_final->SetRhoSeismicResolution(rho_filtered_final);
      if (n_intervals_ > 1)
        blocked_log_final->SetIntervalLog(interval_log);
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    blocked_logs_output.insert(std::pair<std::string, BlockedLogsCommon *>(well_name, blocked_log_final));

  }
}

void CravaResult::CopyWellLog(std::vector<double>       & log_new,
                              const std::vector<double> & log_old)
{
  for (size_t i = 0; i < log_old.size(); i++) {
    log_new.push_back(log_old[i]);
  }

  InterpolateMissing(log_new);
}

void CravaResult::InterpolateMissing(std::vector<double> & well_log)
{
  //Interpolate missing values

  if (well_log[0] == RMISSING || well_log[well_log.size()-1] == RMISSING)
    ExtrapolateLog(well_log);

  for (size_t i = 1; i < well_log.size()-1; i++) {
    if (well_log[i] == RMISSING) {
      int first_m = static_cast<int>(i);
      int j = static_cast<int>(i) + 1;
      while (well_log[j] == RMISSING) { //Find length of missing values
        j++;
      }
      int last_m = j;

      int start = first_m - 1; //value
      int end   = last_m + 1;  //value

      for (int k = first_m; k < last_m + 1; k++) {
        float rel = static_cast<float>(k-start)/static_cast<float>(end-start);
        well_log[k] = rel * well_log[start] + (1 - rel) * well_log[end];
      }
    }
  }
}

void CravaResult::ExtrapolateLog(std::vector<double> & well_log)
{
  if (well_log[0] == RMISSING) {
    int i = 1;
    while (well_log[i] == RMISSING)
      i++;
    for (int j = 0; j < i; j++)
      well_log[j] = well_log[i];
  }

  if (well_log[well_log.size() - 1] == RMISSING) {
    int i = static_cast<int>(well_log.size()) - 2;
    while (well_log[i] == RMISSING)
      i--;
    for (size_t j = i + 1; j < well_log.size(); j++)
      well_log[j] = well_log[i];
  }

}

void CravaResult::ResampleLog(std::vector<double>                                            & final_log,
                              std::vector<std::vector<double> >                              & old_log_interval, //vector intervals used
                              const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals, //vector over all intervals
                              std::vector<int>                                               & intervals_used,
                              MultiIntervalGrid                                              * multi_interval_grid,
                              const BlockedLogsCommon                                        * blocked_log_final,
                              std::string                                                      well_name,
                              float                                                            res_fac,
                              std::vector<int>                                               & interval_log)
{

  //We only resample from the intervals used, i.e. intervals that contains this well.
  int n_intervals_used = static_cast<int>(intervals_used.size());

  std::vector<std::vector<double> > interval_logs_fine(n_intervals_used); //resample to, per interval used
  std::vector<std::vector<double> > z_pos_resampled_fine(n_intervals_used); //Resample z-blocked log, it is used when we combine traces.

  //Get logs per interval and resample to a fine resolution
  for (int i = 0; i < n_intervals_used; i++) {
    int i_interval_global = intervals_used[i];

    int nz_interval           = static_cast<int>(old_log_interval[i].size() * res_fac);
    std::vector<double> z_pos = blocked_logs_intervals[i_interval_global].find(well_name)->second->GetZposBlocked();

    interval_logs_fine[i].resize(nz_interval);
    z_pos_resampled_fine[i].resize(nz_interval);

    ResampleTrace(old_log_interval[i], interval_logs_fine[i], res_fac); //Interpolate missing
    ResampleTrace(z_pos, z_pos_resampled_fine[i], res_fac);
  }

  CombineTraces(final_log, blocked_log_final, multi_interval_grid, interval_logs_fine, z_pos_resampled_fine, interval_log, intervals_used);
}

void CravaResult::ResampleTrace(std::vector<double> & old_trace, //not const, it is changed
                                std::vector<double> & new_trace,
                                const float           res_fac)
{

  int nz_old = static_cast<int>(old_trace.size());
  int nz_new = static_cast<int>(new_trace.size());

  //Remove trend from trace
  size_t n_trace     = old_trace.size();
  double trend_first = old_trace[0];
  double trend_last  = old_trace[n_trace - 1];
  double trend_inc   = (trend_last - trend_first) / (n_trace - 1);
  for (size_t k_trace = 0; k_trace < old_trace.size(); k_trace++) {
    old_trace[k_trace] -= trend_first + k_trace * trend_inc;
  }

  int nt = nz_old;
  int mt = nz_new;

  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  int cnt = nt/2 + 1;
  int rnt = 2*cnt;
  int cmt = mt/2 + 1;
  int rmt = 2*cmt;

  fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
  fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

  CommonData::ResampleTrace(old_trace,
                            fftplan1,
                            fftplan2,
                            rAmpData,
                            rAmpFine,
                            nt,
                            cnt,
                            rnt,
                            cmt,
                            rmt);

  //Add trend
  trend_inc = (trend_last - trend_first) / (res_fac * (nz_old - 1));
  for (int k = 0; k < nz_new; k++) {
    new_trace[k] = rAmpFine[k] + (trend_first + k*trend_inc);
  }

  fftw_free(rAmpData);
  fftw_free(rAmpFine);
  fftwnd_destroy_plan(fftplan1);
  fftwnd_destroy_plan(fftplan2);

}

void CravaResult::CombineTraces(std::vector<double>                     & final_log,
                                const BlockedLogsCommon                 * blocked_log_final,
                                MultiIntervalGrid                       * multi_interval_grid,
                                const std::vector<std::vector<double> > & interval_logs_fine, //Vector over intervals used
                                const std::vector<std::vector<double> > & z_pos_resampled,
                                std::vector<int>                        & interval_log,
                                std::vector<int>                        & intervals_used)
{
  int nz          = static_cast<int>(final_log.size());
  int n_intervals_used = static_cast<int>(interval_logs_fine.size());

  const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();

  for (int k = 0; k < nz; k++) {

    bool two_intervals = false;

    double global_x = blocked_log_final->GetXposBlocked()[k];
    double global_y = blocked_log_final->GetYposBlocked()[k];
    double global_z = blocked_log_final->GetZposBlocked()[k];

    //We only send in and resample the intervals that contain this well
    //i_interval is the index over all used intervals, and i_interval_global is the corresponding index of all intervals
    int i_interval = 0;
    int i_interval_global = 0;
    for (int i = 0; i < n_intervals_used; i++) {
      i_interval = i;
      i_interval_global = intervals_used[i];

      Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval_global);

      if (interval_simbox->IsInside(global_x, global_y, global_z))
        break;
    }
    if (i_interval < (n_intervals_used-1)) { //Also check if it hits the next interval, unless it is the last one.
      Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval_global+1);

      if (interval_simbox->IsInside(global_x, global_y, global_z))
        two_intervals = true;
    }

    double value        = 0.0;
    int interval_index  = i_interval;
    int interval_index_global = i_interval_global;

    if (two_intervals == true) {
      //Use erosion priorities to select between the two intervals
      if (erosion_priorities[i_interval_global] < erosion_priorities[i_interval_global+1]) {
        interval_index = i_interval;
        interval_index_global = i_interval_global;
      }
      else {
        interval_index = i_interval + 1;
        interval_index_global = i_interval_global + 1;

      }
    }

    value = GetResampledTraceValue(interval_logs_fine[interval_index], z_pos_resampled[interval_index], global_z);

    final_log[k] = value;

    interval_log[k] = interval_index_global;
  }
}


void CravaResult::WriteResults(ModelSettings           * model_settings,
                               CommonData              * common_data,
                               SeismicParametersHolder & seismic_parameters) //For crava-writing
{
  //Results are combined to one grid in CombineResults first
  //Wavelets are written out both in commonData and Wavelet1d/3d.cpp (and possibly modelavodynamic if match energies)
  //Estimation model: WriteEstimationResults

  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  const Simbox & simbox            = common_data->GetOutputSimbox();
  int output_grids_elastic         = model_settings->getOutputGridsElastic();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  //Logging different options for writing segy grids. Sets offset from topsurface if it is missing
  if ((model_settings->getOutputGridFormat() & IO::SEGY) > 0)
    LogAndSetSegyOffsetIfNeeded(model_settings, simbox);

  //Update depth_mapping with resampled post_vp
  if ((time_depth_mapping != NULL && time_depth_mapping->getSimbox() == NULL) || (common_data->GetVelocityFromInversion() && time_depth_mapping->getSimbox()->getnz() != simbox.getnz())) {
    LogKit::LogFormatted(LogKit::Low,"\nUsing Vp velocity field from the combined vp grid to map between time and depth grids.\n");

    ExpTransf(post_vp_);
    int format = model_settings->getOutputGridFormat();
    if (model_settings->getWriteAsciiSurfaces() && !(format & IO::ASCII))
      format += IO::ASCII;

    time_depth_mapping->setMappingFromVelocity(post_vp_, &simbox, format);
    LogTransf(post_vp_);
  }

  //Write blocked wells
  if ((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Blocked Logs...");
    WriteBlockedWells(blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

    //Write blocked background logs (CRA-544). Logs that are blocked to extended background model (extended simbox with correlation direction).
    //Do not write if multiple intervals is used
    if (n_intervals_ == 1 && bg_blocked_logs_.size() > 0)
      WriteBlockedWells(bg_blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

    LogKit::LogFormatted(LogKit::Low,"ok\n");
  }
  //if ((model_settings->getWellOutputFlag() & IO::WELLS) > 0) {
  //  WriteWells(common_data->GetWells(), model_settings);
  //}

  if (model_settings->getWritePrediction() && !model_settings->getForwardModeling() && !model_settings->getEstimationMode() && output_grids_elastic > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Prediction Grids\n");

    //From computePostMeanResidAndFFTCov()
    ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_, post_vs_, post_rho_,
                                      output_grids_elastic, -1, false);

    if (write_crava_) {
      std::string file_name_vp  = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Vp");
      std::string file_name_vs  = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Vs");
      std::string file_name_rho = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Rho");
      seismic_parameters.GetPostVp()->writeCravaFile(file_name_vp,   &simbox);
      seismic_parameters.GetPostVs()->writeCravaFile(file_name_vs,   &simbox);
      seismic_parameters.GetPostRho()->writeCravaFile(file_name_rho, &simbox);
    }

    //From doPredictionKriging
    if (model_settings->getKrigingParameter() > 0) {
      ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_,
                                        output_grids_elastic, -1, true);

      if (write_crava_) {
        std::string file_name_vp  = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Vp_Kriged");
        std::string file_name_vs  = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Vs_Kriged");
        std::string file_name_rho = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixPredictions() + "Rho_Kriged");
        seismic_parameters.GetPostVpKriged()->writeCravaFile(file_name_vp,   &simbox);
        seismic_parameters.GetPostVsKriged()->writeCravaFile(file_name_vs,   &simbox);
        seismic_parameters.GetPostRhoKriged()->writeCravaFile(file_name_rho, &simbox);
      }
    }

    //From CKrigingAdmin::KrigAll
    if (model_settings->getDebugFlag()) {
      ParameterOutput::WriteFile(model_settings, block_grid_, "BlockGrid", IO::PathToInversionResults(), &simbox);

      if (write_crava_) {
        std::string file_name = IO::makeFullFileName(IO::PathToInversionResults(), "BlockGrid");
        seismic_parameters.GetBlockGrid()->writeCravaFile(file_name, &simbox);
      }
    }
  }

  //Write seismic data. Resample from CommonData to output_simbox. If CRAVA-format, we resample with padding.
  if (model_settings->getForwardModeling() == false &&
      ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0
      || (model_settings->getOutputGridsSeismic() & IO::RESIDUAL) > 0 )) {

        WriteSeismicData(model_settings, common_data, simbox, time_depth_mapping);

  }

  //Write Background models
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Background Grids\n");
    WriteGridPackage(model_settings,
                     &simbox,
                     background_vp_,
                     background_vs_,
                     background_rho_,
                     time_depth_mapping,
                     IO::PrefixBackground(),
                     IO::PathToBackground(),
                     true);

    if (write_crava_) {
      std::string file_name_vp  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vp");
      std::string file_name_vs  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vs");
      std::string file_name_rho = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Rho");

      ExpTransf(background_vp_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vp, &simbox);

      ExpTransf(background_vs_intervals_[0]);
      background_vs_intervals_[0]->writeCravaFile(file_name_vs, &simbox);

      ExpTransf(background_rho_intervals_[0]);
      background_rho_intervals_[0]->writeCravaFile(file_name_rho, &simbox);
    }
  }

  //Write background vertical trend grids
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0) {
    const NRLib::Grid2D<std::vector<double> > & vertical_trends = common_data->GetBackgroundVerticalTrends();
    if (vertical_trends.GetNI() > 0) {
      LogKit::LogFormatted(LogKit::Low,"\nWrite Background Vertical Trend Grids\n");

      std::string prefix = IO::PrefixBackground()+IO::PrefixTrend();
      WriteGridPackage(model_settings,
                       &simbox,
                       background_trend_vp_,
                       background_trend_vs_,
                       background_trend_rho_,
                       time_depth_mapping,
                       prefix,
                       IO::PathToBackground(),
                       true);
    }
  }

  //Write correlations and post variances. Write per interval.
  if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Correlations\n");

    for (int i = 0; i < n_intervals_; i++) {
      std::string interval_name = common_data->GetMultipleIntervalGrid()->GetIntervalName(i);
      WriteFilePostVariances(post_var0_[i], post_cov_vp00_[i], post_cov_vs00_[i], post_cov_rho00_[i], interval_name);
    }

    std::string file_name_vp    = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp";
    std::string file_name_vs    = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs";
    std::string file_name_rho   = IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho";
    std::string file_name_vpvs  = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs";
    std::string file_name_vprho = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho";
    std::string file_name_vsrho = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho";

    ParameterOutput::WriteFile(model_settings, cov_vp_,        file_name_vp,    IO::PathToCorrelations(), &simbox, false, "Posterior covariance for Vp");
    ParameterOutput::WriteFile(model_settings, cov_vs_,        file_name_vs,    IO::PathToCorrelations(), &simbox, false, "Posterior covariance for Vs");
    ParameterOutput::WriteFile(model_settings, cov_rho_,       file_name_rho,   IO::PathToCorrelations(), &simbox, false, "Posterior covariance for density");
    ParameterOutput::WriteFile(model_settings, cr_cov_vp_vs_,  file_name_vpvs,  IO::PathToCorrelations(), &simbox, false, "Posterior cross-covariance for (Vp,Vs)");
    ParameterOutput::WriteFile(model_settings, cr_cov_vp_rho_, file_name_vprho, IO::PathToCorrelations(), &simbox, false, "Posterior cross-covariance for (Vp,density)");
    ParameterOutput::WriteFile(model_settings, cr_cov_vs_rho_, file_name_vsrho, IO::PathToCorrelations(), &simbox, false, "Posterior cross-covariance for (Vs,density)");

    if (write_crava_) {
      file_name_vp    = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp");
      file_name_vs    = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs");
      file_name_rho   = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho");
      file_name_vpvs  = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs");
      file_name_vprho = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho");
      file_name_vsrho = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho");

      seismic_parameters.GetCovVp()->writeCravaFile(file_name_vp,         &simbox);
      seismic_parameters.GetCovVs()->writeCravaFile(file_name_vs,         &simbox);
      seismic_parameters.GetCovRho()->writeCravaFile(file_name_rho,       &simbox);
      seismic_parameters.GetCrCovVpVs()->writeCravaFile(file_name_vpvs,   &simbox);
      seismic_parameters.GetCrCovVpRho()->writeCravaFile(file_name_vprho, &simbox);
      seismic_parameters.GetCrCovVsRho()->writeCravaFile(file_name_vsrho, &simbox);
    }
  }


  if (model_settings->getKrigingParameter() > 0 && model_settings->getWritePrediction()) {
    KrigingData3D kd(common_data->GetBlockedLogs(), 1); // 1 = full resolution logs
    std::string base_name = "Raw_" + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string file_name = IO::makeFullFileName(IO::PathToInversionResults(), base_name);
    kd.writeToFile(file_name);
  }

  //FaciesProb: prob, undefined, qualitygrid, lhcube
  if (model_settings->getEstimateFaciesProb()) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Facies Probability Grids\n");
    std::vector<std::string> facies_names = common_data->GetFaciesNames();
    int n_facies = static_cast<int>(facies_names.size());

    std::string base_name = IO::PrefixFaciesProbability();
    if (model_settings->getFaciesProbRelative()) {
      if (model_settings->getFaciesProbFromRockPhysics())
        base_name += "Rock_Physics_";
    }
    else if (model_settings->getFaciesProbRelative() == false) {
      base_name += "Absolute_";
      if (model_settings->getFaciesProbFromRockPhysics())
        base_name += "Rock_Physics_";
    }

    if (model_settings->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF) {
      for (int i = 0; i < n_facies; i++) {
        std::string file_name = base_name +"With_Undef_"+ facies_names[i];
        ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, facies_prob_[i], file_name, "");

        if (write_crava_) {
          std::string file_name_crava = IO::makeFullFileName(IO::PathToInversionResults(), file_name);
          seismic_parameters.GetFaciesProb()[i]->writeCravaFile(file_name_crava, &simbox);
        }
      }
      std::string file_name = base_name + "Undef";
      ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, facies_prob_undef_, file_name, "");

      if (write_crava_) {
        std::string file_name_crava = IO::makeFullFileName(IO::PathToInversionResults(), file_name);
        seismic_parameters.GetFaciesProbUndefined()->writeCravaFile(file_name_crava, &simbox);
      }

    }
    if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
      for (int i = 0; i < n_facies; i++) {
        std::string file_name = base_name + facies_names[i];
        ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, facies_prob_geo_[i], file_name, "");

        if (write_crava_) {
          std::string file_name_crava = IO::makeFullFileName(IO::PathToInversionResults(), file_name);
          seismic_parameters.GetFaciesProbGeomodel()[i]->writeCravaFile(file_name_crava, &simbox);
        }
      }
    }
    if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID) {
      std::string file_name = "Seismic_Quality_Grid";
      ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, quality_grid_, file_name, "");

      if (write_crava_) {
        std::string file_name_crava = IO::makeFullFileName(IO::PathToInversionResults(), file_name);
        seismic_parameters.GetQualityGrid()->writeCravaFile(file_name_crava, &simbox);
      }
    }
    if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
      for (int i = 0; i < n_facies; i++) {
        std::string file_name = IO::PrefixLikelihood() + facies_names[i];
        ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, lh_cubes_[i], file_name, "");

        if (write_crava_) {
          file_name = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixLikelihood() + facies_names[i]);
          seismic_parameters.GetLHCube()[i]->writeCravaFile(file_name, &simbox);
        }

      }
    }
  }

  //Simulations
  if (model_settings->getNumberOfSimulations() > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Simulation Grids\n");
    int n_simulations = static_cast<int>(simulations_seed0_.size());
    bool kriging      = model_settings->getKrigingParameter() > 0;
    for (int i = 0; i < n_simulations; i++) {
      ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, simulations_seed0_[i], simulations_seed1_[i], simulations_seed2_[i],
                                        model_settings->getOutputGridsElastic(), i, kriging);

      if (write_crava_) {
        std::string prefix = IO::PrefixSimulations();
        std::string suffix;
        if (kriging)
          suffix += "_Kriged_"+NRLib::ToString(i+1);
        else
          suffix += "_"+NRLib::ToString(i+1);

        std::string file_name_vp  = IO::makeFullFileName(IO::PathToInversionResults(), prefix + "Vp" + suffix);
        std::string file_name_vs  = IO::makeFullFileName(IO::PathToInversionResults(), prefix + "Vs" + suffix);
        std::string file_name_rho = IO::makeFullFileName(IO::PathToInversionResults(), prefix + "Rho" + suffix);
        seismic_parameters.GetSimulationSeed0(i)->writeCravaFile(file_name_vp,  &simbox);
        seismic_parameters.GetSimulationSeed1(i)->writeCravaFile(file_name_vs,  &simbox);
        seismic_parameters.GetSimulationSeed2(i)->writeCravaFile(file_name_rho, &simbox);
      }
    }
  }

  //Synthetic seismic
  if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {

    int n_angles = model_settings->getNumberOfAngles(0); //Only write synthetic seismic for the first vintage
    std::vector<float> angles = model_settings->getAngle(0);

    for (int i = 0; i < n_angles; i++) {

      float theta       = angles[i];
      float theta_deg   = static_cast<float>((theta*180.0/NRLib::Pi));
      std::string angle = NRLib::ToString(theta_deg, 1);

      std::string sgri_label = " Synthetic seismic for incidence angle "+angle;
      std::string file_name  = IO::PrefixSyntheticSeismicData() + angle;

      if (((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) || (model_settings->getForwardModeling() == true)) {
        if (i == 0)
          LogKit::LogFormatted(LogKit::Low,"\nWrite Synthetic Seismic\n");
        ParameterOutput::WriteFile(model_settings, synt_seismic_data_[i], file_name, IO::PathToSeismicData(), &simbox, true, sgri_label, time_depth_mapping, false, theta_deg);
      }
    }
  }

  //Trend Cubes
  if ((model_settings->getOutputGridsOther() & IO::TREND_CUBES) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Trend Cubes\n");

    const std::vector<std::string>  & trend_cube_parameters = model_settings->getTrendCubeParameters();

    for (size_t i = 0; i < trend_cubes_.size(); i++) {
      std::string file_name = IO::PrefixTrendCubes() + trend_cube_parameters[i];
      ParameterOutput::WriteFile(model_settings, trend_cubes_[i], file_name, IO::PathToRockPhysics(), &simbox, false, "trend cube", time_depth_mapping);
    }
  }

  Timings::setTimeWriteResults(wall,cpu);
}


void CravaResult::WriteEstimationResults(ModelSettings * model_settings,
                                         CommonData    * common_data)
{
  double wall=0.0, cpu=0.0;
  TimeKit::getTime(wall,cpu);

  //No CombineResults have been run, so much is unset, but we should be able to live with it.
  const Simbox & simbox            = common_data->GetOutputSimbox();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  blocked_logs_ = common_data->GetBlockedLogsOutput();

  n_intervals_ = common_data->GetMultipleIntervalGrid()->GetNIntervals();
  if (n_intervals_ == 1 && ((model_settings->getOutputGridFormat() & IO::CRAVA) > 0))
    write_crava_ = true;

  //Create missing map
  NRLib::Grid2D<bool> * missing_map = CreateMissingGrid(simbox);

  //Logging different options for writing segy grids. Sets offset from topsurface if it is missing
  if ((model_settings->getOutputGridFormat() & IO::SEGY) > 0)
    LogAndSetSegyOffsetIfNeeded(model_settings, simbox);

  //Estimation model: All estimated parameters are written to file, regardless of output settings
  if (((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) &&
     (model_settings->getEstimateBackground() || model_settings->getEstimateWaveletNoise())) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Blocked Logs...");

    //Fill in seismic if needed.
    if (model_settings->getEstimateWaveletNoise()) {
      //ComputeSeismicImpedance from avoinversion.cpp
      reflection_matrix_ = common_data->GetReflectionMatrixTimeLapse(0);
      wavelets_          = common_data->GetWavelet(0);

      //Resample wavelet to output simbox, and generate seismic logs.
      int n_wavelets = static_cast<int>(wavelets_.size());
      for (int i = 0; i < n_wavelets; i++) {
        wavelets_[i]->resample(static_cast<float>(simbox.getdz()),
                               simbox.getnz(),
                               simbox.GetNZpad());

        if (common_data->GetShiftGrid(0, i) != NULL)
          wavelets_[i]->setShiftGrid(new Grid2D(*common_data->GetShiftGrid(0, i)));
        if (common_data->GetGainGrid(0, i) != NULL)
          wavelets_[i]->setGainGrid(new Grid2D(*common_data->GetGainGrid(0, i)));
      }

      CommonData::GenerateSyntheticSeismicLogs(wavelets_, blocked_logs_, reflection_matrix_, &simbox);

      //We need synt_seis from well wavelets, but it needs to be based on output simbox
      //Estimate a temp wavelet, which adds well_synt_seismic_data to blocked logs
      bool wavelet_estimated = false;
      for (size_t i = 0; i < wavelets_.size(); i++) {
        if (model_settings->getEstimateWavelet(0)[i] == true)
          wavelet_estimated = true;
      }
      if (wavelet_estimated == true)
        GenerateWellOptSyntSeis(model_settings, common_data, blocked_logs_, simbox, reflection_matrix_);

      std::vector<SeismicStorage *> seismic_data = common_data->GetSeismicDataTimeLapse(0);
      int n_angles = static_cast<int>(seismic_data.size());
      ModelAVODynamic::AddSeismicLogsFromStorage(blocked_logs_,
                                                 seismic_data,
                                                 simbox,
                                                 n_angles);
    }
    WriteBlockedWells(common_data->GetBlockedLogsOutput(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());
    LogKit::LogFormatted(LogKit::Low,"ok");
  }

  if (model_settings->getEstimateBackground() == true) {

    if (n_intervals_ > 1) {
      //Resample background grids to one final grid
      //We don't have seismicParametersHolder (with FFT-grids), only CommonData (With NRLib-grids)
      //Have to transform from NRLib-Grid to FFTGrid to be used in CombineResult.
      LogKit::LogFormatted(LogKit::Low,"\nCombine Background Grids");
      MultiIntervalGrid      * multi_interval_grid = common_data->GetMultipleIntervalGrid();

      int nx = simbox.getnx();
      int ny = simbox.getny();
      int nz = simbox.getnz();

      std::vector<FFTGrid *> dummy_fft_grids;

      std::vector<NRLib::Grid<float> *> background_vp_intervals(n_intervals_);
      std::vector<NRLib::Grid<float> *> background_vs_intervals(n_intervals_);
      std::vector<NRLib::Grid<float> *> background_rho_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        background_vp_intervals[i]  = common_data->GetBackgroundVpInterval(i);
        background_vs_intervals[i]  = common_data->GetBackgroundVsInterval(i);
        background_rho_intervals[i] = common_data->GetBackgroundRhoInterval(i);

      }
      background_vp_  = new StormContGrid(simbox, nx, ny, nz);
      background_vs_  = new StormContGrid(simbox, nx, ny, nz);
      background_rho_ = new StormContGrid(simbox, nx, ny, nz);
      std::vector<StormContGrid> zone_prob_grid(multi_interval_grid->GetNIntervals());
      LogKit::LogFormatted(LogKit::Low,"\n Setting up zone probabilities");
      for(size_t i=0; i<zone_prob_grid.size();i++)
        zone_prob_grid[i] = NRLib::StormContGrid(simbox, nx, ny, nz);
      multi_interval_grid->FindZoneProbGrid(zone_prob_grid);

      LogKit::LogFormatted(LogKit::Low,"\n Vp");
      CombineResult(background_vp_,  dummy_fft_grids,  multi_interval_grid, zone_prob_grid, background_vp_intervals, missing_map, model_settings->getFilterMultizoneModel(), model_settings->getMaxHzBackground());
      LogKit::LogFormatted(LogKit::Low,"\n Vs");
      CombineResult(background_vs_,  dummy_fft_grids,  multi_interval_grid, zone_prob_grid, background_vs_intervals, missing_map, model_settings->getFilterMultizoneModel(),model_settings->getMaxHzBackground());
      LogKit::LogFormatted(LogKit::Low,"\n Rho");
      CombineResult(background_rho_, dummy_fft_grids, multi_interval_grid, zone_prob_grid, background_rho_intervals, missing_map, model_settings->getFilterMultizoneModel(), model_settings->getMaxHzBackground());

      for (int i = 0; i < n_intervals_; i++) {
        common_data->ReleaseBackgroundGrids(i, 0);
        common_data->ReleaseBackgroundGrids(i, 1);
        common_data->ReleaseBackgroundGrids(i, 2);
      }

      //Background trends
      if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0) {
        const NRLib::Grid2D<std::vector<double> > & vertical_trends = common_data->GetBackgroundVerticalTrends();
        if (vertical_trends.GetNI() > 0) {

          background_trend_vp_  = new StormContGrid(simbox, nx, ny, nz);
          background_trend_vs_  = new StormContGrid(simbox, nx, ny, nz);
          background_trend_rho_ = new StormContGrid(simbox, nx, ny, nz);

          CombineVerticalTrends(multi_interval_grid,
                                common_data,
                                vertical_trends,
                                zone_prob_grid,
                                background_trend_vp_,
                                background_trend_vs_,
                                background_trend_rho_,
                                missing_map);
        }
      }
    }
    else {
      std::vector<NRLib::Grid<float> *> bg = common_data->GetBackgroundParametersInterval(0);
      background_vp_  = new NRLib::StormContGrid(simbox, *(bg[0]));
      background_vs_  = new NRLib::StormContGrid(simbox, *(bg[1]));
      background_rho_ = new NRLib::StormContGrid(simbox, *(bg[2]));

      common_data->ReleaseBackgroundGrids(0, 0);
      common_data->ReleaseBackgroundGrids(0, 1);
      common_data->ReleaseBackgroundGrids(0, 2);

      //Background trends
      if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0) {

        int nx = simbox.getnx();
        int ny = simbox.getny();
        int nz = simbox.getnz();

        background_trend_vp_  = new StormContGrid(simbox, nx, ny, nz);
        background_trend_vs_  = new StormContGrid(simbox, nx, ny, nz);
        background_trend_rho_ = new StormContGrid(simbox, nx, ny, nz);

        const NRLib::Grid2D<std::vector<double> > & vertical_trends = common_data->GetBackgroundVerticalTrends();

        FFTGrid * trend_grid_vp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nx, ny, nz, false);
        Background::FillInVerticalTrend(trend_grid_vp, vertical_trends(0,0));
        CreateStormGrid(*background_trend_vp_, trend_grid_vp);

        FFTGrid * trend_grid_vs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nx, ny, nz, false);
        Background::FillInVerticalTrend(trend_grid_vs, vertical_trends(0,1));
        CreateStormGrid(*background_trend_vs_, trend_grid_vs);

        FFTGrid * trend_grid_rho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nx, ny, nz, false);
        Background::FillInVerticalTrend(trend_grid_rho, vertical_trends(0,2));
        CreateStormGrid(*background_trend_rho_, trend_grid_rho);

      }
    }

    LogKit::LogFormatted(LogKit::Low,"\nWrite Background Grids\n");
    WriteGridPackage(model_settings,
                     &simbox,
                     background_vp_,
                     background_vs_,
                     background_rho_,
                     time_depth_mapping,
                     IO::PrefixBackground(),
                     IO::PathToBackground(),
                     true);

    if (write_crava_) {
      std::string file_name_vp  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vp");
      std::string file_name_vs  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vs");
      std::string file_name_rho = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Rho");

      ExpTransf(background_vp_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vp, &simbox);

      ExpTransf(background_vs_intervals_[0]);
      background_vs_intervals_[0]->writeCravaFile(file_name_vs, &simbox);

      ExpTransf(background_rho_intervals_[0]);
      background_rho_intervals_[0]->writeCravaFile(file_name_rho, &simbox);
    }

    if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0) {
      LogKit::LogFormatted(LogKit::Low,"\nWrite Background Vertical Trend Grids\n");

      std::string prefix = IO::PrefixBackground()+IO::PrefixTrend();
      WriteGridPackage(model_settings,
                        &simbox,
                        background_trend_vp_,
                        background_trend_vs_,
                        background_trend_rho_,
                        time_depth_mapping,
                        prefix,
                        IO::PathToBackground(),
                        true);
    }
  }

  //Write original seismic data. Resample from CommonData to output_simbox. If CRAVA-format, we resample with padding.
  if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0) {
    WriteSeismicData(model_settings, common_data, simbox, time_depth_mapping);
  }

  Timings::setTimeWriteResults(wall,cpu);

}

void CravaResult::WriteSeismicData(ModelSettings * model_settings,
                                   CommonData    * common_data,
                                   const Simbox  & simbox,
                                   GridMapping   * time_depth_mapping)
{
  //Separate functions, used both in inversion mode (WriteResults) and in estimation mode (WriteEstimationResults)
  int n_timelapses = model_settings->getNumberOfTimeLapses();
  LogKit::LogFormatted(LogKit::Low,"\nWrite Seismic Data\n");

  for (int i = 0; i < n_timelapses; i ++) {

    int n_angles              = model_settings->getNumberOfAngles(i);
    std::vector<float> angles = model_settings->getAngle(i);

    for (int j = 0; j < n_angles; j++) {
      std::string angle           = NRLib::ToString(angles[j]*(180/M_PI), 1);
      std::string file_name_orig  = IO::PrefixOriginalSeismicData() + angle;
      std::string sgri_label      = std::string("Original seismic data for angle stack ") + angle;

      int seismic_type = common_data->GetSeismicDataTimeLapse(i)[j]->GetSeismicType();

      FFTGrid * fft_grid_resampled = NULL;
      bool delete_fft_grid         = true;
      if (seismic_type == 3) {
        LogKit::LogFormatted(LogKit::Warning," Cannot write " + file_name_orig + ". Not allowed to write original seismic to file when input seismic is on CRAVA-format.\n");
      }
      else {
        SegY * segy                           = NULL;
        StormContGrid * storm                 = NULL;
        NRLib::Grid<float> * nrlib_grid       = new NRLib::Grid<float>();
        bool is_segy                          = false;
        bool is_storm                         = false;
        NRLib::Grid2D<bool> * dead_traces_map = new NRLib::Grid2D<bool>();
        if (seismic_type == 0) {//SEGY
          segy = common_data->GetSeismicDataTimeLapse(i)[j]->GetSegY();
          is_segy = true;
        }
        else if (seismic_type == 1 || seismic_type == 2) {//STORM/SGRI
          storm = common_data->GetSeismicDataTimeLapse(i)[j]->GetStorm();
          is_storm = true;
        }
        if (write_crava_)
          fft_grid_resampled = new FFTGrid(simbox.getnx(), simbox.getny(), simbox.getnz(), simbox.GetNXpad(), simbox.GetNYpad(), simbox.GetNZpad());
        else
          fft_grid_resampled = new FFTGrid(simbox.getnx(), simbox.getny(), simbox.getnz(), simbox.getnx(), simbox.getny(), simbox.getnz());
        fft_grid_resampled->createRealGrid();

        int missing_traces_simbox  = 0;
        int missing_traces_padding = 0;
        int dead_traces_simbox     = 0;

        common_data->FillInData(nrlib_grid,
                                fft_grid_resampled,
                                &simbox,
                                storm,
                                segy,
                                model_settings->getSmoothLength(),
                                missing_traces_simbox,
                                missing_traces_padding,
                                dead_traces_simbox,
                                dead_traces_map,
                                FFTGrid::DATA,
                                false,
                                is_segy,
                                is_storm,
                                true);

        delete nrlib_grid;

        if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0 && write_crava_) {
          std::string file_name_crava = IO::makeFullFileName(IO::PathToSeismicData(), IO::PrefixOriginalSeismicData() + angle);
          fft_grid_resampled->writeCravaFile(file_name_crava, &simbox);
        }

        //Setting up storm grid for writing
        if (((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0) || ((model_settings->getOutputGridsSeismic() & IO::RESIDUAL) > 0)) {
          StormContGrid * seismic_storm = CreateStormGrid(simbox, fft_grid_resampled, delete_fft_grid); // deletes fft_grid_resampled

          //Real seismic gives value at cell base, synthetic at cell top. Shift real seismic.
          for (int k = static_cast<int>(seismic_storm->GetNK())-1; k > 0; k--) {
            for (size_t jj = 0; jj < seismic_storm->GetNJ(); jj++) {
              for (size_t ii = 0; ii < seismic_storm->GetNI(); ii++) {
                (*seismic_storm)(ii,jj,k) = (*seismic_storm)(ii,jj,k-1);
              }
            }
          }


          float theta_deg   = static_cast<float>((angles[j]*180.0/NRLib::Pi));
          if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0)
            ParameterOutput::WriteFile(model_settings, seismic_storm, file_name_orig, IO::PathToSeismicData(), &simbox, true, sgri_label, time_depth_mapping, false, theta_deg);

          if (model_settings->getEstimationMode() == false) {
            if ((i==0) && ((model_settings->getOutputGridsSeismic() & IO::RESIDUAL) > 0)) { //residuals only for first vintage.

              StormContGrid residual(*(synt_seismic_data_[j]));
              for (size_t k=0;k<seismic_storm->GetNK();k++) {
                for (size_t jj=0;jj<seismic_storm->GetNJ();jj++) {
                  for (size_t ii=0;ii<seismic_storm->GetNI();ii++) {
                    residual(ii,jj,k) = (*seismic_storm)(ii,jj,k)-residual(ii,jj,k);
                  }
                }
              }

              sgri_label = "Residual computed from synthetic seismic for incidence angle "+angle;
              std::string file_name  = IO::PrefixResiduals() + angle;

              ParameterOutput::WriteFile(model_settings, &residual, file_name, IO::PathToSeismicData(), &simbox, true, sgri_label, time_depth_mapping, false, theta_deg);
            }
          }
        }
      }
    }
  }

}

void CravaResult::WriteFilePriorCorrT(fftw_real   * prior_corr_T,
                                      const int   & nzp,
                                      const float & dt,
                                      std::string   interval_name) const
{
  if (interval_name != "")
    interval_name = "_"+interval_name;

  // This is the cyclic and filtered version of CorrT
  std::string base_name = IO::PrefixPrior() + IO::FileTemporalCorr() + interval_name + IO::SuffixGeneralData();
  std::string file_name = IO::makeFullFileName(IO::PathToCorrelations(), base_name);
  std::ofstream file;
  NRLib::OpenWrite(file, file_name);
  file << std::fixed
       << std::right
       << std::setprecision(6)
       << dt << "\n";
  for (int i = 0; i < nzp; i++) {
    file << std::setw(9) << prior_corr_T[i] << "\n";
  }
  file.close();
}

void CravaResult::WriteFilePostVariances(const NRLib::Matrix      & post_var0,
                                         const std::vector<float> & post_cov_vp00,
                                         const std::vector<float> & post_cov_vs00,
                                         const std::vector<float> & post_cov_rho00,
                                         std::string                interval_name) const
{
  if (interval_name != "")
    interval_name = "_"+interval_name;

  std::string base_name = IO::PrefixPosterior() + IO::FileParameterCov() + interval_name + IO::SuffixGeneralData();
  std::string file_name = IO::makeFullFileName(IO::PathToCorrelations(), base_name);

  std::ofstream file;
  NRLib::OpenWrite(file, file_name);
  file << std::fixed;
  file << std::right;
  file << std::setprecision(6);
  for (int i=0 ; i<3 ; i++) {
    for (int j=0 ; j<3 ; j++) {
      file << std::setw(10) << post_var0(i,j) << " ";
    }
    file << "\n";
  }
  file.close();

  std::string base_name1 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vp" +IO::SuffixGeneralData();
  std::string base_name2 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Vs" +IO::SuffixGeneralData();
  std::string base_name3 = IO::PrefixPosterior() + IO::PrefixTemporalCorr()+"Rho"+IO::SuffixGeneralData();
  WriteFilePostCorrT(post_cov_vp00,  IO::PathToCorrelations(), base_name1);
  WriteFilePostCorrT(post_cov_vs00,  IO::PathToCorrelations(), base_name2);
  WriteFilePostCorrT(post_cov_rho00, IO::PathToCorrelations(), base_name3);
}

void CravaResult::WriteFilePostCorrT(const std::vector<float> & post_cov,
                                     const std::string        & sub_dir,
                                     const std::string        & base_name) const
{
  std::string fileName = IO::makeFullFileName(sub_dir,base_name);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed;
  file << std::setprecision(6);
  file << std::right;
  float c0 = 1.0f/post_cov[0];

  for (int k=0 ; k < static_cast<int>(post_cov.size()) ; k++)
    file << std::setw(9) << post_cov[k]*c0 << "\n";

  file.close();
}

void CravaResult::WriteBlockedWells(const std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                    const ModelSettings                              * model_settings,
                                    std::vector<std::string>                           facies_name,
                                    std::vector<int>                                   facies_label)
{
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    blocked_log->WriteWell(model_settings->getWellFormatFlag(),
                           model_settings->getMaxHzBackground(),
                           model_settings->getMaxHzSeismic(),
                           facies_name,
                           facies_label);

  }
}


void CravaResult::WriteWells(const std::vector<NRLib::Well> & wells,
                             const ModelSettings            * model_settings)
{
  assert(wells.size() >= 0);        //Placeholder code. Do these wells have any function anymore?
  assert(model_settings != NULL);
/*
  for (size_t i = 0; i < wells.size(); i++) {
    wells[i].WriteWell(model_settings->getWellFormatFlag(),
                       model_settings->getMaxHzBackground(),
                       model_settings->getMaxHzSeismic());

  }
  */
}

StormContGrid *
CravaResult::CreateStormGrid(const Simbox & simbox,
                             FFTGrid      * fft_grid,
                             bool           delete_fft)
{
  if (fft_grid != NULL) {
    int nx = fft_grid->getNx();
    int ny = fft_grid->getNy();
    int nz = fft_grid->getNz();

    StormContGrid * storm = new StormContGrid(simbox, nx, ny, nz);

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          float value = fft_grid->getRealValue(i, j, k);
          storm->SetValue(i, j, k, value);
        }
      }
    }

    if (delete_fft)
      delete fft_grid;

    return(storm);
  }

  StormContGrid * storm = new StormContGrid();
  return(storm);

}

void CravaResult::CreateStormGrid(StormContGrid & grid_new,
                                  FFTGrid       * fft_grid,
                                  bool            allow_delete)
{
  int nx = static_cast<int>(grid_new.GetNI());
  int ny = static_cast<int>(grid_new.GetNJ());
  int nz = static_cast<int>(grid_new.GetNK());

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        float value = fft_grid->getRealValue(i, j, k);
        grid_new.SetValue(i, j, k, value);
      }
    }
  }

  if (write_crava_ == false && allow_delete == true)
    delete fft_grid;
}

void CravaResult::WriteGridPackage(const ModelSettings     * model_settings,
                                   const Simbox            * simbox,
                                   StormContGrid           * grid_vp,
                                   StormContGrid           * grid_vs,
                                   StormContGrid           * grid_rho,
                                   GridMapping             * depth_mapping,
                                   const std::string       & prefix,
                                   const std::string       & path,
                                   bool                      exp_transf)
{
  if (depth_mapping != NULL && depth_mapping->getSimbox() == NULL) {

    depth_mapping->setMappingFromVelocity(grid_vp, simbox, model_settings->getOutputGridFormat());
  }

  std::string file_name_vp  = prefix + "Vp" ;
  std::string file_name_vs  = prefix + "Vs" ;
  std::string file_name_rho = prefix + "Rho";

  if (exp_transf) {
    ExpTransf(grid_vp);
    ExpTransf(grid_vs);
    ExpTransf(grid_rho);
  }

  ParameterOutput::WriteFile(model_settings, grid_vp,  file_name_vp,  path, simbox, false, "NO_LABEL", depth_mapping);
  ParameterOutput::WriteFile(model_settings, grid_vs,  file_name_vs,  path, simbox, false, "NO_LABEL", depth_mapping);
  ParameterOutput::WriteFile(model_settings, grid_rho, file_name_rho, path, simbox, false, "NO_LABEL", depth_mapping);

}


void CravaResult::ExpTransf(StormContGrid * grid)
{
  float value = 0.0f;
  for (size_t i = 0; i < grid->GetNI(); i++) {
    for (size_t j = 0; j < grid->GetNJ(); j++) {
      for (size_t k = 0; k < grid->GetNK(); k++) {

        value = grid->GetValue(i, j, k);

        if (value != RMISSING) {
          value = exp(value);
          grid->SetValue(i, j, k, value);
        }

      }
    }
  }
}

void CravaResult::LogTransf(StormContGrid * grid)
{
  float value = 0.0f;
  for (size_t i = 0; i < grid->GetNI(); i++) {
    for (size_t j = 0; j < grid->GetNJ(); j++) {
      for (size_t k = 0; k < grid->GetNK(); k++) {

        value = grid->GetValue(i, j, k);

        if (value != RMISSING) {
          if (value < 0.0) {
            grid->SetValue(i, j, k, 0.0);
          }
          else {
            value = log(value);
            grid->SetValue(i, j, k, value);
          }
        }

      }
    }
  }
}

void CravaResult::ExpTransf(FFTGrid * grid)
{
  float value = 0.0f;
  for (int i = 0; i < grid->getNx(); i++) {
    for (int j = 0; j < grid->getNy(); j++) {
      for (int k = 0; k < grid->getNz(); k++) {

        value = grid->getRealValue(i, j, k);

        if (value != RMISSING) {
          value = exp(value);
          grid->setRealValue(i, j, k, value);
        }

      }
    }
  }
}

void CravaResult::LogTransf(FFTGrid * grid)
{
  float value = 0.0f;
  for (int i = 0; i < grid->getNx(); i++) {
    for (int j = 0; j < grid->getNy(); j++) {
      for (int k = 0; k < grid->getNz(); k++) {

        value = grid->getRealValue(i, j, k);

        if (value == RMISSING || value < 0.0) {
          grid->setRealValue(i, j, k, 0);
        }
        else {
          value = log(value);
          grid->setRealValue(i, j, k, value);
        }
      }
    }
  }
}

void CravaResult::ComputeSyntSeismic(const ModelSettings          * model_settings,
                                     const Simbox                 * simbox,
                                     std::vector<Wavelet *>       & wavelets,
                                     StormContGrid                * vp,
                                     StormContGrid                * vs,
                                     StormContGrid                * rho,
                                     std::vector<StormContGrid *> & synt_seismic_data)
{
  //LogKit::WriteHeader("Compute Synthetic Seismic and Residuals");

  int nx = static_cast<int>(vp->GetNI());
  int ny = static_cast<int>(vp->GetNJ());
  int nz = static_cast<int>(vp->GetNK());

  int nzp = simbox->GetNZpad();

  std::vector<float> angles = model_settings->getAngle(0); //Synt seismic only for first vintage
  int n_theta = static_cast<int>(angles.size());

  for (int l = 0; l < n_theta; l++) {
    StormContGrid * imp = ComputeSeismicImpedance(vp, vs, rho, reflection_matrix_, l);

    for (int i = 0; i <nx; i++) {
      for (int j = 0; j < ny; j++) {
        Wavelet1D impVec(0, nz, nzp);
        //impVec.setupAsVector();
        int k;
        for (k = 0; k < nz; k++) {
          float value = imp->GetValue(i, j, k);
          impVec.setRAmp(value, k);
        }
        //Tapering:
        float fac = 1.0f/static_cast<float>(nzp-nz-1);
        for (; k < nzp; k++) {
          float value = fac*((k-nz)*impVec.getRAmp(0)+(nzp-k-1)*impVec.getRAmp(nz-1));
          impVec.setRAmp(value, k);
        }
        Wavelet1D resultVec(&impVec, Wavelet::FIRSTORDERBACKWARDDIFF);
        resultVec.fft1DInPlace();

        Wavelet1D * localWavelet = wavelets[l]->createLocalWavelet1D(i,j);

        float sf = static_cast<float>(simbox->getRelThick(i, j))*wavelets[l]->getLocalStretch(i,j);

        for (k = 0; k <(nzp/2 +1); k++) {
          fftw_complex r = resultVec.getCAmp(k);
          fftw_complex w = localWavelet->getCAmp(k,static_cast<float>(sf));// returns complex conjugate
          fftw_complex s;
          s.re = r.re*w.re+r.im*w.im; //Use complex conjugate of w
          s.im = -r.re*w.im+r.im*w.re;
          resultVec.setCAmp(s,k);
        }
        delete localWavelet;

        resultVec.invFFT1DInPlace();
        for (k = 0; k < nz; k++) { //nzp
          float value = resultVec.getRAmp(k);
          imp->SetValue(i, j, k, value);
        }
      }
    }

    synt_seismic_data.push_back(imp);

  }
}

StormContGrid *
CravaResult::ComputeSeismicImpedance(StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     StormContGrid       * rho,
                                     const NRLib::Matrix & reflection_matrix,
                                     int                   angle) const
{
  int nx = static_cast<int>(vp->GetNI());
  int ny = static_cast<int>(vp->GetNJ());
  int nz = static_cast<int>(vp->GetNK());

  StormContGrid * impedance = new StormContGrid(nx, ny, nz);

  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        float imp = 0;
        imp += vp->GetValue(i, j, k)*static_cast<float>(reflection_matrix(angle,0));
        imp += vs->GetValue(i, j, k)*static_cast<float>(reflection_matrix(angle,1));
        imp += rho->GetValue(i, j, k)*static_cast<float>(reflection_matrix(angle,2));

        impedance->SetValue(i, j, k, imp);
      }
    }
  }

  return(impedance);
}

void CravaResult::AddBlockedLogs(const std::map<std::string, BlockedLogsCommon *> & blocked_logs)
{
  std::map<std::string, BlockedLogsCommon *> new_blocked_logs;

  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    BlockedLogsCommon * new_blocked_log = new BlockedLogsCommon(*blocked_log); //copy

    std::string name = blocked_log->GetWellName();

    new_blocked_logs.insert(std::pair<std::string, BlockedLogsCommon *>(name, new_blocked_log));
  }

  blocked_logs_intervals_.push_back(new_blocked_logs);
}

void CravaResult::GenerateWellOptSyntSeis(ModelSettings                              * model_settings,
                                          CommonData                                 * common_data,
                                          std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                          const Simbox                               & simbox,
                                          const NRLib::Matrix                        & reflection_matrix)
{
  int n_angles = model_settings->getNumberOfAngles(0);
  int nzp = simbox.GetNZpad();
  int nz  = simbox.getnz();
  std::vector<std::vector<Wavelet *> > wavelets(blocked_wells.size()); //First index is well, seoond is angle.
  for (int i = 0; i < n_angles; i++) {
    std::vector<Wavelet1D *> angle_wavelets = common_data->GetWellWavelets(i);
    for (size_t w=0;w<blocked_wells.size();w++) {
      if (i==0)
        wavelets[w].resize(n_angles);
      if (angle_wavelets.size() == 0 || angle_wavelets[w] == NULL)
        wavelets[w][i] = NULL;
      else {
        Wavelet1D * wavelet = new Wavelet1D(angle_wavelets[w]);
        wavelet->resample(static_cast<float>(simbox.getdz()),
                          simbox.getnz(),
                          simbox.GetNZpad());
        wavelets[w][i] = wavelet;
      }
    }
  }

  size_t w=0;
  //Need to set angles in blocked_log before wavelet estimation
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    BlockedLogsCommon * blocked_log = it->second;
    blocked_log->SetNAngles(n_angles);

    if (blocked_log->GetIsDeviated() == false)
      blocked_log->GenerateSyntheticSeismic(reflection_matrix, wavelets[w], nz, nzp, &simbox, true);
    w++;
  }

  for (int i = 0; i < n_angles; i++) {
    for (size_t ww=0;ww<blocked_wells.size();ww++) {
      delete wavelets[ww][i];
    }
  }
}


void
CravaResult::CreateDownscalingPlans(const std::vector<int>    & nzp,
                                    int                         scale,
                                    std::vector<rfftwnd_plan> & small_plans,
                                    std::vector<rfftwnd_plan> & big_plans)
{
  small_plans.resize(nzp.size());
  big_plans.resize(nzp.size());
  for(size_t zone=0;zone<small_plans.size();zone++) {
    int nt = nzp[zone];;
    int mt = nt*scale;
    small_plans[zone] = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
    big_plans[zone]   = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
  }
}


void
CravaResult::AddPadding(std::vector<float> & trace,
                        int                  nzp)
{
  int old_size = static_cast<int>(trace.size());
  trace.resize(nzp);
  for(int i = old_size;i<nzp;i++) { //Tapering
    float t  = static_cast<float>(i-old_size+1)/static_cast<float>(nzp-old_size+1);
    trace[i] = (1-t)*trace[old_size-1]+t*trace[0];
  }
}

void
CravaResult::DownscaleTrace(const std::vector<float> & trace_in,
                            std::vector<float>       & trace_out,
                            int                        scale,
                            int                        prepad_size,
                            const rfftwnd_plan       & small_plan,
                            const rfftwnd_plan       & big_plan)
{
  //Assumes trace_in is already padded.
  int nt = static_cast<int>(trace_in.size());
  int mt = nt*scale;

  int cnt = nt/2 + 1;
  int rnt = 2*cnt;
  int cmt = mt/2 + 1;
  int rmt = 2*cmt;

  fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
  fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

  CommonData::ResampleTrace(trace_in,
                            small_plan,
                            big_plan,
                            rAmpData,
                            rAmpFine,
                            nt,
                            cnt,
                            rnt,
                            cmt,
                            rmt);

  //Add trend
  size_t out_len = prepad_size*scale-(scale-1); //Points beyond this are contaminated by padding.
  trace_out.resize(out_len);
  for (size_t k = 0; k < out_len; k++) {
    trace_out[k] = rAmpFine[k];
  }

  fftw_free(rAmpData);
  fftw_free(rAmpFine);
}

NRLib::Grid2D<bool> *
CravaResult::CreateMissingGrid(const Simbox & simbox)
{

  int nx = simbox.getnx();
  int ny = simbox.getny();
  NRLib::Grid2D<bool> * missing_map = new NRLib::Grid2D<bool>(nx, ny, false);

  bool has_missing = false;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      if (simbox.getTop(i,j) == RMISSING || simbox.getBot(i,j) == RMISSING) {
        (*missing_map)(i,j) = true;
        has_missing = true;
      }
    }
  }
  if (has_missing == false)
    missing_map = NULL;

  return missing_map;
}

void CravaResult::SetMissingInGrid(StormContGrid       & grid,
                                   NRLib::Grid2D<bool> * missing_map)
{
  for (size_t i = 0; i < grid.GetNI(); i++) {
    for (size_t j = 0; j < grid.GetNJ(); j++) {

      if ((*missing_map)(i,j) == true) {
        for (size_t k = 0; k < grid.GetNK(); k++)
          grid.SetValue(i,j,k,RMISSING);
      }
    }
  }
}


void CravaResult::LogAndSetSegyOffsetIfNeeded(ModelSettings * model_settings,
                                              const Simbox  & simbox)
{

  double simbox_dz = floor(simbox.getdz());
  double segy_dz   = floor(model_settings->getSegyDz());

  //If simbox dz and segy dz match we write output segy matching input segy. This check can be unstable (?) as we go from dz to nz back to dz
  // so if we know that they should match we have set getUseInputSegyDzForOutputSegy variable
  if ((simbox_dz == segy_dz || model_settings->getUseInputSegyDzForOutputSegy() == true) && model_settings->getMatchOutputInputSegy() == true) {
    LogKit::LogFormatted(LogKit::Low, "\n\nThe output segy grid will be written out matching the input seismic segy cubes (dz = " + NRLib::ToString(floor(segy_dz)) + ").\n");
  }
  else if (simbox_dz != segy_dz && model_settings->getSegyDz() != RMISSING && model_settings->getMatchOutputInputSegy() == true) {
    LogKit::LogFormatted(LogKit::Low, "\n\nWarning: The input segy dz (" + NRLib::ToString(model_settings->getSegyDz()) + ") does not match the output dz ("
                          + NRLib::ToString(simbox.getdz()) +"). The output segy grid can therefore not be matched with input segy grid. The output dz will be used.\n");
  }
  //Offset given in modelfile under grid-output, we then override the matching of segy grid
  else if (model_settings->getMatchOutputInputSegy() == false && model_settings->getOutputOffset() != RMISSING) {
    LogKit::LogFormatted(LogKit::Low, "\n\nSegy grids will be written out with offset " + NRLib::ToString(model_settings->getOutputOffset()) + ".\n");
  }
  //No segy cubes and offset not given in modelfile. Take offset from topsurface.
  else if (model_settings->getOutputOffset() == RMISSING) {
    float offset = static_cast<float>(floor(simbox.getTopZMin()));
    model_settings->setOutputOffset(static_cast<float>(floor(simbox.getTopZMin())));
    LogKit::LogFormatted(LogKit::High, "\n\nSegy grids will be written out with offset " + NRLib::ToString(offset) + ", taken from the top of the top surface.\n");
    model_settings->setOutputOffset(offset);
  }

}
