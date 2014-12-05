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

  MultiIntervalGrid * multi_interval_grid     = common_data->GetMultipleIntervalGrid();
  Simbox & output_simbox                      = common_data->GetOutputSimbox();
  n_intervals_                                = multi_interval_grid->GetNIntervals();

  std::vector<NRLib::Grid<float> *> dummy_grids;

  int nx           = output_simbox.getnx();
  int ny           = output_simbox.getny();
  int nz_output    = output_simbox.getnz();

  LogKit::LogFormatted(LogKit::Low,"\nCombine Blocked Logs...");
  blocked_logs_ = common_data->GetBlockedLogsOutput(); //Logs blocked to output_simbox
  CombineBlockedLogs(blocked_logs_, blocked_logs_intervals_, multi_interval_grid, common_data, &output_simbox); //Combine and resample logs create during inversion
  LogKit::LogFormatted(LogKit::Low,"ok");

  //If we want to write grids on CRAVA-format we do not delete the fft-grids in CombineResult so we can write the grids with padding in WriteResults
  if (n_intervals_ == 1 && ((model_settings->getOutputGridFormat() & IO::CRAVA) > 0))
    write_crava_ = true;

  //Output on CRAVA-format requires padding and is only used with single zone. Since they already are on fft-format with padding, we write them out here
  //before they are transformed to storm-grids and deleted.
  //WriteCravaGrids(model_settings, common_data, output_simbox, seismic_parameters_intervals[0]);

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
    CombineResult(post_vp_,  post_vp_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(post_vs_,  post_vs_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(post_rho_, post_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);

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
      CombineResult(post_vp_kriged_,  post_vp_kriged_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
      LogKit::LogFormatted(LogKit::Low,"\n Vs kriged ");
      CombineResult(post_vs_kriged_,  post_vs_kriged_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
      LogKit::LogFormatted(LogKit::Low,"\n Rho kriged ");
      CombineResult(post_rho_kriged_, post_rho_kriged_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
    }
  }

  //Background models
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nCombine Background Model");
    background_vp_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    LogKit::LogFormatted(LogKit::Low,"\n Vp ");
    CombineResult(background_vp_,  background_vp_intervals_,  multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(background_vs_,  background_vs_intervals_,  multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(background_rho_, background_rho_intervals_, multi_interval_grid, zone_prob_grid, dummy_grids);
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

  //Covariance grids
  if (!model_settings->getForwardModeling()) {
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
    CombineResult(cov_vp_,        cov_vp_intervals,        multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Vs ");
    CombineResult(cov_vs_,        cov_vs_intervals,        multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n Rho ");
    CombineResult(cov_rho_,       cov_rho_intervals,       multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n VpVs ");
    CombineResult(cr_cov_vp_vs_,  cr_cov_vp_vs_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n VpRho ");
    CombineResult(cr_cov_vp_rho_, cr_cov_vp_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
    LogKit::LogFormatted(LogKit::Low,"\n VsRho ");
    CombineResult(cr_cov_vs_rho_, cr_cov_vs_rho_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
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
      CombineResult(facies_prob_[j],  facies_prob_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
    }

    //Undef
    facies_prob_undef_ = new StormContGrid(output_simbox, nx, ny, nz_output);
    std::vector<FFTGrid *> facies_prob_intervals_undef(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      facies_prob_intervals_undef[i] = seismic_parameters_intervals[i].GetFaciesProbUndefined();
    }
    LogKit::LogFormatted(LogKit::Low,"\n Undef ");
    CombineResult(facies_prob_undef_,  facies_prob_intervals_undef,  multi_interval_grid, zone_prob_grid, dummy_grids);

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
      CombineResult(facies_prob_geo_[j],  facies_prob_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
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
      CombineResult(lh_cubes_[j], lh_cubes_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
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
    CombineResult(quality_grid_, quality_grid_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
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
      CombineResult(simulations_seed0_[j], simulations_seed0_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
      LogKit::LogFormatted(LogKit::Low,"\n seed1 " + NRLib::ToString(j) + " ");
      CombineResult(simulations_seed1_[j], simulations_seed1_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
      LogKit::LogFormatted(LogKit::Low,"\n seed2 " + NRLib::ToString(j) + " ");
      CombineResult(simulations_seed2_[j], simulations_seed2_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);

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
    CombineResult(block_grid_, block_grid_intervals, multi_interval_grid, zone_prob_grid, dummy_grids);
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
  GenerateSyntheticSeismicLogs(wavelets_, blocked_logs_, reflection_matrix_, output_simbox);

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
      CombineResult(trend_cubes_[i],  trend_cubes_intervals,  multi_interval_grid, zone_prob_grid, dummy_grids);
    }
  }

  //Delete grids from seismicparamtersholder
  if (model_settings->getEstimationMode() || !model_settings->getForwardModeling()) {
    for (int i = 0; i < n_intervals_; i++)
      seismic_parameters_intervals[i].releaseExpGrids();
  }

}

void CravaResult::CombineResult(StormContGrid                    *& final_grid,
                                std::vector<FFTGrid *>            & interval_grids,
                                MultiIntervalGrid                 * multi_interval_grid,
                                const std::vector<StormContGrid>  & zone_probability,
                                std::vector<NRLib::Grid<float> *> & interval_grids_nrlib) //Optional, send in an empty vector if FFTGrids are used
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

  //Resample
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      std::vector<std::vector<float> > new_traces(n_intervals_); //Finely interpolated values

      //Resample each trace to new nz
      for (int i_interval = 0; i_interval < n_intervals_; i_interval++) {

        Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

        std::vector<float> old_trace;
        if (use_nrlib_grids == false) {
          FFTGrid * interval_grid  = interval_grids[i_interval];
          old_trace = interval_grid->getRealTrace(i, j); //old_trace is changed below.
        }
        else {
          old_trace = GetNRLibGridTrace(interval_grids_nrlib[i_interval], i, j);
        }

        int nz_old       = interval_simbox->getnz();

        int scale = 10; //How densely to sample "fine" values.
        new_traces[i_interval].resize(nz_old*scale);

        DownscaleTrace(old_trace, new_traces[i_interval], scale);

      } //n_intervals


      //Combine vectors for each interval to one trace in stormgrid
      for (int k = 0; k < nz; k++) {
        double global_x = 0.0;
        double global_y = 0.0;
        double global_z = 0.0;

        double value = 0;
        final_grid->FindCenterOfCell(i, j, k, global_x, global_y, global_z);
        for (int i_interval = 0; i_interval < n_intervals_; i_interval++) {
          if(zone_probability[i_interval](i,j,k) > 0) {
            Simbox * z_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);
            double dummy1, dummy2, rel_index;
            z_simbox->getInterpolationIndexes(global_x, global_y, global_z, dummy1, dummy2, rel_index);
            rel_index -= 0.5; //First half grid cell is outside interpolation vector.
            rel_index /= static_cast<double>(z_simbox->getnz()-1);
            if(rel_index < 0)
              rel_index = 0;
            else if(rel_index > 1)
              rel_index = 1;

            int index = static_cast<int>(floor(0.5+rel_index*(new_traces[i_interval].size()-1))); //0 to first item, 1 to last item.
            value += zone_probability[i_interval](i,j,k)*new_traces[i_interval][index];
          }
        }
        (*final_grid)(i,j,k) = static_cast<float>(value);
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
                                     CommonData                                                     * common_data,
                                     const Simbox                                                   * output_simbox)
{
  //Resample blocked logs to output_simbox and combine if multiple intervals
  //blocked_logs_output are blocked to output_simbox in commondata
  //Need to resample blocked logs that are added after commondata

  int n_intervals = static_cast<int>(blocked_logs_intervals.size());
  float res_fac   = 10.0; //Degree of refinement, must be integer.

  //Do not resample if there is only one interval, and if output_simbox and has the same resolution as the interval_simbox
  if (n_intervals_ == 1 && output_simbox->getnz() == multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
    blocked_logs_output = blocked_logs_intervals[0];
    return;
  }

  //Set seismic data directly instead of resampling from intervals
  std::vector<SeismicStorage *> seismic_data = common_data->GetSeismicDataTimeLapse(0);
  int n_angles = static_cast<int>(seismic_data.size());
  ModelAVODynamic::AddSeismicLogsFromStorage(blocked_logs_,
                                             seismic_data,
                                             *output_simbox,
                                             n_angles);

  //Loop over wells
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_output.begin(); it != blocked_logs_output.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_output.find(it->first);

    BlockedLogsCommon * blocked_log_final = iter->second;

    std::string well_name = blocked_log_final->GetWellName();
    int n_blocks          = blocked_log_final->GetNumberOfBlocks();

    LogKit::LogFormatted(LogKit::Low,"\n  "+well_name);

    //Predicted Logs
    bool got_predicted = blocked_logs_intervals[0].find(well_name)->second->GetContLogsPredicted().size() > 0;
    if (got_predicted) {
      LogKit::LogFormatted(LogKit::Low,"\n    predicted... ");

      std::vector<double> vp_predicted_final(n_blocks);
      std::vector<double> vs_predicted_final(n_blocks);
      std::vector<double> rho_predicted_final(n_blocks);
      std::vector<std::vector<double> > vp_predicted_intervals(n_intervals);
      std::vector<std::vector<double> > vs_predicted_intervals(n_intervals);
      std::vector<std::vector<double> > rho_predicted_intervals(n_intervals);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        CopyWellLog(vp_predicted_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpPredicted());
        CopyWellLog(vs_predicted_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVsPredicted());
        CopyWellLog(rho_predicted_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoPredicted());
      }

      ResampleLog(vp_predicted_final,  vp_predicted_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(vs_predicted_final,  vs_predicted_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_predicted_final, rho_predicted_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

      blocked_log_final->SetVpPredicted(vp_predicted_final);
      blocked_log_final->SetVsPredicted(vs_predicted_final);
      blocked_log_final->SetRhoPredicted(rho_predicted_final);
      LogKit::LogFormatted(LogKit::Low,"ok");
    }


    //Facies prob
    bool got_facies_prob = blocked_logs_intervals[0].find(well_name)->second->GetNFaciesProb() > 0;
    if (got_facies_prob) {
       LogKit::LogFormatted(LogKit::Low,"\n    facies prob... ");

      int n_faices = blocked_logs_intervals[0].find(well_name)->second->GetNFaciesProb();

      for (int j = 0; j < n_faices; j++) {

        std::vector<double> facies_prob_final(n_blocks);
        std::vector<std::vector<double> > facies_prob_intervals(n_intervals_);

        //Get well logs, missing values are interpolated
        for (int i = 0; i < n_intervals_; i++) {
          CopyWellLog(facies_prob_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRealSeismicData(j));
        }

        ResampleLog(facies_prob_final, facies_prob_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

        blocked_log_final->SetFaciesProb(j, facies_prob_final);
      }
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    //Cpp
    bool got_cpp = blocked_logs_intervals[0].find(well_name)->second->GetNCpp() > 0;
    if (got_cpp) {
      LogKit::LogFormatted(LogKit::Low,"\n    cpp... ");

      int n_angles = blocked_logs_intervals[0].find(well_name)->second->GetNCpp();

      for (int j = 0; j < n_angles; j++) {

        std::vector<double> cpp_final(n_blocks);
        std::vector<std::vector<double> > cpp_intervals(n_intervals_);

        //Get well logs, missing values are interpolated
        for (int i = 0; i < n_intervals_; i++) {
          CopyWellLog(cpp_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetCpp(j));
        }

        ResampleLog(cpp_final, cpp_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

        blocked_log_final->SetCpp(j, cpp_final);
      }
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    //ForFacies logs (vp, rho)
    bool got_vp_rho_fac_log = blocked_logs_intervals[0].find(well_name)->second->GetVpFaciesFiltered().size() > 0;
    if (got_vp_rho_fac_log) {
      LogKit::LogFormatted(LogKit::Low,"\n    vp rho fac log... ");
      std::vector<double> vp_for_facies_final(n_blocks);
      std::vector<double> rho_for_facies_final(n_blocks);
      std::vector<std::vector<double> > vp_for_facies_intervals(n_intervals_);
      std::vector<std::vector<double> > rho_for_facies_intervals(n_intervals_);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        CopyWellLog(vp_for_facies_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpFaciesFiltered());
        CopyWellLog(rho_for_facies_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoFaciesFiltered());
      }

      ResampleLog(vp_for_facies_final,  vp_for_facies_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_for_facies_final, rho_for_facies_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

      blocked_log_final->SetVpFaciesFiltered(vp_for_facies_final);
      blocked_log_final->SetRhoFaciesFiltered(rho_for_facies_final);
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    //Filtered logs
    bool got_filtered_logs = blocked_logs_intervals[0].find(well_name)->second->GetContLogsSeismicRes().size() > 0;
    if (got_filtered_logs) {
      LogKit::LogFormatted(LogKit::Low,"\n    filtered logs... ");
      std::vector<double> vp_filtered_final(n_blocks);
      std::vector<double> vs_filtered_final(n_blocks);
      std::vector<double> rho_filtered_final(n_blocks);
      std::vector<std::vector<double> > vp_filtered_intervals(n_intervals);
      std::vector<std::vector<double> > vs_filtered_intervals(n_intervals);
      std::vector<std::vector<double> > rho_filtered_intervals(n_intervals);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        CopyWellLog(vp_filtered_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpSeismicResolution());
        CopyWellLog(vs_filtered_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVsSeismicResolution());
        CopyWellLog(rho_filtered_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoSeismicResolution());
      }

      ResampleLog(vp_filtered_final,  vp_filtered_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(vs_filtered_final,  vs_filtered_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_filtered_final, rho_filtered_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

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
      LogKit::LogFormatted(LogKit::Low,"ok");
    }

    blocked_logs_output.insert(std::pair<std::string, BlockedLogsCommon *>(well_name, blocked_log_final));

  }
}

void CravaResult::CopyWellLog(std::vector<double>       & log_new,
                              const std::vector<double> & log_old)
{
  //for (int i = start; i <= end; i++) {
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
                              std::vector<std::vector<double> >                              & old_log_interval, //vector interval
                              const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                              MultiIntervalGrid                                              * multi_interval_grid,
                              const BlockedLogsCommon                                        * blocked_log_final,
                              std::string                                                      well_name,
                              float                                                            res_fac)
{
  std::vector<std::vector<double> > interval_logs_fine(n_intervals_); //resample to, per interval
  std::vector<std::vector<double> > z_pos_resampled_fine(n_intervals_); //Resample z-blocked log, it is used when we combine traces.

  //Get logs per interval and resample to a fine resolution
  for (int i_interval = 0; i_interval < n_intervals_; i_interval++) {

    int nz_interval           = static_cast<int>(old_log_interval[i_interval].size() * res_fac);
    std::vector<double> z_pos = blocked_logs_intervals[i_interval].find(well_name)->second->GetZposBlocked();

    interval_logs_fine[i_interval].resize(nz_interval);
    z_pos_resampled_fine[i_interval].resize(nz_interval);

    ResampleTrace(old_log_interval[i_interval], interval_logs_fine[i_interval], res_fac); //Interpolate missing
    ResampleTrace(z_pos, z_pos_resampled_fine[i_interval], res_fac);

  }

  //CombineTraces(final_log, blocked_log_final, blocked_logs_intervals, well_name, multi_interval_grid, interval_logs_fine, z_pos_resampled_fine);
  CombineTraces(final_log, blocked_log_final, multi_interval_grid, interval_logs_fine, z_pos_resampled_fine);
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
                                const std::vector<std::vector<double> > & interval_logs_fine,
                                const std::vector<std::vector<double> > & z_pos_resampled)
{
  int nz          = static_cast<int>(final_log.size());
  int n_intervals = static_cast<int>(interval_logs_fine.size());

  const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();

  for (int k = 0; k < nz; k++) {

    bool two_intervals = false;

    double global_x = blocked_log_final->GetXposBlocked()[k];
    double global_y = blocked_log_final->GetYposBlocked()[k];
    double global_z = blocked_log_final->GetZposBlocked()[k];

    int i_interval = 0;
    for (i_interval = 0; i_interval < n_intervals; i_interval++) {
      Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

      if (interval_simbox->IsInside(global_x, global_y, global_z))
        break;
    }
    if (i_interval < (n_intervals-1)) { //Also check if it hits the next interval, unless it is the last one.
      Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval+1);

      if (interval_simbox->IsInside(global_x, global_y, global_z))
        two_intervals = true;
    }

    double value        = 0.0;
    int interval_index  = i_interval;

    if (two_intervals == true) {
      //Use erorsion priorities to select between the two intervals
      if (erosion_priorities[i_interval] < erosion_priorities[i_interval+1])
        interval_index = i_interval;
      else
        interval_index = i_interval+1;
    }

    //value = GetResampledTraceValue(interval_logs_fine[i_interval], z_pos_resampled[i_interval], dz_resampled, top, global_z, dz_final);
    value = GetResampledTraceValue(interval_logs_fine[i_interval], z_pos_resampled[i_interval], global_z);

    final_log[k] = value;
  }
}


void CravaResult::WriteResults(ModelSettings           * model_settings,
                               CommonData              * common_data,
                               SeismicParametersHolder & seismic_parameters) //For crava-writing
{
  //Results are combined to one grid in CombineResults first
  //Wavelets are written out both in commonData and Wavelet1d/3d.cpp (and possibly modelavodynamic if match energies)
  //Estimation model: WriteEstimationResults

  const Simbox & simbox            = common_data->GetOutputSimbox();
  int output_grids_elastic         = model_settings->getOutputGridsElastic();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  //Write blocked wells
  if ((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Blocked Logs...");
    WriteBlockedWells(blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

    //Write blocked background logs (CRA-544). Logs that are blocked to extended background model (extended simbox with correlation direction).
    //Do not write if multiple intervals is used
    if (n_intervals_ == 1)
      WriteBlockedWells(bg_blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

    LogKit::LogFormatted(LogKit::Low,"ok\n");
  }
  //if ((model_settings->getWellOutputFlag() & IO::WELLS) > 0) {
  //  WriteWells(common_data->GetWells(), model_settings);
  //}

  if (model_settings->getWritePrediction() && !model_settings->getForwardModeling() && !model_settings->getEstimationMode()) {
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
      || (model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0 )) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Seismic Data\n");

    int n_timelapses = model_settings->getNumberOfTimeLapses();

    for (int i = 0; i < n_timelapses; i ++) {

      int n_angles              = model_settings->getNumberOfAngles(i);
      std::vector<float> angles = model_settings->getAngle(i);
      std::vector<float> offset = model_settings->getLocalSegyOffset(i);

      for (int j = 0; j < n_angles; j++) {
        std::string angle           = NRLib::ToString(angles[j]*(180/M_PI), 1);
        std::string file_name_orig  = IO::PrefixOriginalSeismicData() + angle;
        std::string sgri_label      = std::string("Original seismic data for angle stack ") + angle;
        if (offset.size() > 0 && offset[j] < 0)
          offset[j] = model_settings->getSegyOffset(i);

        std::string angle_synt     = NRLib::ToString(j);
        std::string file_name_synt = IO::makeFullFileName(IO::PathToSeismicData(), IO::FileTemporarySeismic()+angle_synt);

        int seismic_type = common_data->GetSeismicDataTimeLapse(i)[j]->GetSeismicType();
        int missing_traces_simbox  = 0;
        int missing_traces_padding = 0;
        int dead_traces_simbox     = 0;

        FFTGrid * fft_grid_resampled = NULL;
        bool delete_fft_grid         = true;
        if (seismic_type == 3) {
          fft_grid_resampled = common_data->GetSeismicDataTimeLapse(i)[j]->GetFFTGrid();
          delete_fft_grid = false;
        }
        else {
          SegY * segy                     = NULL;
          StormContGrid * storm           = NULL;
          NRLib::Grid<float> * nrlib_grid = new NRLib::Grid<float>();
          bool is_segy                    = false;
          bool is_storm                   = false;
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

          common_data->FillInData(nrlib_grid,
                                  fft_grid_resampled,
                                  &simbox,
                                  storm,
                                  segy,
                                  model_settings->getSmoothLength(),
                                  missing_traces_simbox,
                                  missing_traces_padding,
                                  dead_traces_simbox,
                                  FFTGrid::DATA,
                                  false,
                                  is_segy,
                                  is_storm,
                                  true);

          delete nrlib_grid;
        }

        if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0 && write_crava_) {
          std::string file_name_crava = IO::makeFullFileName(IO::PathToSeismicData(), IO::PrefixOriginalSeismicData() + angle);
          fft_grid_resampled->writeCravaFile(file_name_crava, &simbox);
        }

        StormContGrid * seismic_storm = CreateStormGrid(simbox, fft_grid_resampled); // deletes fft_grid_resampled

        //Real seismic gives value at cell base, synthetic at cell top. Shift real seismic.
        for (int k = static_cast<int>(seismic_storm->GetNK())-1; k > 0; k--) {
          for (size_t j = 0; j < seismic_storm->GetNJ(); j++) {
            for (size_t i = 0; i < seismic_storm->GetNI(); i++) {
              (*seismic_storm)(i,j,k) = (*seismic_storm)(i,j,k-1);
            }
          }
        }

        if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0)
          ParameterOutput::WriteFile(model_settings, seismic_storm, file_name_orig, IO::PathToSeismicData(), &simbox, true, sgri_label, offset[j], time_depth_mapping);

        if ((i==0) && (model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0) { //residuals only for first vintage.
          StormContGrid residual(*(synt_seismic_data_[j]));
          for (size_t k=0;k<seismic_storm->GetNK();k++) {
            for (size_t j=0;j<seismic_storm->GetNJ();j++) {
              for (size_t i=0;i<seismic_storm->GetNI();i++) {
                residual(i,j,k) = (*seismic_storm)(i,j,k)-residual(i,j,k);
              }
            }
          }

          sgri_label = "Residual computed from synthetic seismic for incidence angle "+angle;
          std::string file_name  = IO::PrefixSyntheticResiduals() + angle;

          ParameterOutput::WriteFile(model_settings, &residual, file_name, IO::PathToSeismicData(), &simbox, true, sgri_label);
        }
      }
    }
  }

  //Write Background models
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Background Grids\n");
    WriteBackgrounds(model_settings,
                      &simbox,
                      background_vp_,
                      background_vs_,
                      background_rho_,
                      time_depth_mapping,
                      *model_settings->getTraceHeaderFormatOutput());

    if (write_crava_) {
      std::string file_name_vp  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vp");
      std::string file_name_vs  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vs");
      std::string file_name_rho = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Rho");

      ExpTransf(background_vp_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vp, &simbox);

      ExpTransf(background_vs_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vs, &simbox);

      ExpTransf(background_rho_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_rho, &simbox);
    }
  }

  //Write correlations and post variances. Write per interval.
  for (int i = 0; i < n_intervals_; i++) {

    std::string interval_name = common_data->GetMultipleIntervalGrid()->GetIntervalName(i);

    //if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
    //  WriteFilePriorCorrT(corr_T_[i], simbox.GetNZpad(), dt, interval_name);

    //delete corr_T_[i];
    //fftw_free(corr_T_[i]);

    //if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0) {
    //  WriteFilePriorCorrT(corr_T_filtered_[i], simbox.GetNZpad(), dt, interval_name); // No zeros in the middle
    //  //delete [] corr_T_filtered_[i];
    //}

    if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
      LogKit::LogFormatted(LogKit::Low,"\nWrite Correlations\n");
      WriteFilePostVariances(post_var0_[i], post_cov_vp00_[i], post_cov_vs00_[i], post_cov_rho00_[i], interval_name);
      WriteFilePostCovGrids(model_settings, simbox, interval_name);

      if (write_crava_) {
        std::string file_name_vp    = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp");
        std::string file_name_vs    = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs");
        std::string file_name_rho   = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho");
        std::string file_name_vpvs  = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs");
        std::string file_name_vprho = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho");
        std::string file_name_vsrho = IO::makeFullFileName(IO::PathToCorrelations(), IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho");

        seismic_parameters.GetCovVp()->writeCravaFile(file_name_vp,         &simbox);
        seismic_parameters.GetCovVs()->writeCravaFile(file_name_vs,         &simbox);
        seismic_parameters.GetCovRho()->writeCravaFile(file_name_rho,       &simbox);
        seismic_parameters.GetCrCovVpVs()->writeCravaFile(file_name_vpvs,   &simbox);
        seismic_parameters.GetCrCovVpRho()->writeCravaFile(file_name_vprho, &simbox);
        seismic_parameters.GetCrCovVsRho()->writeCravaFile(file_name_vsrho, &simbox);
      }
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
          std::string file_name = IO::makeFullFileName(IO::PathToInversionResults(), IO::PrefixLikelihood() + facies_names[i]);
          seismic_parameters.GetLHCube()[i]->writeCravaFile(file_name, &simbox);
        }

      }
    }
  }

  //Simulations
  if (model_settings->getNumberOfSimulations() > 0) {
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
        seismic_parameters.GetSimulationSeed0(i)->writeCravaFile(file_name_vs,  &simbox);
        seismic_parameters.GetSimulationSeed0(i)->writeCravaFile(file_name_rho, &simbox);
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

      if (((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) ||
        (model_settings->getForwardModeling() == true))
        ParameterOutput::WriteFile(model_settings, synt_seismic_data_[i], file_name, IO::PathToSeismicData(), &simbox, true, sgri_label);
    }
  }

  //Trend Cubes
  if (model_settings->getOutputGridsOther() && IO::TREND_CUBES > 0) {

    const std::vector<std::string>  & trend_cube_parameters = model_settings->getTrendCubeParameters();

    for (size_t i = 0; i < trend_cubes_.size(); i++) {
      std::string file_name = IO::PrefixTrendCubes() + "_" + trend_cube_parameters[i];
      ParameterOutput::WriteFile(model_settings, trend_cubes_[i], file_name, IO::PathToSeismicData(), &simbox, "trend cube");
    }
  }
}


void CravaResult::WriteEstimationResults(ModelSettings * model_settings,
                                         CommonData    * common_data)
{
  //No CombineResults have been run, so much is unset, but we should be able to live with it.
  const Simbox & simbox            = common_data->GetOutputSimbox();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  blocked_logs_ = common_data->GetBlockedLogsOutput();

  n_intervals_ = common_data->GetMultipleIntervalGrid()->GetNIntervals();
  if (n_intervals_ == 1 && ((model_settings->getOutputGridFormat() & IO::CRAVA) > 0))
    write_crava_ = true;

  //Estimation model: All estimated parameters are written to file, regardless of output settings
  if (((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) &&
     (model_settings->getEstimateBackground() || model_settings->getEstimateWaveletNoise())) {
    LogKit::LogFormatted(LogKit::Low,"\nWrite Blocked Logs...");

    //Fill in seismic if needed.
    if (model_settings->getEstimateWaveletNoise()) {
      //ComputeSeismicImpedance from avoinversion.cpp
      reflection_matrix_ = common_data->GetReflectionMatrixTimeLapse(0);
      wavelets_          = common_data->GetWavelet(0);

      std::vector<bool> estimate_wavelet = model_settings->getEstimateWavelet(0);
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

      GenerateSyntheticSeismicLogs(wavelets_, blocked_logs_, reflection_matrix_, simbox);

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
      CombineResult(background_vp_,  dummy_fft_grids,  multi_interval_grid, zone_prob_grid, background_vp_intervals);
      LogKit::LogFormatted(LogKit::Low,"\n Vs");
      CombineResult(background_vs_,  dummy_fft_grids,  multi_interval_grid, zone_prob_grid, background_vs_intervals);
      LogKit::LogFormatted(LogKit::Low,"\n Rho");
      CombineResult(background_rho_, dummy_fft_grids, multi_interval_grid, zone_prob_grid, background_vs_intervals);

      for (int i = 0; i < n_intervals_; i++) {
        common_data->ReleaseBackgroundGrids(i, 0);
        common_data->ReleaseBackgroundGrids(i, 1);
        common_data->ReleaseBackgroundGrids(i, 2);
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
    }

    LogKit::LogFormatted(LogKit::Low,"\nWrite Background Grids\n");
    WriteBackgrounds(model_settings,
                     &simbox,
                     background_vp_,
                     background_vs_,
                     background_rho_,
                     time_depth_mapping,
                     *model_settings->getTraceHeaderFormatOutput());

    if (write_crava_) {
      std::string file_name_vp  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vp");
      std::string file_name_vs  = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Vs");
      std::string file_name_rho = IO::makeFullFileName(IO::PathToBackground(), IO::PrefixBackground() + "Rho");

      ExpTransf(background_vp_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vp, &simbox);

      ExpTransf(background_vs_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_vs, &simbox);

      ExpTransf(background_rho_intervals_[0]);
      background_vp_intervals_[0]->writeCravaFile(file_name_rho, &simbox);
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

void CravaResult::WriteFilePostCovGrids(const ModelSettings * model_settings,
                                        const Simbox        & simbox,
                                        std::string           interval_name) const
{
  if (interval_name != "")
    interval_name = "_" + interval_name;

  std::string file_name;
  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp" + interval_name;
  ParameterOutput::WriteFile(model_settings, cov_vp_, file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for Vp");

  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs" + interval_name;
  ParameterOutput::WriteFile(model_settings, cov_vs_, file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for Vs");

  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho" + interval_name;
  ParameterOutput::WriteFile(model_settings, cov_rho_, file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for density");

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs" + interval_name;
  ParameterOutput::WriteFile(model_settings, cr_cov_vp_vs_, file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vp,Vs)");

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho" + interval_name;
  ParameterOutput::WriteFile(model_settings, cr_cov_vp_rho_, file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vp,density)");

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho" + interval_name;
  ParameterOutput::WriteFile(model_settings, cr_cov_vs_rho_, file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vs,density)");
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

//StormContGrid *
//CravaResult::CreateStormGrid(const Simbox       & simbox,
//                             NRLib::Grid<float> * grid)
//{
//  if (grid != NULL) {
//    int nx = grid->GetNI();
//    int ny = grid->GetNJ();
//    int nz = grid->GetNK();
//
//    StormContGrid * storm = new StormContGrid(simbox, nx, ny, nz);
//
//    for (int i = 0; i < nx; i++) {
//      for (int j = 0; j < ny; j++) {
//        for (int k = 0; k < nz; k++) {
//          float value = grid->GetValue(i, j, k);
//          storm->SetValue(i, j, k, value);
//        }
//      }
//    }
//
//    //Not needed anymore
//    delete grid;
//
//    return(storm);
//  }
//
//  StormContGrid * storm = new StormContGrid();
//  return(storm);
//
//}

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

void CravaResult::WriteBackgrounds(const ModelSettings     * model_settings,
                                   const Simbox            * simbox,
                                   StormContGrid           * background_vp,
                                   StormContGrid           * background_vs,
                                   StormContGrid           * background_rho,
                                   GridMapping             * depth_mapping,
                                   const TraceHeaderFormat & thf)
{
  if (depth_mapping != NULL && depth_mapping->getSimbox() == NULL) {
    int output_format = model_settings->getOutputGridFormat();
    if (model_settings->getWriteAsciiSurfaces() && !(output_format & IO::ASCII))
      output_format += IO::ASCII;

    depth_mapping->setMappingFromVelocity(background_vp, simbox, model_settings->getOutputGridFormat());
  }

  std::string file_name_vp  = IO::PrefixBackground() + "Vp" ;
  std::string file_name_vs  = IO::PrefixBackground() + "Vs" ;
  std::string file_name_rho = IO::PrefixBackground() + "Rho";

  ExpTransf(background_vp_);
  ParameterOutput::WriteFile(model_settings, background_vp, file_name_vp, IO::PathToBackground(), simbox, false, "NO_LABEL", 0, depth_mapping, thf);

  ExpTransf(background_vs_);
  ParameterOutput::WriteFile(model_settings, background_vs, file_name_vs, IO::PathToBackground(), simbox, false, "NO_LABEL", 0, depth_mapping, thf);

  ExpTransf(background_rho_);
  ParameterOutput::WriteFile(model_settings, background_rho, file_name_rho, IO::PathToBackground(), simbox, false, "NO_LABEL", 0, depth_mapping, thf);

  //
  // For debugging: write cubes not in ASCII, with padding, and with flat top.
  //
  //backModel_[0]->writeStormFile(fileName1, simbox, true, false, true, true);
  //backModel_[1]->writeStormFile(fileName2, simbox, true, false, true, true);
  //backModel_[2]->writeStormFile(fileName3, simbox, true, false, true, true);
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

        for (int k = 0; k <(nzp/2 +1); k++) {
          fftw_complex r = resultVec.getCAmp(k);
          fftw_complex w = localWavelet->getCAmp(k,static_cast<float>(sf));// returns complex conjugate
          fftw_complex s;
          s.re = r.re*w.re+r.im*w.im; //Use complex conjugate of w
          s.im = -r.re*w.im+r.im*w.re;
          resultVec.setCAmp(s,k);
        }
        delete localWavelet;

        resultVec.invFFT1DInPlace();
        for (int k = 0; k < nz; k++) { //nzp
          float value = resultVec.getRAmp(k);
          imp->SetValue(i, j, k, value);
        }
      }
    }

    if (((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) ||
      (model_settings->getForwardModeling() == true))
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

void CravaResult::GenerateSyntheticSeismicLogs(std::vector<Wavelet *>                     & wavelet,
                                               std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                               const NRLib::Matrix                        & reflection_matrix,
                                               const Simbox                               & simbox)
{
  int nzp = simbox.GetNZpad();
  int nz  = simbox.getnz();

  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log->GetIsDeviated() == false)
      blocked_log->GenerateSyntheticSeismic(reflection_matrix, wavelet, nz, nzp, &simbox);
  }
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
    for (size_t w=0;w<blocked_wells.size();w++) {
      delete wavelets[w][i];
    }
  }
}

 void
 CravaResult::DownscaleTrace(const std::vector<float> & trace_in,
                             std::vector<float>       & trace_out,
                             int                         scale)
{
  std::vector<float> input_pad = trace_in;
  size_t pad = 101; //Want odd length vector in fft, better for interpolation.
  if(trace_in.size() % 2 == 1)
    pad = 100;
  size_t n_data = trace_in.size();
  input_pad.resize(n_data+pad);
  double step = 1.0/static_cast<float>(pad);
  for(size_t i=n_data; i<input_pad.size();i++) {
    double t = (i-trace_in.size())*step;
    input_pad[i] = static_cast<float>((1-t)*trace_in[n_data-1]+t*trace_in[0]);
  }

  int nt = static_cast<int>(input_pad.size());
  int mt = nt*scale;
  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  int cnt = nt/2 + 1;
  int rnt = 2*cnt;
  int cmt = mt/2 + 1;
  int rmt = 2*cmt;

  fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
  fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

  CommonData::ResampleTrace(input_pad,
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
  size_t out_len = trace_in.size()*scale-(scale-1); //Points beyond this are contaminated by padding.
  trace_out.resize(out_len);
  for (size_t k = 0; k < out_len; k++) {
    trace_out[k] = rAmpFine[k];
  }

  fftw_free(rAmpData);
  fftw_free(rAmpFine);
  fftwnd_destroy_plan(fftplan1);
  fftwnd_destroy_plan(fftplan2);
}
