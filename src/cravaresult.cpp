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
//reflection_matrix_(NULL),
n_intervals_(0)
{
}

CravaResult::~CravaResult()
{
  //
  for (size_t i = 0; i < corr_T_.size(); i++){
    if (corr_T_[i] != NULL){
      delete [] corr_T_[i];
    }
  }
  for (size_t i = 0; i < corr_T_filtered_.size(); i++){
    if (corr_T_filtered_[i] != NULL)
      delete [] corr_T_filtered_[i];
  }
    /*
  for (size_t i = 0; i < background_vp_intervals_.size(); i++){
    if (background_vp_intervals_[i] != NULL)
      delete background_vp_intervals_[i];
  }
  for (size_t i = 0; i < background_vs_intervals_.size(); i++){
    if (background_vs_intervals_[i] != NULL)
      delete background_vs_intervals_[i];
  }
  for (size_t i = 0; i < background_rho_intervals_.size(); i++){
    if (background_rho_intervals_[i] != NULL)
      delete background_rho_intervals_[i];
  }

  if (cov_vp_ != NULL){
    delete cov_vp_;
    cov_vp_ = NULL;
  }
  if (cov_vs_ != NULL){
    delete cov_vs_;
    cov_vs_ = NULL;
  }
  if (cov_rho_ != NULL){
    delete cov_rho_;
    cov_rho_ = NULL;
  }
  if (cr_cov_vp_vs_ != NULL){
    delete cr_cov_vp_vs_;
    cr_cov_vp_vs_ = NULL;
  }
  if (cr_cov_vp_rho_ != NULL){
    delete cr_cov_vp_rho_;
    cr_cov_vp_rho_ = NULL;
  }
  if (cr_cov_vs_rho_ != NULL){
    delete cr_cov_vs_rho_;
    cr_cov_vs_rho_ = NULL;
  }
  if (post_vp_ != NULL){
    delete post_vp_;
    post_vp_ = NULL;
  }
  if (post_vs_ != NULL){
    delete post_vs_;
    post_vs_ = NULL;
  }
  if (post_rho_ != NULL){
    delete post_rho_;
    post_rho_ = NULL;
  }
  if (post_vp_kriged_ != NULL){
    delete post_vp_kriged_;
    post_vp_kriged_ = NULL;
  }
  if (post_vs_kriged_ != NULL){
    delete post_vs_kriged_;
    post_vs_kriged_ = NULL;
  }
  if (post_rho_kriged_ != NULL){
    delete post_rho_kriged_;
    post_rho_kriged_ = NULL;
  }
  if (background_vp_ != NULL){
    delete background_vp_;
    background_vp_ = NULL;
  }
  if (background_vs_ != NULL){
    delete background_vs_;
    background_vs_ = NULL;
  }
  if (background_rho_ != NULL){
    delete background_rho_;
    background_rho_ = NULL;
  }
  for (size_t i = 0; i < simulations_seed0_.size(); i++){
    delete simulations_seed0_[i];
  }
  for (size_t i = 0; i < simulations_seed1_.size(); i++){
    delete simulations_seed1_[i];
  }
  for (size_t i = 0; i < simulations_seed2_.size(); i++){
    delete simulations_seed2_[i];
  }
  for (size_t i = 0; i < synt_seismic_data_.size(); i++){
    delete synt_seismic_data_[i];
  }
  for (size_t i = 0; i < synt_residuals_.size(); i++){
    delete synt_residuals_[i];
  }
  if (block_grid_ != NULL){
    delete block_grid_;
    block_grid_ = NULL;
  }
  for (size_t i = 0; i < facies_prob_.size(); i++){
    delete facies_prob_[i];
  }
  if (facies_prob_undef_!= NULL){
    delete facies_prob_undef_;
    facies_prob_undef_ = NULL;
  }
  for (size_t i = 0; i < facies_prob_geo_.size(); i++){
    delete facies_prob_geo_[i];
  }
  for (size_t i = 0; i < lh_cubes_.size(); i++){
    delete lh_cubes_[i];
  }
  if (quality_grid_!= NULL){
    delete quality_grid_;
    quality_grid_ = NULL;
  }
  
  for (size_t i = 0; i < wavelets_.size(); i++){
    delete reflection_matrix_[i];
    //delete wavelets_[i];
  }
  delete reflection_matrix_;
  */


}

void CravaResult::CombineResults(ModelSettings                        * model_settings,
                                 CommonData                           * common_data,
                                 std::vector<SeismicParametersHolder> & seismic_parameters_intervals)
{
  //Combine interval grids to one stormgrid
  //Compaction: Use the finest dz between intervals

  //H-TODO: CombineResult is not done
  //  There is no smoothing between intervals


  MultiIntervalGrid * multi_interval_grid = common_data->GetMultipleIntervalGrid();
  int n_intervals = multi_interval_grid->GetNIntervals();

  const Simbox & simbox = common_data->GetFullInversionSimbox();

  if (n_intervals > 1) {

   const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();
   int nx = simbox.getnx();
   int ny = simbox.getny();

    //Find minimum dz from all traces
    double dz_min = 0.0;
    FindDzMin(simbox, multi_interval_grid, dz_min);

    int nz_max = static_cast<int>((simbox.getBotZMax() - simbox.getTopZMin()) / dz_min);

    //H-TEST
    //nz_max = nz_max + 50;



    if (model_settings->getWritePrediction()) {

      //Post vp, vs and rho from avoinversion computePostMeanResidAndFFTCov()
      post_vp_  = new StormContGrid(simbox, nx, ny, nz_max);
      post_vs_  = new StormContGrid(simbox, nx, ny, nz_max);
      post_rho_ = new StormContGrid(simbox, nx, ny, nz_max);

      std::vector<FFTGrid *> post_vp_intervals(n_intervals);
      std::vector<FFTGrid *> post_vs_intervals(n_intervals);
      std::vector<FFTGrid *> post_rho_intervals(n_intervals);
      for (int i = 0; i < n_intervals; i++) {
        post_vp_intervals[i]  = seismic_parameters_intervals[i].GetPostVp();
        post_vs_intervals[i]  = seismic_parameters_intervals[i].GetPostVs();
        post_rho_intervals[i] = seismic_parameters_intervals[i].GetPostRho();
      }
      CombineResult(post_vp_,  post_vp_intervals,  multi_interval_grid, erosion_priorities, dz_min);
      CombineResult(post_vs_,  post_vs_intervals,  multi_interval_grid, erosion_priorities, dz_min);
      CombineResult(post_rho_, post_rho_intervals, multi_interval_grid, erosion_priorities, dz_min);

      //Post vp, vs and rho from avoinversion doPredictionKriging()
      if (model_settings->getKrigingParameter() > 0) {
        post_vp_kriged_  = new StormContGrid(simbox, nx, ny, nz_max);
        post_vs_kriged_  = new StormContGrid(simbox, nx, ny, nz_max);
        post_rho_kriged_ = new StormContGrid(simbox, nx, ny, nz_max);

        std::vector<FFTGrid *> post_vp_kriged_intervals(n_intervals);
        std::vector<FFTGrid *> post_vs_kriged_intervals(n_intervals);
        std::vector<FFTGrid *> post_rho_kriged_intervals(n_intervals);
        for (int i = 0; i < n_intervals; i++) {
          post_vp_kriged_intervals[i]  = seismic_parameters_intervals[i].GetPostVpKriged();
          post_vs_kriged_intervals[i]  = seismic_parameters_intervals[i].GetPostVsKriged();
          post_rho_kriged_intervals[i] = seismic_parameters_intervals[i].GetPostRhoKriged();
        }
        CombineResult(post_vp_kriged_,  post_vp_kriged_intervals,  multi_interval_grid, erosion_priorities, dz_min);
        CombineResult(post_vs_kriged_,  post_vs_kriged_intervals,  multi_interval_grid, erosion_priorities, dz_min);
        CombineResult(post_rho_kriged_, post_rho_kriged_intervals, multi_interval_grid, erosion_priorities, dz_min);
      }
    }

    //Background models
    if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
      background_vp_  = new StormContGrid(simbox, nx, ny, nz_max);
      background_vs_  = new StormContGrid(simbox, nx, ny, nz_max);
      background_rho_ = new StormContGrid(simbox, nx, ny, nz_max);

      CombineResult(background_vp_,  background_vp_intervals_,  multi_interval_grid, erosion_priorities, dz_min);
      CombineResult(background_vs_,  background_vs_intervals_,  multi_interval_grid, erosion_priorities, dz_min);
      CombineResult(background_rho_, background_rho_intervals_, multi_interval_grid, erosion_priorities, dz_min);
    }

    //Covariance grids
    cov_vp_        = new StormContGrid(simbox, nx, ny, nz_max);
    cov_vs_        = new StormContGrid(simbox, nx, ny, nz_max);
    cov_rho_       = new StormContGrid(simbox, nx, ny, nz_max);
    cr_cov_vp_vs_  = new StormContGrid(simbox, nx, ny, nz_max);
    cr_cov_vp_rho_ = new StormContGrid(simbox, nx, ny, nz_max);
    cr_cov_vs_rho_ = new StormContGrid(simbox, nx, ny, nz_max);

    std::vector<FFTGrid *> cov_vp_intervals(n_intervals);
    std::vector<FFTGrid *> cov_vs_intervals(n_intervals);
    std::vector<FFTGrid *> cov_rho_intervals(n_intervals);
    std::vector<FFTGrid *> cr_cov_vp_vs_intervals(n_intervals);
    std::vector<FFTGrid *> cr_cov_vp_rho_intervals(n_intervals);
    std::vector<FFTGrid *> cr_cov_vs_rho_intervals(n_intervals);

    for (int i = 0; i < n_intervals; i++) {
      cov_vp_intervals[i]        = seismic_parameters_intervals[i].GetCovVp();
      cov_vs_intervals[i]        = seismic_parameters_intervals[i].GetCovVs();
      cov_rho_intervals[i]       = seismic_parameters_intervals[i].GetCovRho();
      cr_cov_vp_vs_intervals[i]  = seismic_parameters_intervals[i].GetCrCovVpVs();
      cr_cov_vp_rho_intervals[i] = seismic_parameters_intervals[i].GetCrCovVpRho();
      cr_cov_vs_rho_intervals[i] = seismic_parameters_intervals[i].GetCrCovVsRho();
    }
    CombineResult(cov_vp_,        cov_vp_intervals,        multi_interval_grid, erosion_priorities, dz_min);
    CombineResult(cov_vs_,        cov_vs_intervals,        multi_interval_grid, erosion_priorities, dz_min);
    CombineResult(cov_rho_,       cov_rho_intervals,       multi_interval_grid, erosion_priorities, dz_min);
    CombineResult(cr_cov_vp_vs_,  cr_cov_vp_vs_intervals,  multi_interval_grid, erosion_priorities, dz_min);
    CombineResult(cr_cov_vp_rho_, cr_cov_vp_rho_intervals, multi_interval_grid, erosion_priorities, dz_min);
    CombineResult(cr_cov_vs_rho_, cr_cov_vs_rho_intervals, multi_interval_grid, erosion_priorities, dz_min);

    //Facies prob
    if (model_settings->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF) {
      int n_facies = seismic_parameters_intervals[0].GetFaciesProb().size();
      facies_prob_.resize(n_facies);

      for (int j = 0; j < n_facies; j++) {
        facies_prob_[j] = new StormContGrid(simbox, nx, ny, nz_max);

        std::vector<FFTGrid *> facies_prob_intervals(n_intervals);
        for (int i = 0; i < n_intervals; i++) {
          facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProb()[j];
        }
        CombineResult(facies_prob_[j],  facies_prob_intervals,  multi_interval_grid, erosion_priorities, dz_min);
      }

      //Undef
      std::vector<FFTGrid *> facies_prob_intervals_undef(n_intervals);
      for (int i = 0; i < n_intervals; i++) {
        facies_prob_intervals_undef[i] = seismic_parameters_intervals[i].GetFaciesProbUndefined();
      }
      CombineResult(facies_prob_undef_,  facies_prob_intervals_undef,  multi_interval_grid, erosion_priorities, dz_min);

    }
    if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
      int n_facies = seismic_parameters_intervals[0].GetFaciesProbGeomodel().size();
      facies_prob_.resize(n_facies);

      for (int j = 0; j < n_facies; j++) {
        facies_prob_geo_[j] = new StormContGrid(simbox, nx, ny, nz_max);

        std::vector<FFTGrid *> facies_prob_intervals(n_intervals);
        for (int i = 0; i < n_intervals; i++) {
          facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProbGeomodel()[j];
        }
        CombineResult(facies_prob_geo_[j],  facies_prob_intervals,  multi_interval_grid, erosion_priorities, dz_min);
      }
    }

    //LH Cube
    if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
      int n_facies = seismic_parameters_intervals[0].GetLHCube().size();
      lh_cubes_.resize(n_facies);

      for (int j = 0; j < n_facies; j++) {
        lh_cubes_[j] = new StormContGrid(simbox, nx, ny, nz_max);

        std::vector<FFTGrid *> lh_cubes_intervals(n_intervals);
        for (int i = 0; i < n_intervals; i++) {
          lh_cubes_intervals[i] = seismic_parameters_intervals[i].GetLHCube()[j];
        }
        CombineResult(lh_cubes_[j], lh_cubes_intervals, multi_interval_grid, erosion_priorities, dz_min);

      }
    }

    //Quality grid
    if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID) {
      quality_grid_ = new StormContGrid(simbox, nx, ny, nz_max);

      std::vector<FFTGrid *> quality_grid_intervals(n_intervals);
      for (int i = 0; i < n_intervals; i++) {
        quality_grid_intervals[i] = seismic_parameters_intervals[i].GetQualityGrid();
      }
      CombineResult(quality_grid_, quality_grid_intervals, multi_interval_grid, erosion_priorities, dz_min);
    }

    //Simulation grids
    if (model_settings->getNumberOfSimulations() > 0) {
      int n_simulations = seismic_parameters_intervals[0].GetSimulationsSeed0().size();

      for (int j = 0; j < n_simulations; j++) {
        simulations_seed0_.resize(n_simulations);
        simulations_seed1_.resize(n_simulations);
        simulations_seed2_.resize(n_simulations);

        std::vector<FFTGrid *> simulations_seed0_intervals(n_intervals);
        std::vector<FFTGrid *> simulations_seed1_intervals(n_intervals);
        std::vector<FFTGrid *> simulations_seed2_intervals(n_intervals);

        for (int i = 0; i < n_intervals; i++) {
          simulations_seed0_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed0(j);
          simulations_seed1_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed1(j);
          simulations_seed2_intervals[i] = seismic_parameters_intervals[i].GetSimulationSeed2(j);
        }

        CombineResult(simulations_seed0_[j], simulations_seed0_intervals, multi_interval_grid, erosion_priorities, dz_min);
        CombineResult(simulations_seed1_[j], simulations_seed1_intervals, multi_interval_grid, erosion_priorities, dz_min);
        CombineResult(simulations_seed2_[j], simulations_seed2_intervals, multi_interval_grid, erosion_priorities, dz_min);

      }
    }

    //Block grid
    if (model_settings->getWritePrediction() && model_settings->getDebugFlag()) {
      block_grid_ = new StormContGrid(simbox, nx, ny, nz_max);

      std::vector<FFTGrid *> block_grid_intervals(n_intervals);
      for (int i = 0; i < n_intervals; i++) {
        block_grid_intervals[i] = seismic_parameters_intervals[i].GetBlockGrid();
      }
      CombineResult(block_grid_, block_grid_intervals, multi_interval_grid, erosion_priorities, dz_min);
    }

    //Correlations and post variances
    //Do not combine, store and write per interval
    for (int i = 0; i < n_intervals; i++) {
      corr_T_.push_back(seismic_parameters_intervals[i].GetCorrT());
      if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
        corr_T_filtered_.push_back(seismic_parameters_intervals[i].GetCorrTFiltered());

      post_var0_.push_back(seismic_parameters_intervals[i].GetPostVar0());
      post_cov_vp00_.push_back(seismic_parameters_intervals[i].GetPostCovVp00());
      post_cov_vs00_.push_back(seismic_parameters_intervals[i].GetPostCovVs00());
      post_cov_rho00_.push_back(seismic_parameters_intervals[i].GetPostCovRho00());
    }


  }
  else {

    //Final grids are stored as StormContGrids
    // Results in seismic_parameters are stored as fftgrid

    //CreateStormGrid() creates a storm grid from fft grid, and deletes the fft_grid


    std::string interval_name = multi_interval_grid->GetIntervalName(0);

    post_vp_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVp());
    post_vs_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVs());
    post_rho_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostRho());

    if (model_settings->getForwardModeling()) {
      post_vp_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetMeanVp());
      post_vs_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetMeanVs());
      post_rho_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetMeanRho());
    }

    post_vp_kriged_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVpKriged());
    post_vs_kriged_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVsKriged());
    post_rho_kriged_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostRhoKriged());

    background_vp_  = CreateStormGrid(simbox, background_vp_intervals_[0]);
    background_vs_  = CreateStormGrid(simbox, background_vs_intervals_[0]);
    background_rho_ = CreateStormGrid(simbox, background_rho_intervals_[0]);

    corr_T_.push_back(seismic_parameters_intervals[0].GetCorrT());
    if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      corr_T_filtered_.push_back(seismic_parameters_intervals[0].GetCorrTFiltered());

    post_var0_.push_back(seismic_parameters_intervals[0].GetPostVar0());
    post_cov_vp00_.push_back(seismic_parameters_intervals[0].GetPostCovVp00());
    post_cov_vs00_.push_back(seismic_parameters_intervals[0].GetPostCovVs00());
    post_cov_rho00_.push_back(seismic_parameters_intervals[0].GetPostCovRho00());

    cov_vp_        = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCovVp());
    cov_vs_        = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCovVs());
    cov_rho_       = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCovRho());
    cr_cov_vp_vs_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCrCovVpVs());
    cr_cov_vp_rho_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCrCovVpRho());
    cr_cov_vs_rho_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetCrCovVsRho());

    int n_simulations = seismic_parameters_intervals[0].GetSimulationsSeed0().size();
    simulations_seed0_.resize(n_simulations);
    simulations_seed1_.resize(n_simulations);
    simulations_seed2_.resize(n_simulations);
    for (int i = 0; i < n_simulations; i++) {
      simulations_seed0_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetSimulationsSeed0()[i]);
      simulations_seed1_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetSimulationsSeed1()[i]);
      simulations_seed2_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetSimulationsSeed2()[i]);
    }

    block_grid_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetBlockGrid());

    int n_facies = seismic_parameters_intervals[0].GetFaciesProb().size();
    facies_prob_.resize(n_facies);
    for (int i = 0; i < n_facies; i++) {
      facies_prob_[i]    = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetFaciesProb()[i]);
    }
    facies_prob_undef_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetFaciesProbUndefined());

    n_facies = seismic_parameters_intervals[0].GetFaciesProbGeomodel().size();
    facies_prob_geo_.resize(n_facies);
    for (int i = 0; i < n_facies; i++) {
      facies_prob_geo_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetFaciesProbGeomodel()[i]);
    }


    int n_lh_cubes = seismic_parameters_intervals[0].GetLHCube().size();
    lh_cubes_.resize(n_lh_cubes);
    for (int i = 0; i < n_lh_cubes; i++) {
      lh_cubes_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetLHCube()[i]);
    }

    quality_grid_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetQualityGrid());

  }


  //Compute Synt seismic

  //H ComputeSeismicImpedance needs reflection matrix, but in avoinversion the reflection matrix is per interval
  reflection_matrix_ = common_data->GetReflectionMatrixTimeLapse(0);
  wavelets_          = common_data->GetWavelet(0);

  if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {
    if (model_settings->getKrigingParameter() > 0)
      ComputeSyntSeismic(model_settings, &simbox, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_);
    else
      ComputeSyntSeismic(model_settings, &simbox, post_vp_, post_vs_, post_rho_);
  }


  //seismic_parameters_intervals[0].releaseGrids();

}

void CravaResult::CombineResult(StormContGrid         *& final_grid,
                                std::vector<FFTGrid *> & interval_grids,
                                MultiIntervalGrid      * multi_interval_grid,
                                const std::vector<int> & erosion_priorities,
                                double                   dz_min)
{
  int n_intervals = multi_interval_grid->GetNIntervals();
  int nx = final_grid->GetNI();
  int ny = final_grid->GetNJ();
  int nz = final_grid->GetNK();

  bool writing = false;

  //Resample
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      std::vector<std::vector<float> > new_traces(n_intervals);

      //Resample each trace to new nz
      for (int i_interval = 0; i_interval < n_intervals; i_interval++) {

        Simbox * interval_simbox       = multi_interval_grid->GetIntervalSimbox(i_interval);
        FFTGrid * interval_grid        = interval_grids[i_interval];
        //FFTGrid * cov_vp_interval = seismic_parameters_intervals[i_interval].GetPostVp();
        std::vector<float> old_trace = interval_grid->getRealTrace2(i, j); //old_trace is changed below.

        double top_value = interval_simbox->getTop(i,j);
        double bot_value = interval_simbox->getBot(i,j);

        int nz_old       = interval_simbox->getnz();
        //double dz_old    = (bot_value - top_value) / nz_old;

        //Resample to new nz based on minimum dz from all traces and intervals
        int nz_new = static_cast<int>((bot_value - top_value) / dz_min);

        //H-TEST For resampling
        //nz_new = nz_new + 50;

        new_traces[i_interval].resize(nz_new);

        //H-Writing
        if (writing) {
          std::string file_name = "combine/old_trace";
          std::ofstream file;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
                << std::setprecision(6);
          for (size_t k = 0; k < old_trace.size(); k++)
            file << old_trace[k] << " ";
          file.close();
        }

        //Resampling
        //FFT(trace) -> pad with zeros until new nz -> IFFT back
        //std::vector<float> trace_new_simple(nz_new);
        //ResampleSimple(trace_new_simple, old_trace); //H Not working properly. Oscillates at the beginning and end.

        //Alternative: Resampling of traces from CommonData::FillInData
        //Remove trend -> fft -> pad with zeroes -> resample -> ifft -> add trend

        //Remove trend from trace
        size_t n_trace    = old_trace.size();
        float trend_first = old_trace[0];
        float trend_last  = old_trace[n_trace - 1];
        float trend_inc   = (trend_last - trend_first) / (n_trace - 1);
        for (size_t k_trace = 0; k_trace < old_trace.size(); k_trace++) {
          old_trace[k_trace] -= trend_first + k_trace * trend_inc;
        }

        //int nt = CommonData::FindClosestFactorableNumber(nz_new);
        //int mt = 4*nt;

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
                                  cnt,
                                  rnt,
                                  cmt,
                                  rmt);

        //Add trend
        trend_inc = (trend_last - trend_first) / (nz_new - 1);
        for (int k = 0; k < nz_new; k++) {
          new_traces[i_interval][k] = rAmpFine[k] + (trend_first + k*trend_inc);
        }

        //H-Writing REMOVE
        if (writing) {
          //file_name = "combine/new_trace_simple";
          //NRLib::OpenWrite(file, file_name);
          //file << std::fixed
          //      << std::setprecision(6);
          //for (size_t k = 0; k < trace_new_simple.size(); k++)
          //  file << trace_new_simple[k] << " ";
          //file.close();

          std::string file_name = "combine/old_trace_without_trend";
          std::ofstream file;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
                << std::setprecision(6);
          for (size_t k = 0; k < old_trace.size(); k++)
            file << old_trace[k] << " ";
          file.close();

          file_name = "combine/rAmpFine";
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
                << std::setprecision(6);
          for (int k = 0; k < rmt; k++)
            file << rAmpFine[k] << " ";
          file.close();

          file_name = "combine/new_trace_resampled";
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
                << std::setprecision(6);
          for (size_t k = 0; k < new_traces[i_interval].size(); k++)
            file << new_traces[i_interval][k] << " ";
          file.close();
        }


      } //n_intervals

      //Combine vectors for each interval to one trace in stormgrid
      std::vector<int> interval_indexes(n_intervals, 0);

      for (int k = 0; k < nz; k++) {

        bool two_intervals = false;

        double tmp_x = 0.0;
        double tmp_y = 0.0;
        double tmp_z = 0.0;

        final_grid->FindCenterOfCell(i, j, k, tmp_x, tmp_y, tmp_z);

        int i_interval = 0;
        for (i_interval = 0; i_interval < n_intervals; i_interval++) {
          Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

          if (interval_simbox->IsInside(tmp_x, tmp_y, tmp_z))
            break;
        }
        if (i_interval < (n_intervals-1)) { //Also check if it hits the next interval, unless it is the last one.
          Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval+1);

          if (interval_simbox->IsInside(tmp_x, tmp_y, tmp_z))
            two_intervals = true;
        }

        float value = 0.0f;
        int index   = 0;
        if (two_intervals == true) {
          //Use erorsion priorities to select between the two intervals

          if (erosion_priorities[i_interval] < erosion_priorities[i_interval+1]) {
            index = interval_indexes[i_interval];
            value = new_traces[i_interval][index];
          }
          else {
            index = interval_indexes[i_interval+1];
            value = new_traces[i_interval+1][index];
          }

          interval_indexes[i_interval]++;
          interval_indexes[i_interval+1]++;
        }
        else {

          int index = interval_indexes[i_interval];
          value = new_traces[i_interval][index];
          interval_indexes[i_interval]++;
        }

        final_grid->SetValue(i, j, k, value);
      }

    } //ny
  } //nx
}


void CravaResult::WriteResults(ModelSettings * model_settings,
                               CommonData    * common_data)
{
  //Move writing rutines from modelavodynamic and avoinversion here

  //Results are combined to one grid in CombineResults first
  const Simbox & simbox = common_data->GetFullInversionSimbox();

  float dt = static_cast<float>(simbox.getdz());

  //CRA-544: Include blocked background logs when writing blocked logs to file

  //Wavelets are written out both in commonData and modelavodynamic


  int output_grids_elastic = model_settings->getOutputGridsElastic();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  if (model_settings->getEstimationMode()) { //Estimation model: All estimated parameters are written to file, regardless of output settings

     WriteBlockedWells(common_data->GetBlockedLogs(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

    WriteWells(common_data->GetWells(), model_settings);

    WriteBackgrounds(model_settings,
                      &simbox,
                      time_depth_mapping,
                      *model_settings->getTraceHeaderFormatOutput());


    //H More?


  }
  else {

    if (model_settings->getWritePrediction()) {

      //From computePostMeanResidAndFFTCov()
      ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_, post_vs_, post_rho_,
                                       output_grids_elastic, -1, false);

      if (model_settings->getKrigingParameter() > 0) {

        //From doPredictionKriging
        ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_,
                                         output_grids_elastic, -1, true);
      }

      //From CKrigingAdmin::KrigAll
      if (model_settings->getDebugFlag())
        ParameterOutput::WriteFile(model_settings, block_grid_, "BlockGrid", IO::PathToInversionResults(), &simbox);

    }

    //Write Background models
    if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
      WriteBackgrounds(model_settings,
                       &simbox,
                       time_depth_mapping,
                       *model_settings->getTraceHeaderFormatOutput());
    }

    //Write correlations and post variances. Write per interval.
    for (int i = 0; i < n_intervals_; i++) {

      std::string interval_name = common_data->GetMultipleIntervalGrid()->GetIntervalName(i);

      if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
        WriteFilePriorCorrT(corr_T_[i], simbox.GetNZpad(), dt, interval_name);

      //delete corr_T_[i];
      //fftw_free(corr_T_[i]);

      if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0) {
        WriteFilePriorCorrT(corr_T_filtered_[i], simbox.GetNZpad(), dt, interval_name); // No zeros in the middle
        //delete [] corr_T_filtered_[i];
      }

      if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
        WriteFilePostVariances(post_var0_[i], post_cov_vp00_[i], post_cov_vs00_[i], post_cov_rho00_[i], interval_name);
        WriteFilePostCovGrids(model_settings, simbox, interval_name);
      }
    }

    //Write blocked wells
    if ((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
      WriteBlockedWells(common_data->GetBlockedLogs(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());
    }
    if ((model_settings->getWellOutputFlag() & IO::WELLS) > 0) {
      WriteWells(common_data->GetWells(), model_settings);
    }


    if (model_settings->getWritePrediction() || model_settings->getKrigingParameter() > 0) {
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
        }
        std::string file_name = base_name + "Undef";
        ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, facies_prob_undef_, file_name, "");
      }
      if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
        for (int i = 0; i < n_facies; i++) {
          std::string file_name = base_name + facies_names[i];
          ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, facies_prob_geo_[i], file_name, "");
        }
      }
      if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID) {
        std::string file_name = "Seismic_Quality_Grid";
        ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, quality_grid_, file_name, "");
      }
      if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
        for (int i = 0; i < n_facies; i++) {
          std::string file_name = IO::PrefixLikelihood() + facies_names[i];
          ParameterOutput::WriteToFile(&simbox, time_depth_mapping, model_settings, lh_cubes_[i], file_name, "");
        }
      }
    }

    //Simulations
    if (model_settings->getNumberOfSimulations() > 0) {
      int n_simulations = simulations_seed0_.size();
      bool kriging = model_settings->getKrigingParameter() > 0;
      for (int i = 0; i < n_simulations; i++) {
        ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, simulations_seed0_[i], simulations_seed1_[i], simulations_seed2_[i],
                                          model_settings->getOutputGridsElastic(), n_simulations, kriging);
      }
    }

    //Synthetic seismic
    if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {

      int n_angles = model_settings->getNumberOfAngles(0); //Only write synthetic seismic for the first vintage
      std::vector<float> angles = model_settings->getAngle(0);

      for (int i = 0; i < n_angles; i++) {

        //float theta       = common_data->GetSeismicDataTimeLapse(0)[i].GetAngle();
        float theta       = angles[i];
        float theta_deg   = static_cast<float>((theta*180.0/NRLib::Pi));
        std::string angle = NRLib::ToString(theta_deg, 1);

        std::string sgri_label = " Synthetic seismic for incidence angle "+angle;
        std::string file_name  = IO::PrefixSyntheticSeismicData() + angle;

        if(((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) ||
          (model_settings->getForwardModeling() == true))
          ParameterOutput::WriteFile(model_settings, synt_seismic_data_[i], file_name, IO::PathToSeismicData(), &simbox, sgri_label);

        sgri_label = "Residual computed from synthetic seismic for incidence angle "+angle;
        file_name  = IO::PrefixSyntheticResiduals() + angle;

        if((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0)
          ParameterOutput::WriteFile(model_settings, synt_residuals_[i], file_name, IO::PathToSeismicData(), &simbox, sgri_label);
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
  for(int i=0 ; i<3 ; i++) {
    for(int j=0 ; j<3 ; j++) {
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

  for(int k=0 ; k < static_cast<int>(post_cov.size()) ; k++)
    file << std::setw(9) << post_cov[k]*c0 << "\n";

  file.close();
}

void CravaResult::WriteFilePostCovGrids(const ModelSettings * model_settings,
                                        const Simbox        & simbox,
                                        std::string           interval_name) const
{
  if (interval_name != "")
    interval_name = "_"+interval_name;

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
                             const ModelSettings      * model_settings)
{
  for (size_t i = 0; i < wells.size(); i++) {
    wells[i].WriteWell(model_settings->getWellFormatFlag(),
                       model_settings->getMaxHzBackground(),
                       model_settings->getMaxHzSeismic());

  }
}

void CravaResult::FindDzMin(const Simbox      & full_inversion_simbox,
                            MultiIntervalGrid * multi_interval_grid,
                            double            & dz)
{
  double min_dz   = std::numeric_limits<double>::infinity();
  int nx          = full_inversion_simbox.getnx();
  int ny          = full_inversion_simbox.getny();
  int n_intervals = multi_interval_grid->GetNIntervals();

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      double min_dz_trace = std::numeric_limits<double>::infinity();

      //Find minumum dz for each interval
      for (int i_interval = 0; i_interval < n_intervals; i_interval++) {

        Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

        double top_value = interval_simbox->getTop(i,j);
        double bot_value = interval_simbox->getBot(i,j);
        int nz = interval_simbox->getnz();

        double min_dz_interval = (bot_value - top_value) / nz;

        if (min_dz_interval < min_dz_trace)
          min_dz_trace = min_dz_interval;

      }

      if (min_dz_trace < min_dz)
        min_dz = min_dz_trace;

    }
  }

  dz = min_dz;
}

//void CravaResult::ResampleSimple(std::vector<float>       & new_trace, //H REMOVE
//                                 const std::vector<float> & old_trace)
//{
//
//  int nt = old_trace.size();
//  int mt = new_trace.size();
//
//  int cnt = nt/2 + 1;
//  int rnt = 2*cnt;
//  int cmt = mt/2 + 1;
//  int rmt = 2*cmt;
//
//  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
//  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);
//
//  fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
//  fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));
//
//  fftw_complex * cAmpData = reinterpret_cast<fftw_complex*>(rAmpData);
//  fftw_complex * cAmpFine = reinterpret_cast<fftw_complex*>(rAmpFine);
//
//  //
//  // Fill vector to be FFT'ed
//  //
//  int n_data = static_cast<int>(old_trace.size());
//
//  for (int i = 0 ; i < n_data ; i++) {
//    rAmpData[i] = old_trace[i];
//  }
//
//  // Transform to Fourier domain
//  //
//  rfftwnd_one_real_to_complex(fftplan1, rAmpData, cAmpData);
//
//  //
//  // Fill fine-sampled grid
//  //
//  for (int i = 0; i < cnt; i++) {
//    cAmpFine[i].re = cAmpData[i].re;
//    cAmpFine[i].im = cAmpData[i].im;
//  }
//  // Pad with zeros (cmt is always greater than cnt)
//  for (int i = cnt; i < cmt; i++) {
//    cAmpFine[i].re = 0.0f;
//    cAmpFine[i].im = 0.0f;
//  }
//
//  //
//  // Fine-sampled grid: Fourier --> Time
//  //
//  rfftwnd_one_complex_to_real(fftplan2, cAmpFine, rAmpFine);
//
//  //
//  // Scale and fill grid_trace
//  //
//  float scale = 1/static_cast<float>(rnt);
//  for (int i = 0; i < rmt; i++) {
//    rAmpFine[i] = scale*rAmpFine[i];
//  }
//
//  for (int i = 0; i < mt; i++) {
//    new_trace[i] = rAmpFine[i];
//  }
//
//}


StormContGrid *
CravaResult::CreateStormGrid(const Simbox & simbox,
                             FFTGrid      * fft_grid)
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

    //Not needed anymore
    delete fft_grid;

    return(storm);
  }

  StormContGrid * storm = new StormContGrid();
  return(storm);

}

void CravaResult::WriteBackgrounds(const ModelSettings     * model_settings,
                                   const Simbox            * simbox,
                                   GridMapping             * depth_mapping,
                                   const TraceHeaderFormat & thf)
{
  if(depth_mapping != NULL && depth_mapping->getSimbox() == NULL) {
    depth_mapping->setMappingFromVelocity(background_vp_, simbox, model_settings->getOutputGridFormat());
  }

  std::string file_name_vp  = IO::PrefixBackground() + "Vp" ;
  std::string file_name_vs  = IO::PrefixBackground() + "Vs" ;
  std::string file_name_rho = IO::PrefixBackground() + "Rho";

  ExpTransf(background_vp_);
  ParameterOutput::WriteFile(model_settings, background_vp_, file_name_vp, IO::PathToBackground(), simbox, "NO_LABEL", 0, depth_mapping, thf);
  delete background_vp_;

  ExpTransf(background_vs_);
  ParameterOutput::WriteFile(model_settings, background_vs_, file_name_vs, IO::PathToBackground(), simbox, "NO_LABEL", 0, depth_mapping, thf);
  delete background_vs_;

  ExpTransf(background_rho_);
  ParameterOutput::WriteFile(model_settings, background_rho_, file_name_rho, IO::PathToBackground(), simbox, "NO_LABEL", 0, depth_mapping, thf);
  delete background_rho_;

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

void CravaResult::ComputeSyntSeismic(const ModelSettings * model_settings,
                                     const Simbox * simbox,
                                     StormContGrid * vp,
                                     StormContGrid * vs,
                                     StormContGrid * rho)
{
  LogKit::WriteHeader("Compute Synthetic Seismic and Residuals");

  int nx = vp->GetNI();
  int ny = vp->GetNJ();
  int nz = vp->GetNK();

  int nxp = simbox->GetNXpad();
  int nyp = simbox->GetNYpad();
  int nzp = simbox->GetNZpad(); //H A new nz is created in combine-grid. nzp > nz?

  std::vector<float> angles = model_settings->getAngle(0); //Synt seismic only for first vintage
  int n_theta = angles.size();

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
        Wavelet1D resultVec(&impVec, Wavelet::FIRSTORDERFORWARDDIFF);
        resultVec.fft1DInPlace();

        Wavelet1D * localWavelet = wavelets_[l]->createLocalWavelet1D(i,j);

        float sf = static_cast<float>(simbox->getRelThick(i, j))*wavelets_[l]->getLocalStretch(i,j);

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

    float theta_deg        = static_cast<float>(angles[l]*180.0/NRLib::Pi);
    std::string angle      = NRLib::ToString(theta_deg, 1);
    std::string sgri_label = " Synthetic seismic for incidence angle "+angle;
    std::string file_name  = IO::PrefixSyntheticSeismicData() + angle;
    if (((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_SEISMIC_DATA) > 0) ||
      (model_settings->getForwardModeling() == true))
      synt_seismic_data_.push_back(imp);
      //imp->writeFile(file_name, IO::PathToSeismicData(), simbox, sgri_label);

    if ((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0) {
      //FFTGrid * imp_residual = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);
      StormContGrid * imp_residual = new StormContGrid(nx, ny, nz);

      FFTGrid seis(nx, ny, nz, nxp, nyp, nzp);

      std::string file_name = IO::makeFullFileName(IO::PathToSeismicData(), IO::FileTemporarySeismic()+NRLib::ToString(l)+IO::SuffixCrava());
      std::string err_text;
      seis.readCravaFile(file_name, err_text);
      if (err_text == "") {
        seis.setAccessMode(FFTGrid::RANDOMACCESS);
        for (int k = 0; k < nz; k++) {
          for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
              float residual = seis.getRealValue(i, j, k) - imp->GetValue(i, j, k);
              imp_residual->SetValue(i, j, k, residual);
            }
          }
        }
        sgri_label = "Residual computed from synthetic seismic for incidence angle " + angle;
        file_name = IO::PrefixSyntheticResiduals() + angle;
        synt_residuals_.push_back(imp_residual);
        //imp->writeFile(file_name, IO::PathToSeismicData(), simbox, sgri_label);
      }
      else {
        err_text += "\nFailed to read temporary stored seismic data.\n";
        LogKit::LogMessage(LogKit::Error, err_text);
      }
    }
  }
}

StormContGrid *
CravaResult::ComputeSeismicImpedance(StormContGrid * vp,
                                     StormContGrid * vs,
                                     StormContGrid * rho,
                                     const NRLib::Matrix & reflection_matrix,
                                     int             angle) const
{
  int nx = vp->GetNI();
  int ny = vp->GetNJ();
  int nz = vp->GetNK();

  StormContGrid * impedance = new StormContGrid(nx, ny, nz);

  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        float imp = 0;
        imp += vp->GetValue(i, j, k)*reflection_matrix(angle,0);
        imp += vs->GetValue(i, j, k)*reflection_matrix(angle,1);
        imp += rho->GetValue(i, j, k)*reflection_matrix(angle,2);

        impedance->SetValue(i, j, k, imp);
      }
    }
  }

  return(impedance);
}