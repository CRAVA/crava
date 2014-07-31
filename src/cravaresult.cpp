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
n_intervals_(1)
{
}

CravaResult::~CravaResult()
{
  //
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
  for (size_t i = 0; i < synt_residuals_.size(); i++) {
    delete synt_residuals_[i];
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
    if (background_vp_intervals_[i] != NULL)
      delete background_vp_intervals_[i];
  }
  for (size_t i = 0; i < background_vs_intervals_.size(); i++) {
    if (background_vs_intervals_[i] != NULL)
      delete background_vs_intervals_[i];
  }
  for (size_t i = 0; i < background_rho_intervals_.size(); i++) {
    if (background_rho_intervals_[i] != NULL)
      delete background_rho_intervals_[i];
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

  /*
  for (size_t i = 0; i < wavelets_.size(); i++) {
    if (wavelets_[i] != NULL) {
      delete wavelets_[i];
      wavelets_[i] = NULL;
    }
  }
  */



}

void CravaResult::CombineResults(ModelSettings                        * model_settings,
                                 CommonData                           * common_data,
                                 std::vector<SeismicParametersHolder> & seismic_parameters_intervals)
{
  //Combine interval grids to one stormgrid
  //Compaction: Use the finest dz between intervals

  //Final grids are stored as StormContGrids
  // Results in seismic_parameters are stored as fftgrid

  //H CreateStormGrid() (inside CombineResult) creates a storm grid from fft grid, and deletes the fft_grid

  //H-TODO: CombineResult is not done
  //  There is no smoothing between intervals

  MultiIntervalGrid * multi_interval_grid     = common_data->GetMultipleIntervalGrid();
  Simbox & output_simbox                      = common_data->GetOutputSimbox();
  n_intervals_                                = multi_interval_grid->GetNIntervals();
  const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();

  int nx           = output_simbox.getnx();
  int ny           = output_simbox.getny();
  int nz_output    = output_simbox.getnz();
  double dz_output = output_simbox.getdz();

  blocked_logs_ = common_data->GetBlockedLogsOutput(); //Logs blocked to output_simbox
  CombineBlockedLogs(blocked_logs_, blocked_logs_intervals_, multi_interval_grid, &output_simbox); //Combine and resample logs create during inversion

  if (model_settings->getWritePrediction()) {

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
    CombineResult(post_vp_,  post_vp_intervals,  multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(post_vs_,  post_vs_intervals,  multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(post_rho_, post_rho_intervals, multi_interval_grid, erosion_priorities, dz_output);

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
      CombineResult(post_vp_kriged_,  post_vp_kriged_intervals,  multi_interval_grid, erosion_priorities, dz_output);
      CombineResult(post_vs_kriged_,  post_vs_kriged_intervals,  multi_interval_grid, erosion_priorities, dz_output);
      CombineResult(post_rho_kriged_, post_rho_kriged_intervals, multi_interval_grid, erosion_priorities, dz_output);
    }
  }

  //Background models
  if ((model_settings->getOutputGridsElastic() & IO::BACKGROUND) > 0) {
    background_vp_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_vs_  = new StormContGrid(output_simbox, nx, ny, nz_output);
    background_rho_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    CombineResult(background_vp_,  background_vp_intervals_,  multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(background_vs_,  background_vs_intervals_,  multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(background_rho_, background_rho_intervals_, multi_interval_grid, erosion_priorities, dz_output);
  }

  //Covariance grids
  if (!model_settings->getForwardModeling()) {
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
    CombineResult(cov_vp_,        cov_vp_intervals,        multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(cov_vs_,        cov_vs_intervals,        multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(cov_rho_,       cov_rho_intervals,       multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(cr_cov_vp_vs_,  cr_cov_vp_vs_intervals,  multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(cr_cov_vp_rho_, cr_cov_vp_rho_intervals, multi_interval_grid, erosion_priorities, dz_output);
    CombineResult(cr_cov_vs_rho_, cr_cov_vs_rho_intervals, multi_interval_grid, erosion_priorities, dz_output);
  }

  //Facies prob
  if (model_settings->getOutputGridsOther() & IO::FACIESPROB_WITH_UNDEF) {
    int n_facies = seismic_parameters_intervals[0].GetFaciesProb().size();
    facies_prob_.resize(n_facies);

    for (int j = 0; j < n_facies; j++) {
      facies_prob_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> facies_prob_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProb()[j];
      }
      CombineResult(facies_prob_[j],  facies_prob_intervals,  multi_interval_grid, erosion_priorities, dz_output);
    }

    //Undef
    facies_prob_undef_ = new StormContGrid(output_simbox, nx, ny, nz_output);
    std::vector<FFTGrid *> facies_prob_intervals_undef(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      facies_prob_intervals_undef[i] = seismic_parameters_intervals[i].GetFaciesProbUndefined();
    }
    CombineResult(facies_prob_undef_,  facies_prob_intervals_undef,  multi_interval_grid, erosion_priorities, dz_output);

  }
  if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
    int n_facies = seismic_parameters_intervals[0].GetFaciesProbGeomodel().size();
    facies_prob_.resize(n_facies);

    for (int j = 0; j < n_facies; j++) {
      facies_prob_geo_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> facies_prob_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        facies_prob_intervals[i] = seismic_parameters_intervals[i].GetFaciesProbGeomodel()[j];
      }
      CombineResult(facies_prob_geo_[j],  facies_prob_intervals,  multi_interval_grid, erosion_priorities, dz_output);
    }
  }

  //LH Cube
  if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
    int n_facies = seismic_parameters_intervals[0].GetLHCube().size();
    lh_cubes_.resize(n_facies);

    for (int j = 0; j < n_facies; j++) {
      lh_cubes_[j] = new StormContGrid(output_simbox, nx, ny, nz_output);

      std::vector<FFTGrid *> lh_cubes_intervals(n_intervals_);
      for (int i = 0; i < n_intervals_; i++) {
        lh_cubes_intervals[i] = seismic_parameters_intervals[i].GetLHCube()[j];
      }
      CombineResult(lh_cubes_[j], lh_cubes_intervals, multi_interval_grid, erosion_priorities, dz_output);

    }
  }

  //Quality grid
  if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID) {
    quality_grid_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    std::vector<FFTGrid *> quality_grid_intervals(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      quality_grid_intervals[i] = seismic_parameters_intervals[i].GetQualityGrid();
    }
    CombineResult(quality_grid_, quality_grid_intervals, multi_interval_grid, erosion_priorities, dz_output);
  }

  //Simulation grids
  if (model_settings->getNumberOfSimulations() > 0) {
    int n_simulations = seismic_parameters_intervals[0].GetSimulationsSeed0().size();
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

      CombineResult(simulations_seed0_[j], simulations_seed0_intervals, multi_interval_grid, erosion_priorities, dz_output);
      CombineResult(simulations_seed1_[j], simulations_seed1_intervals, multi_interval_grid, erosion_priorities, dz_output);
      CombineResult(simulations_seed2_[j], simulations_seed2_intervals, multi_interval_grid, erosion_priorities, dz_output);

    }
  }

  //Block grid
  if (model_settings->getWritePrediction() && model_settings->getDebugFlag()) {
    block_grid_ = new StormContGrid(output_simbox, nx, ny, nz_output);

    std::vector<FFTGrid *> block_grid_intervals(n_intervals_);
    for (int i = 0; i < n_intervals_; i++) {
      block_grid_intervals[i] = seismic_parameters_intervals[i].GetBlockGrid();
    }
    CombineResult(block_grid_, block_grid_intervals, multi_interval_grid, erosion_priorities, dz_output);
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

  if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {
    if (model_settings->getKrigingParameter() > 0)
      ComputeSyntSeismic(model_settings, &output_simbox, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_);
    else
      ComputeSyntSeismic(model_settings, &output_simbox, post_vp_, post_vs_, post_rho_);
  }

  //From CommonData::WaveletHandling
  GenerateSyntheticSeismicLogs(wavelets_, blocked_logs_, reflection_matrix_, output_simbox);

  //We need synt_seis from well wavelets, but it needs to be based on simbox_final/blocked_logs_final
  //Estimate a temp wavelet, which adds well_synt_seismic_data to blocked logs
  bool wavelet_estimated = false;
  for (int i = 0; i < model_settings->getNumberOfAngles(0); i++) {
    if (model_settings->getEstimateWavelet(0)[i] == true)
      wavelet_estimated = true;
  }
  if (wavelet_estimated == true)
    GenerateWellOptSyntSeis(model_settings, common_data, blocked_logs_, output_simbox, reflection_matrix_);

}

void CravaResult::CombineResult(StormContGrid         *& final_grid,
                                std::vector<FFTGrid *> & interval_grids,
                                MultiIntervalGrid      * multi_interval_grid,
                                const std::vector<int> & erosion_priorities,
                                double                   dz_min)
{
  //int n_intervals = multi_interval_grid->GetNIntervals();
  int nx = final_grid->GetNI();
  int ny = final_grid->GetNJ();
  int nz = final_grid->GetNK();

  bool writing  = false; //debugging
  float res_fac = 10.0; //Degree of refinement, must be integer.

  //If output simbox has the same size as the result grid there is no need to resample
  if (n_intervals_ == 1 && static_cast<int>(final_grid->GetNK()) == multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
    CreateStormGrid(*final_grid, interval_grids[0]); //deletes interval_grids[0]
    return;
  }

  //Resample
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      std::vector<std::vector<float> > new_traces(n_intervals_);

      //Resample each trace to new nz
      for (int i_interval = 0; i_interval < n_intervals_; i_interval++) {

        Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);
        FFTGrid * interval_grid  = interval_grids[i_interval];

        std::vector<float> old_trace = interval_grid->getRealTrace(i, j); //old_trace is changed below.

        double top_value = interval_simbox->getTop(i,j);
        double bot_value = interval_simbox->getBot(i,j);
        int nz_old       = interval_simbox->getnz();

        //Resample to new nz based on minimum dz from all traces and intervals
        int nz_new = static_cast<int>( ((bot_value - top_value) / dz_min) * res_fac);

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

        //FFT(trace) -> pad with zeros until new nz -> IFFT back
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

        //int nt = CommonData::FindClosestFactorableNumber(nz_old);
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
        trend_inc = (trend_last - trend_first) / (res_fac * (nz_old - 1));
        for (int k = 0; k < nz_new; k++) {
          new_traces[i_interval][k] = rAmpFine[k] + (trend_first + k*trend_inc);
        }

        fftw_free(rAmpData);
        fftw_free(rAmpFine);
        fftwnd_destroy_plan(fftplan1);
        fftwnd_destroy_plan(fftplan2);

        //H-Writing REMOVE
        if (writing) {

        //  NRLib::OpenWrite(file, file_name);
        //  file << std::fixed
        //        << std::setprecision(6);
        //  for (size_t k = 0; k < old_trace.size(); k++)
        //    file << old_trace[k] << " ";
        //  file.close();

        //  file_name = "combine/rAmpFine";
        //  NRLib::OpenWrite(file, file_name);
        //  file << std::fixed
        //        << std::setprecision(6);
        //  for (int k = 0; k < rmt; k++)
        //    file << rAmpFine[k] << " ";
        //  file.close();

          std::string file_name = "combine/new_trace_resampled_long";
          std::ofstream file;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
                << std::setprecision(6);
          for (size_t k = 0; k < new_traces[i_interval].size(); k++)
            file << new_traces[i_interval][k] << " ";
          file.close();
        }

      } //n_intervals

      //Combine vectors for each interval to one trace in stormgrid
      for (int k = 0; k < nz; k++) {

        bool two_intervals = false;

        double global_x = 0.0;
        double global_y = 0.0;
        double global_z = 0.0;

        final_grid->FindCenterOfCell(i, j, k, global_x, global_y, global_z);
        double dz_final = (final_grid->GetBotSurface().GetZ(global_x, global_y) - final_grid->GetTopSurface().GetZ(global_x, global_y)) / final_grid->GetNK();

        int i_interval = 0;
        for (i_interval = 0; i_interval < n_intervals_; i_interval++) {
          Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

          if (interval_simbox->IsInside(global_x, global_y, global_z))
            break;
        }
        if (i_interval < (n_intervals_-1)) { //Also check if it hits the next interval, unless it is the last one.
          Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval+1);

          if (interval_simbox->IsInside(global_x, global_y, global_z))
            two_intervals = true;
        }

        float value = 0.0f;
        if (two_intervals == true) {

          //Use erorsion priorities to select between the two intervals
          if (erosion_priorities[i_interval] < erosion_priorities[i_interval+1]) {
            value = GetResampledTraceValue(new_traces[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval), global_x, global_y, global_z, dz_final);
          }
          else {
            value = GetResampledTraceValue(new_traces[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval+1), global_x, global_y, global_z, dz_final);
          }
        }
        else {
          value = GetResampledTraceValue(new_traces[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval), global_x, global_y, global_z, dz_final);
        }

        final_grid->SetValue(i, j, k, value);
      }

      if (writing) {
        std::vector<float> new_trace_resampled(nz);

        for (int k = 0; k < nz; k++) {
          new_trace_resampled[k] = final_grid->GetValue(i, j, k);
        }

        std::string file_name = "combine/new_trace_resampled";
        std::ofstream file;
        NRLib::OpenWrite(file, file_name);
        file << std::fixed
              << std::setprecision(6);
        for (size_t k = 0; k < new_trace_resampled.size(); k++)
          file << new_trace_resampled[k] << " ";
        file.close();
      }

    } //ny
  } //nx

  for (size_t i = 0; i < interval_grids.size(); i++) {
    delete interval_grids[i];
    interval_grids[i] = NULL;
  }

}

float CravaResult::GetResampledTraceValue(const std::vector<float> & resampled_trace,
                                          const Simbox             & interval_simbox,
                                          const double             & global_x,
                                          const double             & global_y,
                                          const double             & global_z, //center of cell
                                          const double             & dz_final)
{
  //double top = interval_simbox.getTop(global_x, global_y);
  //double bot = interval_simbox.getBot(global_x, global_y);
  double top = interval_simbox.GetTopErodedSurface(global_x, global_y);
  double bot = interval_simbox.GetBotErodedSurface(global_x, global_y);

  int nz_resampled    = resampled_trace.size();
  double dz_resampled = (bot - top) / nz_resampled;

  double global_z_top = global_z - 0.5*dz_final;

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
                                           const Simbox              & interval_simbox,
                                           const double              & global_x,
                                           const double              & global_y,
                                           const double              & global_z, //center of cell
                                           const double              & dz_final)
{
  //double top = interval_simbox.getTop(global_x, global_y);
  //double bot = interval_simbox.getBot(global_x, global_y);
  double top = interval_simbox.GetTopErodedSurface(global_x, global_y);
  double bot = interval_simbox.GetBotErodedSurface(global_x, global_y);

  int nz_resampled    = resampled_trace.size();
  double dz_resampled = (bot - top) / nz_resampled;

  double global_z_top = global_z - 0.5*dz_final; //Use top of cell

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

  double value = resampled_trace[index];
  return(value);
}

void CravaResult::CombineBlockedLogs(std::map<std::string, BlockedLogsCommon *>                     & blocked_logs_output, //blocked to ouput_simbox in CommonData
                                     const std::vector<std::map<std::string, BlockedLogsCommon *> > & blocked_logs_intervals,
                                     MultiIntervalGrid                                              * multi_interval_grid,
                                     const Simbox                                                   * output_simbox)
{
  //Resample blocked logs to output_simbox and combine if multi_interval
  //Blocked_logs_output are blocked to output_simbox in commondata
  //Need to resampled blocked logs that are added after commondata

  int n_intervals = blocked_logs_intervals.size();
  //int nz          = output_simbox->getnz();
  float res_fac   = 10.0; //Degree of refinement, must be integer.

  //Do not resample if there is only one interval, and if output_simbox and has the same resolution as the interval_simbox
  if (n_intervals_ == 1 && output_simbox->getnz() == multi_interval_grid->GetIntervalSimbox(0)->getnz()) {
    blocked_logs_output = blocked_logs_intervals[0];
    return;
  }

  //Loop over wells
  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs_output.begin(); it != blocked_logs_output.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_logs_output.find(it->first);

    BlockedLogsCommon * blocked_log_final = iter->second;

    std::string well_name = blocked_log_final->GetWellName();
    int n_blocks          = blocked_log_final->GetNumberOfBlocks();

    //Predicted Logs
    bool got_predicted = blocked_logs_intervals[0].find(well_name)->second->GetVpPredicted().size() > 0;
    if (got_predicted) {
      std::vector<double> vp_predicted_final(n_blocks);
      std::vector<double> vs_predicted_final(n_blocks);
      std::vector<double> rho_predicted_final(n_blocks);
      std::vector<std::vector<double> > vp_predicted_intervals(n_intervals);
      std::vector<std::vector<double> > vs_predicted_intervals(n_intervals);
      std::vector<std::vector<double> > rho_predicted_intervals(n_intervals);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
        int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();
        GetWellLogContributed(vp_predicted_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpPredicted(),  first_B, last_B);
        GetWellLogContributed(vs_predicted_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVsPredicted(),  first_B, last_B);
        GetWellLogContributed(rho_predicted_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoPredicted(), first_B, last_B);
      }

      ResampleLog(vp_predicted_final,  vp_predicted_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(vs_predicted_final,  vs_predicted_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_predicted_final, rho_predicted_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

      blocked_log_final->SetVpPredicted(vp_predicted_final);
      blocked_log_final->SetVsPredicted(vs_predicted_final);
      blocked_log_final->SetRhoPredicted(rho_predicted_final);
    }

    //Real Seismic Data
    bool got_real_seismic_data = blocked_logs_intervals[0].find(well_name)->second->GetNRealSeismicData() > 0;
    if (got_real_seismic_data) {

      int n_angles = blocked_logs_intervals[0].find(well_name)->second->GetNRealSeismicData();

      for (int j = 0; j < n_angles; j++) {

        std::vector<double> real_seismic_data_final(n_blocks);
        std::vector<std::vector<double> > real_seismic_data_intervals(n_intervals_);

        //Get well logs, missing values are interpolated
        for (int i = 0; i < n_intervals_; i++) {
          int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
          int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();
          GetWellLogContributed(real_seismic_data_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRealSeismicData(j), first_B, last_B);
        }

        ResampleLog(real_seismic_data_final,  real_seismic_data_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

        blocked_log_final->SetRealSeismicData(j, real_seismic_data_final);
      }
    }

    //Facies prob
    bool got_facies_prob = blocked_logs_intervals[0].find(well_name)->second->GetNFaciesProb() > 0;
    if (got_facies_prob) {

      int n_faices = blocked_logs_intervals[0].find(well_name)->second->GetNFaciesProb();

      for (int j = 0; j < n_faices; j++) {

        std::vector<double> facies_prob_final(n_blocks);
        std::vector<std::vector<double> > facies_prob_intervals(n_intervals_);

        //Get well logs, missing values are interpolated
        for (int i = 0; i < n_intervals_; i++) {
          int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
          int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();
          GetWellLogContributed(facies_prob_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRealSeismicData(j), first_B, last_B);
        }

        ResampleLog(facies_prob_final, facies_prob_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

        blocked_log_final->SetFaciesProb(j, facies_prob_final);
      }
    }

    //Cpp
    bool got_cpp = blocked_logs_intervals[0].find(well_name)->second->GetNCpp() > 0;
    if (got_cpp) {

      int n_angles = blocked_logs_intervals[0].find(well_name)->second->GetNCpp();

      for (int j = 0; j < n_angles; j++) {

        std::vector<double> cpp_final(n_blocks);
        std::vector<std::vector<double> > cpp_intervals(n_intervals_);

        //Get well logs, missing values are interpolated
        for (int i = 0; i < n_intervals_; i++) {
          int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
          int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();
          GetWellLogContributed(cpp_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetCpp(j), first_B, last_B);
        }

        ResampleLog(cpp_final, cpp_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

        blocked_log_final->SetCpp(j, cpp_final);
      }
    }

    //ForFacies logs (vp, rho)
    bool got_vp_rho_fac_log = blocked_logs_intervals[0].find(well_name)->second->GetVpForFacies().size() > 0;
    if (got_vp_rho_fac_log) {

      std::vector<double> vp_for_facies_final(n_blocks);
      std::vector<double> rho_for_facies_final(n_blocks);
      std::vector<std::vector<double> > vp_for_facies_intervals(n_intervals_);
      std::vector<std::vector<double> > rho_for_facies_intervals(n_intervals_);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
        int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();

        GetWellLogContributed(vp_for_facies_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpForFacies(),  first_B, last_B);
        GetWellLogContributed(rho_for_facies_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoForFacies(), first_B, last_B);
      }

      ResampleLog(vp_for_facies_final,  vp_for_facies_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_for_facies_final, rho_for_facies_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

      blocked_log_final->SetVpForFacies(vp_for_facies_final);
      blocked_log_final->SetRhoForFacies(rho_for_facies_final);
    }

    //Filtered logs
    bool got_filtered_logs = blocked_logs_intervals[0].find(well_name)->second->GetVpSeismicResolution().size() > 0;
    if (got_filtered_logs) {
      std::vector<double> vp_filtered_final(n_blocks);
      std::vector<double> vs_filtered_final(n_blocks);
      std::vector<double> rho_filtered_final(n_blocks);
      std::vector<std::vector<double> > vp_filtered_intervals(n_intervals);
      std::vector<std::vector<double> > vs_filtered_intervals(n_intervals);
      std::vector<std::vector<double> > rho_filtered_intervals(n_intervals);

      //Get well logs, missing values are interpolated
      for (int i = 0; i < n_intervals_; i++) {
        int first_B = blocked_logs_intervals[i].find(well_name)->second->GetFirstB();
        int last_B  = blocked_logs_intervals[i].find(well_name)->second->GetLastB();
        GetWellLogContributed(vp_filtered_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVpSeismicResolution(),  first_B, last_B);
        GetWellLogContributed(vs_filtered_intervals[i],  blocked_logs_intervals[i].find(well_name)->second->GetVsSeismicResolution(),  first_B, last_B);
        GetWellLogContributed(rho_filtered_intervals[i], blocked_logs_intervals[i].find(well_name)->second->GetRhoSeismicResolution(), first_B, last_B);
      }

      ResampleLog(vp_filtered_final,  vp_filtered_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(vs_filtered_final,  vs_filtered_intervals,  blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);
      ResampleLog(rho_filtered_final, rho_filtered_intervals, blocked_logs_intervals, multi_interval_grid, blocked_log_final, well_name, res_fac);

      blocked_log_final->SetVpPredicted(vp_filtered_final);
      blocked_log_final->SetVsPredicted(vs_filtered_final);
      blocked_log_final->SetRhoPredicted(rho_filtered_final);
    }

    blocked_logs_output.insert(std::pair<std::string, BlockedLogsCommon *>(well_name, blocked_log_final));

  }
}

void CravaResult::GetWellLogContributed(std::vector<double>       & log_new,
                                        const std::vector<double> & log_old,
                                        int                         first_B,
                                        int                         last_B)
{
  for (int i = first_B; i < last_B + 1; i++) {
    log_new.push_back(log_old[i]);
  }

  InterpolateMissing(log_new);
}

void CravaResult::InterpolateMissing(std::vector<double> & well_log)
{
  //Interpolate out missing values

  //Assume first and last value ok.
  assert(well_log[0] != RMISSING);
  assert(well_log[well_log.size()-1] != RMISSING);

  for (size_t i = 1; i < well_log.size()-1; i++) {
    if (well_log[i] == RMISSING) {
      int first_m = i;
      int j = i + 1;
      while (well_log[j] == RMISSING) { //Find length of missing values
        j++;
      }
      int last_m = j;

      int start = first_m - 1; //value
      int end   = last_m + 1; //value

      for (int k = first_m; k < last_m + 1; k++) {
        float rel = static_cast<float>(k-start)/static_cast<float>(end-start);
        well_log[k] = rel * well_log[start] + (1 - rel) * well_log[end];
      }
    }
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

  //Get logs per interval and resample to a fine resolution
  for (int i_interval = 0; i_interval < n_intervals_; i_interval++) {

    BlockedLogsCommon * blocked_log = blocked_logs_intervals[i_interval].find(well_name)->second;
    int nz_interval                 = blocked_log->GetNumberOfBlocks() * static_cast<int>(res_fac);

    interval_logs_fine[i_interval].resize(nz_interval);
    ResampleTrace(old_log_interval[i_interval], interval_logs_fine[i_interval],  res_fac); //Interpolate missing
  }

  CombineTraces(final_log, blocked_log_final, multi_interval_grid, interval_logs_fine);
}

void CravaResult::ResampleTrace(std::vector<double> & old_trace, //not const, it is changed
                                std::vector<double> & new_trace,
                                const float           res_fac)
{

  int nz_old = old_trace.size();
  int nz_new = new_trace.size();

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
                                const std::vector<std::vector<double> > & interval_logs_fine)
{
  int nz          = final_log.size();
  int n_intervals = interval_logs_fine.size();

  const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();

  for (int k = 0; k < nz; k++) {

    bool two_intervals = false;

    double global_x = blocked_log_final->GetXposBlocked()[k];
    double global_y = blocked_log_final->GetYposBlocked()[k];
    double global_z = blocked_log_final->GetZposBlocked()[k];

    double dz_final = blocked_log_final->GetDz();

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

    double value = 0.0;
    if (two_intervals == true) {

      //Use erorsion priorities to select between the two intervals
      if (erosion_priorities[i_interval] < erosion_priorities[i_interval+1]) {
        value = GetResampledTraceValue(interval_logs_fine[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval), global_x, global_y, global_z, dz_final);
      }
      else {
        value = GetResampledTraceValue(interval_logs_fine[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval+1), global_x, global_y, global_z, dz_final);
      }
    }
    else {
      value = GetResampledTraceValue(interval_logs_fine[i_interval], *multi_interval_grid->GetIntervalSimbox(i_interval), global_x, global_y, global_z, dz_final);
    }

    final_log[k] = value;
  }
}


void CravaResult::WriteResults(ModelSettings * model_settings,
                               CommonData    * common_data)
{
  //Move writing rutines from modelavodynamic and avoinversion here

  //Results are combined to one grid in CombineResults first
  //const Simbox & simbox = common_data->GetFullInversionSimbox();
  const Simbox & simbox = common_data->GetOutputSimbox();

  //Wavelets are written out both in commonData and Wavelet1d/3d.cpp (and possibly modelavodynamic if match energies)

  int output_grids_elastic = model_settings->getOutputGridsElastic();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  if (model_settings->getEstimationMode()) { //Estimation model: All estimated parameters are written to file, regardless of output settings

    WriteBlockedWells(blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

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
      if (!model_settings->getForwardModeling()) {
        ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_, post_vs_, post_rho_,
                                         output_grids_elastic, -1, false);
      }

      if (model_settings->getKrigingParameter() > 0) {

        //From doPredictionKriging
        ParameterOutput::WriteParameters(&simbox, time_depth_mapping, model_settings, post_vp_kriged_, post_vs_kriged_, post_rho_kriged_,
                                         output_grids_elastic, -1, true);
      }

      //From CKrigingAdmin::KrigAll
      if (model_settings->getDebugFlag())
        ParameterOutput::WriteFile(model_settings, block_grid_, "BlockGrid", IO::PathToInversionResults(), &simbox);

    }

    //Write seismic data
    if( (model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0
        || (model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0 ) {

      int n_timelapses = model_settings->getNumberOfTimeLapses();

      for (int i = 0; i < n_timelapses; i ++) {

        int n_angles              = model_settings->getNumberOfAngles(i);
        std::vector<float> angles = model_settings->getAngle(i);
        std::vector<float> offset = model_settings->getLocalSegyOffset(i);

        for (int j = 0; j < n_angles; j++) {
          std::string angle           = NRLib::ToString(angles[j]*(180/M_PI), 1);
          std::string file_name_orig  = IO::PrefixOriginalSeismicData() + angle;
          std::string sgri_label      = std::string("Original seismic data for angle stack ") + angle;
          if (offset[j] < 0)
            offset[j] = model_settings->getSegyOffset(i);

          std::string angle_synt     = NRLib::ToString(j);
          std::string file_name_synt = IO::makeFullFileName(IO::PathToSeismicData(), IO::FileTemporarySeismic()+angle_synt);

          int seismic_type = common_data->GetSeismicDataTimeLapse(i)[j].GetSeismicType();

          if (seismic_type == 0) { //SEGY

            //Transforms segy to storm to use in ParamterOutput::WriteFile
            SegY * segy = common_data->GetSeismicDataTimeLapse(i)[j].GetSegY();

            NRLib::Grid<float> * nrlib_grid = new NRLib::Grid<float>();
            FFTGrid * fft_grid_tmp          = new FFTGrid(simbox.getnx(), simbox.getny(), simbox.getnz(), simbox.getnx(), simbox.getny(), simbox.getnz());

            StormContGrid * storm_tmp  = NULL;
            int missing_traces_simbox  = 0;
            int missing_traces_padding = 0;
            int dead_traces_simbox     = 0;

            common_data->FillInData(nrlib_grid,
                                    fft_grid_tmp,
                                    &simbox,
                                    storm_tmp,
                                    segy,
                                    model_settings->getSmoothLength(),
                                    missing_traces_simbox,
                                    missing_traces_padding,
                                    dead_traces_simbox,
                                    FFTGrid::DATA);

            //From FFTGrid to StormGrid
            StormContGrid * seismic_storm = CreateStormGrid(simbox, fft_grid_tmp);

            if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0)
              ParameterOutput::WriteFile(model_settings, seismic_storm, file_name_orig, IO::PathToSeismicData(), &simbox, sgri_label);// sgri_label, &simbox, "NO_LABEL", offset[j], time_depth_mapping);

            if ((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0)
              seismic_storm->WriteCravaFile(file_name_synt, simbox.getIL0(), simbox.getXL0(), simbox.getILStepX(), simbox.getILStepY(), simbox.getXLStepX(), simbox.getXLStepY());

            if (fft_grid_tmp != NULL)
              delete fft_grid_tmp;
            if (storm_tmp != NULL)
              delete storm_tmp;
            if (segy != NULL)
              delete segy;
            if (seismic_storm != NULL)
              delete seismic_storm;

          }
          else if (seismic_type == 1 || seismic_type == 2) { //STORM/SGRI
            StormContGrid * storm = common_data->GetSeismicDataTimeLapse(i)[j].GetStorm();

            if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0)
              ParameterOutput::WriteFile(model_settings, storm, file_name_orig, IO::PathToSeismicData(), &simbox, sgri_label);// sgri_label, &simbox, "NO_LABEL", offset[j], time_depth_mapping);
            if ((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0)
              storm->WriteCravaFile(file_name_synt, simbox.getIL0(), simbox.getXL0(), simbox.getILStepX(), simbox.getILStepY(), simbox.getXLStepX(), simbox.getXLStepY());

            if (storm != NULL)
              delete storm;
          }
          else { //FFTGrid

            FFTGrid * fft_grid = common_data->GetSeismicDataTimeLapse(i)[j].GetFFTGrid();

            if ((model_settings->getOutputGridsSeismic() & IO::ORIGINAL_SEISMIC_DATA) > 0)
              fft_grid->writeFile(file_name_orig, IO::PathToSeismicData(), &simbox, sgri_label, offset[i], time_depth_mapping, *model_settings->getTraceHeaderFormatOutput());
            if ((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0)
              fft_grid->writeCravaFile(file_name_synt, &simbox);

            if (fft_grid != NULL)
              delete fft_grid;
          }

        }
      }
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

      //if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      //  WriteFilePriorCorrT(corr_T_[i], simbox.GetNZpad(), dt, interval_name);

      //delete corr_T_[i];
      //fftw_free(corr_T_[i]);

      //if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0) {
      //  WriteFilePriorCorrT(corr_T_filtered_[i], simbox.GetNZpad(), dt, interval_name); // No zeros in the middle
      //  //delete [] corr_T_filtered_[i];
      //}

      if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
        WriteFilePostVariances(post_var0_[i], post_cov_vp00_[i], post_cov_vs00_[i], post_cov_rho00_[i], interval_name);
        WriteFilePostCovGrids(model_settings, simbox, interval_name);
      }
    }

    //Write blocked wells
    if ((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
      WriteBlockedWells(blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());

      //Write blocked background logs (CRA-544). Logs that are blocked to extended background model (extended simbox with correlation direction).
      //Do not write if multiple intervals is used
      if (n_intervals_ == 1)
        WriteBlockedWells(bg_blocked_logs_, model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());
    }

    //if ((model_settings->getWellOutputFlag() & IO::WELLS) > 0) {
    //  WriteWells(common_data->GetWells(), model_settings);
    //}

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
                                          model_settings->getOutputGridsElastic(), i, kriging);
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
                                  FFTGrid       * fft_grid)
{
  int nx = grid_new.GetNI();
  int ny = grid_new.GetNJ();
  int nz = grid_new.GetNK();

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        float value = fft_grid->getRealValue(i, j, k);
        grid_new.SetValue(i, j, k, value);
      }
    }
  }

  delete fft_grid;
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
  //delete background_vp_;

  ExpTransf(background_vs_);
  ParameterOutput::WriteFile(model_settings, background_vs_, file_name_vs, IO::PathToBackground(), simbox, "NO_LABEL", 0, depth_mapping, thf);
  //delete background_vs_;

  ExpTransf(background_rho_);
  ParameterOutput::WriteFile(model_settings, background_rho_, file_name_rho, IO::PathToBackground(), simbox, "NO_LABEL", 0, depth_mapping, thf);
  //delete background_rho_;

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
                                     const Simbox        * simbox,
                                     StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     StormContGrid       * rho)
{
  LogKit::WriteHeader("Compute Synthetic Seismic and Residuals");

  int nx = vp->GetNI();
  int ny = vp->GetNJ();
  int nz = vp->GetNK();

  int nxp = simbox->GetNXpad();
  int nyp = simbox->GetNYpad();
  int nzp = simbox->GetNZpad();

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
      //delete imp_residual;
      //imp_residual = NULL;
    }
    //delete imp;
    //imp = NULL;
  }
}

StormContGrid *
CravaResult::ComputeSeismicImpedance(StormContGrid       * vp,
                                     StormContGrid       * vs,
                                     StormContGrid       * rho,
                                     const NRLib::Matrix & reflection_matrix,
                                     int                   angle) const
{
  int nx = vp->GetNI();
  int ny = vp->GetNJ();
  int nz = vp->GetNK();

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

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_logs.begin(); it != blocked_logs.end(); it++) {
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

  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    if (blocked_log->GetIsDeviated() == false)
      blocked_log->GenerateSyntheticSeismic(reflection_matrix, wavelet, nz, nzp, &simbox);
  }
}

//void CravaResult::SetWellSyntheticSeismic(const std::vector<Wavelet *>                          & wavelet,
//                                          std::map<std::string, BlockedLogsCommon *>            & blocked_wells,
//                                          const std::vector<std::vector<std::vector<double> > > & synt_seis, //vector(angle) vector(wells)
//                                          const Simbox                                          & simbox,
//                                          const std::vector<bool>                               & wavelet_estimated)
//{
//  int nz  = simbox.getnz();
//
//  int w = 0;
//  for (std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
//    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
//    BlockedLogsCommon * blocked_log = iter->second;
//
//    int n_angles = wavelet.size();
//    const std::vector<int> & ipos = blocked_log->GetIposVector();
//    const std::vector<int> & jpos = blocked_log->GetJposVector();
//    float z0 = static_cast<float> (blocked_log->GetZposBlocked()[0]);
//
//    blocked_log->SetNAngles(n_angles);
//    for (int i = 0; i < n_angles; i++) {
//      float dz_well = static_cast<float>(simbox.getRelThick(ipos[0], jpos[0])) * wavelet[i]->getDz();
//
//      if (wavelet_estimated[i] == true && synt_seis[i][w].size() > 0)
//        blocked_log->SetLogFromVerticalTrend(synt_seis[i][w], blocked_log->GetContLogsSeismicRes(), blocked_log->GetActualSyntSeismicData(), blocked_log->GetWellSyntSeismicData(),
//                                             z0, dz_well, nz, "WELL_SYNTHETIC_SEISMIC", i);
//    }
//
//    w++;
//  }
//}

void CravaResult::GenerateWellOptSyntSeis(ModelSettings                              * model_settings,
                                          CommonData                                 * common_data,
                                          std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                          const Simbox                               & simbox,
                                          const NRLib::Matrix                        & reflection_matrix)
{
  int n_angles = model_settings->getNumberOfAngles(0);

  //Need to set angles in blocked_log before wavelet estimation
  for(std::map<std::string, BlockedLogsCommon *>::const_iterator it = blocked_wells.begin(); it != blocked_wells.end(); it++) {
    std::map<std::string, BlockedLogsCommon *>::const_iterator iter = blocked_wells.find(it->first);
    BlockedLogsCommon * blocked_log = iter->second;

    blocked_log->SetNAngles(n_angles);
  }

  std::vector<Surface *> wavelet_estim_interval;
  std::string            err_text_tmp = "";
  int                    error        = 0;

  CommonData::FindWaveletEstimationInterval(common_data->GetWaveletEstIntTop(), common_data->GetWaveletEstIntBot(), wavelet_estim_interval, simbox, err_text_tmp);
  const std::vector<bool> & estimate_wavelet = model_settings->getEstimateWavelet(0);

  for (int i = 0; i < n_angles; i++) {

    if (estimate_wavelet[i] == true) {

      //Adds well_synt_seismic_data to blocked_logs
      Wavelet * wavelet_tmp = new Wavelet1D(&simbox,
                                            &common_data->GetSeismicDataTimeLapse(0)[i],
                                            blocked_wells,
                                            wavelet_estim_interval,
                                            model_settings,
                                            reflection_matrix,
                                            i,
                                            error,
                                            err_text_tmp,
                                            false);

      delete wavelet_tmp;
    }
  }
}
