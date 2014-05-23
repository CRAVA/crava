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
//#include "src/definitions.h"

CravaResult::CravaResult()
{
}

CravaResult::~CravaResult()
{
}

void CravaResult::CombineResults(ModelSettings                        * model_settings,
                                 CommonData                           * common_data,
                                 std::vector<SeismicParametersHolder> & seismic_parameters_intervals)
{
  //Combine interval grids to one grid
  //Store final grid as stormgrid
  //Compaction: Use the finest dz between intervals


  MultiIntervalGrid * multi_interval_grid = common_data->GetMultipleIntervalGrid();
  int n_intervals = multi_interval_grid->GetNIntervals();

  if (n_intervals > 0) { //H Fix

   const Simbox & full_inversion_simbox        = common_data->GetFullInversionSimbox();
   const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();
   int nx = full_inversion_simbox.getnx();
   int ny = full_inversion_simbox.getny();

    //Find minimum dz from all traces
    double dz_min = 0.0;
    FindDzMin(full_inversion_simbox, multi_interval_grid, dz_min);

    int nz_max = static_cast<int>((full_inversion_simbox.getBotZMax() - full_inversion_simbox.getTopZMin()) / dz_min);

    //H-TEST
    nz_max = nz_max + 50;



    if (model_settings->getWritePrediction()) {

      //Post vp, vs and rho from avoinversion computePostMeanResidAndFFTCov()
      post_vp_  = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);
      post_vs_  = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);
      post_rho_ = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);

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
        post_vp_kriged_  = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);
        post_vs_kriged_  = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);
        post_rho_kriged_ = new StormContGrid(full_inversion_simbox, nx, ny, nz_max);

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

  }
  else {

    //Final grids are stored as StormContGrids
    // Results in seismic_parameters are stored as fftgrid

    //CreateStormGrid() creates a storm grid from fft grid, and deletes the fft_grid


    std::string interval_name = multi_interval_grid->GetIntervalName(0);
    const Simbox & simbox = common_data->GetFullInversionSimbox();

    corr_T_          = seismic_parameters_intervals[0].GetCorrT();
    corr_T_filtered_ = seismic_parameters_intervals[0].GetCorrTFiltered();

    post_vp_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVp());
    post_vs_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVs());
    post_rho_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostRho());

    post_vp_kriged_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVpKriged());
    post_vs_kriged_  = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostVsKriged());
    post_rho_kriged_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetPostRhoKriged());

    background_vp_  = CreateStormGrid(simbox, background_vp_intervals_.find(interval_name)->second);
    background_vs_  = CreateStormGrid(simbox, background_vs_intervals_.find(interval_name)->second);
    background_rho_ = CreateStormGrid(simbox, background_rho_intervals_.find(interval_name)->second);

    post_var0_      = seismic_parameters_intervals[0].GetPostVar0();
    post_cov_vp00_  = seismic_parameters_intervals[0].GetPostCovVp00();
    post_cov_vs00_  = seismic_parameters_intervals[0].GetPostCovVs00();
    post_cov_rho00_ = seismic_parameters_intervals[0].GetPostCovRho00();

    cov_vp_        = seismic_parameters_intervals[0].GetCovVp();
    cov_vs_        = seismic_parameters_intervals[0].GetCovVs();
    cov_rho_       = seismic_parameters_intervals[0].GetCovRho();
    cr_cov_vp_vs_  = seismic_parameters_intervals[0].GetCrCovVpVs();
    cr_cov_vp_rho_ = seismic_parameters_intervals[0].GetCrCovVpRho();
    cr_cov_vs_rho_ = seismic_parameters_intervals[0].GetCrCovVsRho();

    simulations_seed0_ = seismic_parameters_intervals[0].GetSimulationsSeed0();
    simulations_seed1_ = seismic_parameters_intervals[0].GetSimulationsSeed1();
    simulations_seed2_ = seismic_parameters_intervals[0].GetSimulationsSeed2();

    //H-TODO Remove this, computeSyntSeismic will be moved until after post_vp_ grids are made
    int n_angles = seismic_parameters_intervals[0].GetSyntSeismicData().size();
    synt_seismic_data_.resize(n_angles);
    synt_residuals_.resize(n_angles);
    for (int i = 0; i < n_angles; i++) {
      synt_seismic_data_[i] = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetSyntSeismicData()[i]);
      synt_residuals_[i]    = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetSyntResiduals()[i]);
    }

    block_grid_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetBlockGrid());

    int n_facies = seismic_parameters_intervals[0].GetFaciesProb().size();
    facies_prob_.resize(n_facies);
    for (int i = 0; i < n_facies; i++) {
      facies_prob_[i]    = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetFaciesProb()[i]);
      facies_prob_undef_ = CreateStormGrid(simbox, seismic_parameters_intervals[0].GetFaciesProbUndefined());
    }

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



  //Compute Synt seismic:
  //H-TODO Move from avoinversion and use the combined Vp, Vs, and Rho grids.









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
        std::vector<float> old_trace = interval_grid->getRealTrace2(i, j); //H old_trace is changed below.

        double top_value = interval_simbox->getTop(i,j);
        double bot_value = interval_simbox->getBot(i,j);

        int nz_old       = interval_simbox->getnz();
        //double dz_old    = (bot_value - top_value) / nz_old;

        //Resample to new nz based on minimum dz from all traces and intervals
        int nz_new = static_cast<int>((bot_value - top_value) / dz_min);

        //H-TEST For resampling
        nz_new = nz_new + 50;

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

        //Alternative as in resampling of traces from CommonData::FillInData
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
                               CommonData    * common_data,
                               const Simbox  & simbox)
{
  //Move writing rutines from modelavodynamic and avoinversion here

  //Results are combined to one grid in CombineResults first

  float dt = static_cast<float>(simbox.getdz());
  //int nz   = simbox.getnz();
  int nzp  = simbox.GetNZpad();

  int output_grids_elastic = model_settings->getOutputGridsElastic();
  //bool file_grid           = model_settings->getFileGrid();
  GridMapping * time_depth_mapping = common_data->GetTimeDepthMapping();

  if (model_settings->getEstimationMode()) { //Estimation model: All estimated parameters are written to file, regardless of output settings
    WriteBlockedWells(common_data->GetBlockedLogs(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());



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

    if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      WriteFilePriorCorrT(corr_T_, nzp, dt);

    fftw_free(corr_T_);

    if ((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0) {
      WriteFilePriorCorrT(corr_T_filtered_, nzp, dt); // No zeros in the middle
      delete [] corr_T_filtered_;
    }

    if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
      WriteFilePostVariances(post_var0_, post_cov_vp00_, post_cov_vs00_, post_cov_rho00_);
      WriteFilePostCovGrids(simbox);
    }

    if ((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
      WriteBlockedWells(common_data->GetBlockedLogs(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());
    }

    if (model_settings->getWritePrediction() || model_settings->getKrigingParameter() > 0) {
      KrigingData3D kd(common_data->GetBlockedLogs(), 1); // 1 = full resolution logs
      std::string base_name = "Raw_" + IO::PrefixKrigingData() + IO::SuffixGeneralData();
      std::string file_name = IO::makeFullFileName(IO::PathToInversionResults(), base_name);
      kd.writeToFile(file_name);
    }
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

  //Synthetic seismic is created if Forward modelleing or getGenerateSeismicAfterInv()
  if (model_settings->getGenerateSeismicAfterInv() || model_settings->getForwardModeling()) {

    int n_angles = model_settings->getNumberOfAngles(0); //Only write synthetic seismic for the first vintage

    for (int i = 0; i < n_angles; i++) {

      float theta       = common_data->GetSeismicDataTimeLapse(0)[i].GetAngle();
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

void CravaResult::WriteFilePriorCorrT(fftw_real   * prior_corr_T,
                                      const int   & nzp,
                                      const float & dt) const
{
  // This is the cyclic and filtered version of CorrT
  std::string base_name = IO::PrefixPrior() + IO::FileTemporalCorr() + IO::SuffixGeneralData();
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
                                         const std::vector<float> & post_cov_rho00) const
{
  std::string base_name = IO::PrefixPosterior() + IO::FileParameterCov() + IO::SuffixGeneralData();
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

void CravaResult::WriteFilePostCovGrids(const Simbox & simbox) const
{
  std::string file_name;
  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vp";
  cov_vp_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cov_vp_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for Vp");
  cov_vp_ ->endAccess();

  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Vs";
  cov_vs_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cov_vs_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for Vs");
  cov_vs_ ->endAccess();

  file_name = IO::PrefixPosterior() + IO::PrefixCovariance() + "Rho";
  cov_rho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cov_rho_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior covariance for density");
  cov_rho_ ->endAccess();

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpVs";
  cr_cov_vp_vs_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cr_cov_vp_vs_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vp,Vs)");
  cr_cov_vp_vs_ ->endAccess();

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VpRho";
  cr_cov_vp_rho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cr_cov_vp_rho_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vp,density)");
  cr_cov_vp_rho_ ->endAccess();

  file_name = IO::PrefixPosterior() + IO::PrefixCrossCovariance() + "VsRho";
  cr_cov_vs_rho_ ->setAccessMode(FFTGrid::RANDOMACCESS);
  cr_cov_vs_rho_ ->writeFile(file_name, IO::PathToCorrelations(), &simbox, "Posterior cross-covariance for (Vs,density)");
  cr_cov_vs_rho_ ->endAccess();
}

void CravaResult::WriteBlockedWells(const std::map<std::string, BlockedLogsCommon *> & blocked_wells,
                                    const ModelSettings                        * model_settings,
                                    std::vector<std::string>                     facies_name,
                                    std::vector<int>                             facies_label)
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

void CravaResult::ResampleSimple(std::vector<float>       & new_trace, //H REMOVE
                                 const std::vector<float> & old_trace)
{

  int nt = old_trace.size();
  int mt = new_trace.size();

  int cnt = nt/2 + 1;
  int rnt = 2*cnt;
  int cmt = mt/2 + 1;
  int rmt = 2*cmt;

  rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
  fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

  fftw_complex * cAmpData = reinterpret_cast<fftw_complex*>(rAmpData);
  fftw_complex * cAmpFine = reinterpret_cast<fftw_complex*>(rAmpFine);

  //
  // Fill vector to be FFT'ed
  //
  int n_data = static_cast<int>(old_trace.size());

  for (int i = 0 ; i < n_data ; i++) {
    rAmpData[i] = old_trace[i];
  }

  // Transform to Fourier domain
  //
  rfftwnd_one_real_to_complex(fftplan1, rAmpData, cAmpData);

  //
  // Fill fine-sampled grid
  //
  for (int i = 0; i < cnt; i++) {
    cAmpFine[i].re = cAmpData[i].re;
    cAmpFine[i].im = cAmpData[i].im;
  }
  // Pad with zeros (cmt is always greater than cnt)
  for (int i = cnt; i < cmt; i++) {
    cAmpFine[i].re = 0.0f;
    cAmpFine[i].im = 0.0f;
  }

  //
  // Fine-sampled grid: Fourier --> Time
  //
  rfftwnd_one_complex_to_real(fftplan2, cAmpFine, rAmpFine);

  //
  // Scale and fill grid_trace
  //
  float scale = 1/static_cast<float>(rnt);
  for (int i = 0; i < rmt; i++) {
    rAmpFine[i] = scale*rAmpFine[i];
  }

  for (int i = 0; i < mt; i++) {
    new_trace[i] = rAmpFine[i];
  }


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