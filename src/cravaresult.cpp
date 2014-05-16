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
//#include "src/commondata.h"
#include "src/definitions.h"

CravaResult::CravaResult()
{
}

CravaResult::~CravaResult()
{
}

void CravaResult::AddParameters(SeismicParametersHolder & seismic_parameters,
                                std::string             & interval_name)
{

  //post_vp_intervals_[interval_name]  = new FFTGrid(*seismic_parameters.GetPostVp());
  //post_vs_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVs());
  //post_rho_intervals_[interval_name] = seismic_parameters.GetPostRho();

  //if (seismic_parameters.GetPostVpKriging() != NULL)
  //  post_vp_kriging_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVpKriging());
  //if (seismic_parameters.GetPostVsKriging() != NULL)
  //  post_vs_kriging_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVsKriging());
  //if (seismic_parameters.GetPostRhoKriging() != NULL)
  //  post_rho_kriging_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetPostRhoKriging());

  //corr_T_intervals_[interval_name]          = seismic_parameters.GetCorrT();
  //corr_T_filtered_intervals_[interval_name] = seismic_parameters.GetCorrTFiltered();

  //post_var0_intervals_[interval_name]      = seismic_parameters.GetPostVar0();
  //post_cov_vp00_intervals_[interval_name]  = seismic_parameters.GetPostCovVp00();
  //post_cov_vs00_intervals_[interval_name]  = seismic_parameters.GetPostCovVs00();
  //post_cov_rho00_intervals_[interval_name] = seismic_parameters.GetPostCovRho00();

  //cov_vp_intervals_[interval_name]        = new FFTGrid(seismic_parameters.GetCovVp());
  //cov_vs_intervals_[interval_name]        = new FFTGrid(seismic_parameters.GetCovVs());
  //cov_rho_intervals_[interval_name]       = new FFTGrid(seismic_parameters.GetCovRho());
  //cr_cov_vp_vs_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetCrCovVpVs());
  //cr_cov_vp_rho_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetCrCovVpRho());
  //cr_cov_vs_rho_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetCrCovVsRho());

  //simulations_seed0_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed0();
  //simulations_seed1_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed1();
  //simulations_seed2_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed2();

  //synt_seismic_data_intervals_[interval_name] = seismic_parameters.GetSyntSeismicData();
  //synt_residuals_intervals_[interval_name]    = seismic_parameters.GetSyntResiduals();

  //if (seismic_parameters.GetBlockGrid() != NULL)
  //  block_grid_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetBlockGrid());

  //facies_prob_intervals_[interval_name] = seismic_parameters.GetFaciesProb();
  //if (seismic_parameters.GetFaciesProbUndefined() != NULL)
  //  facies_prob_undef_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetFaciesProbUndefined());

  //facies_prob_geo_intervals_[interval_name] = seismic_parameters.GetFaciesProb();

  //lh_cubes_intervals_[interval_name] = seismic_parameters.GetLHCube();

  //if (seismic_parameters.GetQualityGrid() != NULL)
  //  quality_grid_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetQualityGrid());




}

void CravaResult::CombineResults(ModelSettings                        * model_settings,
                                 CommonData                           * common_data,
                                 std::vector<SeismicParametersHolder> & seismic_parameters_intervals)
{
  //Combine interval grids to one grid
  MultiIntervalGrid * multi_interval_grid = common_data->GetMultipleIntervalGrid();

  //int n_intervals = post_vp_intervals_.size();
  int n_intervals = multi_interval_grid->GetNIntervals();

  if (n_intervals > 0) {

    const Simbox & full_inversion_simbox = common_data->GetFullInversionSimbox();

    int nx  = full_inversion_simbox.getnx();
    int ny  = full_inversion_simbox.getny();
    //int nz  = full_inversion_simbox.getnz(); //1
    int nxp = full_inversion_simbox.GetNXpad();
    int nyp = full_inversion_simbox.GetNYpad();
    int nzp = full_inversion_simbox.GetNZpad();

    //FFTGrid * post_vp_final = new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

    //Store final grid as stormgrid
    //Compaction: Use the finest dz between intervals

    //Find maximum nz from all traces (based on interval with lowest dz)
    int nz_new    = 0;
    double dz_new = 0;

    FindNz(full_inversion_simbox, multi_interval_grid, nz_new, dz_new);


    //int erosion_priority_top_surface                                = model_settings->getErosionPriorityTopSurface();
    //const std::map<std::string,int> erosion_priority_base_surfaces  = model_settings->getErosionPriorityBaseSurfaces();

    const std::vector<int> & erosion_priorities = multi_interval_grid->GetErosionPriorities();

    //Resample
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {

        //Find minumum dz for each interval
        for (int i_interval = 0; i_interval < n_intervals; i_interval++) {

          //Use erosion priorities to decide overlapping intervals (?)

          Simbox * interval_simbox  = multi_interval_grid->GetIntervalSimbox(i_interval);
          FFTGrid * cov_vp_interval = seismic_parameters_intervals[i_interval].GetPostVp();

          std::vector<float> & trace = cov_vp_interval->getRealTrace2(i, j);

          double top_value = interval_simbox->getTop(i,j);
          double bot_value = interval_simbox->getBot(i,j);

          int nz_old       = interval_simbox->getnz();
          double dz_old    = (bot_value - top_value) / nz_old;

          //Resample to new nz based on minimum dz from all traces and intervals
          int nz_new = (bot_value - top_value) / dz_new;

          //H-TEST For resampling
          nz_new = nz_new + 50;

          std::vector<float> trace_new(nz_new);

          //From CommonData FillInData
          int nt = CommonData::FindClosestFactorableNumber(nz_new);

          int mt = 4*nt;

          rfftwnd_plan fftplan1 = rfftwnd_create_plan(1, &nt, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
          rfftwnd_plan fftplan2 = rfftwnd_create_plan(1, &mt, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

          int cnt = nt/2 + 1;
          int rnt = 2*cnt;
          int cmt = mt/2 + 1;
          int rmt = 2*cmt;

          fftw_real * rAmpData = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rnt));
          fftw_real * rAmpFine = static_cast<fftw_real*>(fftw_malloc(sizeof(float)*rmt));

          fftw_complex * cAmpData = reinterpret_cast<fftw_complex*>(rAmpData);
          fftw_complex * cAmpFine = reinterpret_cast<fftw_complex*>(rAmpFine);

          int n_data = static_cast<int>(trace.size());

          for (int i = 0 ; i < n_data ; i++) {
            rAmpData[i] = trace[i];
          }

          // Pad with zeros
          for (int i = n_data; i < rnt; i++) { //H-CHECK Pad before transforming?
            rAmpData[i] = 0.0f;
          }

          //
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
          std::vector<float> test(rmt);

          float scale = 1/static_cast<float>(rnt);
          for (int i = 0 ; i < rmt ; i++) {
            rAmpFine[i] = scale*rAmpFine[i];

            test[i] = rAmpFine[i];
          }



          //Interplolate
          double dz_fine = dz_old;
          double dz_grid = dz_new;
          int n_fine = rmt;
          //float z0_shift = 0.0f;
          int refk = 0;

          double x  = 0.0;
          double y  = 0.0;
          double z0 = 0.0;
          interval_simbox->getCoord(i, j, 0, x, y, z0);  // Get lateral position and z-start (z0)
          float  xf = static_cast<float>(x);
          float  yf = static_cast<float>(y);

          double z_min = 0.0;
          interval_simbox->getZCoord(0, xf, yf, z_min); //H FFTGrid doesnt have top/bot surfaces. Correct to use simbox?
          float z0_data = static_cast<float>(z_min);

          float z0_grid = static_cast<float>(z0);

          float z0_shift    = z0_grid - z0_data; //= 0 for fftgrid without geometry since simbox is used to getZCoord z0_data
          float inv_dz_fine = 1.0f/dz_fine;

          int n_grid = static_cast<int>(trace_new.size());

          for (int k = 0; k < n_grid; k++) {

            //GetZSimboxIndex(k, nz, nzp);

            int nz = n_grid;
            nzp = nz;

            int refk;
            if (k < (nz + nzp)/2)
              refk=k;
            else
              refk=k-nzp;


            float dl = (z0_shift + static_cast<float>(refk)*dz_grid)*inv_dz_fine;
            int   l1 = static_cast<int>(floor(dl));
            int   l2 = static_cast<int>(ceil(dl));

            if (l2 < 0 || l1 > n_fine - 1) {
              trace_new[k] = 0.0f;
            }
            else {
              if (l1 < 0) {
                trace_new[k] = rAmpFine[l2];
              }
              else if (l2 > n_fine - 1) {
                trace_new[k] = rAmpFine[l1];
              }
              else if (l1 == l2) {
                trace_new[k] = rAmpFine[l1];
              }
              else {
                float w1 = ceil(dl) - dl;
                float w2 = dl - floor(dl);
                trace_new[k] = w1*rAmpFine[l1] + w2*rAmpFine[l2];
              }
            }
          }

          //H-TEMP Writing
          std::string file_name = "old_trace";

          std::ofstream file;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
               << std::setprecision(2);

          for (int k = 0; k < trace.size(); k++)
            file << trace[k] << " ";

          file.close();


          file_name = "new_trace";
          //std::ofstream file_2;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
               << std::setprecision(2);

          for (int k = 0; k < trace_new.size(); k++)
            file << trace_new[k] << " ";

          file.close();

          file_name = "rAmpFine";
          //std::ofstream file_2;
          NRLib::OpenWrite(file, file_name);
          file << std::fixed
               << std::setprecision(2);

          for (int k = 0; k < rmt; k++)
            file << rAmpFine[k] << " ";

          file.close();

          int b = 0;

        }

      }
    }








  }
  else {

    std::string interval_name = multi_interval_grid->GetIntervalName(0);

    corr_T_          = seismic_parameters_intervals[0].GetCorrT();
    corr_T_filtered_ = seismic_parameters_intervals[0].GetCorrTFiltered();

    post_vp_  = seismic_parameters_intervals[0].GetPostVp();
    post_vs_  = seismic_parameters_intervals[0].GetPostVs();
    post_rho_ = seismic_parameters_intervals[0].GetPostRho();

    post_vp_kriging_  = seismic_parameters_intervals[0].GetPostVpKriging();
    post_vs_kriging_  = seismic_parameters_intervals[0].GetPostVsKriging();
    post_rho_kriging_ = seismic_parameters_intervals[0].GetPostRhoKriging();

    background_vp_  = background_vp_intervals_.find(interval_name)->second;
    background_vs_  = background_vs_intervals_.find(interval_name)->second;
    background_rho_ = background_rho_intervals_.find(interval_name)->second;

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

    synt_seismic_data_ = seismic_parameters_intervals[0].GetSyntSeismicData();
    synt_residuals_    = seismic_parameters_intervals[0].GetSyntResiduals();

    block_grid_ = seismic_parameters_intervals[0].GetBlockGrid();

    facies_prob_       = seismic_parameters_intervals[0].GetFaciesProb();
    facies_prob_undef_ = seismic_parameters_intervals[0].GetFaciesProbUndefined();

    facies_prob_geo_ = seismic_parameters_intervals[0].GetFaciesProbGeomodel();

    lh_cubes_ = seismic_parameters_intervals[0].GetLHCube();

    quality_grid_ = seismic_parameters_intervals[0].GetQualityGrid();

  }



  //Compute Synt seismic:
  //H-TODO Move from avoinversion and use the combined Vp, Vs, and Rho grids.













  //seismic_parameters_intervals[0].releaseGrids();






  //corr_T_          = corr_T_intervals_.find(interval_name)->second;
  //corr_T_filtered_ = corr_T_filtered_intervals_.find(interval_name)->second;

  //post_vp_  = post_vp_intervals_.find(interval_name)->second;
  //post_vs_  = post_vs_intervals_.find(interval_name)->second;
  //post_rho_ = post_rho_intervals_.find(interval_name)->second;

  //post_vp_kriging_ = post_vp_kriging_intervals_.find(interval_name)->second;
  //post_vs_kriging_ = post_vs_kriging_intervals_.find(interval_name)->second;
  //post_rho_kriging_ = post_rho_kriging_intervals_.find(interval_name)->second;

  //background_vp_ = background_vp_intervals_.find(interval_name)->second;
  //background_vs_ = background_vs_intervals_.find(interval_name)->second;
  //background_rho_ = background_rho_intervals_.find(interval_name)->second;

  //post_var0_ = post_var0_intervals_.find(interval_name)->second;
  //post_cov_vp00_ = post_cov_vp00_intervals_.find(interval_name)->second;
  //post_cov_vs00_ = post_cov_vs00_intervals_.find(interval_name)->second;
  //post_cov_rho00_ = post_cov_rho00_intervals_.find(interval_name)->second;

  //cov_vp_ = cov_vp_intervals_.find(interval_name)->second;
  //cov_vs_ = cov_vs_intervals_.find(interval_name)->second;
  //cov_rho_ = cov_rho_intervals_.find(interval_name)->second;
  //cr_cov_vp_vs_ = cr_cov_vp_vs_intervals_.find(interval_name)->second;
  //cr_cov_vp_rho_ = cr_cov_vp_rho_intervals_.find(interval_name)->second;
  //cr_cov_vs_rho_ = cr_cov_vs_rho_intervals_.find(interval_name)->second;

  //simulations_seed0_ = simulations_seed0_intervals_.find(interval_name)->second;
  //simulations_seed1_ = simulations_seed1_intervals_.find(interval_name)->second;
  //simulations_seed2_ = simulations_seed2_intervals_.find(interval_name)->second;

  //synt_seismic_data_ = synt_seismic_data_intervals_.find(interval_name)->second;
  //synt_residuals_ = synt_residuals_intervals_.find(interval_name)->second;

  //if (block_grid_intervals_.find(interval_name) != block_grid_intervals_.end())
  //  block_grid_ = block_grid_intervals_.find(interval_name)->second;

  //if (facies_prob_intervals_.find(interval_name) != facies_prob_intervals_.end())
  //  facies_prob_ = facies_prob_intervals_.find(interval_name)->second;
  //if (facies_prob_undef_intervals_.find(interval_name) != facies_prob_undef_intervals_.end())
  //  facies_prob_undef_ = facies_prob_undef_intervals_.find(interval_name)->second;

  //if (facies_prob_geo_intervals_.find(interval_name) != facies_prob_geo_intervals_.end())
  //  facies_prob_geo_ = facies_prob_geo_intervals_.find(interval_name)->second;

  //if (lh_cubes_intervals_.find(interval_name) != lh_cubes_intervals_.end())
  //  lh_cubes_ = lh_cubes_intervals_.find(interval_name)->second;

  //if (quality_grid_intervals_.find(interval_name) != quality_grid_intervals_.end())
  //  quality_grid_ = quality_grid_intervals_.find(interval_name)->second;

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

  if (model_settings->getEstimationMode()) { //Estimation model: All estimated parameters are written to file, regardless of output settings
    WriteBlockedWells(common_data->GetBlockedLogs(), model_settings, common_data->GetFaciesNames(), common_data->GetFaciesNr());



  }
  else {


    if (model_settings->getWritePrediction()) {
    //From computePostMeanResidAndFFTCov()
    //  ParameterOutput::writeParameters(simbox_, modelGeneral_, modelSettings_, postVp_, postVs_, postRho_,
    //                                   outputGridsElastic_, fileGrid_, -1, false);

    //From doPredictionKriging
    //  ParameterOutput::writeParameters(simbox_, modelGeneral_, modelSettings_, postVpKriging_, postVsKriging_, postRhoKriging_,
    //                                   outputGridsElastic_, fileGrid_, -1, true);

    //From CKrigingAdmin::KrigAll
      if (model_settings->getDebugFlag())
        block_grid_->writeFile("BlockGrid", IO::PathToInversionResults(), &simbox);
    }

    if((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0)
      WriteFilePriorCorrT(corr_T_, nzp, dt);

    fftw_free(corr_T_);

    if((model_settings->getOtherOutputFlag() & IO::PRIORCORRELATIONS) > 0) {
      WriteFilePriorCorrT(corr_T_filtered_, nzp, dt);     // No zeros in the middle
      delete [] corr_T_filtered_;
    }

    if (model_settings->getOutputGridsOther() & IO::CORRELATION) {
      WriteFilePostVariances(post_var0_, post_cov_vp00_, post_cov_vs00_, post_cov_rho00_);
      WriteFilePostCovGrids(simbox);
    }

    if((model_settings->getWellOutputFlag() & IO::BLOCKED_WELLS) > 0) {
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
        //ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, facies_prob_[i], fileName,"");
      }
      std::string file_name = base_name + "Undef";
      //ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, facies_prob_undef_, fileName,"");
    }
    if (model_settings->getOutputGridsOther() & IO::FACIESPROB) {
      for (int i = 0; i < n_facies; i++) {
        std::string file_name = base_name + facies_names[i];
        //ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, facies_prob_geo_[i], fileName,"");
      }
    }
    if (model_settings->getOutputGridsOther() & IO::SEISMIC_QUALITY_GRID)
      int a = 0;//ParameterOutput::writeToFile(simbox, modelGeneral, modelSettings, quality_grid_, fileName, "");
    if ((model_settings->getOutputGridsOther() & IO::FACIES_LIKELIHOOD) > 0) {
      for (int i = 0; i < n_facies; i++) {
        std::string file_name = IO::PrefixLikelihood() + facies_names[i];

        //ParameterOutput::writeToFile(simbox_, modelGeneral_, modelSettings_, lh_cubes_[i], fileName,"");
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
        synt_seismic_data_[i]->writeFile(file_name, IO::PathToSeismicData(), &simbox, sgri_label);

      sgri_label = "Residual computed from synthetic seismic for incidence angle "+angle;
      file_name  = IO::PrefixSyntheticResiduals() + angle;

      if((model_settings->getOutputGridsSeismic() & IO::SYNTHETIC_RESIDUAL) > 0)
        synt_residuals_[i]->writeFile(file_name, IO::PathToSeismicData(), &simbox, sgri_label);

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

void CravaResult::FindNz(const Simbox      & full_inversion_simbox,
                         MultiIntervalGrid * multi_interval_grid,
                         int               & nz,
                         double            & dz)
{
  int max_nz      = 0;
  double min_dz   = std::numeric_limits<double>::infinity();
  int nx          = full_inversion_simbox.getnx();
  int ny          = full_inversion_simbox.getny();
  int n_intervals = multi_interval_grid->GetNIntervals();

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      double min_dz_trace = 0.0;

      //Find minumum dz for each interval
      for (int i_interval = 0; i_interval < n_intervals; i_interval++) {

        Simbox * interval_simbox = multi_interval_grid->GetIntervalSimbox(i_interval);

        double top_value = interval_simbox->getTop(i,j);
        double bot_value = interval_simbox->getBot(i,j);
        int nz = interval_simbox->getnz();

        double min_dz_interval = (bot_value - top_value) / nz;

        if (i_interval == 0)
          min_dz_trace = min_dz_interval;
        else if (min_dz_interval < min_dz_trace)
          min_dz_trace = min_dz_interval;

      }

      double top_value = full_inversion_simbox.getTop(i, j);
      double bot_value = full_inversion_simbox.getBot(i, j);

      int max_nz_trace = (bot_value - top_value) / min_dz_trace;

      if (max_nz_trace > max_nz)
        max_nz = max_nz_trace;

      if (min_dz_trace < min_dz)
        min_dz = min_dz_trace;

    }
  }

  nz = max_nz;
  dz = min_dz;
}