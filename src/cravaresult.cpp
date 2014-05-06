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

  post_vp_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVp());
  post_vs_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVs());
  post_rho_intervals_[interval_name] = seismic_parameters.GetPostRho();

  post_vp_kriging_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVpKriging());
  post_vs_kriging_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetPostVsKriging());
  post_rho_kriging_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetPostRhoKriging());

  corr_T_intervals_[interval_name]          = seismic_parameters.GetCorrT();
  corr_T_filtered_intervals_[interval_name] = seismic_parameters.GetCorrTFiltered();

  post_var0_intervals_[interval_name]      = seismic_parameters.GetPostVar0();
  post_cov_vp00_intervals_[interval_name]  = seismic_parameters.GetPostCovVp00();
  post_cov_vs00_intervals_[interval_name]  = seismic_parameters.GetPostCovVs00();
  post_cov_rho00_intervals_[interval_name] = seismic_parameters.GetPostCovRho00();

  cov_vp_intervals_[interval_name]        = new FFTGrid(seismic_parameters.GetCovVp());
  cov_vs_intervals_[interval_name]        = new FFTGrid(seismic_parameters.GetCovVs());
  cov_rho_intervals_[interval_name]       = new FFTGrid(seismic_parameters.GetCovRho());
  cr_cov_vp_vs_intervals_[interval_name]  = new FFTGrid(seismic_parameters.GetCrCovVpVs());
  cr_cov_vp_rho_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetCrCovVpRho());
  cr_cov_vs_rho_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetCrCovVsRho());

  simulations_seed0_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed0();
  simulations_seed1_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed1();
  simulations_seed2_intervals_[interval_name] = seismic_parameters.GetSimulationsSeed2();

  synt_seismic_data_intervals_[interval_name] = seismic_parameters.GetSyntSeismicData();
  synt_residuals_intervals_[interval_name]    = seismic_parameters.GetSyntResiduals();

  if (seismic_parameters.GetBlockGrid() != NULL)
    block_grid_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetBlockGrid());


  facies_prob_intervals_[interval_name] = seismic_parameters.GetFaciesProb();
  facies_prob_undef_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetFaciesProbUndefined());

  facies_prob_geo_intervals_[interval_name] = seismic_parameters.GetFaciesProb();

  lh_cubes_intervals_[interval_name] = seismic_parameters.GetLHCube();

  quality_grid_intervals_[interval_name] = new FFTGrid(seismic_parameters.GetQualityGrid());




}

void CravaResult::CombineResults(MultiIntervalGrid * multi_interval_grid)
{
  //Combine interval grids to one grid

  int n_intervals = post_vp_intervals_.size();

  //H-TEMP Testing for 1 interval;
  std::string interval_name = multi_interval_grid->GetIntervalName(0);

  corr_T_          = corr_T_intervals_.find(interval_name)->second;
  corr_T_filtered_ = corr_T_filtered_intervals_.find(interval_name)->second;

  post_vp_  = post_vp_intervals_.find(interval_name)->second;
  post_vs_  = post_vs_intervals_.find(interval_name)->second;
  post_rho_ = post_rho_intervals_.find(interval_name)->second;

  post_vp_kriging_ = post_vp_kriging_intervals_.find(interval_name)->second;
  post_vs_kriging_ = post_vs_kriging_intervals_.find(interval_name)->second;
  post_rho_kriging_ = post_rho_kriging_intervals_.find(interval_name)->second;

  background_vp_ = background_vp_intervals_.find(interval_name)->second;
  background_vs_ = background_vs_intervals_.find(interval_name)->second;
  background_rho_ = background_rho_intervals_.find(interval_name)->second;

  post_var0_ = post_var0_intervals_.find(interval_name)->second;
  post_cov_vp00_ = post_cov_vp00_intervals_.find(interval_name)->second;
  post_cov_vs00_ = post_cov_vs00_intervals_.find(interval_name)->second;
  post_cov_rho00_ = post_cov_rho00_intervals_.find(interval_name)->second;

  cov_vp_ = cov_vp_intervals_.find(interval_name)->second;
  cov_vs_ = cov_vs_intervals_.find(interval_name)->second;
  cov_rho_ = cov_rho_intervals_.find(interval_name)->second;
  cr_cov_vp_vs_ = cr_cov_vp_vs_intervals_.find(interval_name)->second;
  cr_cov_vp_rho_ = cr_cov_vp_rho_intervals_.find(interval_name)->second;
  cr_cov_vs_rho_ = cr_cov_vs_rho_intervals_.find(interval_name)->second;

  simulations_seed0_ = simulations_seed0_intervals_.find(interval_name)->second;
  simulations_seed1_ = simulations_seed1_intervals_.find(interval_name)->second;
  simulations_seed2_ = simulations_seed2_intervals_.find(interval_name)->second;

  synt_seismic_data_ = synt_seismic_data_intervals_.find(interval_name)->second;
  synt_residuals_ = synt_residuals_intervals_.find(interval_name)->second;

  block_grid_ = block_grid_intervals_.find(interval_name)->second;

  facies_prob_ = facies_prob_intervals_.find(interval_name)->second;
  facies_prob_undef_ = facies_prob_undef_intervals_.find(interval_name)->second;

  facies_prob_geo_ = facies_prob_geo_intervals_.find(interval_name)->second;

  lh_cubes_ = lh_cubes_intervals_.find(interval_name)->second;

  quality_grid_ = quality_grid_intervals_.find(interval_name)->second;


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