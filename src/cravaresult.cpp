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
//#include "src/fftgrid.h"

CravaResult::CravaResult()
{
}

CravaResult::~CravaResult()
{
}

void CravaResult::CombineResults(MultiIntervalGrid * multi_interval_grid)
{
  //Combine interval grids to one grid

  int n_intervals = post_vp_intervals_.size();

  //Testing for 1 interval;
  std::string interval_name = multi_interval_grid->GetIntervalName(0);

  corr_T_          = corr_T_intervals_.find(interval_name)->second;
  corr_T_filtered_ = corr_T_filtered_intervals_.find(interval_name)->second;

  post_vp_  = post_vp_intervals_.find(interval_name)->second;
  post_vs_  = post_vs_intervals_.find(interval_name)->second;
  post_rho_ = post_rho_intervals_.find(interval_name)->second;



}

void CravaResult::WriteResults(ModelSettings * model_settings,
                               const Simbox  & simbox)
{
  //Move writing rutines from modelavodynamic and avoinversion here

  //Results are combined to one grid in CombineResults first

  float dt = static_cast<float>(simbox.getdz());
  int nz   = simbox.getnz();
  int nzp  = simbox.GetNZpad();

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