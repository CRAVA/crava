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


