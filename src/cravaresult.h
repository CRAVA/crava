/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RESULT_H
#define RESULT_H

#include <math.h>
#include <string>

//#include "nrlib/segy/segy.hpp"
#include "src/definitions.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/modelsettings.h"

#include "src/multiintervalgrid.h"

class FFTGrid;
class Simbox;

class CravaResult
{
public:
  CravaResult();

  ~CravaResult();

  void CombineResults(MultiIntervalGrid * multi_interval_grid);

  void WriteResults(ModelSettings * model_settings,
                    const Simbox  & simbox);

  void WriteFilePriorCorrT(fftw_real   * prior_corr_T,
                           const int   & nzp,
                           const float & dt) const;

  //GET FUNCTIONS

  //SET FUNCTIONS

  void AddPostVp(std::string interval_name, FFTGrid * vp)     { post_vp_intervals_[interval_name]  = new FFTGrid(vp)  ;}
  void AddPostVs(std::string interval_name, FFTGrid * vs)     { post_vs_intervals_[interval_name]  = new FFTGrid(vs)  ;}
  void AddPostRho(std::string interval_name, FFTGrid * rho)   { post_rho_intervals_[interval_name] = new FFTGrid(rho) ;}

  void AddCorrT(std::string interval_name, fftw_real * corr_T)              { corr_T_intervals_[interval_name]          = corr_T          ;}
  void AddCorrTFiltered(std::string interval_name, float * corr_T_filtered) { corr_T_filtered_intervals_[interval_name] = corr_T_filtered ;}

private:

  std::map<std::string, FFTGrid *> post_vp_intervals_;
  std::map<std::string, FFTGrid *> post_vs_intervals_;
  std::map<std::string, FFTGrid *> post_rho_intervals_;

  std::map<std::string, fftw_real *> corr_T_intervals_;
  std::map<std::string, fftw_real *> corr_T_filtered_intervals_;

  FFTGrid * post_vp_;
  FFTGrid * post_vs_;
  FFTGrid * post_rho_;

  fftw_real * corr_T_;
  float     * corr_T_filtered_;

  int n_intervals_;

};

#endif
