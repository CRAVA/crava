/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <limits.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <typeinfo>
#include <algorithm>

#ifdef PARALLEL
#include <omp.h>
#endif

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelsettings.h"
#include "src/simbox.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/inputfiles.h"
#include "src/io.h"
#include "src/tasklist.h"
#include "src/timeline.h"
#include "src/state4d.h"
#include "src/seismicparametersholder.h"
#include "src/surfacefrompoints.h"
#include "src/parameteroutput.h"

#include "lib/utils.h"
#include "lib/random.h"
#include "lib/timekit.hpp"
#include "nrlib/iotools/fileio.hpp"
#include "nrlib/iotools/stringtools.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/surface/surfaceio.hpp"
#include "nrlib/surface/surface.hpp"
#include "nrlib/iotools/logkit.hpp"

#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsrock.h"


ModelGeneral::ModelGeneral(ModelSettings           *& model_settings,
                           const InputFiles         * input_files,
                           SeismicParametersHolder  & seismic_parameters,
                           CommonData               * common_data,
                           int                        i_interval)
                           :do_4D_inversion_(model_settings->getDo4DInversion()),
                            do_4D_rock_physics_vnversion_(model_settings->getDo4DRockPhysicsInversion())
{
  random_gen_              = NULL;
  time_line_               = NULL;
  time_depth_mapping_      = NULL;
  velocity_from_inversion_ = false;
  {
    simbox_ = common_data->GetMultipleIntervalGrid()->GetIntervalSimbox(i_interval);

    if (input_files->getSeedFile() == "")
      random_gen_ = new RandomGen(model_settings->getSeed());
    else
      random_gen_ = new RandomGen(input_files->getSeedFile().c_str());

    interval_name_ = model_settings->getIntervalName(i_interval);

    //
    // FORWARD MODELLING
    //
    if (model_settings->getForwardModeling() == true) {
    //
    }
    else {
      //
      // INVERSION/ESTIMATION
      //

      //Facies-names
      facies_names_  = common_data->GetFaciesNames();
      facies_labels_ = common_data->GetFaciesNr();

      //Priorfacies
      if (common_data->GetPriorFacies().size() > 0) {
        prior_facies_ = common_data->GetPriorFaciesInterval(i_interval);
      }
      if (model_settings->getIsPriorFaciesProbGiven()==ModelSettings::FACIES_FROM_CUBES) {
        prior_facies_prob_cubes_.resize(3);
        prior_facies_prob_cubes_[0] = new FFTGrid(common_data->GetPriorFaciesProbCube(i_interval, 0), simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
        prior_facies_prob_cubes_[1] = new FFTGrid(common_data->GetPriorFaciesProbCube(i_interval, 1), simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
        prior_facies_prob_cubes_[2] = new FFTGrid(common_data->GetPriorFaciesProbCube(i_interval, 2), simbox_->GetNXpad(), simbox_->GetNYpad(), simbox_->GetNZpad());
      }

      //TimeDepthMapping if intervals isn't used.
      if (model_settings->getDoDepthConversion()) {
        time_depth_mapping_      = common_data->GetTimeDepthMapping();
        velocity_from_inversion_ = common_data->GetVelocityFromInversion();
      }

      //Wells blocked to this interval simbox
      blocked_logs_ = common_data->GetBlockedLogsInterval(i_interval);

      if (common_data->GetTrendCubes().size() > 0)
        trend_cubes_ = common_data->GetTrendCube(i_interval);

      rock_distributions_  = common_data->GetDistributionsRock();
      reservoir_variables_ = common_data->GetReservoirVariables();

      //Set up timeline
      time_line_ = common_data->GetTimeLine();

      if (model_settings->getDo4DInversion()) {

        //setFaciesNamesFromRockPhysics();
        NRLib::Vector initial_mean(6);
        NRLib::Matrix initial_cov(6,6);

        SetupState4D(seismic_parameters, simbox_, state4d_, initial_mean, initial_cov);

        time_evolution_ = TimeEvolution(10000, *time_line_, rock_distributions_.begin()->second); //NBNB OK 10000->1000 for speed during testing
        time_evolution_.SetInitialMean(initial_mean);
        time_evolution_.SetInitialCov(initial_cov);
      }
    }
  }

}


ModelGeneral::~ModelGeneral(void)
{
  if (time_depth_mapping_!=NULL)
    delete time_depth_mapping_;

  delete random_gen_;
}

std::map<std::string, DistributionsRock *>
ModelGeneral::GetRockDistributionTime0() const
{
  std::map<std::string, DistributionsRock *> rock_dist_t0;

  for (std::map<std::string, std::vector<DistributionsRock *> >::const_iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    std::string name = it->first;
    std::vector<DistributionsRock *> rock_dist = it->second;
    rock_dist_t0[name] = rock_dist[0];
  }

  return(rock_dist_t0);
}

FFTGrid*
ModelGeneral::CreateFFTGrid(int nx, int ny, int nz, int nxp, int nyp, int nzp, bool fileGrid)
{
  FFTGrid* fftGrid;

  if (fileGrid)
    fftGrid =  new FFTFileGrid(nx, ny, nz, nxp, nyp, nzp);
  else
    fftGrid =  new FFTGrid(nx, ny, nz, nxp, nyp, nzp);

  return(fftGrid);
}

void
ModelGeneral::CopyCorrelationsTo4DState(SeismicParametersHolder & seismicParameters,
                                        State4D                 & state4d)
{
  // Allocates the static sigma grids: 6 grids.

  // Static sigma
  FFTGrid * vp_vp_stat;
  FFTGrid * vp_vs_stat;
  FFTGrid * vp_rho_stat;
  FFTGrid * vs_vs_stat;
  FFTGrid * vs_rho_stat;
  FFTGrid * rho_rho_stat;

  // Copying grids for sigma static
  vp_vp_stat   = new FFTGrid( seismicParameters.GetCovVp());
  vp_vs_stat   = new FFTGrid( seismicParameters.GetCrCovVpVs());
  vp_rho_stat  = new FFTGrid( seismicParameters.GetCrCovVpRho());
  vs_vs_stat   = new FFTGrid( seismicParameters.GetCovVs());
  vs_rho_stat  = new FFTGrid( seismicParameters.GetCrCovVsRho());
  rho_rho_stat = new FFTGrid( seismicParameters.GetCovRho());

  state4d.setStaticSigma(vp_vp_stat, vp_vs_stat, vp_rho_stat, vs_vs_stat, vs_rho_stat, rho_rho_stat);
}

void
ModelGeneral::SetupState4D(SeismicParametersHolder & seismicParameters,
                           const Simbox            * simbox,
                           State4D                 & state4d,
                           NRLib::Vector           & initialMean,
                           NRLib::Matrix           & initialCov)
{
  //Earlier a background was made from rockphysics3d, which was copied to seismicParameters and state4d.
  //New version: use background created in CommonData
  state4d.setStaticMu(seismicParameters.GetMeanVp(), seismicParameters.GetMeanVs(), seismicParameters.GetMeanRho());

  CopyCorrelationsTo4DState(seismicParameters, state4d_);

  const int nx    = simbox->getnx();
  const int ny    = simbox->getny();
  const int nz    = simbox->getnz();
  const int nxPad = simbox->GetNXpad();
  const int nyPad = simbox->GetNYpad();
  const int nzPad = simbox->GetNZpad();

  Complete4DBackground(nx, ny, nz, nxPad, nyPad, nzPad, initialMean, initialCov);
}

void
ModelGeneral::Complete4DBackground(const int nx, const int ny, const int nz, const int nxPad, const int nyPad, const int nzPad,NRLib::Vector &initial_mean,NRLib::Matrix &initial_cov)
{
  // Static grids (3 + 6) are set in process4DBackground.
  // Dynamic grids (3 + 6 + 9) are set here.

  FFTGrid * dynamicVp;
  FFTGrid * dynamicVs;
  FFTGrid * dynamicRho;
  FFTGrid * dynamicVpVp;

  FFTGrid *dynamicVpVs;
  FFTGrid *dynamicVpRho;
  FFTGrid *dynamicVsVs;
  FFTGrid *dynamicVsRho;
  FFTGrid *dynamicRhoRho;

  FFTGrid *staticDynamicVpVp;
  FFTGrid *staticDynamicVpVs;
  FFTGrid *staticDynamicVpRho;
  FFTGrid *staticDynamicVsVp;
  FFTGrid *staticDynamicVsVs;
  FFTGrid *staticDynamicVsRho;
  FFTGrid *staticDynamicRhoVp;
  FFTGrid *staticDynamicRhoVs;
  FFTGrid *staticDynamicRhoRho;

  dynamicVp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVp->fillInConstant(0.0);
  dynamicVp->setType(FFTGrid::PARAMETER);
  dynamicVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVs->fillInConstant(0.0);
  dynamicVs->setType(FFTGrid::PARAMETER);
  dynamicRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicRho->fillInConstant(0.0);
  dynamicRho->setType(FFTGrid::PARAMETER);

  state4d_.setDynamicMu(dynamicVp, dynamicVs, dynamicRho);
  initial_mean=state4d_.GetFullMean000();

  dynamicVpVp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpVp->fillInConstant(0.0);
  dynamicVpVp->setType(FFTGrid::COVARIANCE);
  dynamicVpVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpVs->fillInConstant(0.0);
  dynamicVpVs->setType(FFTGrid::COVARIANCE);
  dynamicVpRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVpRho->fillInConstant(0.0);
  dynamicVpRho->setType(FFTGrid::COVARIANCE);
  dynamicVsVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVsVs->fillInConstant(0.0);
  dynamicVsVs->setType(FFTGrid::COVARIANCE);
  dynamicVsRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicVsRho->fillInConstant(0.0);
  dynamicVsRho->setType(FFTGrid::COVARIANCE);
  dynamicRhoRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  dynamicRhoRho->fillInConstant(0.0);
  dynamicRhoRho->setType(FFTGrid::COVARIANCE);

  state4d_.setDynamicSigma(dynamicVpVp, dynamicVpVs, dynamicVpRho,
                                        dynamicVsVs, dynamicVsRho,
                                                     dynamicRhoRho);

  staticDynamicVpVp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpVp->fillInConstant(0.0);
  staticDynamicVpVp->setType(FFTGrid::COVARIANCE);
  staticDynamicVpVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpVs->fillInConstant(0.0);
  staticDynamicVpVs->setType(FFTGrid::COVARIANCE);
  staticDynamicVpRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVpRho->fillInConstant(0.0);
  staticDynamicVpRho->setType(FFTGrid::COVARIANCE);
  staticDynamicVsVp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsVp->fillInConstant(0.0);
  staticDynamicVsVp->setType(FFTGrid::COVARIANCE);
  staticDynamicVsVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsVs->fillInConstant(0.0);
  staticDynamicVsVs->setType(FFTGrid::COVARIANCE);
  staticDynamicVsRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicVsRho->fillInConstant(0.0);
  staticDynamicVsRho->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoVp = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoVp->fillInConstant(0.0);
  staticDynamicRhoVp->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoVs = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoVs->fillInConstant(0.0);
  staticDynamicRhoVs->setType(FFTGrid::COVARIANCE);
  staticDynamicRhoRho = ModelGeneral::CreateFFTGrid(nx, ny, nz, nxPad, nyPad, nzPad, false);
  staticDynamicRhoRho->fillInConstant(0.0);
  staticDynamicRhoRho->setType(FFTGrid::COVARIANCE);

  state4d_.setStaticDynamicSigma(staticDynamicVpVp,  staticDynamicVpVs,  staticDynamicVpRho,
                                 staticDynamicVsVp,  staticDynamicVsVs,  staticDynamicVsRho,
                                 staticDynamicRhoVp, staticDynamicRhoVs, staticDynamicRhoRho);

  initial_cov=state4d_.GetFullCov();

  state4d_.FFT();
  state4d_.setRelativeGridBase(nx, ny, nz, nxPad, nyPad, nzPad);
}

void
ModelGeneral::AdvanceTime(int time_step, SeismicParametersHolder & seismicParameters,ModelSettings* model_settings)
{
  bool debug=false;
  if (debug) Dump4Dparameters(model_settings, "_prior", time_step);  // note this prior should be equal to
                                                                    // next_prior in previous step
  if (debug) DumpSeismicParameters(model_settings,"_posterior", time_step,seismicParameters);
  state4d_.split(seismicParameters);
  if (debug) Dump4Dparameters(model_settings, "_posterior", time_step);
  state4d_.evolve(time_step, time_evolution_); //NBNB grad I grad J
  //if (debug) dump4Dparameters(model_settings, "_next_prior", time_step+1);
  state4d_.merge(seismicParameters);
  if (debug) DumpSeismicParameters(model_settings,"_next_prior", time_step+1,seismicParameters);
  seismicParameters.invFFTAllGrids(); //merge gives FFT-transformed version, need the standard for now.
}

void
ModelGeneral::LastUpdateOfStaticAndDynamicParts(SeismicParametersHolder &  seismicParameters,ModelSettings* model_settings)
{
  bool debug=true;
  int time_step=time_evolution_.GetNTimSteps()-1;
  if (debug) DumpSeismicParameters(model_settings,"_posterior", time_step,seismicParameters);

  state4d_.split(seismicParameters);
  Dump4Dparameters(model_settings, "_posterior", time_step);

}

void
ModelGeneral::mergeState4D(SeismicParametersHolder &  seismicParameters)
{
  state4d_.merge(seismicParameters);
}


void
ModelGeneral::updateState4D(SeismicParametersHolder &  seismicParameters)
{
  state4d_.split(seismicParameters);
}

void
ModelGeneral::updateState4DWithSingleParameter(FFTGrid * EPost,
                                               FFTGrid * CovPost,
                                               int       parameterNumber)
{
  state4d_.updateWithSingleParameter(EPost,
                                     CovPost,
                                     parameterNumber);
}

void
ModelGeneral::updateState4DMu(FFTGrid * mu_vp_static,
                              FFTGrid * mu_vs_static,
                              FFTGrid * mu_rho_static,
                              FFTGrid * mu_vp_dynamic,
                              FFTGrid * mu_vs_dynamic,
                              FFTGrid * mu_rho_dynamic)
{
  state4d_.updateStaticMu(mu_vp_static, mu_vs_static, mu_rho_static);

  state4d_.updateDynamicMu(mu_vp_dynamic, mu_vs_dynamic, mu_rho_dynamic);
}

bool
ModelGeneral::Do4DRockPhysicsInversion(ModelSettings* model_settings)
{
  std::vector<FFTGrid*> predictions = state4d_.doRockPhysicsInversion(*time_line_, rock_distributions_.begin()->second,  time_evolution_);
  int nParamOut = static_cast<int>(predictions.size());

  std::vector<std::string> labels(nParamOut);

  int i=0;

  for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator it = reservoir_variables_.begin(); it != reservoir_variables_.end(); it++)
  {
    labels[i] = it->first;
    i++;
  }

  std::string  outPre =  "mu_";

  for (i=0;i<nParamOut;i++)
  {
     std::string fileName;
     fileName= outPre + labels[i];

     WriteToFile(simbox_, time_depth_mapping_, model_settings, predictions[i] , fileName, labels[i]);
  }

  for (size_t ii = 0; ii < predictions.size(); ii++)
    delete predictions[ii];

  return 0;
}


void
ModelGeneral::DumpSeismicParameters(ModelSettings* model_settings, std::string identifyer, int timestep,SeismicParametersHolder &  current_state)
{

  std::string  label;
  std::string fileName;
  std::stringstream tag;
  bool transformHere=false;

  if (current_state.GetMeanVp()->getIsTransformed())
  {
    transformHere=true;
    current_state.invFFTAllGrids();
  }

  // write mu current
  tag.str(std::string());tag.clear();label = "mean_vp_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings,  current_state.GetMeanVp() , fileName,  tag.str(),true);

   ///*
  tag.str(std::string());tag.clear();label = "mean_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings,  current_state.GetMeanVs(), fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetMeanRho() , fileName, tag.str() ,true);
  // */
  // write sigma current
  tag.str(std::string());tag.clear();label = "cov_vp_vp_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCovVp() , fileName,  tag.str(),true);

  ///*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCrCovVpVs() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCrCovVpRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCovVs() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_current_step_"; tag << label << timestep << identifyer ; fileName= tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCrCovVsRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCovRho() , fileName,  tag.str(),true);
  // */
  if (transformHere)
    current_state.FFTAllGrids();
}

void
ModelGeneral::Dump4Dparameters(const ModelSettings* model_settings, std::string identifyer, int timestep, bool print_padding)
{
  state4d_.iFFT();

  std::string  outPath =  "";
  std::string  label;
  std::string fileName;
  std::stringstream tag;

  if (timestep<0)
    return;

  // write mu static
  tag.str(std::string());tag.clear();label = "mean_vp_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVpStatic() , fileName,  tag.str(),print_padding);
 // /*
  tag.str(std::string());tag.clear();label = "mean_vs_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVsStatic() , fileName, tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "mean_rho_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuRhoStatic() , fileName, tag.str(),print_padding);
  // */
  // write mu dynamic
  tag.str(std::string());tag.clear();label = "mean_vp_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVpDynamic() , fileName, tag.str(),print_padding);
  // /*
  tag.str(std::string());tag.clear();label = "mean_vs_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVsDynamic() , fileName, tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "mean_rho_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuRhoDynamic() , fileName,  tag.str(),print_padding);
  // */

  // write sigma static - static
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpStaticStatic() , fileName,  tag.str(),print_padding);

  ///*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVsStaticStatic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpRhoStaticStatic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsVsStaticStatic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsRhoStaticStatic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovRhoRhoStaticStatic() , fileName,  tag.str(),print_padding);
     // */
  // write sigma dynamic - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpDynamicDynamic() , fileName,  tag.str(),print_padding);

  ///*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVsDynamicDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpRhoDynamicDynamic(), fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsVsDynamicDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsRhoDynamicDynamic(), fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovRhoRhoDynamicDynamic() , fileName,  tag.str(),print_padding);
  // */
  // write sigma static - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpStaticDynamic() , fileName,  tag.str(),print_padding);

  // /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVsStaticDynamic() , fileName, tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpRhoStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsVpStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsVsStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVsRhoStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_rho_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovRhoVpStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_rho_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovRhoVsStaticDynamic() , fileName,  tag.str(),print_padding);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovRhoRhoStaticDynamic() , fileName,  tag.str(),print_padding);
   // */
  state4d_.FFT();
}

void
ModelGeneral::WriteToFile(const Simbox        * simbox,
                          GridMapping         * time_depth_mapping,
                          const ModelSettings * model_settings,
                          FFTGrid             * grid,
                          const std::string   & file_name,
                          const std::string   & sgri_label,
                          bool                  padding)
{
  //GridMapping * timeCutMapping;//   = modelGeneral->getTimeCutMapping(); //Included in the new simbox format.
  float seismic_start_time  = 0.0; //Hack for Sebastian, was: model->getModelSettings()->getSegyOffset();
  TraceHeaderFormat *format = model_settings->getTraceHeaderFormatOutput();

  grid->writeFile(file_name,
                  IO::PathToInversionResults(),
                  simbox,
                  sgri_label,
                  seismic_start_time,
                  time_depth_mapping,
                  *format,
                  padding);
}
