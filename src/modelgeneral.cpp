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
//#include "src/parameteroutput.h"

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
      //if (model_settings->GetMultipleIntervalSetting() == true)
      //  multi_interval_ = true;

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

  for (std::map<std::string, std::vector<DistributionsRock *> >::iterator it = rock_distributions_.begin(); it != rock_distributions_.end(); it++) {
    std::vector<DistributionsRock *> rock = it->second;
    for (size_t i=0; i<rock.size(); i++)
      delete rock[i];
  }

  for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator it = reservoir_variables_.begin(); it != reservoir_variables_.end(); it++) {
    std::vector<DistributionWithTrend *> variable = it->second;
    for (size_t i=0; i<variable.size(); i++)
      delete variable[i];
  }

  // Erik N: Time line is deleted in common data
  //if (time_line_ != NULL)
  //  delete time_line_;

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

//void
//ModelGeneral::CalculateCovariancesFromRockPhysics(const std::vector<DistributionsRock *>           & rock_distribution,
//                                                  const std::vector<float>                         & probability,
//                                                  NRLib::Grid2D<double>                            & param_corr,
//                                                  std::string                                      & errTxt)
//{
//
//  LogKit::LogFormatted(LogKit::Low,"\nGenerating covariances from rock physics\n");
//
//  bool has_trend = false;
//  for (size_t i=0; i<rock_distribution.size(); i++) {
//    std::vector<bool> rock_has_trend = rock_distribution[i]->HasTrend();
//
//    for (int j=0; j<2; j++) {
//      if (rock_has_trend[j] == true)
//        has_trend = true;
//    }
//  }
//
//  if (has_trend == true) {
//
//    std::vector<int> trend_cube_size = trend_cubes_.GetSizeTrendCubes();
//
//    const int nx = trend_cube_size[0];
//    const int ny = trend_cube_size[1];
//    const int nz = trend_cube_size[2];
//
//    int modulus = 100;
//
//    float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
//    float nextMonitor = monitorSize;
//    std::cout
//      << "\n  0%       20%       40%       60%       80%      100%"
//      << "\n  |    |    |    |    |    |    |    |    |    |    |  "
//      << "\n  ^";
//
//
//    // Local storage for summed combined variances
//    NRLib::Grid2D<double> sumVariance(3,3,0);
//
//    int n_samples = 0;
//
//    for (int k = 0; k < nz; k++) {
//      for (int j = 0; j < ny; j++) {
//        for (int i = 0; i < nx; i++) {
//
//          if ( ( (i+1)*(j+1)*(k+1) ) % modulus == 0) {
//
//            std::vector<double> trend_position = trend_cubes_.GetTrendPosition(i,j,k);
//
//            NRLib::Grid2D<double> sigma_sum(3,3,0);
//
//            CalculateCovarianceInTrendPosition(rock_distribution,
//                                               probability,
//                                               trend_position,
//                                               sigma_sum);
//
//
//            for (size_t a=0; a<3; a++){
//              for (size_t b=0; b<3; b++)
//                sumVariance(a,b) += sigma_sum(a,b);
//            }
//
//            n_samples++;
//          }
//        }
//      }
//
//      // Log progress
//      if (k+1 >= static_cast<int>(nextMonitor) && k < nz) {
//        nextMonitor += monitorSize;
//        std::cout << "^";
//        fflush(stdout);
//      }
//    }
//
//    if (n_samples > 0) {
//      for (int i=0; i<3; i++) {
//        for (int j=0; j<3; j++)
//          param_corr(i,j) = sumVariance(i,j)/n_samples;
//      }
//    }
//
//    else
//      errTxt += "Could not build a covariance structure from rock physics.\n";
//  }
//
//  else {
//    std::vector<double> trend_position(2, 0.0);
//
//    NRLib::Grid2D<double> sigma_sum(3,3,0);
//
//    CalculateCovarianceInTrendPosition(rock_distribution,
//                                       probability,
//                                       trend_position,
//                                       sigma_sum);
//
//
//    for (int i=0; i<3; i++) {
//      for (int j=0; j<3; j++)
//        param_corr(i,j) = sigma_sum(i,j);
//    }
//
//  }
//}

//void
//ModelGeneral::CalculateCovarianceInTrendPosition(const std::vector<DistributionsRock *> & rock_distribution,
//                                                 const std::vector<float>               & probability,
//                                                 const std::vector<double>              & trend_position,
//                                                 NRLib::Grid2D<double>                  & sigma_sum) const
//{
//  int number_of_facies = static_cast<int>(rock_distribution.size());
//
//  std::vector<std::vector<double> > expectation_m(number_of_facies);
//
//  for (int f = 0; f < number_of_facies; f++)
//    expectation_m[f] = rock_distribution[f]->GetLogExpectation(trend_position);
//
//  // Sum up for all facies: probability for a facies multiplied with the expectations of (vp, vs, rho) given the facies
//  std::vector<float> expectations(3, 0);
//  for (int f = 0; f < number_of_facies; f++){
//    expectations[0] += static_cast<float>(expectation_m[f][0] * probability[f]);
//    expectations[1] += static_cast<float>(expectation_m[f][1] * probability[f]);
//    expectations[2] += static_cast<float>(expectation_m[f][2] * probability[f]);
//  }
//
//  // Compute combined variance for all facies in the given grid cell.
//  // Calculation of variance for a set of pdfs with given probability in the rock physics model:
//  //
//  // Var(X) = E([X-E(X)]^2) = E([X - E(X|facies)]^2) + E([E(X|facies) -E(X)]^2) = E(Var(X|facies)) + Var(E(X|facies))
//  //        = Sum_{over all facies} (probability of facies * variance given facies) + sum_{over all facies} probability of facies * (expected value given facies - EX)^2,
//  // where EX is the sum of probability of a facies multiplied with expectation of \mu given facies
//  //
//  // For all facies: Summing up expected value of variances and variance of expected values
//
//  for (int f = 0; f < number_of_facies; f++) {
//    NRLib::Grid2D<double> sigma = rock_distribution[f]->GetLogCovariance(trend_position);
//
//    // For all elements in the 3x3 matrix of the combined variance
//    for (size_t a=0; a<3; a++) {
//      for (size_t b=0; b<3; b++) {
//        double sigma_weigth = probability[f] * (sigma(a,b) + (expectation_m[f][a] - expectations[a])*(expectation_m[f][b] - expectations[b]));
//        sigma_sum(a,b)     += sigma_weigth;
//      }
//    }
//  }
//}

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

//void ModelGeneral::ValidateCorrelationMatrix(float               ** C,
//                                             const ModelSettings *  model_settings,
//                                             std::string         &  errTxt)
//{
//  float minAlpha = model_settings->getVarVpMin();
//  float maxAlpha = model_settings->getVarVpMax();
//  float minBeta  = model_settings->getVarVsMin();
//  float maxBeta  = model_settings->getVarVsMax();
//  float minRho   = model_settings->getVarRhoMin();
//  float maxRho   = model_settings->getVarRhoMax();
//
//  float C00      = C[0][0];
//  float C11      = C[1][1];
//  float C22      = C[2][2];
//  float C01      = C[0][1];
//  float C10      = C[1][0];
//  float C02      = C[0][2];
//  float C20      = C[2][0];
//  float C12      = C[1][2];
//  float C21      = C[2][1];
//
//  if (C00 < minAlpha || C00 > maxAlpha) {
//    errTxt += "The prior Vp variance is outside valid range:\n";
//    errTxt += "  Given value   : " + NRLib::ToString(C00) + "\n";
//    errTxt += "  Minimum value : " + NRLib::ToString(minAlpha) + "\n";
//    errTxt += "  Maximum value : " + NRLib::ToString(maxAlpha) + "\n";
//  }
//  if (C11 < minBeta || C11 > maxBeta) {
//    errTxt += "The prior Vs variance is outside valid range:\n";
//    errTxt += "  Given value   : " + NRLib::ToString(C11) + "\n";
//    errTxt += "  Minimum value : " + NRLib::ToString(minBeta) + "\n";
//    errTxt += "  Maximum value : " + NRLib::ToString(maxBeta) + "\n";
//  }
//  if (C22 < minRho || C22 > maxRho) {
//    errTxt += "The prior density variance is outside valid range:\n";
//    errTxt += "  Given value   : " + NRLib::ToString(C22) + "\n";
//    errTxt += "  Minimum value : " + NRLib::ToString(minRho) + "\n";
//    errTxt += "  Maximum value : " + NRLib::ToString(maxRho) + "\n";
//  }
//
//  float corr01 = C01/(std::sqrt(C00)*std::sqrt(C11));
//  float corr02 = C02/(std::sqrt(C00)*std::sqrt(C22));
//  float corr12 = C12/(std::sqrt(C11)*std::sqrt(C22));
//
//  if (corr01 < -1.0 || corr01 > 1.0) {
//    errTxt += "The prior Vp-Vs correlation is illegal (" + NRLib::ToString(corr01) + ")\n";
//  }
//  if (corr02 < -1.0 || corr02 > 1.0) {
//    errTxt += "The prior Vp-Rho correlation is illegal (" + NRLib::ToString(corr02) + ")\n";
//  }
//  if (corr12 < -1.0 || corr12 > 1.0) {
//    errTxt += "The prior Vs-Rho correlation is illegal (" + NRLib::ToString(corr12) + ")\n";
//  }
//
//  if (std::abs(C01 - C10) > 0.0f) {
//    errTxt += "The prior covariance matrix is not symmetric in Vp and Vs\n";
//    errTxt += "  Corr(Vp,Vs) : " + NRLib::ToString(C01) + "\n";
//    errTxt += "  Corr(Vs,Vp) : " + NRLib::ToString(C10) + "\n";
//  }
//  if (std::abs(C02 - C20) > 0.0f) {
//    errTxt += "The prior covariance matrix is not symmetric in Vp and Rho\n";
//    errTxt += "  Corr(Vp,Rho) : " + NRLib::ToString(C02) + "\n";
//    errTxt += "  Corr(Rho,Vp) : " + NRLib::ToString(C20) + "\n";
//  }
//  if (std::abs(C12 - C21) > 0.0f) {
//    errTxt += "The prior covariance matrix is not symmetric in Vs and Rho\n";
//    errTxt += "  Corr(Vs,Rho) : " + NRLib::ToString(C12) + "\n";
//    errTxt += "  Corr(Rho,Vs) : " + NRLib::ToString(C21) + "\n";
//  }
//}

//void
//ModelGeneral::EstimateCorrXYFromSeismic(Surface *& corrXY,
//                                        FFTGrid ** seisCube,
//                                        int numberOfAngles)
//{
//  FFTGrid * transf;
//  float   * grid;
//
//  int n = static_cast<int>(corrXY->GetNI()*corrXY->GetNJ());
//  grid = new float[n];
//
//  for (int i=0 ; i<n ; i++)
//    grid[i] = 0.0;
//
//  for (int i=0 ; i<numberOfAngles ; i++)
//  {
//    if (seisCube[i]->isFile())
//      transf = new FFTFileGrid(static_cast<FFTFileGrid *>(seisCube[i])); //move new out of loop? Copy grid instead
//    else
//      transf = new FFTGrid(seisCube[i]); //move new out of loop? Copy grid instead
//
//    transf->setAccessMode(FFTGrid::RANDOMACCESS);
//    transf->fftInPlace();
//    transf->square();
//    transf->invFFTInPlace();
//    transf->collapseAndAdd( grid ); //the result of the collapse (the result for z=0) is is added to grid
//    transf->endAccess();
//    delete transf;
//  }
//  float sill = grid[0];
//  for (int i=0;i<n;i++)
//    (*corrXY)(i) = grid[i]/sill;
//  delete [] grid;
//}

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

bool
ModelGeneral::Do4DRockPhysicsInversion(ModelSettings* model_settings)
{

  std::vector<FFTGrid*> predictions = state4d_.doRockPhysicsInversion(*time_line_, rock_distributions_.begin()->second,  time_evolution_);
  int nParamOut =predictions.size();

  std::vector<std::string> labels(nParamOut);

  int i=0;

  for (std::map<std::string, std::vector<DistributionWithTrend *> >::iterator it = reservoir_variables_.begin(); it != reservoir_variables_.end(); it++)
  {
    labels[i] = it->first;
    i++;
  }

  std::string  outPre =  "mu_";

  for (int i=0;i<nParamOut;i++)
  {
     std::string fileName;
     fileName= outPre + labels[i];

     WriteToFile(simbox_, time_depth_mapping_, model_settings, predictions[i] , fileName, labels[i]);
  }

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

  /*
  tag.str(std::string());tag.clear();label = "mean_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings,  current_state.GetMuBeta(), fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetMuRho() , fileName, tag.str() ,true);
  // */
  // write sigma current
  tag.str(std::string());tag.clear();label = "cov_vp_vp_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, current_state.GetCovVp() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovAlphaBeta() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovAlphaRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCovBeta() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_current_step_"; tag << label << timestep << identifyer ; fileName= tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCrCovBetaRho() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_current_step_"; tag << label << timestep << identifyer ; fileName=  tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, current_state.GetCovRho() , fileName,  tag.str(),true);
  // */
  if (transformHere)
    current_state.FFTAllGrids();
}

void
ModelGeneral::Dump4Dparameters(ModelSettings* model_settings, std::string identifyer, int timestep)
{
  state4d_.iFFT();

  std::string  outPath =  "";
  std::string  label;
  std::string fileName;
  std::stringstream tag;

  // write mu static
  tag.str(std::string());tag.clear();label = "mean_vp_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVpStatic() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "mean_vs_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVsStatic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuRhoStatic() , fileName, tag.str(),true);
  // */
  // write mu dynamic
  tag.str(std::string());tag.clear();label = "mean_vp_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getMuVpDynamic() , fileName, tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "mean_vs_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuVsDynamic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "mean_rho_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getMuRhoDynamic() , fileName,  tag.str(),true);
  // */


  // write sigma static - static
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpStaticStatic() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoStaticStatic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_static_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoStaticStatic() , fileName,  tag.str(),true);
     // */
  // write sigma dynamic - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpDynamicDynamic() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsDynamicDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoDynamicDynamic(), fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsDynamicDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoDynamicDynamic(), fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_dynamic_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoDynamicDynamic() , fileName,  tag.str(),true);
  // */
  // write sigma static - dynamic
  tag.str(std::string());tag.clear();label = "cov_vp_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  WriteToFile(simbox_, time_depth_mapping_, model_settings, state4d_.getCovVpVpStaticDynamic() , fileName,  tag.str(),true);

  /*
  tag.str(std::string());tag.clear();label = "cov_vp_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpVsStaticDynamic() , fileName, tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vp_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVpRhoStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVpStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsVsStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_vs_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovVsRhoStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_vp_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoVpStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_vs_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoVsStaticDynamic() , fileName,  tag.str(),true);
  tag.str(std::string());tag.clear();label = "cov_rho_rho_static_dynamic_step_"; tag << label << timestep << identifyer ; fileName= outPath + tag.str();
  ParameterOutput::writeToFile(simbox_,this, model_settings, state4d_.getCovRhoRhoStaticDynamic() , fileName,  tag.str(),true);
   // */
  state4d_.FFT();
}

//void
//ModelGeneral::MakeCorr2DPositiveDefinite(Surface         * corrXY)
//{
//  int      nxp    = corrXY->GetNI();
//  int      nyp    = corrXY->GetNJ();
//  FFTGrid  helper = FFTGrid(nxp,nyp,1,nxp,nyp,1);
//  helper.createRealGrid();
//  helper.setType(FFTGrid::COVARIANCE);
//
//  for (int i =0;i<nxp;i++)
//    for (int j =0;j<nyp;j++)
//    {
//      float value = float((*corrXY)(i+j*nxp));
//      helper.setRealValue(i,j,0,value);
//    }
//
//  helper.fftInPlace();
//  int cnxp =helper.getCNxp();
//  for (int i =0;i<cnxp;i++)
//    for (int j =0;j<nyp;j++)
//    {
//      fftw_complex value;
//      value=helper.getComplexValue(i,j,0);
//      value.re=std::sqrt(value.re*value.re+value.im*value.im);
//      value.im=0.0f;
//      helper.setComplexValue(i,j,0,value);
//    }
//
//  helper.invFFTInPlace();
//  double scale=1.0/double(helper.getRealValue(0,0,0));
//
//  printf("\nFix in latteral correlation in CRAVA results in a variance increase of %f %% (of 100%%) \n",(scale-1.0)*100);
//
//  for (int i =0;i<nxp;i++)
//    for (int j =0;j<nyp;j++)
//       (*corrXY)(i+j*nxp)=helper.getRealValue(i,j,0)*scale;
//}

void
ModelGeneral::WriteToFile(const Simbox        * simbox,
                          GridMapping         * time_depth_mapping,
                          const ModelSettings * model_settings,
                          FFTGrid             * grid,
                          const std::string   & file_name,
                          const std::string   & sgri_label,
                          bool                  padding)
{
  //GridMapping * timeDepthMapping = modelGeneral->GetTimeDepthMapping();
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