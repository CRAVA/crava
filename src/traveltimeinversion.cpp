/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/traveltimeinversion.h"
#include "src/modeltraveltimestatic.h"
#include "src/modeltraveltimedynamic.h"
#include "src/seismicparametersholder.h"
#include "src/simbox.h"
#include "src/modelgeneral.h"
#include "src/rmstrace.h"
#include "src/kriging2d.h"
#include "src/krigingdata2d.h"
#include "src/covgrid2d.h"
#include "src/definitions.h"
#include "src/gridmapping.h"
#include "lib/lib_matr.h"
#include "nrlib/flens/nrlib_flens.hpp"

#include "lib/timekit.hpp"


TravelTimeInversion::TravelTimeInversion(ModelGeneral            * modelGeneral,
                                         ModelTravelTimeStatic   * modelTravelTimeStatic,
                                         ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                                         SeismicParametersHolder & seismicParameters)
{
  LogKit::WriteHeader("Building Stochastic Travel Time Inversion Model");
  time_t time_start;
  time_t time_end;
  time(&time_start);

  int  this_time_lapse    = modelTravelTimeDynamic->getThisTimeLapse();
  bool horizon_data_given = modelTravelTimeDynamic->getHorizonDataGiven();
  bool rms_data_given     = modelTravelTimeDynamic->getRMSDataGiven();

  if (horizon_data_given == true && this_time_lapse > 0)
  {
    doHorizonInversion(modelGeneral,
                       modelTravelTimeStatic,
                       modelTravelTimeDynamic,
                       seismicParameters);
  }

  if (rms_data_given == true){
     doRMSInversion(modelGeneral,
                     modelTravelTimeStatic,
                     modelTravelTimeDynamic,
                     seismicParameters);
  }

  time(&time_end);
  LogKit::LogFormatted(LogKit::DebugLow, "\nTime elapsed :  %d\n", time_end - time_start);

}
//-----------------------------------------------------------------------------------------//

TravelTimeInversion::~TravelTimeInversion()
{
}


//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::doHorizonInversion(ModelGeneral            * modelGeneral,
                                        ModelTravelTimeStatic   * modelTravelTimeStatic,
                                        ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                                        SeismicParametersHolder & seismicParameters) const
{
  std::string text = "\nInverting push down data:";
  LogKit::LogFormatted(LogKit::Low, text);

  const State4D & state_4D = modelGeneral->getState4D();

  FFTGrid * mu_log_vp_dynamic  = new FFTGrid(state_4D.getMuVpDynamic());
  FFTGrid * cov_log_vp_dynamic = new FFTGrid(state_4D.getCovVpVpDynamicDynamic());

  mu_log_vp_dynamic ->invFFTInPlace();
  cov_log_vp_dynamic->invFFTInPlace();

  std::vector<double>   cov_log_vp   = getCovLogVp(cov_log_vp_dynamic);
  NRLib::Grid2D<double> Sigma_log_vp = generateSigmaModel(cov_log_vp);

  const std::vector<Surface>     initial_horizons      = modelTravelTimeStatic->getInitialHorizons();
  const std::vector<std::string> initial_horizon_names = modelTravelTimeStatic->getInitialHorizonNames();
  const std::vector<Surface>     push_down_horizons    = modelTravelTimeDynamic->getPushDownHorizons();
  const std::vector<std::string> push_down_names       = modelTravelTimeDynamic->getHorizonNames();
  const std::vector<double>      standard_deviation    = modelTravelTimeDynamic->getHorizonStandardDeviation();

  std::vector<Surface> sorted_initial_horizons = sortHorizons(initial_horizons,
                                                              push_down_horizons,
                                                              initial_horizon_names,
                                                              push_down_names);
  const Simbox * timeSimboxPrev  = modelGeneral->getTimeSimbox();

  const Surface top_simbox  = dynamic_cast<const Surface &> (timeSimboxPrev->GetTopSurface());
  const Surface base_simbox = dynamic_cast<const Surface &> (timeSimboxPrev->GetBotSurface());
  FFTGrid* relativeVelocityPrev = modelGeneral->getRelativeVelocity();

  int nx      = timeSimboxPrev->getnx();
  int ny      = timeSimboxPrev->getny();
  int nzp     = mu_log_vp_dynamic->getNzp();

  float monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  std::vector<KrigingData2D> mu_log_vp_post(nzp);
  NRLib::Matrix cov_mat_post(nzp, nzp);//
   for (int j = 0; j < nzp; j++)
    for (int k = 0; k < nzp; k++)
      cov_mat_post(j, k) = 0.0;

  int n_traces = 0;
  // inversion for parameters
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {

      std::vector<double>   mu_post;
      NRLib::Grid2D<double> Sigma_post;


      do1DHorizonInversion(mu_log_vp_dynamic,
                           Sigma_log_vp,
                           timeSimboxPrev,
                           relativeVelocityPrev,
                           sorted_initial_horizons,
                           push_down_horizons,
                           standard_deviation,
                           top_simbox,
                           base_simbox,
                           i,
                           j,
                           mu_post,
                           Sigma_post,
                           false);


      if (mu_post.size() > 0) {

        setExpectation(i,
                       j,
                       mu_post,
                       mu_log_vp_post);

        addCovarianceMat(Sigma_post, cov_mat_post);

        n_traces++;

      }
      if (i*ny + j  +1 >= static_cast<int>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }

  for (int i = 0; i < nzp; i++)
    for (int j = 0; j < nzp; j++)
    cov_mat_post(i,j) /= n_traces;

  float corrGradI;
  float corrGradJ;
  modelGeneral->getCorrGradIJ(corrGradI, corrGradJ);

  FFTGrid          * stationary_d          = NULL;
  FFTGrid          * stationary_covariance = NULL;
  std::vector<int>   observation_filter;

  generateStationaryDistribution(timeSimboxPrev,
                                 mu_log_vp_post,
                                 cov_log_vp,
                                 cov_mat_post,
                                 n_traces,
                                 modelTravelTimeDynamic->getErrorCorrXY(),
                                 corrGradI,
                                 corrGradJ,
                                 mu_log_vp_dynamic,
                                 stationary_d,
                                 stationary_covariance,
                                 observation_filter);


  FFTGrid * post_mu_log_vp = NULL;
  cov_log_vp_dynamic->fftInPlace();
  mu_log_vp_dynamic->fftInPlace();
  stationary_covariance->fftInPlace();
  stationary_d->fftInPlace();

  calculateExpectation(observation_filter,
                            mu_log_vp_dynamic,
                            cov_log_vp_dynamic,
                            stationary_d,
                            stationary_covariance,
                            post_mu_log_vp);

  FFTGrid * post_cov_log_vp = NULL;
  calculateLogVpCovariance(observation_filter,
                           cov_log_vp_dynamic,
                           stationary_covariance,
                           post_cov_log_vp);

  modelGeneral->updateState4DWithSingleParameter(post_mu_log_vp,
                                                 post_cov_log_vp,
                                                 3);

   mu_log_vp_dynamic->invFFTInPlace();
   post_mu_log_vp   ->invFFTInPlace(); // out
   post_cov_log_vp  ->invFFTInPlace(); // out

 // The resampling is allways done from push down data only (OK)
 // Generate new simbox, and resample expectation grids in State4D
 // inversion for time reference
  text = "\nComputing alignment:";
  LogKit::LogFormatted(LogKit::Low, text);
  std::vector<KrigingData2D> relativeVelocityNew(nzp);

   for (int j = 0; j < nzp; j++)
    for (int k = 0; k < nzp; k++)
      cov_mat_post(j, k) = 0.0;

  monitorSize = std::max(1.0f, static_cast<float>(nx * ny) * 0.02f);
  nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  n_traces = 0;
  double meanVpRelative=0.0;
  // inversion for parameters
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {

      std::vector<double>   mu_post;
      NRLib::Grid2D<double> Sigma_post;
      NRLib::Grid2D<double> Sigma_log_post;

      do1DHorizonInversion(relativeVelocityPrev,
                           Sigma_log_vp,
                           timeSimboxPrev,
                           relativeVelocityPrev,
                           sorted_initial_horizons,
                           push_down_horizons,
                           standard_deviation,
                           top_simbox,
                           base_simbox,
                           i,
                           j,
                           mu_post,
                           Sigma_log_post,
                           true);


      if (mu_post.size() > 0) {

        setExpectation(i,
                       j,
                       mu_post,
                       relativeVelocityNew);

        getSigmaFromLogSigma(mu_post,Sigma_log_post,Sigma_post);
        addCovarianceMat(Sigma_post, cov_mat_post);
        double tracesum =0.0;
        for(int k=0;k<nzp;k++)
          tracesum +=mu_post[k];
        meanVpRelative+= tracesum/nzp;
        n_traces++;

      }
      if (i*ny + j  +1 >= static_cast<int>(nextMonitor)) {
        nextMonitor += monitorSize;
        std::cout << "^";
        fflush(stdout);
      }
    }
  }
  meanVpRelative/=n_traces;
  std::vector<double>  cov_vp(nzp);

  for (int i = 0; i < nzp; i++){
    cov_vp[i]= (exp(cov_log_vp[i])-1)*meanVpRelative*meanVpRelative;

    for (int j = 0; j < nzp; j++)
      cov_mat_post(i,j) /= n_traces;
  }
  modelGeneral->getCorrGradIJ(corrGradI, corrGradJ);


  generateStationaryDistribution(timeSimboxPrev,
                                 relativeVelocityNew, // not Log transformed
                                 cov_vp,    // not Log transformed
                                 cov_mat_post,  // not Log transformed
                                 n_traces,
                                 modelTravelTimeDynamic->getErrorCorrXY(),
                                 corrGradI,
                                 corrGradJ,
                                 relativeVelocityPrev,  // not Log transformed
                                 stationary_d,          // not Log transformed
                                 stationary_covariance, // not Log transformed
                                 observation_filter);   // not Log transformed


  FFTGrid * relativeVelocityGridNew = NULL;

  relativeVelocityPrev->fftInPlace();
  stationary_covariance->fftInPlace();
  stationary_d->fftInPlace();
  cov_log_vp_dynamic->invFFTInPlace();
  FFTGrid* cov_vp_rel = getCovFunkFromLogCovFunk(cov_log_vp_dynamic,meanVpRelative);
  cov_vp_rel->fftInPlace();

  calculateExpectation(observation_filter,
                            relativeVelocityPrev,
                            cov_vp_rel,
                            stationary_d,
                            stationary_covariance,
                            relativeVelocityGridNew);


    relativeVelocityGridNew->invFFTInPlace();
    relativeVelocityPrev->invFFTInPlace();
    relativeVelocityGridNew->writeAsciiFile("relativeVelocityGridNew.ascii");
    relativeVelocityPrev->writeAsciiFile("relativeVelocityPrev.ascii");

    NRLib::Grid<double> divided_grid = calculateRelativeVelocityUpdate(relativeVelocityGridNew,relativeVelocityPrev);
    // divided_grid = (Vp_current/Vp_previous) in the previous timeframe


    Simbox * new_simbox = NULL;
    std::string errTxt  = "";
    generateNewSimbox(divided_grid,
                      modelTravelTimeStatic->getLzLimit(),
                      timeSimboxPrev,
                      new_simbox,
                      errTxt);

    //NRLib::Grid<double> distance;
    //calculateDistanceGrid(new_simbox,
    //                      divided_grid,
    //                      distance);

    // Resample grid (defined below) tells where in the old simbox we should look for a value for the given cell in the new simbox;
    // uses centerpoint of cell.
    NRLib::Grid<double> resample_grid;
    generateResampleGrid(divided_grid,
                         timeSimboxPrev,
                         new_simbox,
                         resample_grid);


    FFTGrid * mu_log_vp_static   = new FFTGrid(state_4D.getMuVpStatic());
    FFTGrid * mu_log_vs_static   = new FFTGrid(state_4D.getMuVsStatic());
    FFTGrid * mu_log_vs_dynamic  = new FFTGrid(state_4D.getMuVsDynamic());
    FFTGrid * mu_log_rho_static  = new FFTGrid(state_4D.getMuRhoStatic());
    FFTGrid * mu_log_rho_dynamic = new FFTGrid(state_4D.getMuRhoDynamic());

    mu_log_vp_static  ->invFFTInPlace();
    mu_log_vs_static  ->invFFTInPlace();
    mu_log_rho_static ->invFFTInPlace();
    mu_log_vs_dynamic ->invFFTInPlace();
    mu_log_rho_dynamic->invFFTInPlace();

    resampleState4D(resample_grid,
                    timeSimboxPrev,
                    mu_log_vp_static,
                    mu_log_vs_static,
                    mu_log_rho_static,
                    mu_log_vp_dynamic,
                    mu_log_vs_dynamic,
                    mu_log_rho_dynamic);

    mu_log_vp_static  ->fftInPlace();
    mu_log_vs_static  ->fftInPlace();
    mu_log_rho_static ->fftInPlace();
    mu_log_vp_dynamic ->fftInPlace();
    mu_log_vs_dynamic ->fftInPlace();
    mu_log_rho_dynamic->fftInPlace();

    modelGeneral->updateState4DMu(mu_log_vp_static,
                                  mu_log_vs_static,
                                  mu_log_rho_static,
                                  mu_log_vp_dynamic,
                                  mu_log_vs_dynamic,
                                  mu_log_rho_dynamic);


    FFTGrid * relativeVelocity20_2 = generateResampleAveragePreserve(relativeVelocityGridNew,   //  v2v0  contains Vp_2(t1)/Vp_0(t1)    in the previous timeframe (t1)
                                                     relativeVelocityPrev,   //  v1v0_1  contains Vp_1(t1)/Vp_0(t1)  in the previous timeframe (t1)
                                                     timeSimboxPrev,  //  simbox in the previous timeframe (t1)
                                                     new_simbox);//  simbox in the current timeframe (t2)

    modelGeneral->updateState4DAllignment(relativeVelocity20_2);

    modelGeneral->setTimeSimbox(new_simbox);

    delete new_simbox;
    delete mu_log_vp_static;
    delete mu_log_vs_static;
    delete mu_log_vs_dynamic;
    delete mu_log_rho_static;
    delete mu_log_rho_dynamic;
    delete relativeVelocityGridNew;
    delete relativeVelocity20_2;

  modelGeneral->mergeState4D(seismicParameters);

  seismicParameters.invFFTAllGrids();

  delete stationary_d;
  delete stationary_covariance;
  delete post_mu_log_vp;
  delete post_cov_log_vp;
  delete mu_log_vp_dynamic;
  delete cov_log_vp_dynamic;


  delete cov_vp_rel;
}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::doRMSInversion(ModelGeneral            * modelGeneral,
                                    ModelTravelTimeStatic   * modelTravelTimeStatic,
                                    ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                                    SeismicParametersHolder & seismicParameters) const
{
  std::string text = "\nInverting RMS data:";
  LogKit::LogFormatted(LogKit::Low, text);

  FFTGrid * mu_log_vp_grid                 = seismicParameters.GetMuAlpha();
  FFTGrid * cov_log_vp_grid                = seismicParameters.GetCovAlpha();

  const Simbox * timeSimbox                = modelGeneral->getTimeSimbox();
  const Simbox * simbox_above              = modelTravelTimeStatic ->getSimboxAbove();
  const Simbox * simbox_below              = modelTravelTimeDynamic->getSimboxBelow();
  const std::vector<RMSTrace *> rms_traces = modelTravelTimeDynamic->getRMSTraces();
  const double mu_vp_top                   = modelTravelTimeStatic ->getMeanVpTop();
  const double mu_vp_base                  = modelTravelTimeStatic ->getMeanVpBase();
  const double var_vp_above                = modelTravelTimeStatic ->getVarVpAbove();
  const double var_vp_below                = modelTravelTimeStatic ->getVarVpBelow();
  const double range_above                 = modelTravelTimeStatic ->getRangeAbove();
  const double range_below                 = modelTravelTimeStatic ->getRangeBelow();
  const double standard_deviation          = modelTravelTimeDynamic->getRMSStandardDeviation();
  const int    this_time_lapse             = modelTravelTimeDynamic->getThisTimeLapse();

  double dtAbove                           = simbox_above->getMinDz();
  double dtBelow                           = simbox_below->getMinDz();

  int n_above                              = simbox_above->getnz();
  int n_below                              = simbox_below->getnz();
  int n_pad_above                          = FFTGrid::findClosestFactorableNumber(n_above + static_cast<int>(std::ceil(range_above / dtAbove)));
  int n_pad_below                          = FFTGrid::findClosestFactorableNumber(n_below + static_cast<int>(std::ceil(range_below / dtBelow)));
  int n_pad_model                          = mu_log_vp_grid->getNzp();

  int n_rms_traces                         = static_cast<int>(rms_traces.size());

  float corrGradI;
  float corrGradJ;
  modelGeneral->getCorrGradIJ(corrGradI, corrGradJ);

  const Surface * errorCorrXY = modelTravelTimeDynamic->getErrorCorrXY();

  double  dt_max_above    = simbox_above->getdz();
  double  dt_max_below    = simbox_below->getdz();
  Vario * variogram_above = new GenExpVario(1, static_cast<float>(range_above));
  Vario * variogram_below = new GenExpVario(1, static_cast<float>(range_below));

  NRLib::Grid2D<double> Sigma_vp_above = generateSigmaVp(dt_max_above, n_pad_above, var_vp_above, variogram_above);
  NRLib::Grid2D<double> Sigma_vp_below = generateSigmaVp(dt_max_below, n_pad_below, var_vp_below, variogram_below);

  delete variogram_above;
  delete variogram_below;

  FFTGrid * mu_log_vp_above    = NULL;
  FFTGrid * Sigma_log_vp_above = NULL;
  generateMuSigmaLogVpAbove(n_above,
                            n_pad_above,
                            mu_vp_top,
                            Sigma_vp_above,
                            modelGeneral->getPriorCorrXY(),
                            corrGradI,
                            corrGradJ,
                            mu_log_vp_grid,
                            mu_log_vp_above,
                            Sigma_log_vp_above);

  std::vector<double> cov_log_vp_above = getCovLogVp(Sigma_log_vp_above);
  std::vector<double> cov_log_vp_model = getCovLogVp(cov_log_vp_grid);

  float monitorSize = std::max(1.0f, static_cast<float>(n_rms_traces) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  //std::vector<double>        cov_circulant_above(n_pad_above, 0);  // Will only be used at the first time lapse to calculate depth at top reservoir
  std::vector<KrigingData2D> mu_log_vp_post_above(n_pad_above);    // Will only be used at the first time lapse to calculate depth at top reservoir
  NRLib::Matrix cov_mat_log_vp_above_post(n_pad_above, n_pad_above);// Will only be used at the first time lapse to calculate depth at top reservoir
   for (int j = 0; j < n_pad_above; j++)
    for (int k = 0; k < n_pad_above; k++)
      cov_mat_log_vp_above_post(j, k) = 0.0;


 // std::vector<double>        cov_circulant_model(n_pad_model, 0);
  std::vector<KrigingData2D> mu_log_vp_post_model(n_pad_model);
  NRLib::Matrix cov_mat_log_vp_model_post(n_pad_model, n_pad_model);
  for (int j = 0; j < n_pad_model; j++)
    for (int k = 0; k < n_pad_model; k++)
      cov_mat_log_vp_model_post(j, k) = 0.0;




  for (int i = 0; i < n_rms_traces; i++) {

    std::vector<double>   mu_post;
    NRLib::Grid2D<double> Sigma_post;

    do1DRMSInversion(mu_vp_base,
                     Sigma_vp_below,
                     standard_deviation,
                     rms_traces[i],
                     mu_log_vp_above,
                     mu_log_vp_grid,
                     cov_log_vp_above,
                     cov_log_vp_model,
                     simbox_above,
                     simbox_below,
                     timeSimbox,
                     mu_post,
                     Sigma_post);

    if (this_time_lapse == 0) {

      std::vector<double> mu_above(n_pad_above);
      for (int j = 0; j < n_pad_above; j++)
        mu_above[j] = mu_post[j];

      NRLib::Grid2D<double> cov_above(n_pad_above, n_pad_above);
      for (int j = 0; j < n_pad_above; j++) {
        for (int k = 0; k < n_pad_above; k++)
          cov_above(j, k) = Sigma_post(j, k);
      }

      setExpectation(rms_traces[i]->getIIndex(),
                     rms_traces[i]->getJIndex(),
                     mu_above,
                     mu_log_vp_post_above);

      addCovarianceMat(cov_above, cov_mat_log_vp_above_post);
    }

    std::vector<double> mu_model(n_pad_model);
    for (int j = 0; j < n_pad_model; j++)
      mu_model[j] = mu_post[j + n_pad_above];

    NRLib::Grid2D<double> cov_model(n_pad_model, n_pad_model);
    for (int j = 0; j < n_pad_model; j++) {
      for (int k = 0; k < n_pad_model; k++)
        cov_model(j, k) = Sigma_post(j + n_pad_above, k + n_pad_above);
    }

    setExpectation(rms_traces[i]->getIIndex(),
                   rms_traces[i]->getJIndex(),
                   mu_model,
                   mu_log_vp_post_model);

    addCovarianceMat(cov_model, cov_mat_log_vp_model_post);

    if (i + 1 >= static_cast<int>(nextMonitor)) {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }

  for (int i = 0; i < n_pad_model; i++)
     for (int j = 0; j < n_pad_model; j++)
       cov_mat_log_vp_model_post(i,j) /= n_rms_traces;


  FFTGrid          * stationary_d          = NULL;
  FFTGrid          * stationary_covariance = NULL;
  std::vector<int>   observation_filter;

  generateStationaryDistribution(timeSimbox,
                                 mu_log_vp_post_model,
                                 cov_log_vp_model,
                                 cov_mat_log_vp_model_post,
                                 n_rms_traces,
                                 errorCorrXY,
                                 corrGradI,
                                 corrGradJ,
                                 mu_log_vp_grid,
                                 stationary_d,
                                 stationary_covariance,
                                 observation_filter);



    seismicParameters.FFTAllGrids();

    calculateFullPosteriorModel(observation_filter,
                                seismicParameters,
                                stationary_d,
                                stationary_covariance);

    if (this_time_lapse == 0 && modelGeneral->getTimeDepthMapping() == NULL) {
      // Calculate time/depth mapping

      for (int i = 0; i < n_pad_above; i++)
        for (int j = 0; j < n_pad_above; j++)
          cov_mat_log_vp_above_post(i,j) /= n_rms_traces;

      FFTGrid          * stationary_d_above          = NULL;
      FFTGrid          * stationary_covariance_above = NULL;
      std::vector<int>   observation_filter_above;



      generateStationaryDistribution(simbox_above,
                                     mu_log_vp_post_above,
                                     cov_log_vp_above,
                                     cov_mat_log_vp_above_post,
                                     n_rms_traces,
                                     errorCorrXY,
                                     corrGradI,
                                     corrGradJ,
                                     mu_log_vp_above,
                                     stationary_d_above,
                                     stationary_covariance_above,
                                     observation_filter_above);

      stationary_covariance_above->fftInPlace();
      mu_log_vp_above->fftInPlace();
      Sigma_log_vp_above         ->fftInPlace();
      stationary_d_above->fftInPlace();

      FFTGrid * post_mu_log_vp_above = NULL;
      calculateExpectation(observation_filter_above,
                                mu_log_vp_above,
                                Sigma_log_vp_above,
                                stationary_d_above,
                                stationary_covariance_above,
                                post_mu_log_vp_above);

      FFTGrid * post_cov_log_vp_above = NULL;

      calculateLogVpCovariance(observation_filter_above,
                               Sigma_log_vp_above,
                               stationary_covariance_above,
                               post_cov_log_vp_above);

      post_mu_log_vp_above ->invFFTInPlace();
      post_cov_log_vp_above->invFFTInPlace();

      GridMapping * time_depth_mapping = NULL;

      mu_log_vp_grid->invFFTInPlace();
      cov_log_vp_grid->invFFTInPlace();

      generateTimeDepthMapping(post_mu_log_vp_above,
                               post_cov_log_vp_above,
                               mu_log_vp_grid,
                               cov_log_vp_grid,
                               modelTravelTimeStatic->getOutputGridFormat(),
                               simbox_above,
                               timeSimbox,
                               time_depth_mapping);

      mu_log_vp_grid->fftInPlace();
      cov_log_vp_grid->fftInPlace();

      delete stationary_d_above;
      delete stationary_covariance_above;
      delete post_mu_log_vp_above;
      delete post_cov_log_vp_above;

      modelGeneral->setTimeDepthMapping(time_depth_mapping); // delete time_depth_mapping in modelGeneral

  // }

    modelGeneral->updateState4D(seismicParameters);
  }

  delete mu_log_vp_above;
  delete Sigma_log_vp_above;
  delete stationary_d;
  delete stationary_covariance;

}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::do1DHorizonInversion(FFTGrid                     * mu_prior,
                                          const NRLib::Grid2D<double> & Sigma_log_vp_prior,
                                          const Simbox                * timeSimbox,
                                          const FFTGrid               * relativeVelocity,
                                          const std::vector<Surface>  & initial_horizons,
                                          const std::vector<Surface>  & push_down_horizons,
                                          const std::vector<double>   & standard_deviation,
                                          const Surface               & top_simbox,
                                          const Surface               & base_simbox,
                                          int                           i_ind,
                                          int                           j_ind,
                                          std::vector<double>         & mu_post_out,
                                          NRLib::Grid2D<double>       & Sigma_post_log_vp,
                                          bool                          velosityForTimeToTimeMap) const
{
  double x, y;
  timeSimbox->getXYCoord(i_ind, j_ind, x, y);

  // Get data in current position
  int                 n_horizons = static_cast<int>(initial_horizons.size());
  std::vector<double> time_P0(n_horizons);
  std::vector<double> push_down(n_horizons);
  double              top  = 0;

  for (int k = 0; k < n_horizons; ++k) {

    double base = 0;

    time_P0[k]   = initial_horizons[k].GetZ(x, y);
    push_down[k] = push_down_horizons[k].GetZ(x, y);

    top  = top_simbox.GetZ(x, y);
    base = base_simbox.GetZ(x, y);
    if (base < time_P0[k])
      time_P0[k] = base;
  }

  double missing_value = initial_horizons[0].GetMissingValue();
  int    n_nonmissing  = 0;

  for (int k = 0; k < n_horizons; k++) {
    if (time_P0[k] != missing_value && push_down[k] != missing_value)
      n_nonmissing++;
  }


  if (n_nonmissing > 0) {// if data in trace
    std::vector<double>   d(n_nonmissing);
    NRLib::Grid2D<double> Sigma_d(n_nonmissing, n_nonmissing, 0);

    int dataInd=0;
    for (int k = 0; k < n_horizons; ++k) {
      if (time_P0[k] != missing_value && push_down[k] != missing_value) {
        d[dataInd]          = time_P0[k] + push_down[k] - top;
        Sigma_d(dataInd,dataInd) = std::pow(standard_deviation[k], 2);
        dataInd++;
      }
      else
      {
        double dummy=0.0; // NBNB OK to check during debug
      }
    }

    double dt  = timeSimbox->getdz(i_ind, j_ind); // NBNB Check this  for original(static) vs current.
    int    nz  = timeSimbox->getnz();
    int    nzp = mu_prior->getNzp();
    std::vector<double> relativeVelocityVec(nzp);
    for(int k=0;k<nzp;k++)
      relativeVelocityVec[k] = relativeVelocity->getRealValue(i_ind,j_ind,k);



    NRLib::Grid2D<double> G = calculateGHorizon(dt,     // from previous timestep
                                                relativeVelocityVec, //from previous timestep (velocity  relative to initial time )
                                                top,        //  common all timesteps
                                                missing_value,
                                                n_nonmissing,
                                                nz,
                                                nzp,
                                                time_P0,   //  time of base syurvey
                                                push_down);//  at current timestep

    std::vector<double> mu_log_vp;
    std::vector<double> mu_vp;

    if(velosityForTimeToTimeMap==false)
    {
      mu_log_vp = generateMuFromGrid(mu_prior, i_ind, j_ind);
    }else
    {
       mu_vp = generateMuFromGrid(mu_prior, i_ind, j_ind);
       getLogMuFromMu(mu_log_vp , mu_vp ,Sigma_log_vp_prior);
    }

    // Transform to Vp^(-1)
    std::vector<double>   mu_vp_minus;
    NRLib::Grid2D<double> Sigma_vp_minus;

    calculateMinusFirstCentralMomentLogNormal(mu_log_vp,
                                              Sigma_log_vp_prior,
                                              mu_vp_minus,
                                              Sigma_vp_minus);
    if(velosityForTimeToTimeMap==true)
    {
      for(int i=0;i<nzp;i++)
        mu_vp_minus[i]=1.0/mu_vp[i];
    }

    std::vector<double>   mu_post;
    NRLib::Grid2D<double> Sigma_post;

    calculatePosteriorModel(d,
                            Sigma_d,
                            mu_vp_minus,
                            Sigma_vp_minus,
                            G,
                            mu_post,
                            Sigma_post);

     /* printout
  int  n_layers=nzp;
  NRLib::Matrix Sigma_m1(n_layers, n_layers);
  NRLib::Vector Mu_m1(n_layers);

  for (int i = 0; i < n_layers; i++) {
    Mu_m1(i)=mu_vp_minus[i];
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_vp_minus(i,j);
  }
  NRLib::WriteVectorToFile("Mu_prior.dat",Mu_m1);
  NRLib::WriteMatrixToFile("Sigma_prior.dat",Sigma_m1);
  for (int i = 0; i < n_layers; i++) {
    Mu_m1(i)=mu_post[i];
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post(i,j);
  }
  NRLib::WriteVectorToFile("Mu_post.dat",Mu_m1);
  NRLib::WriteMatrixToFile("Sigma_post.dat",Sigma_m1);

  int nData =d.size();
  NRLib::Matrix Gmat(nData,n_layers);
  for (int i = 0; i < nData; i++) {

    for (int j = 0; j < n_layers; j++)
      Gmat(i,j) = G(i,j);
  }
   NRLib::WriteMatrixToFile("GHor.dat",Gmat);
  // */

  transformCovarianceToRelativeScale(mu_post,Sigma_post, 0,0,nz,nzp,0, 0);

    /* printout
  for (int i = 0; i < n_layers; i++) {
     Mu_m1(i)=mu_post[i];
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post(i,j);
  }
  NRLib::WriteVectorToFile("Mu_post_trans.dat",Mu_m1);
  NRLib::WriteMatrixToFile("Sigma_post_trans.dat",Sigma_m1);
   // */


  // Transform to log(Vp)

    if(velosityForTimeToTimeMap==false)
    {
      transformVpMinusToLogVp(mu_post,
                            Sigma_post,
                            mu_post_out,
                            Sigma_post_log_vp);
    }else
    {
      std::vector<double>  mu_post_log_vp;
      transformVpMinusToLogVp(mu_post,
                            Sigma_post,
                            mu_post_log_vp,
                            Sigma_post_log_vp);
      mu_post_out.resize(nzp);
      for(int i=0;i<nzp;i++)
         mu_post_out[i]=1.0/mu_post[i];

    }

      /* printout
  for (int i = 0; i < n_layers; i++) {
     Mu_m1(i)=mu_post_out[i];
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post_log_vp(i,j);
  }
  NRLib::WriteVectorToFile("Mu_post_log.dat",Mu_m1);
  NRLib::WriteMatrixToFile("Sigma_post_log.dat",Sigma_m1);
   // */
  }
}



//----------------------------------------------------------------------------------------//
void
TravelTimeInversion::do1DRMSInversion(const double                & mu_vp_base,
                                      const NRLib::Grid2D<double> & Sigma_m_below,
                                      const double                & standard_deviation,
                                      const RMSTrace              * rms_trace,
                                      FFTGrid                     * mu_log_vp_above,
                                      FFTGrid                     * mu_log_vp_model,
                                      const std::vector<double>   & cov_grid_log_vp_above,
                                      const std::vector<double>   & cov_grid_log_vp_model,
                                      const Simbox                * simbox_above,
                                      const Simbox                * simbox_below,
                                      const Simbox                * timeSimbox,
                                      std::vector<double>         & mu_post_log_vp,
                                      NRLib::Grid2D<double>       & Sigma_post_log_vp) const
{

  int i_ind = rms_trace->getIIndex();
  int j_ind = rms_trace->getJIndex();

  double t_top     = timeSimbox->getTop(i_ind, j_ind);
  double t_bot     = timeSimbox->getBot(i_ind, j_ind);

  double dt_above  = simbox_above->getdz(i_ind, j_ind);
  double dt_below  = simbox_below->getdz(i_ind, j_ind);
  double dt_simbox = timeSimbox  ->getdz(i_ind,  j_ind);

  int n_above = simbox_above->getnz();
  int n_model = timeSimbox->getnz();
  int n_below = simbox_below->getnz();

  int n_pad_above = mu_log_vp_above->getNzp();
  int n_pad_model = mu_log_vp_model->getNzp();
  int n_pad_below = static_cast<int>(Sigma_m_below.GetNI());

  const std::vector<double> time = rms_trace->getTime();

  NRLib::Grid2D<double> G = calculateG(time,
                                       t_top,
                                       t_bot,
                                       dt_above,
                                       dt_simbox,
                                       dt_below,
                                       n_above,
                                       n_model,
                                       n_below,
                                       n_pad_above,
                                       n_pad_model,
                                       n_pad_below);

  std::vector<double> mu_log_vp_model_profile = generateMuFromGrid(mu_log_vp_model, i_ind, j_ind);
  std::vector<double> mu_log_vp_above_profile = generateMuFromGrid(mu_log_vp_above, i_ind, j_ind);

  std::vector<double>   mu_m_square;
  NRLib::Grid2D<double> Sigma_m_square;

  calculateMuSigma_mSquare(mu_log_vp_above_profile,
                           mu_log_vp_model_profile,
                           cov_grid_log_vp_above,
                           cov_grid_log_vp_model,
                           mu_vp_base,
                           Sigma_m_below,
                           n_below,
                           mu_m_square,
                           Sigma_m_square);

  std::vector<double> d_square = calculateDSquare(rms_trace->getVelocity());
  biasAdjustDsquare(d_square, standard_deviation * standard_deviation);

  NRLib::Grid2D<double> Sigma_d_square = calculateSigmaDSquare(rms_trace->getVelocity(), standard_deviation);

  std::vector<double>   mu_post;
  NRLib::Grid2D<double> Sigma_post;

  calculatePosteriorModel(d_square,
                          Sigma_d_square,
                          mu_m_square,
                          Sigma_m_square,
                          G,
                          mu_post,
                          Sigma_post);

  /* printout
  int  n_layers=n_pad_above+ n_pad_model+n_pad_below;
  NRLib::Matrix Sigma_m1(n_layers, n_layers);

  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_m_square(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_prior.dat",Sigma_m1);
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_post.dat",Sigma_m1);
  // */
  transformCovarianceToRelativeScale(mu_post,Sigma_post, n_above, n_pad_above,n_model,n_pad_model,n_below, n_pad_below);

  /* printout
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_post_trans.dat",Sigma_m1);
   // */

  transformVpSquareToLogVp(mu_post,
                           Sigma_post,
                           mu_post_log_vp,
                           Sigma_post_log_vp);
  /* printout
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_post_log_vp(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_post_log.dat",Sigma_m1);
   // */

}
//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::calculatePosteriorModel(const std::vector<double>   & d,
                                             const NRLib::Grid2D<double> & Sigma_d,
                                             const std::vector<double>   & mu_m,
                                             const NRLib::Grid2D<double> & Sigma_m,
                                             const NRLib::Grid2D<double> & G,
                                             std::vector<double>         & mu_post,
                                             NRLib::Grid2D<double>       & Sigma_post) const
{
  int n_layers = static_cast<int>(mu_m.size());
  int n_data   = static_cast<int>(d.size());

  NRLib::Vector mu_m1(n_layers);
  for (int i = 0; i < n_layers; i++)
    mu_m1(i) = mu_m[i];

  NRLib::Matrix Sigma_m1(n_layers, n_layers);
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_m1(i,j) = Sigma_m(i,j);
  }

  NRLib::Matrix G1(n_data, n_layers);
  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j < n_layers; j++)
      G1(i,j) = G(i,j);
  }

  NRLib::Matrix G1_transpose(n_layers, n_data);
  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j < n_layers; j++)
      G1_transpose(j,i) = G(i,j);
  }

  NRLib::Vector d1(n_data);
  for (int i = 0; i < n_data; i++)
    d1(i) = d[i];

  NRLib::Matrix Sigma_d1(n_data, n_data);
  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j < n_data; j++)
      Sigma_d1(i,j) = Sigma_d(i,j);
  }

  NRLib::Vector dataMean(n_data);
  NRLib::Vector diff(n_data);
  NRLib::Matrix dataModelCovariance(n_data, n_layers);
  NRLib::Matrix modelDataCovariance(n_layers, n_data);
  NRLib::Matrix dataCovariance(n_data, n_data);

  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j < n_layers; j++)
      dataModelCovariance(i, j) = 0;
  }

  dataMean            = G1 * mu_m1;
  diff                = d1 - dataMean;
  dataModelCovariance = G1 * Sigma_m1;
  modelDataCovariance = Sigma_m1 * G1_transpose;
  dataCovariance      = dataModelCovariance * G1_transpose + Sigma_d1;

  NRLib::SymmetricMatrix data_covariance_inv_sym(n_data);

  for (int i = 0; i < n_data; i++) {
    for (int j = 0; j <= i; j++)
      data_covariance_inv_sym(j,i) = dataCovariance(i,j);
  }

  NRLib::CholeskyInvert(data_covariance_inv_sym);

  NRLib::Matrix data_covariance_inv(n_data, n_data);
  for (int i = 0; i < n_data; i++) {
    for (int j = i; j < n_data; j++) {
      data_covariance_inv(i,j) = data_covariance_inv_sym(i,j);
      data_covariance_inv(j,i) = data_covariance_inv(i,j);
    }
  }

  NRLib::Matrix helpMat(n_layers, n_data);
  NRLib::Vector mu_help(n_layers);
  NRLib::Matrix Sigma_help(n_layers, n_layers);

  NRLib::Vector mu_post1(n_layers);
  NRLib::Matrix Sigma_post1(n_layers, n_layers);

  helpMat     = modelDataCovariance * data_covariance_inv;
  mu_help     = helpMat * diff;
  Sigma_help  = helpMat * dataModelCovariance;

  mu_post1    = mu_m1 + mu_help;
  Sigma_post1 = Sigma_m1 - Sigma_help;

  mu_post.resize(n_layers);
  for (int i = 0; i < n_layers; i++)
    mu_post[i] = mu_post1(i);

  Sigma_post.Resize(n_layers, n_layers);
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_post(i,j) = Sigma_post1(i,j);
  }
  /*
  NRLib::WriteVectorToFile("mu_m1.dat",mu_m1);
  NRLib::WriteVectorToFile("d1.dat",d1);
  NRLib::WriteVectorToFile("dataMean.dat",dataMean);
  NRLib::WriteVectorToFile("mu_post1.dat",mu_post1);
  NRLib::WriteVectorToFile("mu_help.dat", mu_help);
  NRLib::WriteMatrixToFile("Sigma_post.dat",Sigma_post1);
  NRLib::WriteMatrixToFile("G.dat",G1);
  NRLib::WriteMatrixToFile("Sigma_d1.dat",Sigma_d1);
  NRLib::WriteMatrixToFile("G1_transpose.dat",G1_transpose);
  NRLib::WriteMatrixToFile("Sigma_m1.dat",Sigma_m1);
  NRLib::WriteMatrixToFile("helpMat.dat",helpMat);
  NRLib::WriteMatrixToFile("Sigma_help.dat",Sigma_help);
  */
}

//-----------------------------------------------------------------------------------------//
std::vector<double>
TravelTimeInversion::calculateDSquare(const std::vector<double> & d) const
{
  int n = static_cast<int>(d.size());

  std::vector<double> d_square(n);

  for (int i = 0; i < n; i++)
    d_square[i] = std::pow(d[i],2);

  return d_square;
}
//-----------------------------------------------------------------------------------------//


//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::biasAdjustDsquare(std::vector<double> & d2, double variance) const
{
  int n = static_cast<int>(d2.size());

  for (int i = 0; i < n; i++)
    d2[i] -= variance;
}
//----------

NRLib::Grid2D<double>
TravelTimeInversion::calculateGHorizon(double                      dt,  // time increment in  timeframe used. "used time frame"
                                       std::vector<double>         relativeVelocity, // velocity relative to initial timeframe. V_used/V_Initial, reference in "used time frame"
                                       double                      top,
                                       double                      missing_value,
                                       int                         n_nonmissing,
                                       int                         nz,
                                       int                         nzp,
                                       const std::vector<double> & time_P0, // time to horixzon in initial time frame
                                       const std::vector<double> & push_down) const // push down of horizon under consideration (next compared with dt)
{
  int n_horizons = static_cast<int>(push_down.size());


  NRLib::Grid2D<double> G(n_nonmissing, nzp, 0);

  for (int l = 0; l < n_horizons; ++l) {
    int indData=0;
    if (time_P0[l] != missing_value && push_down[l] != missing_value) {
      double time0Start = missing_value;
      double time0End   = top;
      for(int k=0; k<nz;k++)
      {
        time0Start = time0End;
        time0End   += relativeVelocity[k]*dt;
        double p  = std::max(0.0, std::min( (time_P0[l] - time0Start)/(time0End-time0Start)  ,  1.0));
        // returns 1 if time_P0[l] is below cell end,
        // returns 0 if time_P0[l] is above cell start,
        // returns percentage distance from cell top to time_P0[l] otherwise
        G(indData, k)     = p*relativeVelocity[k]*dt;
      }
      indData++;
    }
  }
  // sum over k in G, will give the traveltime to horizon  for time t0, i.e time_P0[l].
  return G;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
TravelTimeInversion::calculateG(const std::vector<double> & rms_time,
                                const double              & t_top,
                                const double              & t_bot,
                                const double              & dt_above,
                                const double              & dt_simbox,
                                const double              & dt_below,
                                const int                 & n_above,
                                const int                 & n_model,
                                const int                 & n_below,
                                const int                 & n_pad_above,
                                const int                 & n_pad_model,
                                const int                 & n_pad_below) const
{
  int n_layers = n_pad_above + n_pad_model + n_pad_below;

  std::vector<double> t(n_layers + 1, 0);
  std::vector<double> dt(n_layers + 1, 0);

  for (int j = 0; j < n_layers + 1; j++) {
    if(j < n_above) {
      t[j]  = j * dt_above;
      dt[j] = dt_above;
    }
    else if(j >= n_pad_above && j < n_pad_above + n_model) {
      t[j]  = t_top + (j - n_pad_above) * dt_simbox;
      dt[j] = dt_simbox;
    }
    else if(j >= n_pad_above + n_pad_model && j < n_pad_above + n_pad_model + n_below) {
      t[j]  = t_bot + (j - n_pad_above - n_pad_model) * dt_below;
      dt[j] = dt_below;
    }
  }

  int n_rms_data = static_cast<int>(rms_time.size());

  NRLib::Grid2D<double> G(n_rms_data, n_layers, 0);

  for (int j = 0; j < n_rms_data; j++) {
    int k=0;
    double prev_base_t = t[0];
    while (rms_time[j] >= t[k] + dt[k] && k < n_layers) { // while rms_time is below the base of the current cell
      G(j,k) = dt[k] / rms_time[j];
      if (t[k] != 0)  // no increment  if this is a padding cell
        prev_base_t = t[k]+dt[k];
      k++;
    }
    if (k < n_layers) { // the rms_time is within the current cell
        G(j,k) = (rms_time[j] - prev_base_t) / rms_time[j];
    }
  }

  return G;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
TravelTimeInversion::calculateSigmaDSquare(const std::vector<double> & rms_velocity,
                                           const double              & standard_deviation) const
{
  int n                 = static_cast<int>(rms_velocity.size());
  const double variance = std::pow(standard_deviation,2);

  NRLib::Grid2D<double> I(n, n, 0);
  for (int i = 0; i < n; i++)
    I(i, i) = 1;

  NRLib::Grid2D<double> Sigma_d_square(n, n, 0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      Sigma_d_square(i, j) = 4 * std::pow(rms_velocity[i], 2) * variance * I(i, j) + 2 * std::pow(variance, 2) * I(i, j);
  }

  return Sigma_d_square;
}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::calculateMuSigma_mSquare(const std::vector<double>   & mu_log_vp_above,
                                              const std::vector<double>   & mu_log_vp_model,
                                              const std::vector<double>   & cov_grid_log_vp_above,
                                              const std::vector<double>   & cov_grid_log_vp_model,
                                              const double                & mu_vp_base,
                                              const NRLib::Grid2D<double> & Sigma_vp_below,
                                              const int                   & n_below,
                                              std::vector<double>         & mu_vp_square,
                                              NRLib::Grid2D<double>       & Sigma_vp_square) const
{

  int n_layers_simbox = static_cast<int>(mu_log_vp_model.size());

  // Model
  NRLib::Grid2D<double> Sigma_log_vp_model = generateSigmaModel(cov_grid_log_vp_model);

  // Above
  NRLib::Grid2D<double> Sigma_log_vp_above = generateSigmaModel(cov_grid_log_vp_above);

  // Below
  int                   n_pad_below    = static_cast<int>(Sigma_vp_below.GetNI());
  std::vector<double>   mu_vp_below    = generateMuVpBelow(std::exp(mu_log_vp_model[n_layers_simbox-1]), mu_vp_base, n_below, n_pad_below);
  std::vector<double>   mu_log_vp_below;
  NRLib::Grid2D<double> Sigma_log_vp_below;
  calculateCentralMomentLogNormalInverse(mu_vp_below, Sigma_vp_below, mu_log_vp_below, Sigma_log_vp_below);

  // Combine
  std::vector<double>   mu_log_m    = generateMuCombined(mu_log_vp_above, mu_log_vp_model, mu_log_vp_below);
  NRLib::Grid2D<double> Sigma_log_m = generateSigmaCombined(Sigma_log_vp_above, Sigma_log_vp_model, Sigma_log_vp_below);

  // Transform to Vp^2
  calculateSecondCentralMomentLogNormal(mu_log_m, Sigma_log_m, mu_vp_square, Sigma_vp_square);

}
//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
TravelTimeInversion::generateSigmaVp(const double & dt,
                                     const int    & n_layers,
                                     const double & var_vp,
                                     const Vario  * variogram) const
{

  // Circulant corrT
  std::vector<double> corrT(n_layers, 0);
  corrT[0] = 1;

  for (int i = 1; i < n_layers / 2; i++) {
    corrT[i] = static_cast<double>(variogram->corr(static_cast<float>(i * dt), 0));
    corrT[n_layers - i] = corrT[i];
  }

  NRLib::Grid2D<double> Sigma_vp = generateSigma(var_vp, corrT);

  return Sigma_vp;
}

//-----------------------------------------------------------------------------------------//
std::vector<double>
TravelTimeInversion::getCovLogVp(FFTGrid * cov_log_vp) const
{
  int n_layers_padding = cov_log_vp->getNzp();

  cov_log_vp->setAccessMode(FFTGrid::READ);

  std::vector<double> cov_grid_log_vp(n_layers_padding, 0);
  for (int j = 0; j < n_layers_padding; j++)
    cov_grid_log_vp[j] = cov_log_vp->getRealValue(0, 0, j, true);

  cov_log_vp->endAccess();

  return cov_grid_log_vp;

}

//-----------------------------------------------------------------------------------------//

std::vector<double>
TravelTimeInversion::generateMuFromGrid(FFTGrid   * muGrid,
                                             const int & i_ind,
                                             const int & j_ind) const
{

  int n_layers_padding = muGrid->getNzp();

  muGrid->setAccessMode(FFTGrid::READ);

  std::vector<double> mu_profile(n_layers_padding, 0);
  for (int j = 0; j < n_layers_padding; j++)
    mu_profile[j] = muGrid->getRealValue(i_ind, j_ind, j, true);

  muGrid->endAccess();

  return mu_profile;

}

//-----------------------------------------------------------------------------------------//

std::vector<double>
TravelTimeInversion::generateMuVpBelow(const double & top_value,
                                       const double & base_value,
                                       const int    & nz,
                                       const int    & nzp) const
{
  std::vector<double> mu_vp = generateMuVp(top_value, base_value, nz);

  for (int i = 0; i < nz; i++)
    mu_vp[i] = mu_vp[i + 1];

  std::vector<double> mu_vp_out(nzp, mu_vp[nz-1]);

  for (int i = 0; i < nz; i++)
    mu_vp_out[i] = mu_vp[i];


  return mu_vp_out;
}

//-----------------------------------------------------------------------------------------//

std::vector<double>
TravelTimeInversion::generateMuVp(const double & top_value,
                                  const double & base_value,
                                  const int    & n_layers) const
{
  std::vector<double> mu_vp(n_layers+1);

  for (int j = 0; j < n_layers + 1; j++)
    mu_vp[j] = top_value + j * (base_value - top_value) / n_layers;

  return mu_vp;
}

//-----------------------------------------------------------------------------------------//

std::vector<double>
TravelTimeInversion::generateMuCombined(const std::vector<double> & mu_above,
                                        const std::vector<double> & mu_model,
                                        const std::vector<double> & mu_below) const
{
  int n_layers_above        = static_cast<int>(mu_above.size());
  int n_layers_below        = static_cast<int>(mu_below.size());
  int n_layers_model        = static_cast<int>(mu_model.size());
  int n_layers              = n_layers_above + n_layers_model + n_layers_below;

  std::vector<double> mu_m(n_layers, RMISSING);

  for (int j = 0; j < n_layers_above; j++)
    mu_m[j] = mu_above[j];

  for (int j = n_layers_above; j < n_layers_above + n_layers_model; j++)
    mu_m[j] = mu_model[j - n_layers_above];

  for (int j = n_layers_above+n_layers_model; j < n_layers; j++)
    mu_m[j] = mu_below[j - n_layers_above - n_layers_model];

  return mu_m;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
TravelTimeInversion::generateSigmaCombined(const NRLib::Grid2D<double> & Sigma_above,
                                           const NRLib::Grid2D<double> & Sigma_model,
                                           const NRLib::Grid2D<double> & Sigma_below) const
{
  int n_layers_above = static_cast<int>(Sigma_above.GetNI());
  int n_layers_below = static_cast<int>(Sigma_below.GetNI());
  int n_layers_model = static_cast<int>(Sigma_model.GetNI());
  int n_layers       = n_layers_above + n_layers_model + n_layers_below;

  NRLib::Grid2D<double> Sigma_m(n_layers, n_layers, 0);

  for (int j = 0; j < n_layers_above; j++) {
    for (int k = j; k < n_layers_above; k++) {
      Sigma_m(j, k) = Sigma_above(j, k);
      Sigma_m(k, j) = Sigma_m(j, k);
    }
  }

  for (int j = n_layers_above; j < n_layers_above+n_layers_model; j++) {
    for (int k = j; k < n_layers_above+n_layers_model; k++) {
      Sigma_m(j, k) = Sigma_model(j - n_layers_above, k - n_layers_above);
      Sigma_m(k, j) = Sigma_m(j, k);
    }
  }

  for (int j = n_layers_above + n_layers_model; j < n_layers; j++) {
    for (int k = j; k < n_layers; k++) {
      Sigma_m(j, k) = Sigma_below(j - n_layers_above - n_layers_model, k - n_layers_above - n_layers_model);
      Sigma_m(k, j) = Sigma_m(j, k);
    }
  }

  return Sigma_m;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
TravelTimeInversion::generateSigma(const double              & var,
                                   const std::vector<double> & corrT) const
{

  int n_layers = static_cast<int>(corrT.size());

  NRLib::Grid2D<double> Sigma(n_layers, n_layers, 0);

  for (int j = 0; j < n_layers; j++) {
    int count = 0;
    for (int k = j; k < n_layers; k++) {
      Sigma(j, k) = var * corrT[count];
      Sigma(k, j) = Sigma(j, k);
      count ++;
    }
  }

  return Sigma;
}

//-----------------------------------------------------------------------------------------//
NRLib::Grid2D<double>
TravelTimeInversion::generateSigmaModel(const std::vector<double> & cov_grid) const
{
  int n = static_cast<int>(cov_grid.size());

  NRLib::Grid2D<double> Sigma_m(n, n);

  for (int j = 0; j < n; j++) {
    int count = 0;
    for (int k = j; k < n; k++) {
      double cov_log_vp = cov_grid[count];
      Sigma_m(j, k) = cov_log_vp;
      Sigma_m(k, j) = Sigma_m(j, k);
      count++;
    }
  }

  return Sigma_m;
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::transformVpSquareToLogVp(const std::vector<double>   & mu_vp_square,
                                              const NRLib::Grid2D<double> & Sigma_vp_square,
                                              std::vector<double>         & mu_log_vp,
                                              NRLib::Grid2D<double>       & Sigma_log_vp) const
{
  int n_layers = static_cast<int>(mu_vp_square.size());

  calculateCentralMomentLogNormalInverse(mu_vp_square, Sigma_vp_square, mu_log_vp, Sigma_log_vp);

  for (int i = 0; i < n_layers; i++)
    mu_log_vp[i] = mu_log_vp[i] / 2;

  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_log_vp(i, j) = Sigma_log_vp(i, j) / 4;
  }
}

void
TravelTimeInversion::getMuFromLogMu(std::vector<double>  & mu_post,const std::vector<double>  & mu_post_log, const NRLib::Grid2D<double>& Sigma_log) const
{
  int n = static_cast<int>(mu_post_log.size());
  mu_post.resize(n);
  for(int i=0;i<n;i++)
  {
    mu_post[i]=exp( mu_post_log[i] + 0.5*Sigma_log(i,i));
  }
}

void
TravelTimeInversion::getLogMuFromMu(std::vector<double>  & mu_post_log,const std::vector<double>   & mu_post, const NRLib::Grid2D<double>& Sigma_log) const
{
  int n = static_cast<int>(mu_post.size());
  mu_post_log.resize(n);

  for(int i=0;i<n;i++)
  {
    mu_post_log[i] = log(mu_post[i])-Sigma_log(i,i)*0.5;
  }
}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::getSigmaFromLogSigma(std::vector<double>   & mu,const NRLib::Grid2D<double> &Sigma_log, NRLib::Grid2D<double> & Sigma) const
{
  int n = static_cast<int>(mu.size());
  Sigma.Resize(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      Sigma(i,j)=(exp(Sigma_log(i,j))-1.0)*mu[i]*mu[j];
  }
}

FFTGrid*
TravelTimeInversion::getCovFunkFromLogCovFunk(FFTGrid*  cov_log,double meanVpRelative) const
{
  bool transformedOninput =cov_log->getIsTransformed();
  if(transformedOninput)
   cov_log->fftInPlace();

  FFTGrid* cov =  new FFTGrid( cov_log);
  int nxp=cov->getNxp();
  int nyp=cov->getNyp();
  int nzp=cov->getNzp();

  for(int i=0;i<nxp;i++)
    for(int j=0;j<nyp;j++)
      for(int k=0;k<nzp;k++)
      {
        double v =   cov_log->getRealValue(i,j,k);
        double value=(exp(v)-1.0)* meanVpRelative*meanVpRelative;
        cov->setRealValue(i,j,k,static_cast<float>(value));
      }
  if(transformedOninput)
    cov_log->invFFTInPlace();

  return cov;
}

void
TravelTimeInversion::transformVpMinusToLogVp(const std::vector<double>   & mu_vp_minus,
                                             const NRLib::Grid2D<double> & Sigma_vp_minus,
                                             std::vector<double>         & mu_log_vp,
                                             NRLib::Grid2D<double>       & Sigma_log_vp) const
{
  int n_layers = static_cast<int>(mu_vp_minus.size());

  calculateCentralMomentLogNormalInverse(mu_vp_minus, Sigma_vp_minus, mu_log_vp, Sigma_log_vp);

  for (int i = 0; i < n_layers; i++)
    mu_log_vp[i] = mu_log_vp[i] / (-1);

  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      Sigma_log_vp(i, j) = Sigma_log_vp(i, j) / 1;
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                     const NRLib::Grid2D<double> & variance_log_vp,
                                                     std::vector<double>         & mu_vp_trans,
                                                     NRLib::Grid2D<double>       & variance_vp_trans) const
{
  int n = static_cast<int>(mu_log_vp.size());

  mu_vp_trans.resize(n);
  for (int i = 0; i < n; i++)
    mu_vp_trans[i] = std::exp(mu_log_vp[i] + 0.5 * variance_log_vp(i, i));

  variance_vp_trans.Resize(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      variance_vp_trans(i, j) = mu_vp_trans[i] * mu_vp_trans[j] * (std::exp(variance_log_vp(i, j)) - 1);
  }

}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                            const NRLib::Grid2D<double> & variance_vp_trans,
                                                            std::vector<double>         & mu_log_vp,
                                                            NRLib::Grid2D<double>       & variance_log_vp) const
{
  int n = static_cast<int>(mu_vp_trans.size());

  variance_log_vp.Resize(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      variance_log_vp(i, j) = std::log(1 + variance_vp_trans(i, j) / (mu_vp_trans[i] * mu_vp_trans[j]));
  }

  mu_log_vp.resize(n);
  for (int i = 0; i < n; i++)
    mu_log_vp[i] = std::log(mu_vp_trans[i]) - 0.5 * variance_log_vp(i, i);

}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::calculateSecondCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                           const NRLib::Grid2D<double> & variance_log_vp,
                                                           std::vector<double>         & mu_vp_square,
                                                           NRLib::Grid2D<double>       & variance_vp_square) const
{
  int n_layers = static_cast<int>(mu_log_vp.size());

  std::vector<double> mu(n_layers);
  for (int i = 0; i < n_layers; i++)
    mu[i] = mu_log_vp[i] * 2;

  NRLib::Grid2D<double> variance(n_layers, n_layers);
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      variance(i, j) = variance_log_vp(i, j) * 4;
  }

  mu_vp_square.resize(n_layers);
  variance_vp_square.Resize(n_layers, n_layers);

  calculateCentralMomentLogNormal(mu, variance, mu_vp_square, variance_vp_square);
}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::calculateMinusFirstCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                               const NRLib::Grid2D<double> & variance_log_vp,
                                                               std::vector<double>         & mu_vp_minus,
                                                               NRLib::Grid2D<double>       & variance_vp_minus) const
{
  // log(Vp) \sim N(\mu,Sigma)
  // Vp^(-1) \sim N(\mu*, Sigma*)

  int n_layers = static_cast<int>(mu_log_vp.size());

  std::vector<double> mu(n_layers);
  for (int i = 0; i < n_layers; i++)
    mu[i] = mu_log_vp[i] * (-1);

  NRLib::Grid2D<double> variance(n_layers, n_layers);
  for (int i = 0; i < n_layers; i++) {
    for (int j = 0; j < n_layers; j++)
      variance(i, j) = variance_log_vp(i, j) * 1;
  }

  mu_vp_minus.resize(n_layers);
  variance_vp_minus.Resize(n_layers, n_layers);

  calculateCentralMomentLogNormal(mu, variance, mu_vp_minus, variance_vp_minus);
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::addCovariance(const NRLib::Grid2D<double> & Sigma_post,
                                   std::vector<double>         & cov_stationary,
                                   int                           n_nopad) const
{

  int nzp = static_cast<int>(Sigma_post.GetNI());

  std::vector<double> circulant = makeCirculantCovariance(Sigma_post, n_nopad);

  for (int i = 0; i < nzp; i++)
    cov_stationary[i] += circulant[i];
}

void
TravelTimeInversion::addCovarianceMat(const NRLib::Grid2D<double> & Sigma_post,
                                       NRLib::Matrix         & cov_cum) const
{

  int nzp = static_cast<int>(Sigma_post.GetNI());

  for (int i = 0; i < nzp; i++)
    for (int j = 0; j < nzp; j++)
      cov_cum(i,j) += Sigma_post(i,j);
}


//-----------------------------------------------------------------------------------------//

std::vector<double>
TravelTimeInversion::makeCirculantCovariance(const NRLib::Matrix & cov,
                                             const int                   & n_nopad) const
{
  int n=cov.numCols();
  NRLib::Grid2D<double> covGrid(n,n , 0);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      covGrid(i,j)=cov(i,j);

  std::vector<double> covT=makeCirculantCovariance( covGrid, n_nopad);
  return covT;
}


std::vector<double>
TravelTimeInversion::makeCirculantCovariance(const NRLib::Grid2D<double> & cov,
                                             const int                   & n_nopad) const
{
  int n = static_cast<int>(cov.GetNI());

  NRLib::Grid2D<int> I(n, n, 0);
  for (int i = 0; i < n_nopad; i++) {
    for (int j = i; j < n; j++)
      I(i, j) = 1;
  }

  std::vector<double> covT(n, 0);
  std::vector<int>    count(n, 0);

  for (int i = 0; i < n; i++) {
    int place = 0;
    for (int j = i; j < n; j++) {
      covT[place]  += cov(i, j) * I(i, j);
      count[place] += I(i, j);
      place++;
    }
  }

  for (int i = 0; i < n; i++)
    covT[i] /= count[i];

  for (int i = 1; i <n / 2; i++) {
    covT[i] = (covT[i] + covT[n - i]) / 2;
    covT[n - i] = covT[i];
  }

  return covT;
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::setExpectation(int                          i_ind,
                                    int                          j_ind,
                                    const std::vector<double>  & post_vp,
                                    std::vector<KrigingData2D> & mu_log_vp_post) const
{
  int nzp = static_cast<int>(mu_log_vp_post.size());
  for (int i = 0; i < nzp; i++)
    mu_log_vp_post[i].addData(i_ind, j_ind, static_cast<float>(post_vp[i]));
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::krigeExpectation3D(const Simbox                * simbox,
                                        std::vector<KrigingData2D>  & kriging_post,
                                        const int                   & nxp,
                                        const int                   & nyp,
                                        FFTGrid                    *& mu_post) const
{
  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  const int    rnxp = 2 * (nx / 2 + 1);

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

  Vario * lateralVario = new GenExpVario(1, 500, 500);
  const CovGrid2D & covGrid2D = Kriging2D::makeCovGrid2D(simbox, lateralVario, 0);

  FFTGrid * bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nx, ny, nz, false); // Grid without padding
  bgGrid->createRealGrid();
  bgGrid->setType(FFTGrid::PARAMETER);

  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);

  bgGrid->setAccessMode(FFTGrid::WRITE);

  for (int k = 0 ; k < nz ; k++) {

    kriging_post[k].findMeanValues();

    int n_data = kriging_post[k].getNumberOfData();

    double mean_value = 0;
    for (int i = 0; i < n_data; i++)
      mean_value += static_cast<double>(kriging_post[k].getData(i) / n_data);

    surface.Assign(mean_value);

    // Kriging of layer
    Kriging2D::krigSurface(surface, kriging_post[k], covGrid2D);

    // Set layer in background model from surface
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < rnxp ; i++) {
        if (i < nx)
          bgGrid->setNextReal(float(surface(i, j)));
        else
          bgGrid->setNextReal(0);
      }
    }
  }

  bgGrid->endAccess();

  delete &covGrid2D;

  const int nzp  = static_cast<int>(kriging_post.size());

  mu_post = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);

  Background::createPaddedParameter(mu_post, bgGrid);

  delete bgGrid;

  delete lateralVario;
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::generateStationaryDistribution(const Simbox                * timeSimbox,
                                                    std::vector<KrigingData2D>  & kriging_post,
                                                    const std::vector<double>   & pri_circulant_cov,
                                                    NRLib::Matrix                 post_cov,
                                                    const int                   & n_traces,
                                                    const Surface  * errorCorrXY,
                                                    const float & corrGradI,
                                                    const float & corrGradJ,
                                                    FFTGrid                     * pri_mu,
                                                    FFTGrid                    *& stationary_observations,
                                                    FFTGrid                    *& stationary_covariance,
                                                    std::vector<int>            & observation_filter) const
{
  assert(!pri_mu->getIsTransformed());
   // Checked
  int nx   = pri_mu->getNx();
  int ny   = pri_mu->getNy();
  int nz   = pri_mu->getNz();
  int nxp  = pri_mu->getNxp();
  int nyp  = pri_mu->getNyp();
  int nzp  = pri_mu->getNzp();
  int cnzp = nzp/2 + 1;
  int rnzp = 2 * cnzp;

  NRLib::Matrix pri_cov(nzp,nzp,0);


  FFTGrid * post_mu = NULL;

  krigeExpectation3D(timeSimbox,
                     kriging_post,
                     nxp,
                     nyp,
                     post_mu);

  /* NBNB print out for debug
  pri_mu->writeAsciiRaw("pri_mu.dat");
  post_mu->writeAsciiRaw("interpolated_mu.dat");

  //  */
  //  pri_mu->calculateStatistics();
  //double scalePri=pri_mu->getAvgReal();
  //post_mu->calculateStatistics();
  //double scalePost=post_mu->getAvgReal();

  // adjust covariance by scale factor because the change of absolute level
  // will adjust the total uncertainty level (which is unfortunate)
  // The argument beeing that changes made to the absolute level
  // will influence variability in all frequencies, this is not wanted.
  // The result is used for the computation of properties of likelihood, which in turn is used further.
  //double scaleFactor = exp(2*(scalePri-scalePost));


  fftw_real    * pri_cov_r = new fftw_real[rnzp];
  fftw_complex * pri_cov_c = reinterpret_cast<fftw_complex*>(pri_cov_r);

  for (int i = 0; i < nzp; i++)
    pri_cov_r[i] = static_cast<float>(pri_circulant_cov[i]);//*scaleFactor; NBNB

  Utils::fft(pri_cov_r, pri_cov_c, nzp);

  for (int i = 0; i < cnzp; i++)
    pri_cov_c[i].im = 0;

  /*
  std::vector<double>  post_circulant_cov= makeCirculantCovariance(post_cov,nz);
  fftw_real    * post_cov_r = new fftw_real[rnzp];
  fftw_complex * post_cov_c = reinterpret_cast<fftw_complex*>(post_cov_r);

  for (int i = 0; i < nzp; i++)
    post_cov_r[i] = static_cast<float>(post_circulant_cov[i]);

  Utils::fft(post_cov_r, post_cov_c, nzp);

  for (int i = 0; i < cnzp; i++)
    post_cov_c[i].im = 0;
    */


  fftw_real    * var_e_r_lo = new fftw_real[rnzp];
  fftw_complex * var_e_c_lo = reinterpret_cast<fftw_complex*>(var_e_r_lo);
  fftw_real    * var_e_r_hi = new fftw_real[rnzp];
  fftw_complex * var_e_c_hi = reinterpret_cast<fftw_complex*>(var_e_r_hi);

  // Find those frequencies where there has been a relevant change
  // calculateFilter(pri_cov_c, post_cov_c, nzp, observation_filter);
  // calculate error for traces where the traveltime data are observed
  // calculateErrorVariance(pri_cov_c, post_cov_c,nzp, var_e_c);

  calculateFilterAndErrorVariance(pri_cov_c,post_cov,nz,var_e_c_lo,var_e_c_hi,observation_filter);

  // compute the observations corresponding to prior vp posterior for each observed profile
  stationary_observations = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  stationary_observations->createRealGrid();
  stationary_observations->setType(FFTGrid::DATA);

  calculateStationaryObservations(pri_cov_c,
                                  var_e_c_lo,
                                  observation_filter,
                                  pri_mu,
                                  post_mu,
                                  stationary_observations);

  // lib_matrDumpVecCpx("var_e_c.dat", var_e_c,cnzp);
  // Multiply var_e_c with  factor below to increase the variance as the RMS data are only observed in some of the traces
  // Note that the information in one profile is distributed to all neighbours, thus there is no "loss" of information
  // just redistribution/reinterpretation.
  float factor = static_cast<float>(nx * ny) / static_cast<float>(n_traces);
  for (int i = 0; i < cnzp; i++)
    var_e_c_hi[i].re *= factor;


  /* NBNB printout for debug
  NRLib::WriteMatrixToFile("postCovMat.dat",post_cov);
  lib_matrDumpVecCpx("pri_cov_c.dat", pri_cov_c, cnzp);
  lib_matrDumpVecCpx("var_e_c_lo.dat", var_e_c_lo,cnzp);
  lib_matrDumpVecCpx("var_e_c_hi.dat", var_e_c_hi,cnzp);
  // */

  stationary_covariance = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  stationary_covariance->createRealGrid();
  stationary_covariance->setType(FFTGrid::COVARIANCE);

  Utils::fftInv(var_e_c_hi, var_e_r_hi, nzp);

  stationary_covariance->fillInParamCorr( errorCorrXY,var_e_r_hi,corrGradI,corrGradJ);


    /* NBNB print out for debug
   stationary_covariance->writeAsciiRaw("StationaryErrorCov.dat");
   stationary_observations->writeAsciiRaw("StationaryObs.dat");
   // */

  delete [] pri_cov_r;
  delete [] var_e_r_lo;
  delete [] var_e_r_hi;

  delete post_mu;
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::calculateStationaryObservations(const fftw_complex  *  pri_cov_c,
                                                     const fftw_complex  *  var_e_c,
                                                     const std::vector<int>  filter_c,
                                                     FFTGrid             *  pri_mu,
                                                     FFTGrid             *  post_mu,
                                                     FFTGrid             *& stat_d) const
{
  // Calculate d = mu_pri + (var_pri+var_e)/conj_var_pri * (mu_post-mu_pri)
   // Checked 11.03 2014 by Odd Kolbjørnsen
  assert(!pri_mu->getIsTransformed());assert(!post_mu->getIsTransformed());
  int nxp = pri_mu->getNxp();
  int nyp = pri_mu->getNyp();
  int nzp = pri_mu->getNzp();
  int cnzp = nzp/2 + 1;

  fftw_complex * add_c        = new fftw_complex[cnzp];
  fftw_complex * conj_pri_cov = new fftw_complex[cnzp];
  fftw_complex * divide_c     = new fftw_complex[cnzp];

  addComplex(pri_cov_c, var_e_c, cnzp, add_c);
  complexConjugate(pri_cov_c, cnzp, conj_pri_cov);
  divideComplex(add_c, conj_pri_cov, cnzp, divide_c);

  delete [] conj_pri_cov;
  delete [] add_c;

  pri_mu ->setAccessMode(FFTGrid::RANDOMACCESS);
  post_mu->setAccessMode(FFTGrid::RANDOMACCESS);
  stat_d->setAccessMode(FFTGrid::RANDOMACCESS);
  fftw_complex* subtract=new fftw_complex[cnzp];
  fftw_complex* multiply=new fftw_complex[cnzp];
  fftw_complex * muPri_c  = new fftw_complex[cnzp];
  fftw_real *    muPri_r  = reinterpret_cast<fftw_real*>(muPri_c);
  fftw_complex * muPost_c = new fftw_complex[cnzp];
  fftw_real *    muPost_r = reinterpret_cast<fftw_real*>(muPost_c);
  fftw_complex * dObs_c   = new fftw_complex[cnzp];
  fftw_real *    dObs_r   = reinterpret_cast<fftw_real*>( dObs_c);

  for (int i = 0; i < nxp; i++) {
    for (int j = 0; j < nyp; j++) {
       for (int k = 0; k < nzp; k++) {
         muPri_r[k]  = pri_mu->getRealValue(i,j,k,true);
         muPost_r[k] = post_mu->getRealValue(i,j,k,true);
       }


       Utils::fft(muPri_r, muPri_c, nzp);
       Utils::fft(muPost_r, muPost_c, nzp);
       for (int k = 0; k < cnzp; k++) {
         muPri_c[k].re/=sqrt(static_cast<double>(nzp));
         muPri_c[k].im/=sqrt(static_cast<double>(nzp));
         muPost_c[k].re/=sqrt(static_cast<double>(nzp));
         muPost_c[k].im/=sqrt(static_cast<double>(nzp));
       }



       subtractComplex(muPost_c, muPri_c, cnzp, subtract);
       multiplyComplex(divide_c, subtract, cnzp, multiply);
       addComplex(muPri_c, multiply, cnzp,dObs_c);
       Utils::fftInv(dObs_c, dObs_r, nzp);
       for (int k = 0; k < nzp; k++) {
         stat_d->setRealValue(i,j,k,static_cast<fftw_real>(dObs_r[k]*sqrt(static_cast<double>(nzp))),true);
      }
    }
  }

  pri_mu ->endAccess();
  post_mu->endAccess();
  stat_d ->endAccess();

  delete [] divide_c;
  delete [] subtract;
  delete [] multiply;
  delete [] muPri_c;
  delete [] muPost_c;
  delete [] dObs_c;
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::calculateErrorVariance(const fftw_complex * pri_cov_c,
                                            const fftw_complex * post_cov_c,
                                            const int          & nzp,
                                            fftw_complex       * var_e_c) const
{
  // Calculate var_e = var_pri *(conj_var_pri - var_pri + var_post) / (var_pri-var_post)
   // Checked 28.November 2013 by Odd Kolbjørnsen
  int cnzp = nzp / 2 + 1;

  fftw_complex * conj_pri_cov  = new fftw_complex[cnzp];
  fftw_complex * diff_prior    = new fftw_complex[cnzp];
  fftw_complex * nominator_sum = new fftw_complex[cnzp];
  fftw_complex * mult_c        = new fftw_complex[cnzp];
  fftw_complex * subtract_c    = new fftw_complex[cnzp];

  complexConjugate(pri_cov_c, cnzp, conj_pri_cov);
  subtractComplex(conj_pri_cov, pri_cov_c, cnzp, diff_prior);
  addComplex(diff_prior, post_cov_c, cnzp, nominator_sum);
  multiplyComplex(pri_cov_c, nominator_sum, cnzp, mult_c);
  subtractComplex(pri_cov_c, post_cov_c, cnzp, subtract_c);
  divideComplex(mult_c, subtract_c, cnzp, var_e_c);

  if(var_e_c[0].re<0)
     LogKit::LogFormatted(LogKit::Low,"\n Warning: TraveltimeInversion has no or little  effect :-( \n");

  for(int i=1; i<cnzp;i++){
    var_e_c[i].im  = 0;
    if(var_e_c[i].re < 0)
      var_e_c[i].re = 0.0; //Keep non negative even though it is not well defined
  }                        //Note: std::vector filter from calculateFilter() is used to avoid problems

  delete [] conj_pri_cov;
  delete [] diff_prior;
  delete [] nominator_sum;
  delete [] mult_c;
  delete [] subtract_c;
}

//-------------------------------------------------------------------------------

void
TravelTimeInversion::calculateFilter(const fftw_complex * pri_cov_c,
                                     const fftw_complex * post_cov_c,
                                     const int          & nzp,
                                     std::vector<int>   & filter) const
{
   // Fudge factor:  factor
   // Rules out cases where the error standard deviation is very large
   // factor=1.04 =>5 times larger than the standard deviation in the parameter, 5 = sqrt(1/0.04)
   // factor=1.01 =>10 times larger than the standard deviation in the parameter, 10 = sqrt(1/0.01)
   // factor=1.0001 =>100 times larger than the standard deviation in the parameter, 100 = sqrt(1/0.0001)

  // Checked 22.Januar 2014 by Odd Kolbjørnsen
  double factor =1.001; // 31.62 times larger
  int cnzp = nzp/2 + 1;

  std::vector<double> abs_pri;
  std::vector<double> abs_post;

  absoulteComplex(pri_cov_c, cnzp, abs_pri);
  absoulteComplex(post_cov_c, cnzp, abs_post);

  filter.resize(cnzp, 0);

  for (int i = 0; i < cnzp; i++) {
    if (abs_pri[i] > abs_post[i] *  factor)
      filter[i] = 1;
  }

 // for (int i = cnzp; i < nzp; i++)
 //   filter[i] = filter[nzp-i];

}

void
TravelTimeInversion::calculateFilterAndErrorVariance(const fftw_complex * pri_cov_c,
                                                                const NRLib::Matrix post_cov,
                                                                const int       &  nz,
                                                                fftw_complex       * var_e_c_lo,
                                                                fftw_complex       * var_e_c_hi,
                                                                std::vector<int>   & filter) const
{  // Fudge factor:  factor
   // Rules out cases where the error standard deviation is very large
   // factor=1.04 =>5 times larger than the standard deviation in the parameter, 5 = sqrt(1/0.04)
   // factor=1.01 =>10 times larger than the standard deviation in the parameter, 10 = sqrt(1/0.01)
   // factor=1.0001 =>100 times larger than the standard deviation in the parameter, 100 = sqrt(1/0.0001)
  double factor =1.01; // 10 times larger

  int nzp=post_cov.numCols();
  int cnzp = nzp/2+1;
  filter.resize(cnzp, 0);
  fftw_complex *  pri_cov_cc= new fftw_complex[cnzp];
  fftw_real * pri_cov_r = reinterpret_cast<fftw_real*>(pri_cov_cc);

  for(int i=0;i<cnzp;i++)
  {
    pri_cov_cc[i].re=pri_cov_c[i].re;
    pri_cov_cc[i].im=pri_cov_c[i].im;
  }

  Utils::fftInv(pri_cov_cc,pri_cov_r,nzp);

  NRLib::Matrix pri_cov(nzp,nzp,0);
  for(int i=0;i<nzp;i++)
    for(int j=0;j<nzp;j++)
      pri_cov(i,j)=pri_cov_r[std::min(std::abs(j-i),std::abs(i-j))];

 // NRLib::WriteMatrixToFile("priCovMat.dat",pri_cov);

  // constant  term
  double priVarRe=0.0;
  double postVarRe=0.0;

  for(int i=0;i<nz;i++) // Note: We use nz This is the interesting part ...
    for(int j=0;j<nz;j++){
      priVarRe  += pri_cov(i,j);
      postVarRe += post_cov(i,j);
    }

  bool adj=false;
  double redRe =  postVarRe/ priVarRe;
  if(redRe < 0.0001) {// robust for numerical errors providing  negative variances.
    redRe = 0.0001;
    adj=true;
  }


  if(redRe*factor < 1.0 ){
    var_e_c_lo[0].re = static_cast<fftw_real>(pri_cov_c[0].re/(1.0/redRe-1.0));
    var_e_c_lo[0].im = 0.0;
    var_e_c_hi[0].re = static_cast<fftw_real>(pri_cov_c[0].re/(1.0/redRe-1.0));
    var_e_c_hi[0].im = 0.0;
    filter[0]     = 1;
  }
  else
  {
     LogKit::LogFormatted(LogKit::Low,"\n Warning: TraveltimeInversion has no or little  effect. \nCheck:\n  1)Error variance in inversion of traveltime data \n  2) Low cut for inversion should be -1.0 when inverting traveltime and RMS \n");
     var_e_c_lo[0].re = 0.0;
     var_e_c_lo[0].im = 0.0;
     var_e_c_hi[0].re = 0.0;
     var_e_c_hi[0].im = 0.0;
     filter[0]     = 0;
  }

  double priVarIm=0.0;
  double postVarIm=0.0;
  double redIm;
  bool haveResolution  = true;


  for(int k=1;k<cnzp;k++)
  {
    if(haveResolution || k < 20){ // checks at least 20 smallest components this is enough in most cases
      priVarRe  = 0.0;
      postVarRe = 0.0;
      priVarIm  = 0.0;
      postVarIm = 0.0;

      for(int i=0;i<nz;i++){ // Note: We use nz This is the interesting part of signal...
        for(int j=0;j<nz;j++){
          double cosi=std::cos(2.0*NRLib::Pi*k*static_cast<double>(i)/static_cast<double>(nzp));
          double cosj=std::cos(2.0*NRLib::Pi*k*static_cast<double>(j)/static_cast<double>(nzp));
          double sini=std::sin(2.0*NRLib::Pi*k*static_cast<double>(i)/static_cast<double>(nzp));
          double sinj=std::sin(2.0*NRLib::Pi*k*static_cast<double>(j)/static_cast<double>(nzp));
          priVarRe  += pri_cov(i,j) * cosi*cosj;
          postVarRe += post_cov(i,j)* cosi*cosj;;
          priVarIm  += pri_cov(i,j) * sini*sinj;;
          postVarIm += post_cov(i,j)* sini*sinj;;
        }
      }

      redRe =  postVarRe/ priVarRe;
      redIm =  postVarIm/ priVarIm;

      if(redRe < 0.0001){ // robust for numerical errors providing  negative/very-small  posterior variances.
        redRe = 0.0001;
        adj=true;
      }
      if(redIm < 0.0001){ // robust for numerical errors providing  negative/very-small posterior variances. 1/10000;
        redIm = 0.0001;
        adj=true;
      }

      if(redRe*factor < 1.0 && redIm*factor < 1.0 )
      {
        var_e_c_lo[k].re = static_cast<fftw_real>(pri_cov_c[k].re/(1.0/std::min(redRe,redIm)-1.0));
        var_e_c_hi[k].re = static_cast<fftw_real>(pri_cov_c[k].re/(1.0/std::max(redRe,redIm)-1.0));
        var_e_c_lo[k].im = 0.0;
        var_e_c_hi[k].im = 0.0;
        filter[k]        = 1;
        haveResolution   = true;
      }
      else{
        var_e_c_lo[k].re = 0.0f;
        var_e_c_lo[k].im = 0.0f;
        var_e_c_hi[k].re = 0.0f;
        var_e_c_hi[k].im = 0.0f;
        filter[k]        = 0;
        haveResolution   = false;
      }
    }
    else
    {
      var_e_c_lo[k].re = 0.0f;
      var_e_c_lo[k].im = 0.0f;
      var_e_c_hi[k].re = 0.0f;
      var_e_c_hi[k].im = 0.0f;
      filter[k]        = 0;
    }
  }

   if(adj){
     LogKit::LogFormatted(LogKit::Low,"\n Warning: TraveltimeInversion has very small error, check scale of error standard deviation\n");
   }

   delete [] pri_cov_cc;
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::multiplyComplex(const fftw_complex * z1,
                                     const fftw_complex * z2,
                                     const int          & n,
                                     fftw_complex       * z) const
{
  for (int i = 0; i < n; i++) {
    z[i].re = z1[i].re * z2[i].re - z1[i].im * z2[i].im;
    z[i].im = z1[i].re * z2[i].im + z1[i].im * z2[i].re;
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::divideComplex(const fftw_complex * z1,
                                   const fftw_complex * z2,
                                   const int          & n,
                                   fftw_complex       * z) const
{
  for (int i = 0; i < n; i++) {
    z[i].re = (z1[i].re * z2[i].re + z1[i].im * z2[i].im) / (std::pow(z2[i].re, 2) + std::pow(z2[i].im, 2));
    z[i].im = (z2[i].re * z1[i].im - z1[i].re * z2[i].im) / (std::pow(z2[i].re, 2) + std::pow(z2[i].im, 2));
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::addComplex(const fftw_complex * z1,
                                const fftw_complex * z2,
                                const int          & n,
                                fftw_complex       * z) const
{
  for (int i = 0; i < n; i++) {
    z[i].re = z1[i].re + z2[i].re;
    z[i].im = z1[i].im + z2[i].im;
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::subtractComplex(const fftw_complex * z1,
                                     const fftw_complex * z2,
                                     const int          & n,
                                     fftw_complex       * z) const
{
  for (int i = 0; i < n; i++) {
    z[i].re = z1[i].re - z2[i].re;
    z[i].im = z1[i].im - z2[i].im;
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::absoulteComplex(const fftw_complex  * z,
                                     const int           & n,
                                     std::vector<double> & abs_z) const
{
  abs_z.resize(n);

  for (int i = 0; i < n; i++)
    abs_z[i] = std::sqrt(std::pow(z[i].re, 2) + std::pow(z[i].im, 2));
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::complexConjugate(const fftw_complex * z,
                                      const int          & n,
                                      fftw_complex       * conj_z) const
{

  for (int i = 0; i < n; i++) {
    conj_z[i].re =   z[i].re;
    conj_z[i].im = - z[i].im;
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateFullPosteriorModel(const std::vector<int>  & observation_filter,
                                                 SeismicParametersHolder & seismic_parameters,
                                                 FFTGrid                 * stationary_observations,
                                                 FFTGrid                 * stationary_observation_covariance) const
{
  FFTGrid * mu_vp  = seismic_parameters.GetMuAlpha();
  FFTGrid * mu_vs  = seismic_parameters.GetMuBeta();
  FFTGrid * mu_rho = seismic_parameters.GetMuRho();

  FFTGrid * cov_vp  = seismic_parameters.GetCovAlpha();
  FFTGrid * cov_vs  = seismic_parameters.GetCovBeta();
  FFTGrid * cov_rho = seismic_parameters.GetCovRho();

  FFTGrid * cr_cov_vp_vs  = seismic_parameters.GetCrCovAlphaBeta();
  FFTGrid * cr_cov_vp_rho = seismic_parameters.GetCrCovAlphaRho();
  FFTGrid * cr_cov_vs_rho = seismic_parameters.GetCrCovBetaRho();

  /*  print out for debug

  mu_vp->invFFTInPlace();
  cov_vp->invFFTInPlace();
  mu_vp->writeAsciiRaw("prior_meanVp.dat");
  cov_vp->writeAsciiRaw("prior_covVp.dat");
  mu_vp->fftInPlace();
  cov_vp->fftInPlace();
  // */

  mu_vp ->setAccessMode(FFTGrid::READANDWRITE);
  mu_vs ->setAccessMode(FFTGrid::READANDWRITE);
  mu_rho->setAccessMode(FFTGrid::READANDWRITE);

  cov_vp ->setAccessMode(FFTGrid::READANDWRITE);
  cov_vs ->setAccessMode(FFTGrid::READANDWRITE);
  cov_rho->setAccessMode(FFTGrid::READANDWRITE);

  cr_cov_vp_vs ->setAccessMode(FFTGrid::READANDWRITE);
  cr_cov_vp_rho->setAccessMode(FFTGrid::READANDWRITE);
  cr_cov_vs_rho->setAccessMode(FFTGrid::READANDWRITE);

  if (stationary_observations->getIsTransformed() == false)
    stationary_observations->fftInPlace();

  if (stationary_observation_covariance->getIsTransformed() == false)
    stationary_observation_covariance->fftInPlace();

  int cnxp = mu_vp->getCNxp();
  int nyp  = mu_vp->getNyp();
  int nzp  = mu_vp->getNzp();

  fftw_complex * mu_m    = new fftw_complex[3];
  fftw_complex * mu_post = new fftw_complex[3];

  fftw_complex data;
  fftw_complex data_variance;
  fftw_complex var_d;
  fftw_complex diff;
  fftw_complex divide;
  fftw_complex mu_help;
  fftw_complex multiply;
  fftw_complex Sigma_help;

  fftw_complex ** Sigma_m = new fftw_complex*[3];
  for (int i = 0; i < 3; i++)
    Sigma_m[i] = new fftw_complex[3];

  fftw_complex ** Sigma_post = new fftw_complex*[3];
  for (int i = 0; i < 3; i++)
    Sigma_post[i] = new fftw_complex[3];

  std::string text = "\nBuilding posterior RMS distribution:";
  LogKit::LogFormatted(LogKit::Low, text);

  float monitorSize = std::max(1.0f, static_cast<float>(nzp) * 0.02f);
  float nextMonitor = monitorSize;
  std::cout
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^";

  for (int k = 0; k < nzp; k++) {
    int indK = std::min(k,nzp-k);
    int filter = observation_filter[indK];

    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {

        data          = stationary_observations          ->getNextComplex();
        data_variance = stationary_observation_covariance->getNextComplex();

        mu_m[0] = mu_vp ->getNextComplex();
        mu_m[1] = mu_vs ->getNextComplex();
        mu_m[2] = mu_rho->getNextComplex();

        seismic_parameters.getNextParameterCovariance2(Sigma_m); // NBNB OK disturbes test suite use line under to check if OK
         //seismic_parameters.getNextParameterCovariance(Sigma_m);  seismic_parameters.getNextParameterCovariance2 should be used allways

        if (filter == 0) {
          // Prior = posterior
          mu_vp ->setNextComplex(mu_m[0]);
          mu_vs ->setNextComplex(mu_m[1]);
          mu_rho->setNextComplex(mu_m[2]);

          cov_vp ->setNextComplex(Sigma_m[0][0]);
          cov_vs ->setNextComplex(Sigma_m[1][1]);
          cov_rho->setNextComplex(Sigma_m[2][2]);

          cr_cov_vp_vs ->setNextComplex(Sigma_m[0][1]);
          cr_cov_vp_rho->setNextComplex(Sigma_m[0][2]);
          cr_cov_vs_rho->setNextComplex(Sigma_m[1][2]);
        }
        else {
          addComplex(&Sigma_m[0][0], &data_variance, 1, &var_d);
          subtractComplex(&data, &mu_m[0], 1, &diff);
          divideComplex(&diff, &var_d, 1, &divide);

          for (int ii = 0; ii < 3; ii++) {
            multiplyComplex(&Sigma_m[ii][0], &divide, 1, &mu_help);
            addComplex(&mu_m[ii], &mu_help, 1, &mu_post[ii]);
          }

          for (int ii = 0; ii < 3; ii++) {
            for (int jj = 0; jj < 3; jj++) {
              multiplyComplex(&Sigma_m[ii][0], &Sigma_m[0][jj], 1, &multiply);
              divideComplex(&multiply, &var_d, 1, &Sigma_help);
              subtractComplex(&Sigma_m[ii][jj], &Sigma_help, 1, &Sigma_post[ii][jj]);
            }
          }

          mu_vp ->setNextComplex(mu_post[0]);
          mu_vs ->setNextComplex(mu_post[1]);
          mu_rho->setNextComplex(mu_post[2]);

          cov_vp ->setNextComplex(Sigma_post[0][0]);
          cov_vs ->setNextComplex(Sigma_post[1][1]);
          cov_rho->setNextComplex(Sigma_post[2][2]);

          cr_cov_vp_vs ->setNextComplex(Sigma_post[0][1]);
          cr_cov_vp_rho->setNextComplex(Sigma_post[0][2]);
          cr_cov_vs_rho->setNextComplex(Sigma_post[1][2]);

        }
      }
    }
    if (k + 1 >= static_cast<int>(nextMonitor)) {
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }

  mu_vp ->endAccess();
  mu_vs ->endAccess();
  mu_rho->endAccess();

  cov_vp ->endAccess();
  cov_vs ->endAccess();
  cov_rho->endAccess();

  cr_cov_vp_vs ->endAccess();
  cr_cov_vp_rho->endAccess();
  cr_cov_vs_rho->endAccess();

  delete [] mu_m;
  delete [] mu_post;

  for (int i = 0; i < 3; ++i)
    delete [] Sigma_m[i];
  delete [] Sigma_m;

  for (int i = 0; i < 3; ++i)
    delete [] Sigma_post[i];
  delete [] Sigma_post;
  /*  print out for debug
  mu_vp->invFFTInPlace();
  cov_vp->invFFTInPlace();
  mu_vp->writeAsciiRaw("post_meanVp.dat");
  cov_vp->writeAsciiRaw("post_covVp.dat");
  mu_vp->fftInPlace();
  cov_vp->fftInPlace();
  // */
}
//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::transformCovarianceToRelativeScale( std::vector<double> mu, NRLib::Grid2D<double> & Sigma, int n1, int n1p,int n2, int n2p,int n3, int n3p) const
{
  int ntot=n1p+n2p+n3p;
  NRLib::Vector localCoefficientOfVariation(ntot);
  NRLib::Vector multiplyer(ntot);
  double meanCoefficientOfVariation1;
  double meanCoefficientOfVariation2;
  double meanCoefficientOfVariation3;

  for(int i=0;i<n1p;i++)
    localCoefficientOfVariation(i)=sqrt(Sigma(i,i))/mu[i];

   meanCoefficientOfVariation1=0.0;
   for(int i=0;i<n1;i++) // Here  we use nz  to find the reduction in the reservoir
     meanCoefficientOfVariation1 += localCoefficientOfVariation(i)/n1;

   for(int i=n1p;i<n1p+n2p;i++)
    localCoefficientOfVariation(i)=sqrt(Sigma(i,i))/mu[i];

   meanCoefficientOfVariation2=0.0;
   for(int i=n1p;i<n1p+n2;i++) // Here  we use nz  to find the reduction in the reservoir
     meanCoefficientOfVariation2 += localCoefficientOfVariation(i)/n2;

   for(int i=n1p+n2p;i<ntot;i++)
    localCoefficientOfVariation(i)=sqrt(Sigma(i,i))/mu[i];

   meanCoefficientOfVariation3=0.0;
   for(int i=n1p+n2p;i<n1p+n2p+n3;i++) // Here  we use nz  to find the reduction in the reservoir
     meanCoefficientOfVariation3 += localCoefficientOfVariation(i)/n3;

   for(int i=0;i<n1p;i++)
     multiplyer(i)=meanCoefficientOfVariation1/localCoefficientOfVariation(i);
   for(int i=n1p;i<n1p+n2p;i++)
     multiplyer(i)=meanCoefficientOfVariation2/localCoefficientOfVariation(i);
   for(int i=n1p+n2p;i<ntot;i++)
     multiplyer(i)=meanCoefficientOfVariation3/localCoefficientOfVariation(i);


  for(int i=0;i<ntot;i++)
    for(int j=0;j<ntot;j++)
    {
      Sigma(i,j) = Sigma(i,j)*multiplyer(i)*multiplyer(j);
    }

}


//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::generateMuSigmaLogVpAbove(const int                    & nz,
                                               const int                    & nzp,
                                               const double                 & mu_vp_top,
                                               const NRLib::Grid2D<double>  & Sigma_vp_above,
                                               const Surface                * parCorrXY,
                                               const float                  & corrGradI,
                                               const float                  & corrGradJ,
                                               FFTGrid                      * mu_log_vp_grid,
                                               FFTGrid                     *& mu_log_vp_above,
                                               FFTGrid                     *& Sigma_log_vp_above) const
{

  const int nx  = mu_log_vp_grid->getNx();
  const int ny  = mu_log_vp_grid->getNy();
  const int nxp = mu_log_vp_grid->getNxp();
  const int nyp = mu_log_vp_grid->getNyp();

  const int    rnxp = 2 * (nxp / 2 + 1);

  FFTGrid * bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nx, ny, nz, false); // Grid without padding
  bgGrid->createRealGrid();
  bgGrid->setType(FFTGrid::PARAMETER);

  mu_log_vp_grid->setAccessMode(FFTGrid::READ);
  bgGrid        ->setAccessMode(FFTGrid::WRITE);

  for (int i = 0; i < rnxp ; i++) {
    for (int j = 0; j < ny; j++) {

      double base_value = mu_log_vp_grid->getRealValue(i, j, 0);

      std::vector<double> mu_vp_above = generateMuVp(mu_vp_top, std::exp(base_value), nz);
      for (int k = 0; k < nz ; k++) {
        if (i < nx)
          bgGrid->setRealValue(i, j, k, static_cast<float>(mu_vp_above[k]));
        else
          bgGrid->setRealValue(i, j, k, 0);
      }
    }
  }

  mu_log_vp_grid->endAccess();
  bgGrid        ->endAccess();


  mu_log_vp_above = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);

  Background::createPaddedParameter(mu_log_vp_above, bgGrid);

  delete bgGrid;

  mu_log_vp_above->setAccessMode(FFTGrid::READANDWRITE);

  std::vector<double> cov_circulant_above(nzp, 0);

  for (int i = 0; i < rnxp ; i++) {
    for (int j = 0; j < nyp; j++) {

      std::vector<double> mu_vp_above(nzp);

      for (int k = 0; k < nzp; k++)
      {
        double a=static_cast<double>(mu_log_vp_above->getRealValue(i, j, k, true));
        if(a<=0) // avoid problems with log of negative values
          a=1.0;
        mu_vp_above[k] =a ;
      }

      std::vector<double>   mu_log_vp_above_trace;
      NRLib::Grid2D<double> Sigma_log_vp_above;
      calculateCentralMomentLogNormalInverse(mu_vp_above, Sigma_vp_above, mu_log_vp_above_trace, Sigma_log_vp_above);

      for (int k = 0; k < nzp; k++)
          mu_log_vp_above->setRealValue(i, j, k, static_cast<float>(mu_log_vp_above_trace[k]), true);

      if (i < nx && j < ny)
        addCovariance(Sigma_log_vp_above, cov_circulant_above, nz);

    }
  }

  mu_log_vp_above->endAccess();

  for (int i = 0; i < nzp; i++)
    cov_circulant_above[i] /= (nx * ny);

  int cnzp = nzp/2 + 1;
  int rnzp = 2 * cnzp;

  fftw_real * cov_circ_r = new fftw_real[rnzp];
  for (int i = 0; i < nzp; i++)
    cov_circ_r[i] = static_cast<float>(cov_circulant_above[i]);

  Sigma_log_vp_above = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  Sigma_log_vp_above->createRealGrid();
  Sigma_log_vp_above->setType(FFTGrid::COVARIANCE);

  Sigma_log_vp_above->fillInParamCorr(parCorrXY, cov_circ_r, corrGradI, corrGradJ);

  delete [] cov_circ_r;

}
//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateExpectation(const std::vector<int>  & observation_filter,
                                               FFTGrid                 * prior_mu,
                                               FFTGrid                 * prior_cov,
                                               FFTGrid                 * stationary_observations,
                                               FFTGrid                 * stationary_observation_covariance,
                                               FFTGrid                *& post_mu) const
{
  int nx   = prior_mu->getNx();
  int ny   = prior_mu->getNy();
  int nz   = prior_mu->getNz();
  int nxp  = prior_mu->getNxp();
  int cnxp = prior_mu->getCNxp();
  int nyp  = prior_mu->getNyp();
  int nzp  = prior_mu->getNzp();

  fftw_complex mu_m;
  fftw_complex Sigma_m;
  fftw_complex mu_post;

  fftw_complex data;
  fftw_complex data_variance;
  fftw_complex var_d;
  fftw_complex diff;
  fftw_complex divide;
  fftw_complex mu_help;

  post_mu = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  post_mu->createComplexGrid();
  post_mu->setType(FFTGrid::PARAMETER);

  prior_mu                         ->setAccessMode(FFTGrid::READ);
  prior_cov                        ->setAccessMode(FFTGrid::READ);
  stationary_observations          ->setAccessMode(FFTGrid::READ);
  stationary_observation_covariance->setAccessMode(FFTGrid::READ);
  post_mu                          ->setAccessMode(FFTGrid::WRITE);

  for (int k = 0; k < nzp; k++) {
    int indK   = std::min(k,nzp-k);
    int filter = observation_filter[indK];

    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {

        data          = stationary_observations          ->getNextComplex();
        data_variance = stationary_observation_covariance->getNextComplex();
        mu_m          = prior_mu                            ->getNextComplex();
        Sigma_m       = prior_cov                           ->getNextComplex();

        if (filter == 0)
          post_mu  ->setNextComplex(mu_m);
        else {
          addComplex(&Sigma_m, &data_variance, 1, &var_d);
          subtractComplex(&data, &mu_m, 1, &diff);
          divideComplex(&diff, &var_d, 1, &divide);

          multiplyComplex(&Sigma_m, &divide, 1, &mu_help);
          addComplex(&mu_m, &mu_help, 1, &mu_post);

          post_mu->setNextComplex(mu_post);

        }
      }
    }
  }

  prior_mu                         ->endAccess();
  prior_cov                        ->endAccess();
  stationary_observations          ->endAccess();
  stationary_observation_covariance->endAccess();
  post_mu                          ->endAccess();

}
//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateLogVpCovariance(const std::vector<int>  & observation_filter,
                                              FFTGrid                 * cov_vp,
                                              FFTGrid                 * stationary_observation_covariance,
                                              FFTGrid                *& post_cov_vp) const
{
  int nx   = cov_vp->getNx();
  int ny   = cov_vp->getNy();
  int nz   = cov_vp->getNz();
  int nxp  = cov_vp->getNxp();
  int cnxp = cov_vp->getCNxp();
  int nyp  = cov_vp->getNyp();
  int nzp  = cov_vp->getNzp();


  fftw_complex Sigma_m;


  fftw_complex data_variance;
  fftw_complex var_d;
  fftw_complex multiply;
  fftw_complex Sigma_help;
  fftw_complex Sigma_post;

  post_cov_vp = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  post_cov_vp->createComplexGrid();
  post_cov_vp->setType(FFTGrid::COVARIANCE);


  cov_vp                           ->setAccessMode(FFTGrid::READ);
  stationary_observation_covariance->setAccessMode(FFTGrid::READ);
  post_cov_vp                      ->setAccessMode(FFTGrid::WRITE);

  for (int k = 0; k < nzp; k++) {

    int indK   = std::min(k,nzp-k);
    int filter = observation_filter[indK];

    for (int j = 0; j < nyp; j++) {
      for (int i = 0; i < cnxp; i++) {
        data_variance = stationary_observation_covariance->getNextComplex();
        Sigma_m       = cov_vp                           ->getNextComplex();

        if (filter == 0)
          post_cov_vp  ->setNextComplex(Sigma_m);
        else {
          addComplex(&Sigma_m, &data_variance, 1, &var_d);
          multiplyComplex(&Sigma_m, &Sigma_m, 1, &multiply);
          divideComplex(&multiply, &var_d, 1, &Sigma_help);
          subtractComplex(&Sigma_m, &Sigma_help, 1, &Sigma_post);

          post_cov_vp->setNextComplex(Sigma_post);

        }
      }
    }
  }

  cov_vp                           ->endAccess();
  stationary_observation_covariance->endAccess();
  post_cov_vp                      ->endAccess();

}
//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateDistanceGrid(const Simbox              * simbox2,   // timesimbox of current timestep.
                                           const NRLib::Grid<double> & divided_grid21, // divided_grid contains v2/v1
                                           NRLib::Grid<double>       & distance) const  // returned value
{ // computes the time increment in the previous timesim box (time=1) which correspond to regular sampling of the current (time=2) .

  int nx  = static_cast<int>(divided_grid21.GetNI());
  int ny  = static_cast<int>(divided_grid21.GetNJ());
  int nz  = static_cast<int>(divided_grid21.GetNK());

  distance.Resize(nx, ny, nz);

  // The distance is calculated from the relation d1*v1 = d2*v2,
  // v1 is velocity of previous= (current-1)  timestep
  // v2 is velocity of current timestep
  // divided_grid contains v2/v1, and d2 is given in the simbox
  double d1;   // time from top to base of a cell in the previous = (current-1)  time grid
  double d2;   // time from top to base of a cell in the current time grid
  double v2v1; // ratio between current and previous velocity. (Not initial)


  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      d2 = simbox2    ->getdz(i, j);

      for (int k = 0; k < nz; k++) {
        v2v1 = divided_grid21(i, j, k); // v2/v1
        d1 = v2v1 * d2;
        distance(i, j, k) = d1;
      }
    }
  }
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateEVpGrid(FFTGrid  * mu_log_vp,
                                      FFTGrid  * cov_log_vp,
                                      FFTGrid *& mu_vp) const
{
  // Calculate E[vp] from log(vp)

  int nx  = mu_log_vp->getNx();
  int ny  = mu_log_vp->getNy();
  int nz  = mu_log_vp->getNz();
  int nxp = mu_log_vp->getNxp();
  int nyp = mu_log_vp->getNyp();
  int nzp = mu_log_vp->getNzp();

  std::vector<double>   cov_log_vp_profile = getCovLogVp(cov_log_vp);
  NRLib::Grid2D<double> Sigma_log_vp       = generateSigmaModel(cov_log_vp_profile);

  mu_vp = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, false);
  mu_vp->createRealGrid();
  mu_vp->setType(FFTGrid::PARAMETER);

  mu_log_vp->setAccessMode(FFTGrid::READ);
  mu_vp    ->setAccessMode(FFTGrid::WRITE);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      std::vector<double>   mu_log_vp_profile(nz);
      std::vector<double>   mu_vp_profile;
      NRLib::Grid2D<double> Sigma_mu_vp;

      for (int k = 0; k < nz; k++)
        mu_log_vp_profile[k]  = mu_log_vp->getRealValue(i, j, k);

      calculateCentralMomentLogNormal(mu_log_vp_profile,
                                      Sigma_log_vp,
                                      mu_vp_profile,
                                      Sigma_mu_vp);

      for (int k = 0; k < nz; k++)
        mu_vp->setRealValue(i, j, k, static_cast<float>(mu_vp_profile[k]));


    }
  }

  mu_log_vp->endAccess();
  mu_vp    ->endAccess();

}

//-----------------------------------------------------------------------------------------//

NRLib::Grid<double>
TravelTimeInversion::calculateDividedGridRMS(FFTGrid * pri_vp,
                                             FFTGrid * post_vp) const
{

  int nx  = pri_vp->getNx();
  int ny  = pri_vp->getNy();
  int nz  = pri_vp->getNz();

  NRLib::Grid<double> divided_grid(nx, ny, nz);

  pri_vp ->setAccessMode(FFTGrid::READ);
  post_vp->setAccessMode(FFTGrid::READ);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++)
        divided_grid(i, j, k) = static_cast<double>(pri_vp->getRealValue(i, j, k) / post_vp->getRealValue(i, j, k));
     }
  }

  pri_vp ->endAccess();
  post_vp->endAccess();

  return divided_grid;
}

//-----------------------------------------------------------------------------------------//
NRLib::Grid<double>
TravelTimeInversion::calculateRelativeVelocityUpdate(FFTGrid *relativeVelocityNew,
                                                     FFTGrid *relativeVelocityOld) const
{
  // divided_grid = vp_current/vp_previous
  int nx  = relativeVelocityNew->getNx();
  int ny  = relativeVelocityNew->getNy();
  int nz  = relativeVelocityNew->getNz();

  NRLib::Grid<double> divided_grid(nx, ny, nz);

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
       for (int k = 0; k < nz; k++)
         divided_grid(i, j, k) = relativeVelocityNew->getRealValue(i,j,k)/relativeVelocityOld->getRealValue(i,j,k);

  return  divided_grid;
}


NRLib::Grid<double>
TravelTimeInversion::calculateDividedGridHorizon(FFTGrid * post_mu_vp,
                                                 FFTGrid * post_cov_mu_vp) const
{
  // In the horizon posterior, E[ln(Vp1/Vp0)|d] and Cov(ln(Vp1/Vp0)|d) have been calculated
  // Of interest in the divided grid is E[(Vp0/Vp1)|d], which can be calculated from post_mu_vp and post_cov_mu_vp

  std::vector<double>   cov_log_vp   = getCovLogVp(post_cov_mu_vp);
  NRLib::Grid2D<double> Sigma_log_vp = generateSigmaModel(cov_log_vp);

  int nx  = post_mu_vp->getNx();
  int ny  = post_mu_vp->getNy();
  int nz  = post_mu_vp->getNz();

  NRLib::Grid<double> divided_grid(nx, ny, nz);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

      std::vector<double> mu_log_vp = generateMuFromGrid(post_mu_vp, i, j);

      // Transform to (Vp0/Vp1)
      std::vector<double>   mu_vp_minus;
      NRLib::Grid2D<double> Sigma_vp_minus;
      calculateMinusFirstCentralMomentLogNormal(mu_log_vp,
                                                Sigma_log_vp,
                                                mu_vp_minus,
                                                Sigma_vp_minus);

      for (int k = 0; k < nz; k++)
        divided_grid(i, j, k) = mu_vp_minus[k];
    }
  }

  return divided_grid;
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::generateNewSimbox(const NRLib::Grid<double>  & vpRatio,
                                       const double               & lz_limit,
                                       const Simbox               * simbox,
                                       Simbox                    *& new_simbox,
                                       std::string                & errTxt) const
{
  Surface base_surface;
  calculateBaseSurface(vpRatio, simbox, base_surface);

  const int nz = simbox->getnz();

  Surface top_surface(dynamic_cast<const Surface &> (simbox->GetTopSurface()));

  new_simbox = new Simbox(simbox);
  new_simbox->setDepth(top_surface, base_surface, nz);
  new_simbox->calculateDz(lz_limit, errTxt);

}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::generateResampleGrid(const NRLib::Grid<double> & v2v1,// v2v1  contains Vp_2(t1)/Vp_1(t1)  in the previous timeframe (t1)
                                          const Simbox              * old_simbox,
                                          const Simbox              * new_simbox,
                                          NRLib::Grid<double>       & resample_grid) const
{
  // Resample grid (computed here) tells where in the old simbox we should look for a value for the given cell in the new simbox;
  // The top of the top cell is zero by definition
  // The center of the cell is top + d2/2
  int nx  = new_simbox->getnx();
  int ny  = new_simbox->getny();
  int nz  = new_simbox->getnz();

  resample_grid.Resize(nx, ny, nz, 0);

  double tk; // is center of cell k
  double d1; // time from top to base of a cell in the old time grid
  double d2; // time from top to base of a cell in the new time grid
  double r;  // Position of tk between t2(l), t2(l+1)
  int    l;  // Largest integer satisfying t2(l) <= tk

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
     // double b1 = old_simbox->getBot(i, j);
     // double b2 = new_simbox->getBot(i,j);
      d1 = old_simbox->getdz(i, j);
      d2 = new_simbox->getdz(i, j);
      std::vector<double> t2(nz + 1, 0); // edges of t2 as a function of t1, i.e. t2[k]=t2(k*d1)
      for (int k = 1; k < nz + 1; ++k)
        t2[k] = t2[k - 1] + d1/v2v1(i, j, k - 1);

      l = 0;
      for (int k = 0; k < nz; k++) {
        tk = (k+0.5) * d2;

        while(tk >= t2[l])
          l++;
        l--;

        r = (tk - t2[l]) / (t2[l+1] - t2[l]);

        resample_grid(i, j, k) =  (l + r) * d1;
      }

    }
  }
}

FFTGrid *
TravelTimeInversion::generateResampleAveragePreserve( FFTGrid        * v2v0_1,   //  v2v0  contains Vp_2(t1)/Vp_0(t1)    in the previous timeframe (t1)
                                                     const FFTGrid        * v1v0_1,   //  v1v0_1  contains Vp_1(t1)/Vp_0(t1)  in the previous timeframe (t1)
                                                     const Simbox         * simbox1,  //  simbox in the previous timeframe (t1)
                                                     const Simbox         * simbox2) const  //  simbox in the current timeframe (t2)

{
  // returns  FFTGrid * v2v0_2 which contains Vp_2(t2)/Vp_0(t2)  in the current timeframe (t2)

  int nx  = simbox1->getnx();
  int ny  = simbox1->getny();
  int nz  = simbox1->getnz();

   FFTGrid * v2v0_2= new  FFTGrid(v2v0_1);


  double d1; // time from top to base of a cell in the old time grid
  double d2; // time from top to base of a cell in the new time grid
  double r;  // Position of tk between t2(l), t2(l+1)

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {

   //   double b1 = simbox1->getBot(i, j);
   //   double b2 = simbox2->getBot(i,j);

      d1 = simbox1->getdz(i, j);
      d2 = simbox2->getdz(i, j);

      std::vector<double> t0_1(nz + 1, 0); // edges of t0 as a function of t1, i.e. t0_1[k]=t0_1(k*d1)
      std::vector<double> t2_1(nz + 1, 0); // edges of t2 as a function of t1, i.e. t2_1[k]=t2_1(k*d1)
      for (int k =0 ; k < nz ; ++k){
        t0_1[k+1] = t0_1[k] + v1v0_1->getRealValue(i,j,k)*d1;
        t2_1[k+1] = t2_1[k] + v1v0_1->getRealValue(i,j,k)/v2v0_1->getRealValue(i,j,k)*d1;
      }
      t2_1[nz]= static_cast<double>(nz)*d2;

      std::vector<double> t1_2(nz + 1, 0); // edges of t1 as a function of t2, i.e. t1_2[k]=t1_2(k*d2)

      int l = 0;// Largest integer satisfying t2(l) <  tkend  but still <= nz, purely mathematically
      for (int k = 0; k < nz; k++) {  // the last tkend = t2[nz] =(nz)*d2;
        double tkend = static_cast<double>(k+1) * d2;

        while(l<=nz && tkend > t2_1[l] )
          l++;
        l--;  int indHi= std::min(nz,l+1);
        r = (tkend - t2_1[l]) / (t2_1[indHi] - t2_1[l]);
        t1_2[k+1] =  (l + r) * d1;
      }

      std::vector<double> t0_2(nz + 1, 0); // edges of t0 as a function of t2, i.e. t0_2[k]=t0_2(k*d2)
      for (int k = 1; k < nz+1; k++) {

        int indLo = static_cast<int>(floor(t1_2[k]/d1));
        int indHi = static_cast<int>( ceil(t1_2[k]/d1));
        double w  = ( t1_2[k]/d1- (indLo));
        indLo     = std::max(0,indLo);
        indHi     = std::min(nz,indHi);

        double valueLo = t0_1[indLo];
        double valueHi = t0_1[indHi];

        double value = valueLo * (1 - w) + valueHi * w;
        t0_2[k]=value;
        float relVel=  static_cast<float>((t0_2[k]-t0_2[k-1]) /d2); //    dt0 /dt2= v2/v0;
        v2v0_2->setRealValue(i,j,k-1,relVel );
      }
    }
  }

  return v2v0_2;
}



//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::resampleFFTGrid(const NRLib::Grid<double> &  resample_grid,
                                     const Simbox              *  previous_simbox,
                                     FFTGrid                   *& grid) const
{
  grid->setAccessMode(FFTGrid::READANDWRITE);

  int nx = previous_simbox->getnx();
  int ny = previous_simbox->getny();
  int nz = previous_simbox->getnz();

  // Resample grid tells where in the old simbox we should look for a value for the given cell in the new simbox;
  // the location is given as a time value, so any useful interpolation may be applied
  // note that the values are considred as the center of a cell. This is why the -0.5 and +0.5 appear
  // ;


  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {

      double dz = previous_simbox->getdz(i, j);


      std::vector<double> profile(nz);
      for (int k = 0; k < nz; ++k) {
         profile[k]=grid->getRealValue(i, j, k);
      }


      for (int k = 0; k < nz; ++k) {
        int indLo = static_cast<int>(floor(resample_grid(i, j, k) / dz - 0.5));
        int indHi = static_cast<int>( ceil(resample_grid(i, j, k) / dz - 0.5));
        double t  = (resample_grid(i, j, k)/dz - (indLo+0.5));
        indLo     = std::max(0,indLo);
        indHi     = std::min(nz-1,indHi);

        double valueLo = profile[indLo];
        double valueHi = profile[indHi];


        double value = valueLo * (1 - t) + valueHi * t;

        grid->setRealValue(i, j, k, static_cast<float>(value));
      }
    }
  }

  grid->endAccess();

}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::resampleState4D(const NRLib::Grid<double> &  resample_grid,
                                     const Simbox              *  previous_simbox,
                                     FFTGrid                   *& mu_vp_static,
                                     FFTGrid                   *& mu_vs_static,
                                     FFTGrid                   *& mu_rho_static,
                                     FFTGrid                   *& mu_vp_dynamic,
                                     FFTGrid                   *& mu_vs_dynamic,
                                     FFTGrid                   *& mu_rho_dynamic) const
{
  std::string text = "\nResampling background...";
  LogKit::LogFormatted(LogKit::Low, text);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_vp_static);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_vp_dynamic);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_vs_static);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_vs_dynamic);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_rho_static);

  resampleFFTGrid(resample_grid,
                  previous_simbox,
                  mu_rho_dynamic);

  LogKit::LogFormatted(LogKit::Low,"done\n\n");

}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::resampleSeismicParameters(const NRLib::Grid<double> & resample_grid,
                                               const Simbox              * old_simbox,
                                               SeismicParametersHolder   & seismic_parameters) const
{
  std::string text = "\nResampling background...";
  LogKit::LogFormatted(LogKit::Low, text);

  FFTGrid * mu_alpha = seismic_parameters.GetMuAlpha();
  FFTGrid * mu_beta  = seismic_parameters.GetMuBeta();
  FFTGrid * mu_rho   = seismic_parameters.GetMuRho();

  resampleFFTGrid(resample_grid,
                  old_simbox,
                  mu_alpha);

  resampleFFTGrid(resample_grid,
                  old_simbox,
                  mu_beta);

  resampleFFTGrid(resample_grid,
                  old_simbox,
                  mu_rho);

  seismic_parameters.setBackgroundParameters(mu_alpha,
                                             mu_beta,
                                             mu_rho);

  LogKit::LogFormatted(LogKit::Low,"done\n\n");
}

//-----------------------------------------------------------------------------------------//

void
TravelTimeInversion::calculateBaseSurface(const NRLib::Grid<double> & ratioVp,
                                          const Simbox              * simbox,
                                          Surface                   & base_surface) const
{

  Surface top_surface(dynamic_cast<const Surface &> (simbox->GetTopSurface()));

  base_surface = top_surface;

  //
  // If a constant time surface has been used, it may have only four grid
  // nodes. To handle this situation we use the grid resolution whenever
  // this is larger than the surface resolution.
  //

  int ni = static_cast<int>(top_surface.GetNI());
  int nj = static_cast<int>(top_surface.GetNJ());

  int maxNx = std::max(simbox->getnx(), ni);
  int maxNy = std::max(simbox->getny(), nj);

  base_surface.Resize(maxNx, maxNy);

  double dx = 0.5 * top_surface.GetDX();
  double dy = 0.5 * top_surface.GetDY();

  int nz = simbox->getnz();

  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni; ++i) {

      double x, y;
      top_surface.GetXY(i, j, x, y);

      int ii, jj;
      simbox->getIndexes(x, y, ii, jj);


      if (ii != IMISSING && jj != IMISSING) {
        double sum = 0.0;
        double dt =simbox->getdz(ii,jj);
        for (int k = 0; k < nz; ++k)
          sum += dt/ratioVp(ii, jj, k);

        base_surface(i, j) = sum;
      }
      else {
        int i1, i2, i3, i4, j1, j2, j3, j4;
        simbox->getIndexes(x + dx, y + dy, i1, j1);
        simbox->getIndexes(x - dx, y - dy, i2, j2);
        simbox->getIndexes(x + dx, y - dy, i3, j3);
        simbox->getIndexes(x - dx, y + dy, i4, j4);

        int n = 0;
        if (i1 != IMISSING && j1 != IMISSING)
          n++;
        if (i2 != IMISSING && j2 != IMISSING)
          n++;
        if (i3 != IMISSING && j3 != IMISSING)
          n++;
        if (i4 != IMISSING && j4 != IMISSING)
          n++;

        if (n == 0)
          base_surface.SetMissing(i, j);

        else {
          double sum = 0.0;
          double dt1 =simbox->getdz(i1,j1);
          double dt2 =simbox->getdz(i2,j2);
          double dt3 =simbox->getdz(i3,j3);
          double dt4 =simbox->getdz(i4,j4);

          for (int k = 0; k < nz; ++k) {
            if (i1 != IMISSING && j1 != IMISSING)
              sum += dt1/ratioVp(i1, j1, k);

            if (i2 != IMISSING && j2 != IMISSING)
              sum += dt2/ratioVp(i2, j2, k);

            if (i3 != IMISSING && j3 != IMISSING)
              sum += dt3/ratioVp(i3, j3, k);

            if (i4 != IMISSING && j4 != IMISSING)
              sum += dt4/ratioVp(i4, j4, k);
          }

          base_surface(i, j) = sum / static_cast<double>(n);
        }
      }
    }
  }

  base_surface.AddNonConform(&top_surface);
}

//-----------------------------------------------------------------------------------------//
void
TravelTimeInversion::generateTimeDepthMapping(FFTGrid       * post_mu_log_vp_above,
                                              FFTGrid       * post_cov_log_vp_above,
                                              FFTGrid       * mu_log_vp_grid,
                                              FFTGrid       * cov_log_vp_grid,
                                              int             output_format,
                                              const Simbox  * simbox_above,
                                              const Simbox  * timeSimbox,
                                              GridMapping  *& grid_depth_mapping) const
{
  LogKit::LogFormatted(LogKit::Low, "\nGenerating depth simbox...");

  FFTGrid * post_mu_vp_above  = NULL;
  calculateEVpGrid(post_mu_log_vp_above,
                   post_cov_log_vp_above,
                   post_mu_vp_above);

  GridMapping * grid_mapping_above = new GridMapping();

  Surface top_surface_above(dynamic_cast<const Surface &> (simbox_above->GetTopSurface()));

  grid_mapping_above->setTopSurface(top_surface_above);

  post_mu_vp_above  ->setAccessMode(FFTGrid::RANDOMACCESS);
  grid_mapping_above->calculateSurfaceFromVelocity(post_mu_vp_above, simbox_above);
  post_mu_vp_above  ->endAccess();

  Surface top_depth_surface(dynamic_cast<const Surface &> (grid_mapping_above->GetBaseSurface()));

  grid_depth_mapping = new GridMapping();

  grid_depth_mapping->setTopSurface(top_depth_surface);

  FFTGrid * post_mu_vp_grid_model = NULL;
  calculateEVpGrid(mu_log_vp_grid,
                   cov_log_vp_grid,
                   post_mu_vp_grid_model);

  bool failed = false;
  std::string errTxt = "";

  post_mu_vp_grid_model->setAccessMode(FFTGrid::RANDOMACCESS);

  grid_depth_mapping->calculateSurfaceFromVelocity(post_mu_vp_grid_model,
                                                   timeSimbox);

  LogKit::LogFormatted(LogKit::Low,"done\n\n");

  grid_depth_mapping->setDepthSimbox(timeSimbox,
                                     timeSimbox->getnz(),
                                     output_format,
                                     failed,
                                     errTxt);

  grid_depth_mapping->makeTimeDepthMapping(post_mu_vp_grid_model,
                                           timeSimbox);

  post_mu_vp_grid_model->endAccess();

  delete grid_mapping_above;
  delete post_mu_vp_above;
  delete post_mu_vp_grid_model;


}

//-----------------------------------------------------------------------------------------//

std::vector<Surface>
TravelTimeInversion::sortHorizons(const std::vector<Surface> & initial_horizons,
                                  const std::vector<Surface> & push_down_horizons,
                                  const std::vector<std::string> & initial_horizon_names,
                                  const std::vector<std::string> & push_down_names) const
{
  int n_push_down = static_cast<int>(push_down_horizons.size());
  int n_initial   = static_cast<int>(initial_horizons.size());

  std::vector<Surface> sorted_initial(n_push_down);

  for (int i = 0; i < n_push_down; i++) {
    for (int j = 0; j < n_initial; j++) {
      if (push_down_names[i] == initial_horizon_names[j]) {
        sorted_initial[i] = initial_horizons[j];
        break;
      }
    }
  }

  return sorted_initial;
}
