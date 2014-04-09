/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef TRAVEL_TIME_INVERSION_H
#define TRAVEL_TIME_INVERSION_H

#include "nrlib/flens/nrlib_flens.hpp"
#include "nrlib/grid/grid2d.hpp"
#include "src/definitions.h"
#include "lib/utils.h"

class ModelGeneral;
class SeismicParametersHolder;
class ModelTravelTimeStatic;
class ModelTravelTimeDynamic;
class Simbox;
class RMSTrace;
class Vario;
class FFTGrid;
class KrigingData2D;
class State4D;
class GridMapping;

class TravelTimeInversion
{
public:
  TravelTimeInversion(ModelGeneral            * modelGeneral,
                      ModelTravelTimeStatic   * modelTravelTimeStatic,
                      ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                      SeismicParametersHolder & seismicParameters);

  ~TravelTimeInversion();

private:
  void                          doHorizonInversion(ModelGeneral            * modelGeneral,
                                                   ModelTravelTimeStatic   * modelTravelTimeStatic,
                                                   ModelTravelTimeDynamic  * modelTravelTimeDynamic) const;

  void                          do1DHorizonInversion(FFTGrid               * mu_prior,
                                                     const NRLib::Grid2D<double> & Sigma_log_vp,
                                                     const Simbox                * timeSimbox,
                                                     const FFTGrid               * relativeVelocityPrev,
                                                     const std::vector<Surface>  & initial_horizons,
                                                     const std::vector<Surface>  & push_down_horizons,
                                                     const std::vector<double>   & standard_deviation,
                                                     const Surface               & top_simbox,
                                                     const Surface               & base_simbox,
                                                     int                           i_ind,
                                                     int                           j_ind,
                                                     std::vector<double>         & mu_post,
                                                     NRLib::Grid2D<double>       & Sigma_post_log_vp,
                                                     bool                          logTransMean) const;  // true if the mean is logarithmic transformed fals if mu_prior is E(V_P)

  void                          doRMSInversion(ModelGeneral            * modelGeneral,
                                               ModelTravelTimeStatic   * modelTravelTimeStatic,
                                               ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                                               SeismicParametersHolder & seismicParameters) const;

  void                          do1DRMSInversion(const double                & mu_vp_base,
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
                                                 NRLib::Grid2D<double>       & Sigma_post_log_vp) const;

  void                          calculatePosteriorModel(const std::vector<double>   & d,
                                                        const NRLib::Grid2D<double> & Sigma_d,
                                                        const std::vector<double>   & mu_m,
                                                        const NRLib::Grid2D<double> & Sigma_m,
                                                        const NRLib::Grid2D<double> & G,
                                                        std::vector<double>         & mu_post,
                                                        NRLib::Grid2D<double>       & Sigma_post) const;

  NRLib::Grid2D<double>         calculateGHorizon(double                      dt,
                                                  std::vector<double>         relativeVelocity,
                                                  double                      top,
                                                  double                      missing_value,
                                                  int                         n_nonmissing,
                                                  int                         nz,
                                                  int                         nzp,
                                                  const std::vector<double> & time_P0,
                                                  const std::vector<double> & push_down) const;

  NRLib::Grid2D<double>         calculateG(const std::vector<double> & rms_time,
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
                                           const int                 & n_pad_below) const;

  std::vector<double>           calculateDSquare(const std::vector<double> & d) const;
  void                          biasAdjustDsquare(std::vector<double> & d2, double variance) const;

  NRLib::Grid2D<double>         calculateSigmaDSquare(const std::vector<double> & rms_velocity,
                                                    const double              & standard_deviation) const;

  void                          calculateMuSigma_mSquare(const std::vector<double>   & mu_log_vp_above,
                                                         const std::vector<double>   & mu_log_vp_model,
                                                         const std::vector<double>   & cov_grid_log_vp_above,
                                                         const std::vector<double>   & cov_grid_log_vp_model,
                                                         const double                & mu_vp_base,
                                                         const NRLib::Grid2D<double> & Sigma_vp_below,
                                                         const int                   & n_below,
                                                         std::vector<double>         & mu_vp_square,
                                                         NRLib::Grid2D<double>       & Sigma_vp_square) const;
  void                          transformCovarianceToRelativeScale(const std::vector<double>     mu,
                                                                         NRLib::Grid2D<double> & Sigma,
                                                                         int n1, int n1p,
                                                                         int n2, int n2p,
                                                                         int n3, int n3p) const;

  void                          generateMuSigmaLogVpAbove(const int                    & nz,
                                                          const int                    & nzp,
                                                          const double                 & mu_vp_top,
                                                          const NRLib::Grid2D<double>  & Sigma_vp_above,
                                                          const Surface                * errorCorrXY,
                                                          const float                  & corrGradI,
                                                          const float                  & corrGradJ,
                                                          FFTGrid                      * mu_log_vp_grid,
                                                          FFTGrid                     *& mu_log_vp_above,
                                                          FFTGrid                     *& Sigma_log_vp_above) const;

  std::vector<double>           generateMuFromGrid(FFTGrid   * muGrid,
                                                        const int & i_ind,
                                                        const int & j_ind) const;

  std::vector<double>           generateMuVpBelow(const double & top_value,
                                                  const double & base_value,
                                                  const int    & nz,
                                                  const int    & nzp) const;

  std::vector<double>           generateMuVp(const double & top_value,
                                             const double & base_value,
                                             const int    & n_layers) const;

  std::vector<double>           generateMuCombined(const std::vector<double> & mu_above,
                                                   const std::vector<double> & mu_model,
                                                   const std::vector<double> & mu_below) const;

  NRLib::Grid2D<double>         generateSigmaCombined(const NRLib::Grid2D<double> & Sigma_above,
                                                      const NRLib::Grid2D<double> & Sigma_model,
                                                      const NRLib::Grid2D<double> & Sigma_below) const;

  NRLib::Grid2D<double>         generateSigma(const double              & var,
                                              const std::vector<double> & corrT) const;

  NRLib::Grid2D<double>         generateSigmaVp(const double & dt,
                                                const int    & n_layers,
                                                const double & var_vp,
                                                const Vario  * variogram) const;

  NRLib::Grid2D<double>         generateSigmaModel(const std::vector<double> & cov_grid) const;

  std::vector<double>           getCovLogVp(FFTGrid * cov_log_vp) const;

  void                          transformVpSquareToLogVp(const std::vector<double>   & mu_vp_square,
                                                         const NRLib::Grid2D<double> & Sigma_vp_square,
                                                         std::vector<double>         & mu_log_vp,
                                                         NRLib::Grid2D<double>       & Sigma_log_vp) const;

  void                          getMuFromLogMu(std::vector<double>   & mu , const std::vector<double>   &mu_log, const NRLib::Grid2D<double> &Sigma_log) const;
  void                          getLogMuFromMu(std::vector<double>   & mu_log,  const std::vector<double>   &mu_, const NRLib::Grid2D<double> &Sigma_log) const;
  void                          getSigmaFromLogSigma(std::vector<double>   & mu,const NRLib::Grid2D<double> &Sigma_log, NRLib::Grid2D<double> & Sigma) const;
  FFTGrid*                      getCovFunkFromLogCovFunk(FFTGrid*  cov_log,double meanVpRelative) const;

  void                          transformVpMinusToLogVp(const std::vector<double>   & mu_vp_minus,
                                                        const NRLib::Grid2D<double> & Sigma_vp_minus,
                                                        std::vector<double>         & mu_log_vp,
                                                        NRLib::Grid2D<double>       & Sigma_log_vp) const;

  void                          calculateCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                const NRLib::Grid2D<double> & variance_log_vp,
                                                                std::vector<double>         & mu_vp_trans,
                                                                NRLib::Grid2D<double>       & variance_vp_trans) const;

  void                          calculateSecondCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                      const NRLib::Grid2D<double> & variance_log_vp,
                                                                      std::vector<double>         & mu_vp_square,
                                                                      NRLib::Grid2D<double>       & variance_vp_square) const;

  void                          calculateMinusFirstCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                          const NRLib::Grid2D<double> & variance_log_vp,
                                                                          std::vector<double>         & mu_vp_minus,
                                                                          NRLib::Grid2D<double>       & variance_vp_minus) const;

  void                          calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                                       const NRLib::Grid2D<double> & variance_vp_trans,
                                                                       std::vector<double>         & mu_log_vp,
                                                                       NRLib::Grid2D<double>       & variance_log_vp) const;

  void                          setExpectation(int                          i_ind,
                                               int                          j_ind,
                                               const std::vector<double>  & post_vp,
                                               std::vector<KrigingData2D> & mu_log_vp_post) const;

  void                          addCovariance(const NRLib::Grid2D<double> & Sigma_post,
                                              std::vector<double>         & cov_stationary,
                                              int                           n_nopad) const;
  void                          addCovarianceMat(const NRLib::Grid2D<double> & Sigma_post,
                                                  NRLib::Matrix         & cov_cum) const;

  std::vector<double>           makeCirculantCovariance(const NRLib::Matrix & cov,
                                                        const int                   & n_nopad) const;
  std::vector<double>           makeCirculantCovariance(const NRLib::Grid2D<double> & cov,
                                                        const int                   & n_nopad) const;

  void                           krigeExpectation3D(const Simbox                * simbox,
                                                   std::vector<KrigingData2D>  & kriging_post,
                                                   const int                   & nxp,
                                                   const int                   & nyp,
                                                   FFTGrid                    *& mu_post) const;

  void                          generateStationaryDistribution(const Simbox                * timeSimbox,
                                                               std::vector<KrigingData2D>  & kriging_post,
                                                               const std::vector<double>   & pri_circulant_cov,
                                                               NRLib::Matrix                 post_cov,
                                                               const int                   & n_rms_traces,
                                                               const Surface   * errorCorrXY,
                                                               const float & corrGradI,
                                                               const float & corrGradJ,
                                                               FFTGrid                     * pri_mu,
                                                               FFTGrid                    *& stationary_observations,
                                                               FFTGrid                    *& stationary_covariance,
                                                               std::vector<int>            & observation_filter) const;

  void                          calculateStationaryObservations(const fftw_complex  *  pri_cov_c,
                                                                const fftw_complex  *  var_e_c,
                                                                const std::vector<int> filter_c,
                                                                FFTGrid             *  pri_mu,
                                                                FFTGrid             *  post_mu,
                                                                FFTGrid             *& stat_d) const;

  void                          calculateErrorVariance(const fftw_complex * pri_cov_c,
                                                       const fftw_complex * post_cov_c,
                                                       const int          & nzp,
                                                       fftw_complex       * var_e_c) const;

  void                          calculateFilter(const fftw_complex * pri_cov_c,
                                                const fftw_complex * post_cov_c,
                                                const int          & nzp,
                                                std::vector<int>   & filter) const;

  void                          calculateFilterAndErrorVariance(const fftw_complex * pri_cov_c,
                                                                const NRLib::Matrix post_cov,
                                                                const int       &  nz,
                                                                fftw_complex       * var_e_c_lo,
                                                                fftw_complex       * var_e_c_hi,
                                                                std::vector<int>   & filter) const;

  void                          multiplyComplex(const fftw_complex * z1,
                                                const fftw_complex * z2,
                                                const int          & n,
                                                fftw_complex       * z) const;

  void                          divideComplex(const fftw_complex * z1,
                                              const fftw_complex * z2,
                                              const int          & n,
                                              fftw_complex       * z) const;

  void                          addComplex(const fftw_complex * z1,
                                           const fftw_complex * z2,
                                           const int          & n,
                                           fftw_complex       * z) const;

  void                          subtractComplex(const fftw_complex * z1,
                                                const fftw_complex * z2,
                                                const int          & n,
                                                fftw_complex       * z) const;

  void                          absoulteComplex(const fftw_complex  * z,
                                                const int           & n,
                                                std::vector<double> & abs_z) const;

  void                          complexConjugate(const fftw_complex * z,
                                                 const int          & n,
                                                 fftw_complex       * conj_z) const;

  void                          calculateFullPosteriorModel(const std::vector<int>  & observation_filter,
                                                            SeismicParametersHolder & seismic_parameters,
                                                            FFTGrid                 * stationary_observations,
                                                            FFTGrid                 * stationary_observation_covariance) const;

  void                          calculateExpectation(const std::vector<int>  & observation_filter,
                                                          FFTGrid                 * prior_mu,
                                                          FFTGrid                 * prior_cov,
                                                          FFTGrid                 * stationary_observations,
                                                          FFTGrid                 * stationary_observation_covariance,
                                                          FFTGrid                *& post_mu) const;

  void                          calculateLogVpCovariance(const std::vector<int>  & observation_filter,
                                                      //   FFTGrid                 * mu_vp,
                                                         FFTGrid                 * cov_vp,
                                                       //  FFTGrid                 * stationary_observations,
                                                         FFTGrid                 * stationary_observation_covariance,
                                                         FFTGrid                *& post_cov_vp) const;

  void                          calculateDistanceGrid(const Simbox              * simbox,
                                                      const NRLib::Grid<double> & divided_grid,
                                                      NRLib::Grid<double>       & distance) const;

  void                          calculateEVpGrid(FFTGrid  * mu_log_vp,
                                                 FFTGrid  * cov_log_vp,
                                                 FFTGrid *& mu_vp) const;
 NRLib::Grid<double>            calculateRelativeVelocityUpdate(FFTGrid *relativeVelocityGridNew,
                                                                FFTGrid *relativeVelocityPrev) const;

 NRLib::Grid<double>            calculateDividedGridRMS(FFTGrid * pri_vp,
                                                        FFTGrid * post_vp) const;

 NRLib::Grid<double>            calculateDividedGridHorizon(FFTGrid * post_mu_vp,
                                                            FFTGrid * post_cov_mu_vp) const;

  void                          generateNewSimbox(const NRLib::Grid<double>  & distance,
                                                  const double               & lz_limit,
                                                  const Simbox               * simbox,
                                                  Simbox                    *& new_simbox,
                                                  std::string                & errTxt) const;

  void                          generateResampleGrid(const NRLib::Grid<double> & v2v1,
                                                     const Simbox              * old_simbox,
                                                     const Simbox              * new_simbox,
                                                     NRLib::Grid<double>       & resample_grid) const;

  FFTGrid *                    generateResampleAveragePreserve( FFTGrid        * v2v0_1,          //  v2v0  contains Vp_2(t1)/Vp_0(t1)    in the previous timeframe (t1)
                                                          const FFTGrid        * v1v0_1,          //  v1v0_1  contains Vp_1(t1)/Vp_0(t1)  in the previous timeframe (t1)
                                                          const Simbox         * simbox1,         //  simbox in the previous timeframe (t1)
                                                          const Simbox         * simbox2) const;  //  simbox in the current timeframe (t2)

  void                          calculateBaseSurface(const NRLib::Grid<double> & distance,
                                                     const Simbox              * simbox,
                                                     Surface                   & base_surface) const;

  void                          resampleState4D(const NRLib::Grid<double> &  resample_grid,
                                                const Simbox              *  previous_simbox,
                                                FFTGrid                   *& mu_vp_static,
                                                FFTGrid                   *& mu_vs_static,
                                                FFTGrid                   *& mu_rho_static,
                                                FFTGrid                   *& mu_vp_dynamic,
                                                FFTGrid                   *& mu_vs_dynamic,
                                                FFTGrid                   *& mu_rho_dynamic) const;

  void                          resampleSeismicParameters(const NRLib::Grid<double> & resample_grid,
                                                          const Simbox              * new_simbox,
                                                          SeismicParametersHolder   & seismic_parameters) const;

  void                          resampleFFTGrid(const NRLib::Grid<double> &  resample_grid,
                                                const Simbox              *  old_simbox,
                                                FFTGrid                   *& grid) const;

  void                          generateTimeDepthMapping(FFTGrid       * post_mu_log_vp_above,
                                                         FFTGrid       * post_cov_log_vp_above,
                                                         FFTGrid       * mu_log_vp_grid,
                                                         FFTGrid       * cov_log_vp_grid,
                                                         int             output_format,
                                                         const Simbox  * simbox_above,
                                                         const Simbox  * timeSimbox,
                                                         GridMapping  *& grid_depth_mapping) const;

  std::vector<Surface>          sortHorizons(const std::vector<Surface> & initial_horizons,
                                             const std::vector<Surface> & push_down_horizons,
                                             const std::vector<std::string> & initial_horizon_names,
                                             const std::vector<std::string> & push_down_names) const;

};

#endif
