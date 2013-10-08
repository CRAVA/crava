/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RMS_INVERSION_H
#define RMS_INVERSION_H

#include "nrlib/grid/grid2d.hpp"
#include "src/definitions.h"
#include "lib/utils.h"

class ModelGeneral;
class SeismicParametersHolder;
class ModelTravelTimeDynamic;
class Simbox;
class RMSTrace;
class Vario;
class FFTGrid;
class KrigingData2D;

class RMSInversion
{
public:
  RMSInversion(const ModelGeneral      * modelGeneral,
               ModelTravelTimeDynamic  * modelTravelTimeDynamic,
               SeismicParametersHolder & seismicParameters);

  ~RMSInversion();

private:
  void                          do1DInversion(const double                & mu_vp_top,
                                              const double                & mu_vp_base,
                                              const NRLib::Grid2D<double> & Sigma_m_above,
                                              const NRLib::Grid2D<double> & Sigma_m_below,
                                              const double                & standard_deviation,
                                              const RMSTrace              * rms_trace,
                                              const FFTGrid               * mu_log_vp,
                                              const std::vector<double>   & cov_grid_log_vp,
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

  NRLib::Grid2D<double>         calculateG(const std::vector<double> & rms_time,
                                           const double              & t_top,
                                           const double              & t_bot,
                                           const double              & dt_above,
                                           const double              & dt_below,
                                           const double              & dt_simbox) const;

  std::vector<double>           calculateDSquare(const std::vector<double> & d) const;

  NRLib::Grid2D<double>         calculateSigmaDSquare(const std::vector<double> & rms_velocity,
                                                    const double              & standard_deviation) const;

  void                          calculateMuSigma_mSquare(const std::vector<double>   & mu_log_vp_model,
                                                         const std::vector<double>   & cov_grid_log_vp,
                                                         const double                & mu_vp_top,
                                                         const double                & mu_vp_base,
                                                         const NRLib::Grid2D<double> & Sigma_vp_above,
                                                         const NRLib::Grid2D<double> & Sigma_vp_below,
                                                         std::vector<double>         & mu_vp_square,
                                                         NRLib::Grid2D<double>       & Sigma_vp_square) const;

  std::vector<double>           generateMuLogVpModel(const FFTGrid * mu_log_vp,
                                                     const int     & i_ind,
                                                     const int     & j_ind) const;

  std::vector<double>           generateMuVpAbove(const double & top_value,
                                                  const double & base_value) const;

  std::vector<double>           generateMuVpBelow(const double & top_value,
                                                  const double & base_value) const;

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

  std::vector<double>           getCovLogVp(const FFTGrid * cov_log_vp) const;

  void                          transformVpToVpSquare(const std::vector<double>   & mu_vp,
                                                      const NRLib::Grid2D<double> & Sigma_vp,
                                                      std::vector<double>         & mu_vp_square,
                                                      NRLib::Grid2D<double>       & Sigma_vp_square) const;

  void                          transformVpSquareToVp(const std::vector<double>   & mu_vp_square,
                                                      const NRLib::Grid2D<double> & Sigma_vp_square,
                                                      std::vector<double>         & mu_vp,
                                                      NRLib::Grid2D<double>       & Sigma_vp) const;

  void                          transformVpSquareToLogVp(const std::vector<double>   & mu_vp_square,
                                                         const NRLib::Grid2D<double> & Sigma_vp_square,
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

  void                          calculateHalfCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                    const NRLib::Grid2D<double> & variance_log_vp,
                                                                    std::vector<double>         & mu_vp,
                                                                    NRLib::Grid2D<double>       & variance_vp) const;

  void                          calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                                       const NRLib::Grid2D<double> & variance_vp_trans,
                                                                       std::vector<double>         & mu_log_vp,
                                                                       NRLib::Grid2D<double>       & variance_log_vp) const;

  void                          setExpectation(const RMSTrace             * rms_trace,
                                               const std::vector<double>  & post_vp,
                                               std::vector<KrigingData2D> & mu_log_vp_post_above,
                                               std::vector<KrigingData2D> & mu_log_vp_post_model,
                                               std::vector<KrigingData2D> & mu_log_vp_post_below) const;

  void                          addCovariance(const int                   & n_rms_traces,
                                              const NRLib::Grid2D<double> & Sigma_post,
                                              std::vector<double>         & cov_stationary_above,
                                              std::vector<double>         & cov_stationary_model,
                                              std::vector<double>         & cov_stationary_below) const;

  std::vector<double>           makeCirculantCovariance(const NRLib::Grid2D<double> & cov,
                                                        const int                   & n_nopad) const;

  void                          writeVector(std::string                 file_name,
                                            const std::vector<double> & vec) const;

  void                          writeMatrix(std::string                   file_name,
                                            const NRLib::Grid2D<double> & grid2d) const;

 void                           krigeExpectation3D(const Simbox               * simbox,
                                                   std::vector<KrigingData2D> & kriging_post,
                                                   FFTGrid                    * mu_post) const;

  void                          generateStationaryDistribution(const std::vector<double> & pri_circulant_cov,
                                                               const std::vector<double> & post_circulant_cov,
                                                               const int                 & n_rms_traces,
                                                               const Surface             * priorCorrXY,
                                                               const float               & corrGradI,
                                                               const float               & corrGradJ,
                                                               FFTGrid                   * pri_mu,
                                                               FFTGrid                   * post_mu,
                                                               FFTGrid                  *& stationary_observations,
                                                               FFTGrid                  *& stationary_covariance) const;

  void                          calculateStationaryObservations(const fftw_complex * pri_cov_c,
                                                                const fftw_complex * var_e_c,
                                                                FFTGrid            * pri_mu,
                                                                FFTGrid            * post_mu,
                                                                FFTGrid            * stat_d) const;

  void                          calculateErrorVariance(const fftw_complex * pri_cov_c,
                                                       const fftw_complex * post_cov_c,
                                                       const int          & nzp,
                                                       fftw_complex       * var_e_c) const;

  void                          calculateFilter(const fftw_complex * pri_cov_c,
                                                const fftw_complex * post_cov_c,
                                                const int          & nzp,
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

  int n_above_;
  int n_below_;
  int n_model_;
  int n_pad_above_;
  int n_pad_below_;
  int n_pad_model_;

};

#endif
