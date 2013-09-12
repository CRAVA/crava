/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef RMS_INVERSION_H
#define RMS_INVERSION_H

#include "nrlib/grid/grid2d.hpp"

class ModelGeneral;
class SeismicParametersHolder;
class ModelTravelTimeDynamic;
class Simbox;
class RMSTrace;
class Vario;
class FFTGrid;

class RMSInversion
{
public:
  RMSInversion(const ModelGeneral      * modelGeneral,
               ModelTravelTimeDynamic  * modelTravelTimeDynamic,
               SeismicParametersHolder & seismicParameters);

  ~RMSInversion();

private:
  void                          do1DInversion(const int                   & n_layers_above,
                                              const int                   & n_layers_below,
                                              const int                   & n_layers_simbox,
                                              const int                   & n_layers_padding,
                                              const Vario                 * variogram_above,
                                              const Vario                 * variogram_below,
                                              const double                & mu_vp_top,
                                              const double                & mu_vp_base,
                                              const double                & var_vp_above,
                                              const double                & var_vp_below,
                                              const double                & max_time,
                                              const RMSTrace              * rms_trace,
                                              const FFTGrid               * mu_log_vp,
                                              const std::vector<double>   & cov_grid_log_vp,
                                              const Simbox                * timeSimbox) const;

  NRLib::Grid2D<double>         calculateG(const std::vector<double> & rms_time,
                                           const double              & t_top,
                                           const double              & t_bot,
                                           const double              & dt_simbox,
                                           const double              & max_time,
                                           const int                 & n_layers_above,
                                           const int                 & n_layers_below,
                                           const int                 & n_layers_simbox,
                                           const int                 & n_layers_padding) const;

  std::vector<double>           generateMuLogVpModel(const FFTGrid * mu_log_vp,
                                                     const int     & i_ind,
                                                     const int     & j_ind) const;

  std::vector<double>           generateMuVpAbove(const double & top_value,
                                                  const double & base_value,
                                                  const int    & n_layers) const;

  std::vector<double>           generateMuVpBelow(const double & top_value,
                                                  const double & base_value,
                                                  const int    & n_layers) const;

  std::vector<double>           generateMuVp(const double & top_value,
                                             const double & base_value,
                                             const int    & n_layers) const;

  std::vector<double>           generateMuCombined(const std::vector<double> & mu_above,
                                                   const std::vector<double> & mu_model,
                                                   const std::vector<double> & mu_below) const;

  NRLib::Grid2D<double>         generateSigmaCombined(const NRLib::Grid2D<double> & Sigma_above,
                                                      const NRLib::Grid2D<double> & Sigma_model,
                                                      const NRLib::Grid2D<double> & Sigma_below) const;

  NRLib::Grid2D<double>         generateSigma(const double             & var,
                                              const std::vector<float> & corrT) const;

  NRLib::Grid2D<double>         generateSigmaVp(const float  & dt,
                                                const int    & n_layers,
                                                const double & var_vp,
                                                const Vario  * variogram) const;

  NRLib::Grid2D<double>         generateSigmaModel(const std::vector<double> & cov_grid) const;

  std::vector<double>           getCovLogVp(const FFTGrid * cov_log_vp) const;

  void                          transformVpToVpSqrt(const std::vector<double>   & mu_vp,
                                                    const NRLib::Grid2D<double> & Sigma_vp,
                                                    std::vector<double>         & mu_vp_sqrt,
                                                    NRLib::Grid2D<double>       & Sigma_vp_sqrt) const;

  void                          calculateCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                const NRLib::Grid2D<double> & variance_log_vp,
                                                                std::vector<double>         & mu_vp_trans,
                                                                NRLib::Grid2D<double>       & variance_vp_trans) const;

  void                          calculateSecondCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                                      const NRLib::Grid2D<double> & variance_log_vp,
                                                                      std::vector<double>         & mu_vp_sqrt,
                                                                      NRLib::Grid2D<double>       & variance_vp_sqrt) const;

  void                          calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                                       const NRLib::Grid2D<double> & variance_vp_trans,
                                                                       std::vector<double>         & mu_log_vp,
                                                                       NRLib::Grid2D<double>       & variance_log_vp) const;

  double                        findMaxTime(const std::vector<RMSTrace *> & rms_traces) const;

};

#endif
