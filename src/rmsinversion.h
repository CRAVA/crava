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
                                              const double                & mean_vp_top,
                                              const double                & mean_vp_base,
                                              const double                & var_vp_above,
                                              const double                & var_vp_below,
                                              const double                & max_time,
                                              const RMSTrace              * rms_trace,
                                              const FFTGrid               * mean_log_vp,
                                              const std::vector<double>   & cov_grid_log_vp,
                                              const Simbox                * timeSimbox) const;

  double                        findMaxTime(const std::vector<RMSTrace *> & rms_traces) const;

  std::vector<double>           getCovLogVp(const FFTGrid * cov_log_vp) const;

  std::vector<double>           calculateMeanVp(const double & top_value,
                                                const double & base_value,
                                                const int    & n_layers) const;

  NRLib::Grid2D<double>         calculateG(const std::vector<double> & rms_time,
                                           const double              & t_top,
                                           const double              & t_bot,
                                           const double              & dt_simbox,
                                           const double              & max_time,
                                           const int                 & n_layers_above,
                                           const int                 & n_layers_below,
                                           const int                 & n_layers_simbox,
                                           const int                 & n_layers_padding) const;

  NRLib::Grid2D<double>         makeSigma_m(const double              & t_top,
                                            const double              & t_bot,
                                            const double              & max_time,
                                            const Vario               * variogram_above,
                                            const Vario               * variogram_below,
                                            const double              & var_vp_above,
                                            const double              & var_vp_below,
                                            const std::vector<double> & cov_grid_log_vp,
                                            const int                 & n_layers_above,
                                            const int                 & n_layers_below) const;


};

#endif
