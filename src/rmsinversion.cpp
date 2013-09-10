/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/rmsinversion.h"
#include "src/modeltraveltimedynamic.h"
#include "src/seismicparametersholder.h"
#include "src/simbox.h"
#include "src/modelgeneral.h"
#include "src/rmstrace.h"

#include "lib/timekit.hpp"

RMSInversion::RMSInversion(const ModelGeneral      * modelGeneral,
                           ModelTravelTimeDynamic  * modelTravelTimeDynamic,
                           SeismicParametersHolder & seismicParameters)
{

  LogKit::WriteHeader("Building Stochastic RMS Inversion Model");

  time_t time_start;
  time_t time_end;
  time(&time_start);

  double wall = 0.0;
  double cpu  = 0.0;
  TimeKit::getTime(wall,cpu);

  FFTGrid * mean_log_vp = seismicParameters.GetMuAlpha();
  FFTGrid * cov_log_vp  = seismicParameters.GetCovBeta();

  Simbox * timeSimbox = modelGeneral->getTimeSimbox();

  const std::vector<RMSTrace *> rms_traces = modelTravelTimeDynamic->getRMSTraces();

  const int n_layers_above  = modelTravelTimeDynamic->getNLayersAbove();
  const int n_layers_below  = modelTravelTimeDynamic->getNLayersBelow();
  const double mean_vp_top  = modelTravelTimeDynamic->getMeanVpTop();
  const double mean_vp_base = modelTravelTimeDynamic->getMeanVpBase();
  const double var_vp_above = modelTravelTimeDynamic->getVarVpAbove();
  const double var_vp_below = modelTravelTimeDynamic->getVarVpBelow();
  const double range_above  = modelTravelTimeDynamic->getRangeAbove();
  const double range_below  = modelTravelTimeDynamic->getRangeBelow();

  int n_rms_traces     = static_cast<int>(rms_traces.size());
  int n_layers_simbox  = timeSimbox->getnz();
  int n_layers_padding = cov_log_vp->getNzp();

  Vario * variogram_above = new GenExpVario(1, static_cast<float>(range_above));
  Vario * variogram_below = new GenExpVario(1, static_cast<float>(range_below));

  std::vector<float> corrT_above(n_layers_above + 1);
  std::vector<float> corrT_below(n_layers_below + 1);

  double max_time = findMaxTime(rms_traces);

  std::vector<double> cov_grid_log_vp = getCovLogVp(cov_log_vp);

  for(int i=0; i<n_rms_traces; i++) {

    do1DInversion(n_layers_above,
                  n_layers_below,
                  n_layers_simbox,
                  n_layers_padding,
                  variogram_above,
                  variogram_below,
                  mean_vp_top,
                  mean_vp_base,
                  var_vp_above,
                  var_vp_below,
                  max_time,
                  rms_traces[i],
                  mean_log_vp,
                  cov_grid_log_vp,
                  timeSimbox);

  }

  time(&time_end);
  LogKit::LogFormatted(LogKit::DebugLow,"\nTime elapsed :  %d\n",time_end-time_start);

  delete variogram_above;
  delete variogram_below;
}

//-----------------------------------------------------------------------------------------//

RMSInversion::~RMSInversion()
{
}

//-----------------------------------------------------------------------------------------//
void
RMSInversion::do1DInversion(const int                   & n_layers_above,
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
                            const Simbox                * timeSimbox) const
{
  const std::vector<double> time = rms_trace->getTime();

  double utmx = rms_trace->getUtmx();
  double utmy = rms_trace->getUtmy();

  int i_ind;
  int j_ind;

  timeSimbox->getIndexes(utmx, utmy, i_ind, j_ind);

  double t_top     = timeSimbox->getTop(i_ind, j_ind);
  double t_bot     = timeSimbox->getBot(i_ind, j_ind);
  double dt_simbox = timeSimbox->getdz(i_ind,  j_ind);

  double log_vp_top_reservoir  = mean_log_vp->getRealValue(i_ind, j_ind, 0);
  double log_vp_base_reservoir = mean_log_vp->getRealValue(i_ind, j_ind, n_layers_simbox-1);

  std::vector<double> mean_vp_sqrt_above = calculateMeanVp(std::log(mean_vp_top), log_vp_top_reservoir, n_layers_above);
  std::vector<double> mean_vp_sqrt_below = calculateMeanVp(log_vp_base_reservoir, std::log(mean_vp_base), n_layers_below);

  NRLib::Grid2D<double> G = calculateG(time,
                                       t_top,
                                       t_bot,
                                       dt_simbox,
                                       max_time,
                                       n_layers_above,
                                       n_layers_below,
                                       n_layers_simbox,
                                       n_layers_padding);

  NRLib::Grid2D<double> sigma_m = makeSigma_m(t_top,
                                              t_bot,
                                              max_time,
                                              variogram_above,
                                              variogram_below,
                                              var_vp_above,
                                              var_vp_below,
                                              cov_grid_log_vp,
                                              n_layers_above,
                                              n_layers_below);


}
//-----------------------------------------------------------------------------------------//
std::vector<double>
RMSInversion::getCovLogVp(const FFTGrid * cov_log_vp) const
{
  int n_layers_padding = cov_log_vp->getNzp();

  std::vector<double> cov_grid_log_vp(n_layers_padding, 0);
  for(int j=0; j<n_layers_padding; j++)
    cov_grid_log_vp[j] = cov_log_vp->getRealValue(0, 0, j);

  return cov_grid_log_vp;

}
//-----------------------------------------------------------------------------------------//
std::vector<double>
RMSInversion::calculateMeanVp(const double & top_value,
                              const double & base_value,
                              const int    & n_layers) const
{
  std::vector<double> mean_vp(n_layers+1);

  for(int j=0; j<n_layers+1; j++)
    mean_vp[j] = top_value + j * (base_value - top_value) / n_layers;

  return mean_vp;
}


//-----------------------------------------------------------------------------------------//

double
RMSInversion::findMaxTime(const std::vector<RMSTrace *> & rms_traces) const
{

  int n_rms_traces = rms_traces.size();

  double max_time = 0;

  for(int i=0; i<n_rms_traces; i++) {
    const std::vector<double> & rms_time = rms_traces[i]->getTime();
    double max = rms_time[rms_time.size()-1];

    if(max > max_time)
      max_time = max;
  }

  return max_time;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
RMSInversion::calculateG(const std::vector<double> & rms_time,
                         const double              & t_top,
                         const double              & t_bot,
                         const double              & dt_simbox,
                         const double              & max_time,
                         const int                 & n_layers_above,
                         const int                 & n_layers_below,
                         const int                 & n_layers_simbox,
                         const int                 & n_layers_padding) const
{
  int n_layers = n_layers_above + n_layers_padding + n_layers_below;

  std::vector<double> t(n_layers + 1, 0);
  std::vector<double> dt(n_layers + 1, 0);

  double dt_above  = static_cast<double>(  t_top             / n_layers_above);
  double dt_below  = static_cast<double>( (max_time - t_bot) / n_layers_below);

  for(int j=0; j<n_layers + 1; j++) {
    if(j < n_layers_above) {
      t[j]  = j * dt_above;
      dt[j] = dt_above;
    }
    else if(j >= n_layers_above && j < n_layers_above + n_layers_simbox) {
      t[j]  = t_top + (j - n_layers_above) * dt_simbox;
      dt[j] = dt_simbox;
    }
    else if(j >= n_layers_above + n_layers_padding) {
      t[j]  = t_bot + (j - n_layers_above - n_layers_padding) * dt_below;
      dt[j] = dt_below;
    }
  }

  int n_rms_data = static_cast<int>(rms_time.size());

  NRLib::Grid2D<double> G(n_rms_data, n_layers, 0);

  for(int j=0; j<n_rms_data; j++) {
    int k=0;
    while(rms_time[j] >= t[k] && k < n_layers) {
      G(j,k) = dt[k] / rms_time[j];
      k++;
    }
    if(k < n_layers)
      G(j,k) = (rms_time[j] - t[k-1]) / rms_time[j];
  }

  return G;
}

//-----------------------------------------------------------------------------------------//
NRLib::Grid2D<double>
RMSInversion::makeSigma_m(const double              & t_top,
                          const double              & t_bot,
                          const double              & max_time,
                          const Vario               * variogram_above,
                          const Vario               * variogram_below,
                          const double              & var_vp_above,
                          const double              & var_vp_below,
                          const std::vector<double> & cov_grid_log_vp,
                          const int                 & n_layers_above,
                          const int                 & n_layers_below) const
{
  int n_layers_padding = static_cast<int>(cov_grid_log_vp.size());
  int n_layers         = n_layers_above + n_layers_padding + n_layers_below;

  float dt_above  = static_cast<float>( t_top             / n_layers_above);
  float dt_below  = static_cast<float>((max_time - t_bot) / n_layers_below);

  std::vector<float> corrT_above(n_layers_above + 1);
  std::vector<float> corrT_below(n_layers_below + 1);

  for(int j=0; j<=n_layers_above; j++)
    corrT_above[j] = variogram_above->corr(j*dt_above, 0);
  for(int j=0; j<=n_layers_below; j++)
    corrT_below[j] = variogram_below->corr(j*dt_below, 0);

  NRLib::Grid2D<double> Sigma_m(n_layers, n_layers, 0);

  for(int j=0; j<n_layers_above; j++) {
    int count = 0;
    for(int k=j; k<n_layers_above; k++) {
      Sigma_m(j,k) = var_vp_above * corrT_above[count];
      Sigma_m(k,j) = Sigma_m(j,k);
      count ++;
    }
  }

  for(int j=n_layers_above; j<n_layers_above+n_layers_padding; j++) {
    int count = 0;
    for(int k=j; k<n_layers_above+n_layers_padding; k++) {
      double cov_log_vp = cov_grid_log_vp[count];
      Sigma_m(j,k) = cov_log_vp;
      Sigma_m(k,j) = cov_log_vp;
      count ++;
    }
  }

  for(int j=n_layers_above+n_layers_padding; j<n_layers; j++) {
    int count = 0;
    for(int k=j; k<n_layers; k++) {
      Sigma_m(j,k) = var_vp_below * corrT_below[count];
      Sigma_m(k,j) = Sigma_m(j,k);
      count ++;
    }
  }

  return Sigma_m;
}


