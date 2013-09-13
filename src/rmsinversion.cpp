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

  FFTGrid * mu_log_vp  = seismicParameters.GetMuAlpha();
  FFTGrid * cov_log_vp = seismicParameters.GetCovBeta();

  Simbox * timeSimbox = modelGeneral->getTimeSimbox();

  const std::vector<RMSTrace *> rms_traces = modelTravelTimeDynamic->getRMSTraces();

  const int n_layers_above        = modelTravelTimeDynamic->getNLayersAbove();
  const int n_layers_below        = modelTravelTimeDynamic->getNLayersBelow();
  const double mu_vp_top          = modelTravelTimeDynamic->getMeanVpTop();
  const double mu_vp_base         = modelTravelTimeDynamic->getMeanVpBase();
  const double var_vp_above       = modelTravelTimeDynamic->getVarVpAbove();
  const double var_vp_below       = modelTravelTimeDynamic->getVarVpBelow();
  const double range_above        = modelTravelTimeDynamic->getRangeAbove();
  const double range_below        = modelTravelTimeDynamic->getRangeBelow();
  const double standard_deviation = modelTravelTimeDynamic->getStandardDeviation();

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
                  mu_vp_top,
                  mu_vp_base,
                  var_vp_above,
                  var_vp_below,
                  max_time,
                  standard_deviation,
                  rms_traces[i],
                  mu_log_vp,
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
                            const double                & mu_vp_top,
                            const double                & mu_vp_base,
                            const double                & var_vp_above,
                            const double                & var_vp_below,
                            const double                & max_time,
                            const double                & standard_deviation,
                            const RMSTrace              * rms_trace,
                            const FFTGrid               * mu_log_vp,
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

  NRLib::Grid2D<double> G = calculateG(time,
                                       t_top,
                                       t_bot,
                                       dt_simbox,
                                       max_time,
                                       n_layers_above,
                                       n_layers_below,
                                       n_layers_simbox,
                                       n_layers_padding);


  NRLib::Grid2D<double> Sigma_d_sqrt = calculateSigmaDSqrt(rms_trace->getVelocity(), standard_deviation);

  //// Calculate mu_m and Sigma_m
  float  dt_above = static_cast<float>( t_top             / n_layers_above);
  float  dt_below = static_cast<float>((max_time - t_bot) / n_layers_below);

  // Model
  std::vector<double>   mu_log_vp_model    = generateMuLogVpModel(mu_log_vp, i_ind, j_ind); //Her kan RMISSING legges til ved behov
  NRLib::Grid2D<double> Sigma_log_vp_model = generateSigmaModel(cov_grid_log_vp); // Har fjernet paddingen fra denne
  std::vector<double>   mu_vp_sqrt_model;
  NRLib::Grid2D<double> Sigma_vp_sqrt_model;

  calculateSecondCentralMomentLogNormal(mu_log_vp_model, Sigma_log_vp_model, mu_vp_sqrt_model, Sigma_vp_sqrt_model);

  // Above
  std::vector<double>   mu_vp_above    = generateMuVpAbove(mu_vp_top, std::exp(mu_log_vp_model[0]), n_layers_above);
  NRLib::Grid2D<double> Sigma_vp_above = generateSigmaVp(dt_above, n_layers_above, var_vp_above, variogram_above);
  std::vector<double>   mu_vp_sqrt_above;
  NRLib::Grid2D<double> Sigma_vp_sqrt_above;

  transformVpToVpSqrt(mu_vp_above, Sigma_vp_above, mu_vp_sqrt_above, Sigma_vp_sqrt_above);

  // Below
  std::vector<double>   mu_vp_below    = generateMuVpBelow(std::exp(mu_log_vp_model[n_layers_simbox-1]), mu_vp_base, n_layers_below);
  NRLib::Grid2D<double> Sigma_vp_below = generateSigmaVp(dt_below, n_layers_below, var_vp_below, variogram_below);
  std::vector<double>   mu_vp_sqrt_below;
  NRLib::Grid2D<double> Sigma_vp_sqrt_below;

  transformVpToVpSqrt(mu_vp_below, Sigma_vp_below, mu_vp_sqrt_below, Sigma_vp_sqrt_below);

  // Combine
  std::vector<double>   mu_m    = generateMuCombined(mu_vp_sqrt_above, mu_vp_sqrt_model, mu_vp_sqrt_below);
  NRLib::Grid2D<double> sigma_m = generateSigmaCombined(Sigma_vp_sqrt_above, Sigma_vp_sqrt_model, Sigma_vp_sqrt_below);

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
RMSInversion::calculateSigmaDSqrt(const std::vector<double> & rms_velocity,
                                  const double              & standard_deviation) const
{
  int n                 = static_cast<int>(rms_velocity.size());
  const double variance = std::sqrt(standard_deviation);

  NRLib::Grid2D<double> I(n, n, 0);

  for(int i=0; i<n; i++)
    I(i,i) = 1;

  NRLib::Grid2D<double> Sigma_d_sqrt(n, n, 0);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      Sigma_d_sqrt(i,j) = 4 * std::sqrt(rms_velocity[i]) * variance * I(i,j) + 2* std::sqrt(variance) * I(i,j);
    }
  }

  return Sigma_d_sqrt;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
RMSInversion::generateSigmaVp(const float  & dt,
                              const int    & n_layers,
                              const double & var_vp,
                              const Vario  * variogram) const
{

  std::vector<float> corrT(n_layers);
  for(int j=0; j<n_layers; j++)
    corrT[j] = variogram->corr(j*dt, 0);

  NRLib::Grid2D<double> Sigma_vp = generateSigma(var_vp, corrT);

  return Sigma_vp;
}

//-----------------------------------------------------------------------------------------//
std::vector<double>
RMSInversion::getCovLogVp(const FFTGrid * cov_log_vp) const
{
  int n_layers_padding = cov_log_vp->getNz(); //cov_log_vp->getNzp();

  std::vector<double> cov_grid_log_vp(n_layers_padding, 0);
  for(int j=0; j<n_layers_padding; j++)
    cov_grid_log_vp[j] = cov_log_vp->getRealValue(0, 0, j);

  return cov_grid_log_vp;

}

//-----------------------------------------------------------------------------------------//

std::vector<double>
RMSInversion::generateMuLogVpModel(const FFTGrid * mu_log_vp,
                                   const int     & i_ind,
                                   const int     & j_ind) const
{
  std::vector<float> mu_log_vp_model_float = mu_log_vp->getRealTrace2(i_ind, j_ind);

  int n = static_cast<int>(mu_log_vp_model_float.size());

  std::vector<double> mu_double(n);

  for(int i=0; i<n; i++)
    mu_double[i] = static_cast<double>(mu_log_vp_model_float[i]);

  return mu_double;

}

//-----------------------------------------------------------------------------------------//

std::vector<double>
RMSInversion::generateMuVpAbove(const double & top_value,
                                const double & base_value,
                                const int    & n_layers) const
{
  std::vector<double> mu_vp = generateMuVp(top_value, base_value, n_layers);

  mu_vp.resize(n_layers);

  return mu_vp;
}

//-----------------------------------------------------------------------------------------//

std::vector<double>
RMSInversion::generateMuVpBelow(const double & top_value,
                                const double & base_value,
                                const int    & n_layers) const
{
  std::vector<double> mu_vp = generateMuVp(top_value, base_value, n_layers);

  for(int i=0; i<n_layers; i++)
    mu_vp[i] = mu_vp[i+1];

  mu_vp.resize(n_layers);

  return mu_vp;
}

//-----------------------------------------------------------------------------------------//

std::vector<double>
RMSInversion::generateMuVp(const double & top_value,
                           const double & base_value,
                           const int    & n_layers) const
{
  std::vector<double> mu_vp(n_layers+1);

  for(int j=0; j<n_layers+1; j++)
    mu_vp[j] = top_value + j * (base_value - top_value) / n_layers;

  return mu_vp;
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

std::vector<double>
RMSInversion::generateMuCombined(const std::vector<double> & mu_above,
                                 const std::vector<double> & mu_model,
                                 const std::vector<double> & mu_below) const
{
  int n_layers_above        = static_cast<int>(mu_above.size());
  int n_layers_below        = static_cast<int>(mu_below.size());
  int n_layers_model        = static_cast<int>(mu_model.size());
  int n_layers              = n_layers_above + n_layers_model + n_layers_below;

  std::vector<double> mu_m(n_layers, RMISSING);

  for(int j=0; j<n_layers_above; j++)
    mu_m[j] = mu_above[j];

  for(int j=n_layers_above; j<n_layers_above+n_layers_model; j++)
    mu_m[j] = mu_model[j - n_layers_above];

  for(int j=n_layers_above+n_layers_model; j<n_layers; j++)
    mu_m[j] = mu_below[j-n_layers_above-n_layers_model];

  return mu_m;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
RMSInversion::generateSigmaCombined(const NRLib::Grid2D<double> & Sigma_above,
                                    const NRLib::Grid2D<double> & Sigma_model,
                                    const NRLib::Grid2D<double> & Sigma_below) const
{
  int n_layers_above = static_cast<int>(Sigma_above.GetNI());
  int n_layers_below = static_cast<int>(Sigma_below.GetNI());
  int n_layers_model = static_cast<int>(Sigma_model.GetNI());
  int n_layers       = n_layers_above + n_layers_model + n_layers_below;

  NRLib::Grid2D<double> Sigma_m(n_layers, n_layers, 0);

  for(int j=0; j<n_layers_above; j++) {
    for(int k=j; k<n_layers_above; k++) {
      Sigma_m(j,k) = Sigma_above(j, k);
      Sigma_m(k,j) = Sigma_m(j,k);
    }
  }

  for(int j=n_layers_above; j<n_layers_above+n_layers_model; j++) {
    for(int k=j; k<n_layers_above+n_layers_model; k++) {
      Sigma_m(j,k) = Sigma_model(j-n_layers_above, k-n_layers_above);
      Sigma_m(k,j) = Sigma_m(j,k);
    }
  }

  for(int j=n_layers_above+n_layers_model; j<n_layers; j++) {
    for(int k=j; k<n_layers; k++) {
      Sigma_m(j,k) = Sigma_below(j-n_layers_above-n_layers_model, k-n_layers_above-n_layers_model);
      Sigma_m(k,j) = Sigma_m(j,k);
    }
  }

  return Sigma_m;
}

//-----------------------------------------------------------------------------------------//

NRLib::Grid2D<double>
RMSInversion::generateSigma(const double             & var,
                            const std::vector<float> & corrT) const
{

  int n_layers = static_cast<int>(corrT.size());

  NRLib::Grid2D<double> Sigma(n_layers, n_layers, 0);

  for(int j=0; j<n_layers; j++) {
    int count = 0;
    for(int k=j; k<n_layers; k++) {
      Sigma(j,k) = var * corrT[count];
      Sigma(k,j) = Sigma(j,k);
      count ++;
    }
  }

  return Sigma;
}

//-----------------------------------------------------------------------------------------//
NRLib::Grid2D<double>
RMSInversion::generateSigmaModel(const std::vector<double> & cov_grid) const
{
  int n = static_cast<int>(cov_grid.size());

  NRLib::Grid2D<double> Sigma_m(n, n);

  for(int j=0; j<n; j++) {
    int count = 0;
    for(int k=j; k<n; k++) {
      double cov_log_vp = cov_grid[count];
      Sigma_m(j,k) = cov_log_vp;
      Sigma_m(k,j) = Sigma_m(j,k);
      count ++;
    }
  }

  return Sigma_m;
}
//-----------------------------------------------------------------------------------------//

void
RMSInversion::transformVpToVpSqrt(const std::vector<double>   & mu_vp,
                                  const NRLib::Grid2D<double> & Sigma_vp,
                                  std::vector<double>         & mu_vp_sqrt,
                                  NRLib::Grid2D<double>       & Sigma_vp_sqrt) const
{
  std::vector<double>   mu_log_vp;
  NRLib::Grid2D<double> Sigma_log_vp;

  calculateCentralMomentLogNormalInverse(mu_vp, Sigma_vp, mu_log_vp, Sigma_log_vp);

  calculateSecondCentralMomentLogNormal(mu_log_vp, Sigma_log_vp, mu_vp_sqrt, Sigma_vp_sqrt);
}

//-----------------------------------------------------------------------------------------//

void
RMSInversion::calculateCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                              const NRLib::Grid2D<double> & variance_log_vp,
                                              std::vector<double>         & mu_vp_trans,
                                              NRLib::Grid2D<double>       & variance_vp_trans) const
{
  int n = static_cast<int>(mu_log_vp.size());

  mu_vp_trans.resize(n);
  variance_vp_trans.Resize(n, n, 0);

  NRLib::Grid2D<double> I(n, n, 0);
  for(int i=0; i<n; i++)
    I(i,i) = 1;


  for(int i=0; i<n; i++)
    mu_vp_trans[i] = std::exp(mu_log_vp[i] + 0.5 * variance_log_vp(i, i));

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      variance_vp_trans(i,j) = mu_vp_trans[i] * mu_vp_trans[j] *(std::exp(variance_log_vp(i,j)) - I(i,j));
  }

}

//-----------------------------------------------------------------------------------------//

void
RMSInversion::calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                     const NRLib::Grid2D<double> & variance_vp_trans,
                                                     std::vector<double>         & mu_log_vp,
                                                     NRLib::Grid2D<double>       & variance_log_vp) const
{
  int n = static_cast<int>(mu_vp_trans.size());

  mu_log_vp.resize(n);
  variance_log_vp.Resize(n, n, 0);

  NRLib::Grid2D<double> I(n, n, 0);
  for(int i=0; i<n; i++)
    I(i,i) = 1;

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      variance_log_vp(i,j) = std::log(I(i,j) + variance_vp_trans(i,j) / (mu_vp_trans[i] * mu_vp_trans[j]));
  }

  for(int i=0; i<n; i++)
    mu_log_vp[i] = std::log(mu_vp_trans[i]) - 0.5 * variance_log_vp(i,i);

}

//-----------------------------------------------------------------------------------------//
void
RMSInversion::calculateSecondCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                    const NRLib::Grid2D<double> & variance_log_vp,
                                                    std::vector<double>         & mu_vp_sqrt,
                                                    NRLib::Grid2D<double>       & variance_vp_sqrt) const
{
  int n_layers = static_cast<int>(mu_log_vp.size());

  mu_vp_sqrt.resize(n_layers);
  variance_vp_sqrt.Resize(n_layers, n_layers);

  std::vector<double>   mu(n_layers);
  NRLib::Grid2D<double> variance(n_layers, n_layers);

  for(int i=0; i<n_layers; i++)
    mu[i] = mu_log_vp[i] * 2;

  for(int i=0; i<n_layers; i++) {
    for(int j=0; j<n_layers; j++)
      variance(i,j) = variance_log_vp(i,j) * 4;
  }

  calculateCentralMomentLogNormal(mu, variance, mu_vp_sqrt, variance_vp_sqrt);
}
