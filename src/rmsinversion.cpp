/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/rmsinversion.h"
#include "src/modeltraveltimedynamic.h"
#include "src/seismicparametersholder.h"
#include "src/simbox.h"
#include "src/modelgeneral.h"
#include "src/rmstrace.h"

#include "nrlib/flens/nrlib_flens.hpp"

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

  std::vector<double>   mu_log_vp_model = generateMuLogVpModel(mu_log_vp, i_ind, j_ind);

  std::vector<double>   mu_m_square;
  NRLib::Grid2D<double> Sigma_m_square;

  calculateMuSigma_mSquare(mu_log_vp_model,
                           cov_grid_log_vp,
                           max_time,
                           n_layers_above,
                           n_layers_below,
                           t_top,
                           t_bot,
                           mu_vp_top,
                           mu_vp_base,
                           var_vp_above,
                           var_vp_below,
                           variogram_above,
                           variogram_below,
                           mu_m_square,
                           Sigma_m_square);

  std::vector<double>   d_square       = calculateDSquare(rms_trace->getVelocity());
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
}
//-----------------------------------------------------------------------------------------//
void
RMSInversion::calculatePosteriorModel(const std::vector<double>   & d,
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
  for(int i=0; i<n_layers; i++)
    mu_m1(i) = mu_m[i];

  NRLib::Matrix Sigma_m1(n_layers, n_layers);
  for(int i=0; i<n_layers; i++) {
    for(int j=0; j<n_layers; j++)
      Sigma_m1(i,j) = Sigma_m(i,j);
  }

  NRLib::Matrix G1(n_data, n_layers);
  for(int i=0; i<n_data; i++) {
    for(int j=0; j<n_layers; j++)
      G1(i,j) = G(i,j);
  }

  NRLib::Matrix G1_transpose(n_layers, n_data);
  for(int i=0; i<n_data; i++) {
    for(int j=0; j<n_layers; j++)
      G1_transpose(j,i) = G(i,j);
  }

  NRLib::Vector d1(n_data);
  for(int i=0; i<n_data; i++)
    d1(i) = d[i];

  NRLib::Matrix Sigma_d1(n_data, n_data);
  for(int i=0; i<n_data; i++) {
    for(int j=0; j<n_data; j++)
      Sigma_d1(i,j) = Sigma_d(i,j);
  }

  NRLib::Vector dataMean            = G1 * mu_m1;
  NRLib::Vector diff                = d1 - dataMean;
  NRLib::Matrix dataModelCovariance = G1 * Sigma_m1;
  NRLib::Matrix modelDataCovariance = Sigma_m1 * G1_transpose;
  NRLib::Matrix dataCovariance      = dataModelCovariance * G1_transpose + Sigma_d1;

  NRLib::SymmetricMatrix data_covariance_inv_sym(n_data);

  for(int i=0; i<n_data; i++) {
    for(int j=0; j<=i; j++)
      data_covariance_inv_sym(j,i) = dataCovariance(i,j);
  }

  NRLib::CholeskyInvert(data_covariance_inv_sym);

  NRLib::Matrix data_covariance_inv(n_data, n_data);
  for(int i=0; i<n_data; i++) {
    for(int j=i; j<n_data; j++) {
      data_covariance_inv(i,j) = data_covariance_inv_sym(i,j);
      data_covariance_inv(j,i) = data_covariance_inv(i,j);
    }
  }

  NRLib::Matrix helpMat     = modelDataCovariance * data_covariance_inv;
  NRLib::Vector mu_help     = helpMat * diff;
  NRLib::Matrix Sigma_help  = helpMat * dataModelCovariance;

  NRLib::Vector mu_post1    = mu_m1 + mu_help;
  NRLib::Matrix Sigma_post1 = Sigma_m1 - Sigma_help;

  mu_post.resize(n_layers);
  for(int i=0; i<n_layers; i++)
    mu_post[i] = mu_post1(i);

  Sigma_post.Resize(n_layers, n_layers);
  for(int i=0; i<n_layers; i++) {
    for(int j=0; j<n_layers; j++)
      Sigma_post(i,j) = Sigma_post1(i,j);
  }

  /*
  NRLib::WriteMatrixToFile("G", G1);

  NRLib::WriteMatrixToFile("Sigma_m", Sigma_m1);

  NRLib::WriteMatrixToFile("Sigma_d", Sigma_d1);

  NRLib::WriteMatrixToFile("Sigma_post", Sigma_post1);

  NRLib::WriteVectorToFile("mu_post", mu_post1);

  NRLib::WriteVectorToFile("mu_m", mu_m1);

  NRLib::WriteVectorToFile("d", d1);
  */
}

//-----------------------------------------------------------------------------------------//
std::vector<double>
RMSInversion::calculateDSquare(const std::vector<double> & d) const
{
  int n = static_cast<int>(d.size());

  std::vector<double> d_square(n);

  for(int i=0; i<n; i++)
    d_square[i] = std::pow(d[i],2);

  return d_square;
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
RMSInversion::calculateSigmaDSquare(const std::vector<double> & rms_velocity,
                                    const double              & standard_deviation) const
{
  int n                 = static_cast<int>(rms_velocity.size());
  const double variance = std::pow(standard_deviation,2);

  NRLib::Grid2D<double> I(n, n, 0);
  for(int i=0; i<n; i++)
    I(i,i) = 1;

  NRLib::Grid2D<double> Sigma_d_square(n, n, 0);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      Sigma_d_square(i,j) = 4 * std::pow(rms_velocity[i],2) * variance * I(i,j) + 2 * std::pow(variance,2) * I(i,j);
    }
  }

  return Sigma_d_square;
}

//-----------------------------------------------------------------------------------------//
void
RMSInversion::calculateMuSigma_mSquare(const std::vector<double> & mu_log_vp_model,
                                       const std::vector<double> & cov_grid_log_vp,
                                       const double              & max_time,
                                       const int                 & n_layers_above,
                                       const int                 & n_layers_below,
                                       const double              & t_top,
                                       const double              & t_base,
                                       const double              & mu_vp_top,
                                       const double              & mu_vp_base,
                                       const double              & var_vp_above,
                                       const double              & var_vp_below,
                                       const Vario               * variogram_above,
                                       const Vario               * variogram_below,
                                       std::vector<double>       & mu_vp_square,
                                       NRLib::Grid2D<double>     & Sigma_vp_square) const
{

  int n_layers_simbox = static_cast<int>(mu_log_vp_model.size());

  float  dt_above = static_cast<float>( t_top              / n_layers_above);
  float  dt_below = static_cast<float>((max_time - t_base) / n_layers_below);

  // Model
  NRLib::Grid2D<double> Sigma_log_vp_model = generateSigmaModel(cov_grid_log_vp);

  // Above
  std::vector<double>   mu_vp_above    = generateMuVpAbove(mu_vp_top, std::exp(mu_log_vp_model[0]), n_layers_above);
  NRLib::Grid2D<double> Sigma_vp_above = generateSigmaVp(dt_above, n_layers_above, var_vp_above, variogram_above);
  std::vector<double>   mu_log_vp_above;
  NRLib::Grid2D<double> Sigma_log_vp_above;

  calculateCentralMomentLogNormalInverse(mu_vp_above, Sigma_vp_above, mu_log_vp_above, Sigma_log_vp_above);

  // Below
  std::vector<double>   mu_vp_below    = generateMuVpBelow(std::exp(mu_log_vp_model[n_layers_simbox-1]), mu_vp_base, n_layers_below);
  NRLib::Grid2D<double> Sigma_vp_below = generateSigmaVp(dt_below, n_layers_below, var_vp_below, variogram_below);
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
  int n_layers_padding = cov_log_vp->getNzp();

  std::vector<double> cov_grid_log_vp(n_layers_padding, 0);
  for(int j=0; j<n_layers_padding; j++)
    cov_grid_log_vp[j] = cov_log_vp->getRealValue(0, 0, j, true);

  return cov_grid_log_vp;

}

//-----------------------------------------------------------------------------------------//

std::vector<double>
RMSInversion::generateMuLogVpModel(const FFTGrid * mu_log_vp,
                                   const int     & i_ind,
                                   const int     & j_ind) const
{

  int n_layers_padding = mu_log_vp->getNzp();

  std::vector<double> mu_grid_log_vp(n_layers_padding, 0);
  for(int j=0; j<n_layers_padding; j++)
    mu_grid_log_vp[j] = mu_log_vp->getRealValue(i_ind, j_ind, j, true);

  return mu_grid_log_vp;

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
RMSInversion::transformVpToVpSquare(const std::vector<double>   & mu_vp,
                                    const NRLib::Grid2D<double> & Sigma_vp,
                                    std::vector<double>         & mu_vp_square,
                                    NRLib::Grid2D<double>       & Sigma_vp_square) const
{
  std::vector<double>   mu_log_vp;
  NRLib::Grid2D<double> Sigma_log_vp;

  calculateCentralMomentLogNormalInverse(mu_vp, Sigma_vp, mu_log_vp, Sigma_log_vp);

  calculateSecondCentralMomentLogNormal(mu_log_vp, Sigma_log_vp, mu_vp_square, Sigma_vp_square);
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
  for(int i=0; i<n; i++)
    mu_vp_trans[i] = std::exp(mu_log_vp[i] + 0.5 * variance_log_vp(i, i));

  variance_vp_trans.Resize(n, n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      variance_vp_trans(i,j) = mu_vp_trans[i] * mu_vp_trans[j] *(std::exp(variance_log_vp(i,j)) - 1);
  }

  /*
  // Write to file
  NRLib::Vector mu(n);
  for(int i=0; i<n; i++)
    mu(i) = mu_log_vp[i];
  NRLib::WriteVectorToFile("mu_log_m", mu);

  NRLib::Matrix var_mat(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      var_mat(i,j) = variance_log_vp(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_log_m", var_mat);

  NRLib::Matrix Sigma_trans(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      Sigma_trans(i,j) = variance_vp_trans(i,j);
  }
  NRLib::WriteMatrixToFile("Sigma_trans_m", Sigma_trans);

  NRLib::Vector mu_trans(n);
  for(int i=0; i<n; i++)
    mu_trans(i) = mu_vp_trans[i];
  NRLib::WriteVectorToFile("mu_trans_m", mu_trans);
  */
}

//-----------------------------------------------------------------------------------------//

void
RMSInversion::calculateCentralMomentLogNormalInverse(const std::vector<double>   & mu_vp_trans,
                                                     const NRLib::Grid2D<double> & variance_vp_trans,
                                                     std::vector<double>         & mu_log_vp,
                                                     NRLib::Grid2D<double>       & variance_log_vp) const
{
  int n = static_cast<int>(mu_vp_trans.size());

  variance_log_vp.Resize(n, n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++)
      variance_log_vp(i,j) = std::log(1+variance_vp_trans(i,j)/(mu_vp_trans[i] * mu_vp_trans[j]));
  }

  mu_log_vp.resize(n);
  for(int i=0; i<n; i++)
    mu_log_vp[i] = std::log(mu_vp_trans[i]) - 0.5 * variance_log_vp(i,i);

}

//-----------------------------------------------------------------------------------------//
void
RMSInversion::calculateSecondCentralMomentLogNormal(const std::vector<double>   & mu_log_vp,
                                                    const NRLib::Grid2D<double> & variance_log_vp,
                                                    std::vector<double>         & mu_vp_square,
                                                    NRLib::Grid2D<double>       & variance_vp_square) const
{
  int n_layers = static_cast<int>(mu_log_vp.size());

  std::vector<double> mu(n_layers);
  for(int i=0; i<n_layers; i++)
    mu[i] = mu_log_vp[i] * 2;

  NRLib::Grid2D<double> variance(n_layers, n_layers);
  for(int i=0; i<n_layers; i++) {
    for(int j=0; j<n_layers; j++)
      variance(i,j) = variance_log_vp(i,j) * 4;
  }

  mu_vp_square.resize(n_layers);
  variance_vp_square.Resize(n_layers, n_layers);

  calculateCentralMomentLogNormal(mu, variance, mu_vp_square, variance_vp_square);
}
