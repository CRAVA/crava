#include <src/timeevolution.h>
#include <assert.h>

#include "nrlib/statistics/statistics.hpp"

#include "src/seismicparametersholder.h"
#include "src/fftgrid.h"
#include "src/timeline.h"
#include "src/correlatedrocksamples.h"
#include "rplib/distributionsrockt0.h"


TimeEvolution::TimeEvolution(int number_of_timesteps,
                             int i_max,
                             TimeLine & time_line,
                             const DistributionsRockT0 * dist_rock,
                             const DistributionsSaturation * dist_sat,
                             const DistributionsGeochemical * dist_geochem)
                             : number_of_timesteps_(number_of_timesteps)
{
  SetUpEvolutionMatrix(evolution_matrix_, i_max, time_line, dist_rock, dist_sat, dist_geochem);
  SetUpDeltaMuTilde(delta_mu_tilde_);
  SetUpDeltaTilde(delta_tilde_);
}

void TimeEvolution::Evolve(int time_step, SeismicParametersHolder &m){

  FFTGrid * mu_alpha_ = m.GetMuAlpha();
  FFTGrid * mu_beta_  = m.GetMuBeta();
  FFTGrid * mu_rho_   = m.GetMuRho();
  FFTGrid * covAlpha_ = m.GetCovAlpha();
  FFTGrid * covBeta_  = m.GetCovBeta();
  FFTGrid * covRho_   = m.GetCovRho();
  FFTGrid * crCovAlphaBeta_ = m.GetCrCovAlphaBeta();
  FFTGrid * crCovAlphaRho_  = m.GetCrCovAlphaRho();
  FFTGrid * crCovBetaRho_   = m.GetCrCovBetaRho();

  assert(mu_alpha_->getIsTransformed());
  assert(mu_beta_->getIsTransformed());
  assert(mu_rho_->getIsTransformed());
  assert(covAlpha_->getIsTransformed());
  assert(covBeta_->getIsTransformed());
  assert(covRho_->getIsTransformed());
  assert(crCovAlphaBeta_->getIsTransformed());
  assert(crCovAlphaRho_->getIsTransformed());
  assert(crCovBetaRho_->getIsTransformed());

  NRLib::Vector mu_real(3);
  NRLib::Vector mu_imag(3);
  NRLib::Vector mu_real_next(3);
  NRLib::Vector mu_imag_next(3);

  NRLib::Matrix sigma_real(3,3);
  NRLib::Matrix sigma_imag(3,3);
  NRLib::Matrix sigma_real_next(3,3);
  NRLib::Matrix sigma_imag_next(3,3);

  fftw_complex return_value;

  int nzp_ = mu_alpha_->getNzp();
  int nyp_ = mu_alpha_->getNyp();
  int nxp_ = mu_alpha_->getNxp();
  int cnxp = nxp_/2+1;

  for(int k = 0; k < nzp_; k++) {
    for(int j = 0; j < nyp_; j++){
      for(int i = 0; i < cnxp; i++)
      {
        // Get values from the FFT-grids
        mu_real = (mu_alpha_->getNextComplex()).re,
          (mu_beta_->getNextComplex()).re,
          (mu_rho_->getNextComplex()).re;

        mu_imag = (mu_alpha_->getNextComplex()).im,
          (mu_beta_->getNextComplex()).im,
          (mu_rho_->getNextComplex()).im;

        //split

        // Evolve values
        mu_real_next = evolution_matrix_[time_step]*mu_real + delta_mu_tilde_[time_step];
        mu_imag_next = evolution_matrix_[time_step]*mu_imag;

        // Set updated values in the FFT-grids
        return_value.re = static_cast<float>(mu_real_next(0));
        return_value.im = static_cast<float>(mu_imag_next(0));
        mu_alpha_->setNextComplex(return_value);

        return_value.re = static_cast<float>(mu_real_next(1));
        return_value.im = static_cast<float>(mu_imag_next(1));
        mu_beta_->setNextComplex(return_value);

        return_value.re = static_cast<float>(mu_real_next(2));
        return_value.im = static_cast<float>(mu_imag_next(2));
        mu_rho_->setNextComplex(return_value);

        //Get values from the FFT-grids
        sigma_real = (covAlpha_->getNextComplex()).re, (crCovAlphaBeta_->getNextComplex()).re, (crCovAlphaRho_->getNextComplex()).re,
          (crCovAlphaBeta_->getNextComplex()).re, (covBeta_->getNextComplex()).re, (crCovBetaRho_->getNextComplex()).re,
          (crCovAlphaRho_->getNextComplex()).re,(crCovBetaRho_->getNextComplex()).re, (covRho_->getNextComplex()).re;

        sigma_imag = (covAlpha_->getNextComplex()).im, (crCovAlphaBeta_->getNextComplex()).im, (crCovAlphaRho_->getNextComplex()).im,
          (crCovAlphaBeta_->getNextComplex()).im, (covBeta_->getNextComplex()).im, (crCovBetaRho_->getNextComplex()).im,
          (crCovAlphaRho_->getNextComplex()).im,(crCovBetaRho_->getNextComplex()).im, (covRho_->getNextComplex()).im;

        // Evolve values
        sigma_real_next = (evolution_matrix_[time_step]*sigma_real);
        sigma_real_next = sigma_real_next*transpose(evolution_matrix_[time_step]) + delta_tilde_[time_step];

        sigma_imag_next = evolution_matrix_[time_step]*sigma_imag;
        sigma_imag_next = sigma_imag_next*transpose(evolution_matrix_[time_step]);

        // Sjekke symmetri etc

        // Set updated values in the FFT-grids
        return_value.re = static_cast<float>(sigma_real_next(0,0));
        return_value.im = static_cast<float>(sigma_imag_next(0,0));
        covAlpha_->setNextComplex(return_value);

        return_value.re = static_cast<float>(sigma_real_next(0,1));
        return_value.im = static_cast<float>(sigma_imag_next(0,1));
        crCovAlphaBeta_->setNextComplex(return_value);

        return_value.re = static_cast<float>(sigma_real_next(0,2));
        return_value.im = static_cast<float>(sigma_imag_next(0,2));
        crCovAlphaRho_->setNextComplex(return_value);

        return_value.re = static_cast<float>(sigma_real_next(1,1));
        return_value.im = static_cast<float>(sigma_imag_next(1,1));
        covBeta_->setNextComplex(return_value);

        return_value.re = static_cast<float>(sigma_real_next(1,2));
        return_value.im = static_cast<float>(sigma_imag_next(1,2));
        crCovBetaRho_->setNextComplex(return_value);

        return_value.re = static_cast<float>(sigma_real_next(2,2));
        return_value.im = static_cast<float>(sigma_imag_next(2,2));
        covRho_->setNextComplex(return_value);

        // Merge

      }
    }
  }
}


void TimeEvolution::SetUpEvolutionMatrix(std::vector< NRLib::Matrix> &evolution_matrix,
                                         int i_max,
                                         TimeLine & time_line,
                                         const DistributionsRockT0 * dist_rock,
                                         const DistributionsSaturation * dist_sat,
                                         const DistributionsGeochemical * dist_geochem)
{
  CorrelatedRockSamples correlated_rock_samples;
  std::vector<std::vector<std::vector<double> > > m_ik = correlated_rock_samples.CreateSamples(i_max, time_line, dist_rock, dist_sat, dist_geochem);

  NRLib::Matrix Ak(3,3);

  int K = number_of_timesteps_;

  NRLib::Vector logVp_k(i_max);
  NRLib::Vector logVs_k(i_max);
  NRLib::Vector logRho_k(i_max);

  NRLib::Vector logVp_km1(i_max);
  NRLib::Vector logVs_km1(i_max);
  NRLib::Vector logRho_km1(i_max);

  //størrelse fra m_ik
  NRLib::Vector E_mk(3);
  NRLib::Vector E_mkm1(3);
  NRLib::Matrix E_mk_mkm1T(3,3);
  NRLib::Matrix E_mkm1_mkm1T(3,3);

  // for all time steps, find A_k
  for (int k = 1; k < K; ++k) {

    // Get m_{k} and m_{k-1} from samples matrix
    // "_km1" is used for m_{k-1}
    for (int i = 0; i < i_max; ++i) {
      logVp_k(i) = m_ik[k][i][0];
      logVs_k(i) = m_ik[k][i][1];
      logRho_k(i) = m_ik[k][i][2];

      logVp_km1(i) = m_ik[k-1][i][0];
      logVs_km1(i) = m_ik[k-1][i][1];
      logRho_km1(i) = m_ik[k-1][i][2];
    }

    // E(m_k)
    E_mk(0) = NRLib::Mean(logVp_k);
    E_mk(1) = NRLib::Mean(logVs_k);
    E_mk(2) = NRLib::Mean(logRho_k);

    // E(m_{k-1})
    E_mkm1(0) = NRLib::Mean(logVp_km1);
    E_mkm1(1) = NRLib::Mean(logVs_km1);
    E_mkm1(2) = NRLib::Mean(logRho_km1);

    // E(m_k, m_{k-1}^T)
    E_mk_mkm1T(0, 0) = NRLib::Cov(logVp_k, logVp_km1);
    E_mk_mkm1T(1, 1) = NRLib::Cov(logVs_k, logVs_km1);
    E_mk_mkm1T(2, 2) = NRLib::Cov(logRho_k, logRho_km1);
    E_mk_mkm1T(0, 1) = NRLib::Cov(logVp_k, logVs_km1);
    E_mk_mkm1T(0, 2) = NRLib::Cov(logVp_k, logRho_km1);
    E_mk_mkm1T(1, 0) = NRLib::Cov(logVs_k, logVp_km1);
    E_mk_mkm1T(1, 2) = NRLib::Cov(logVs_k, logRho_km1);
    E_mk_mkm1T(2, 0) = NRLib::Cov(logRho_k, logVp_km1);
    E_mk_mkm1T(2, 1) = NRLib::Cov(logRho_k, logVs_km1);

    // E(m_{k-1}, m_{k-1}^T)
    E_mkm1_mkm1T(0, 0) = NRLib::Cov(logVp_km1, logVp_km1);
    E_mkm1_mkm1T(1, 1) = NRLib::Cov(logVs_km1, logVs_km1);
    E_mkm1_mkm1T(2, 2) = NRLib::Cov(logRho_km1, logRho_km1);
    E_mkm1_mkm1T(0, 1) = NRLib::Cov(logVp_km1, logVs_km1);
    E_mkm1_mkm1T(0, 2) = NRLib::Cov(logVp_km1, logRho_km1);
    E_mkm1_mkm1T(1, 0) = NRLib::Cov(logVs_km1, logVp_km1);
    E_mkm1_mkm1T(1, 2) = NRLib::Cov(logVs_km1, logRho_km1);
    E_mkm1_mkm1T(2, 0) = NRLib::Cov(logRho_km1, logVp_km1);
    E_mkm1_mkm1T(2, 1) = NRLib::Cov(logRho_km1, logVs_km1);

    // E(m_k m_{k-1}^T) - E(m_k)*E_{k-1}^T = Cov(m_{k}, m_{k-1})
    Ak(0,0) = E_mk_mkm1T(0,0) -  E_mkm1(0)*E_mkm1(0);
    Ak(0,1) = E_mk_mkm1T(0,1) -  E_mkm1(0)*E_mkm1(1);
    Ak(0,2) = E_mk_mkm1T(0,2) -  E_mkm1(0)*E_mkm1(2);
    Ak(1,0) = E_mk_mkm1T(1,0) -  E_mkm1(1)*E_mkm1(0);
    Ak(1,1) = E_mk_mkm1T(1,1) -  E_mkm1(1)*E_mkm1(1);
    Ak(1,2) = E_mk_mkm1T(1,2) -  E_mkm1(1)*E_mkm1(2);
    Ak(2,0) = E_mk_mkm1T(2,0) -  E_mkm1(2)*E_mkm1(0);
    Ak(2,1) = E_mk_mkm1T(2,1) -  E_mkm1(2)*E_mkm1(1);
    Ak(2,2) = E_mk_mkm1T(2,2) -  E_mkm1(2)*E_mkm1(2);

    // A_k = Cov(m_{k}, m_{k-1})\Sigma_{k-1}^{-1}
    NRLib::Matrix B = E_mkm1_mkm1T;  // \Sigma_{k-1}
    NRLib::invert(B);  // \Sigma_{k-1}^{-1}
    Ak = Ak * B;

    evolution_matrix.push_back(Ak);
  }
}

void TimeEvolution::SetUpDeltaMuTilde(std::vector< NRLib::Vector> &delta_mu_tilde)
{
  NRLib::Vector v(3);
  v = 1,1,1;
  for (int i = 1; i<=number_of_timesteps_; i++)
  {
    delta_mu_tilde.push_back(v);
  }
}

void TimeEvolution::SetUpDeltaTilde(std::vector< NRLib::Matrix> &delta_tilde)
{
  NRLib::Matrix m(3,3);
  m = 1,1,1,
    1,1,1,
    1,1,1;
  for (int i = 1; i<=number_of_timesteps_; i++)
  {
    delta_tilde.push_back(m);
  }
}

void TimeEvolution::Split(int fourier_point)
{

}

void TimeEvolution::Merge(int fourier_point)
{

}
