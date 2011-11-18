#include <src/timeevolution.h>
#include <assert.h>

#include "src/seismicparametersholder.h"
#include "src/fftgrid.h"

TimeEvolution::TimeEvolution(int number_of_timesteps)
: number_of_timesteps_(number_of_timesteps)
{
  SetUpEvolutionMatrix(evolution_matrix_);
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
        //?? Still symmetric values??
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


void TimeEvolution::SetUpEvolutionMatrix(std::vector< NRLib::Matrix> &evolution_matrix)
{
  NRLib::Matrix Ak(3,3);
  Ak = 1,1,1,
       1,1,1,
       1,1,1;
  for (int i = 1; i<=number_of_timesteps_; i++)
  {
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
