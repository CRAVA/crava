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
  SetUpEvolutionMatrices(evolution_matrix_, cov_correction_term_, mean_correction_term_, i_max, time_line, dist_rock, dist_sat, dist_geochem);
}

void TimeEvolution::Evolve(int time_step, SeismicParametersHolder &m_combined){

  // Dummy holders for dynamic and static grids
  SeismicParametersHolder m_dynamic, m_static;
  Split(m_combined, m_dynamic, m_static);

  // Holders of FFTGrid pointers
  std::vector<FFTGrid *> mu(6);
  std::vector<FFTGrid *> sigma(12);

  mu[0] = m_static.GetMuAlpha();
  mu[1] = m_static.GetMuBeta();
  mu[2] = m_static.GetMuRho();
  mu[3] = m_dynamic.GetMuAlpha();
  mu[4] = m_dynamic.GetMuBeta();
  mu[5] = m_dynamic.GetMuRho();

  // Note the order of the grids here: The order is used later in set up of symmetric sigma matrix
  sigma[0] = m_static.GetCovAlpha();
  sigma[1] = m_static.GetCrCovAlphaBeta();
  sigma[2] = m_static.GetCrCovAlphaRho();
  sigma[3] = m_static.GetCovBeta();
  sigma[4] = m_static.GetCrCovBetaRho();
  sigma[5] = m_static.GetCovRho();
  sigma[6] = m_dynamic.GetCovAlpha();
  sigma[7] = m_dynamic.GetCrCovAlphaBeta();
  sigma[8] = m_dynamic.GetCrCovAlphaRho();
  sigma[9] = m_dynamic.GetCovBeta();
  sigma[10] = m_dynamic.GetCrCovBetaRho();
  sigma[11] = m_dynamic.GetCovRho();

  // We assume FFT transformed grids
  for(int i = 0; i<6; i++)
  {
    assert(mu[i]->getIsTransformed());
  }

  for(int i = 0; i<12; i++)
  {
    assert(sigma[i]->getIsTransformed());
  }


  NRLib::Vector mu_real(6);
  NRLib::Vector mu_imag(6);
  NRLib::Vector mu_real_next(6);
  NRLib::Vector mu_imag_next(6);

  NRLib::Matrix sigma_real(6,6);
  NRLib::Matrix sigma_imag(6,6);
  NRLib::Matrix sigma_real_next(6,6);
  NRLib::Matrix sigma_imag_next(6,6);

  fftw_complex return_value;

  int nzp_ = mu[0]->getNzp();
  int nyp_ = mu[0]->getNyp();
  int nxp_ = mu[0]->getNxp();
  int cnxp = nxp_/2+1;


  // Iterate through all points in the grid and perform forward transition in time
  for(int k = 0; k < nzp_; k++)
  {
    for(int j = 0; j < nyp_; j++)
    {
      for(int i = 0; i < cnxp; i++)
      {
        //Set up vectors from the FFT grids
        for(int d = 0; d < 6; d++)
        {
          mu_real(d) = (mu[d]->getNextComplex()).re;
          mu_imag(d) = (mu[d]->getNextComplex()).im;
        }

        // Evolve values
        mu_real_next = evolution_matrix_[time_step]*mu_real + mean_correction_term_[time_step];
        mu_imag_next = evolution_matrix_[time_step]*mu_imag;

        // Update values in the FFT-grids
        for(int d = 0; d < 6; d++)
        {
          return_value.re = static_cast<float>(mu_real_next(d));
          return_value.im = static_cast<float>(mu_imag_next(d));
          mu[d]->setNextComplex(return_value);
        }


        //Set up matrices from the FFT-grids. Covariances between static and dynamic parts not included.
        //Note: Here we assume a specific order of the elements in the sigma-vector.
        int counter = 0;
        for(int d1 = 0; d1 < 3; d1++)
        {
          for(int d2 = d1; d2 < 3; d2++)
          {
            sigma_real(d1, d2) = (sigma[counter]->getNextComplex()).re;
            sigma_imag(d1, d2) = (sigma[counter]->getNextComplex()).im;

            sigma_real(d1+3, d2+3) = (sigma[counter+6]->getNextComplex()).re;  // d1+3 and d2+3 due to block structure of matrix
            sigma_imag(d1+3, d2+3) = (sigma[counter+6]->getNextComplex()).im;

            counter++;

            //Enforcing symmetry
            if(d1 != d2){
              sigma_real(d2, d1) = sigma_real(d1, d2);
              sigma_imag(d2 ,d1) = sigma_real(d1, d2);
            }
          }
        }

        // Evolve values
        sigma_real_next = (evolution_matrix_[time_step]*sigma_real);
        sigma_real_next = sigma_real_next*transpose(evolution_matrix_[time_step]) + cov_correction_term_[time_step];

        sigma_imag_next = evolution_matrix_[time_step]*sigma_imag;
        sigma_imag_next = sigma_imag_next*transpose(evolution_matrix_[time_step]);

        //Update values in the FFT-grids. Covariances between static and dynamic parts not included.
        counter = 0;
        for(int d1 = 0; d1 < 3; d1++)
        {
          for(int d2 = d1; d2 < 3; d2++)
          {
            return_value.re = static_cast<float>(sigma_real(d1, d2));
            return_value.im = static_cast<float>(sigma_imag(d1, d2));
            sigma[counter]->setNextComplex(return_value);

            return_value.re = static_cast<float>(sigma_real(d1+3, d2+3));  //d1+3 and d2+3 due to block structure of matrix
            return_value.im = static_cast<float>(sigma_imag(d1+3, d2+3));
            sigma[counter+6]->setNextComplex(return_value);

            counter++;
          }
        }

      }
    }
  }
  Merge(m_combined, m_dynamic, m_static);
}


void TimeEvolution::SetUpEvolutionMatrices(std::vector< NRLib::Matrix> &evolution_matrix,
                                           std::vector< NRLib::Matrix> &cov_correction_term,
                                           std::vector< NRLib::Vector> &mean_correction_term,
                                           int i_max,
                                           TimeLine & time_line,
                                           const DistributionsRockT0 * dist_rock,
                                           const DistributionsSaturation * dist_sat,
                                           const DistributionsGeochemical * dist_geochem)
{
  // A note on variable names in this function:
  // _mk in the variable names refers to the seismic parameter m with subscript k, time index.
  // _mkm1 in the variable names refers to the seismic parameter m with subscript k-1.

  // A_k:        Denotes the time evolution matrix in the expression m_{k} = A_{k}*m_{k-1} + \delta m_{k}
  // delta_k:    Denotes the covariance of the correction term \delta m_{k}, denoted with symbol \delta_k
  // delta_mu_k: Denotes the mean of the correction term \delta m_{k}, denoted with symbol \delta \mu_{k}

  // E_mk:         Denotes the expected value of m_{k},   E(m_{k})
  // E_mkm1:       Denotes the expected value of m_{k-1}, E(m_{k-1})
  // Cov_mk_mk:     Denotes the covariance of m_{k} and m_{k},     Cov(m_{k}, m_{k})
  // Cov_mk_mkm1:   Denotes the covariance of m_{k} and m_{k-1},   Cov(m_{k}, m_{k-1})
  // Cov_mkm1_mkm1: Denotes the covariance of m_{k-1} and m_{k-1}, Cov(m_{k-1}, m_{k-1})

  CorrelatedRockSamples correlated_rock_samples;
  std::vector<std::vector<std::vector<double> > > m_ik = correlated_rock_samples.CreateSamples(i_max, time_line, dist_rock, dist_sat, dist_geochem);

  int dim = static_cast<int>(m_ik[0][0].size());  //This dimension is expected to be equal to 3, in other words we do not expect to have samples splitted into dynamic and static parts.
  int K = number_of_timesteps_;

  // Data structures for dynamic parts of evolution matrix and correction term
  NRLib::Matrix A_k(dim, dim);
  NRLib::Matrix delta_k(dim, dim);
  NRLib::Vector delta_mu_k(dim);

  NRLib::Vector E_mk(dim);
  NRLib::Vector E_mkm1(dim);
  NRLib::Matrix Cov_mk_mk(dim,dim);
  NRLib::Matrix Cov_mk_mkm1(dim,dim);
  NRLib::Matrix Cov_mkm1_mkm1(dim,dim);

  // For conversion of m_ik to NRLib::Vector
  std::vector<NRLib::Vector> m_k(dim);
  std::vector<NRLib::Vector> m_km1(dim);

  // For all time steps, find A_k, delta_k, delta_mu_k
  // The following part applies only to the dynamic part of the evolution matrix and correction term
  for (int k = 1; k < K; ++k)
  {
    // Data spooling: Get m_{k} and m_{k-1} from samples matrix
    NRLib::Vector temp_m_k(i_max);
    NRLib::Vector temp_m_km1(i_max);
    for (int d = 0; d < dim; d++)
    {
      for (int i = 0; i < i_max; ++i) {
        temp_m_k(i)   = m_ik[k][i][d];
        temp_m_km1(i) = m_ik[k-1][i][d];
      }
      m_k[d]   = temp_m_k;
      m_km1[d] = temp_m_km1;
    }

    // Estimation of mean and covariance elements:
    for (int d1 = 0; d1<dim; d1++)
    {
      E_mk(d1)   = NRLib::Mean(m_k[d1]);
      E_mkm1(d1) = NRLib::Mean(m_km1[d1]);

      for (int d2=0; d2<dim; d2++)
      {
        Cov_mk_mk(d1,d2)     = NRLib::Cov(m_k[d1], m_k[d2]);
        Cov_mk_mkm1(d1,d2)   = NRLib::Cov(m_k[d1], m_km1[d2]);
        Cov_mkm1_mkm1(d1,d2) = NRLib::Cov(m_km1[d1], m_km1[d2]);
      }
    }

    NRLib::Matrix SigmaInv = Cov_mkm1_mkm1;
    NRLib::invert(SigmaInv);
    A_k = Cov_mk_mkm1*SigmaInv;               // A_k = Cov(m_{k}, m_{k-1})\Sigma_{k-1}^{-1}

    NRLib::Matrix A_kT = A_k;
    flens::transpose(A_kT);
    delta_k = Cov_mk_mk - Cov_mk_mkm1*A_kT;   // \Sigma_{k} - Cov(m_{k}m m_{k-1})A_{k}^T
    delta_mu_k = E_mk - A_k*E_mkm1;           // \mu_k - A_k\mu_{k-1}

    // Here we expand the evolution matrices to include the static contributors
    NRLib::Matrix A_k_full(6,6);
    NRLib::Matrix delta_k_full(6,6);
    NRLib::Vector delta_mu_k_full(6);

    // Building full matrix consisting of dynamic and static parts
    for(int d1=0; d1<3; d1++)
    {
      A_k_full(d1,d1)     = 1; // Identity matrix for the static block
      delta_k_full(d1,d1) = 0; // 0 contribution for static block

      delta_mu_k_full(d1) = 0; // 0 contribution for static part

      // Dynamic block
      delta_mu_k_full(d1+3) = delta_mu_k(d1);

      for(int d2=0; d2<3; d2++)
      {
        A_k_full(d1+3, d2+3)     = A_k(d1,d2);
        delta_k_full(d1+3, d2+3) = delta_k(d1,d2);
      }
    }

    // Push back of computed matrices and vectors into return parameters. These matrices and vector are "full", meaning of size 6.
    evolution_matrix.push_back(A_k_full);
    cov_correction_term.push_back(delta_k_full);
    mean_correction_term.push_back(delta_mu_k_full);
  }
}

void TimeEvolution::Split(SeismicParametersHolder &m, SeismicParametersHolder &m_dynamic, SeismicParametersHolder &m_static)
{
  // Dummy implementation
  m_dynamic = m;
  m_static = m;
}

void TimeEvolution::Merge(SeismicParametersHolder &m, SeismicParametersHolder &m_dynamic, SeismicParametersHolder &m_static)
{
  //Dummy implementation
  m = m_dynamic;
}
