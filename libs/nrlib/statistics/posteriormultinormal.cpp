// $Id: posteriormultinormal.cpp 1134 2013-01-02 11:56:53Z anner $

// Copyright (c)  2011, Norwegian Computing Center
// All rights reserved.
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// •  Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// •  Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <vector>

#include "posteriormultinormal.hpp"
#include "../math/constants.hpp"


using namespace NRLib;

//
// Calculates conditional mean (exp) and covariance matrix (cov)
//
//

bool NRLib::CondDistrMultiNormal(const std::vector<double> &x_known,
                                 const std::vector<double> &known_values,
                                 const std::vector<double> &x_unknown,
                                 const std::vector<double> &exp_known,
                                 const std::vector<double> &exp_unknown,
                                 const std::vector<double> &std_known,
                                 const std::vector<double> &std_unknown,
                                 const Variogram &vario,
                                 Vector &exp,
                                 SymmetricMatrix &cov)
{

  int n_known = static_cast<int>(x_known.size());
  int n_unknown = static_cast<int>(x_unknown.size());
  assert(known_values.size() == n_known);
  assert(exp_known.size() == n_known);
  assert(exp_unknown.size() == n_unknown);

  exp = Vector(n_unknown);
  cov = SymmetricMatrix(n_unknown);
  if(n_known == 0){ // no observations, return unconditional mean and covariance
    for(int i = 0; i < n_unknown; i++){
      exp(i) = exp_unknown[i];
      for(int j = 0; j < n_unknown; j++){
        cov(i, j) = std_unknown[i]*std_unknown[j]*vario.GetCorr(x_unknown[i]-x_unknown[j]);
      }
    }
      return true;
  }

  NRLib::Vector mu_known(n_known);
  NRLib::Vector mu_unknown(n_unknown);

  for(int i = 0; i < n_known; i++){
    mu_known(i) =  known_values[i] - exp_known[i];
  }
  for(int i = 0; i < n_unknown; i++)
    mu_unknown(i) = exp_unknown[i];

  SymmetricMatrix sigma_11, sigma_22;
  Matrix sigma_12;
  sigma_11 = SymmetricMatrix(n_unknown);
  sigma_22 = SymmetricMatrix(n_known);
  sigma_12 = Matrix(n_unknown, n_known);
  Matrix sigma_21 = Matrix(n_known, n_unknown);
  for(int i = 0; i < n_unknown; i++)
    for(int j = 0; j <= i; j++)
      sigma_11(j,i) = std_unknown[i]*std_unknown[j]*vario.GetCorr(x_unknown[i]-x_unknown[j]);

  for(int i = 0; i < n_known; i++)
    for(int j = 0; j <= i; j++)
      sigma_22(j,i) = std_known[i]*std_known[j]*vario.GetCorr(x_known[i]-x_known[j]);

  for(int i = 0; i < n_known; i++)
    for(int j = 0; j < n_unknown; j++){
      sigma_21(i,j) = std_known[i]*std_unknown[j]*vario.GetCorr(x_known[i]-x_unknown[j]);
      sigma_12(j,i) = std_known[i]*std_unknown[j]*vario.GetCorr(x_known[i]-x_unknown[j]);
    }

  CholeskyInvert(sigma_22);

  NRLib::Vector sigma_22_mu_known = sigma_22*mu_known;
  Matrix sigma_22_sigma_21 = sigma_22*sigma_21;
  exp = mu_unknown + sigma_12*sigma_22_mu_known;
  Matrix help = sigma_12*sigma_22_sigma_21;
  SymmetricMatrix symmetric_help = help.upper();
  cov = sigma_11 - symmetric_help;

  return true;

}

// Calculates posterior distribution in binorma distribution when one of the variables is observed

void NRLib::Posterior1DNormal(double exp_prior, double cov_prior,
                       double exp_obs, double cov_obs,
                       double corr, double obs,
                       double & exp_post, double & cov_post)
{
  exp_post = exp_prior + (corr*sqrt(cov_prior/cov_obs))*(obs - exp_obs);
  cov_post = (1-corr*corr)*cov_prior;

  return;
}

double NRLib::PotentialMultiNormal(const Vector          & x,
                                   const Vector          & mu,
                                   const SymmetricMatrix & sigma)
{
  Matrix e_mat;
  Vector e_vec;
  ComputeEigenVectorsSymmetric(sigma, e_vec, e_mat);

  double det(1.0);
  int dim = sigma.dim();
  for (int i = 0; i < dim; i++)
    det *= e_vec(i);
  assert(det > 0);
  Vector diff = x - mu;
  Vector e_vec_inv = e_vec;
  for(int i = 0; i < dim; i++)
    e_vec_inv(i) = 1.0/e_vec(i);
  double sum(0.0);
  for (int i = 0; i < dim; i++) {
    //double r1 = x(i) - mu(i);
    double r1 = diff(i);
    for (int j = 0; j < dim; j++) {
      double q1 = e_mat(i,j);
      //double s = 1.0 / e_vec(j);
      double r1q1s = r1*q1*e_vec_inv(j);
      for (int k = 0; k < dim; k++) {
        //double r2 = x(k) - mu(k);
        double r2 = diff(k);
        double q2 = e_mat(k,j);
        //sum += q1 * q2 * r1 * r2 * s;
        sum += r1q1s*r2*q2;
      }
    }
  }

  double potential = 0.5 * (dim * NRLib::Log2Pi + std::log(det) + sum);

  return(potential);
}
