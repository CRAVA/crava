// $Id: kriging.cpp 1134 2013-01-02 11:56:53Z anner $

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

#include "kriging.hpp"


using namespace NRLib;

void NRLib::Krig1D(std::vector<double>       &field,
                   const std::vector<bool> &is_known,
                   const std::vector<double> &obs,
                   double                     dx,
                   const Variogram           &vario)
{

  size_t n_obs = obs.size();
  size_t n_unobs = field.size() - n_obs;

  NRLib::Vector residual(static_cast<int> (n_obs));
  int j = 0;
  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true){
      residual(j) = (obs[j]-field[i]);
      j++;
    }

  }


  SymmetricMatrix K;
  NRLib::Vector k_vec(static_cast<int> (n_obs));
  K = SymmetricMatrix(static_cast<int> (n_obs));
  std::vector<double> x_known;
  std::vector<double> x_unknown;
  NRLib::Vector K_res(static_cast<int> (n_obs));

  for(size_t i = 0; i < is_known.size(); i++){
    if(is_known[i] == true)
      x_known.push_back(i*dx);
    else
      x_unknown.push_back(i*dx);
  }


  for(int i = 0; i < static_cast<int> (n_obs); i++)
    for(int j = 0; j <= i; j++)
      K(j,i) = vario.GetCorr(x_known[i]-x_known[j]);

  CholeskyInvert(K);
  K_res = K * residual;
  j = 0;
  size_t kk = 0;
  for(size_t k = 0; k < field.size(); k++){
    if(is_known[k] == false){
    for(int i = 0; i < static_cast<int> (n_obs); i++)
      k_vec(i) = vario.GetCorr(x_known[i]-x_unknown[j]);

    double prod = k_vec * K_res;
    field[k] += prod;
    j++;
    }
    else{
      field[k] = obs[kk];
      kk++;
    }
  }


}

