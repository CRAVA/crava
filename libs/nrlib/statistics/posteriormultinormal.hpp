// $Id: posteriormultinormal.hpp 1092 2012-10-09 14:24:00Z georgsen $

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

#ifndef NRLIB_STATISTICS_POSTERIORMULTINORMAL_HPP
#define NRLIB_STATISTICS_POSTERIORMULTINORMAL_HPP

#include <cstdlib>
#include <vector>


#include "../flens/nrlib_flens.hpp"
#include "../variogram/variogram.hpp"

namespace NRLib {


  bool CondDistrMultiNormal(const std::vector<double> &x_known,
                            const std::vector<double> &known_values,
                            const std::vector<double> &x_unknown,
                            const std::vector<double> &exp_known,
                            const std::vector<double> &exp_unknown,
                            const std::vector<double> &std_known,
                            const std::vector<double> &std_unknown,
                            const Variogram &vario,
                            Vector &exp,
                            SymmetricMatrix &cov);

  void Posterior1DNormal(double exp_prior, double cov_prior,
                       double exp_obs, double cov_obs,
                       double corr, double obs,
                       double & exp_post, double & cov_post);


  double PotentialMultiNormal(const Vector          & x,
                              const Vector          & mu,
                              const SymmetricMatrix & sigma);


}




#endif
