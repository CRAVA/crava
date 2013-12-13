#ifndef TABULATED_H
#define TABULATED_H

#include "rplib/distributionwithtrend.h"
#include "nrlib/grid/grid2d.hpp"

#include <vector>

class DistributionWithTrend;

namespace NRLib {
  template <typename T>
  class Distribution;
}

//
// Class that calculates correlated samples of elastic variables
// being for example (Vp, Vs, density) or (bulk-modulus, shear-modulus, density)
// using independent standard uniform variables

// Have A ~ MultiNormal(0, Sigma)
// with f(A) = C * exp(-1/2*A'*Sigma*A) = C * exp(-1/2*B'*I*B)
// where B ~ MultiNormal(0, I)
// with B = Sigma^(-1/2)*A
// if Sigma is positive definite
//
// Using U ~ Uniform(0,1) independent as quantiles in F^(-1),
// B = F^(-1)(u) ~ MultiNormal(0, I)
// We then have A ~ MultiNormal(0, Sigma)
// using A = Sigma^(1/2)*B.
//
// The correlated variables in A can be mapped back to correlated uniform variables
// using U_corr = CDF(A)
//
// The correlated uniform variables are used as quantiles for getting correlated elastic variables
//


class Tabulated {
public:

  Tabulated();

  Tabulated(std::vector<DistributionWithTrend *> elastic,
            NRLib::Grid2D<double>                correlation_matrix);

  ~Tabulated();

  std::vector<double> GenerateSample(std::vector<double> & u, double s1, double s2);
  std::vector<double> GetQuantileValues(const std::vector<double> & u,
                                        double                      s1,
                                        double                      s2);

private:

  void CalculateSigmaSqrt(const NRLib::Grid2D<double> & Sigma);



  std::vector<DistributionWithTrend *>         elastic_variables_;

  int                                          n_variables_;

  NRLib::Grid2D<double>                        Sigma_sqrt_;

  NRLib::Distribution<double>                * normal_;

};

#endif
