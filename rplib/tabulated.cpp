#include "rplib/tabulated.h"
#include "rplib/distributionwithtrend.h"

#include "nrlib/flens/nrlib_flens.hpp"
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/uniform.hpp"
#include "nrlib/random/normal.hpp"

Tabulated::Tabulated(std::vector<const DistributionWithTrend *> elastic,
                     NRLib::Grid2D<double>                      correlation_matrix)
{
  n_variables_ = static_cast<int>(elastic.size());

  elastic_variables_ = elastic;

  normal_ = new NRLib::Normal();

  CalculateSigmaSqrt(correlation_matrix);

}

Tabulated::~Tabulated()
{
  delete normal_;
  delete uniform_;
}

std::vector<double>
Tabulated::GetQuantileValues(const std::vector<double> u, double s1, double s2)
{

  NRLib::Vector ind_normal_samples(n_variables_);
  for(size_t i=0; i<u.size(); i++) {
    double sample = normal_->Quantile(u[i]);
    ind_normal_samples(i) = sample;
  }

  NRLib::Matrix Sigma_sqrt(n_variables_,n_variables_);

  for(int i=0; i<n_variables_; i++) {
    for(int j=0; j<n_variables_; j++)
      Sigma_sqrt(i,j) = Sigma_sqrt_(i,j);
  }

  NRLib::Vector A(n_variables_);
  A = Sigma_sqrt * ind_normal_samples;

  std::vector<double> correlated_samples(n_variables_);

  for(int i=0; i<n_variables_; i++)
    correlated_samples[i] = A(i);

  std::vector<double> correlated_u(n_variables_);
  for(int i=0; i<n_variables_; i++)
    correlated_u[i] = normal_->Cdf(correlated_samples[i]);

  std::vector<double> correlated_elastic_variables(n_variables_);
  for(int i=0; i<n_variables_; i++)
    correlated_elastic_variables[i] = elastic_variables_[i]->GetQuantileValue(correlated_u[i], s1, s2);


  return(correlated_elastic_variables);

}

std::vector<double>
Tabulated::GenerateSample(std::vector<double> & u, double s1, double s2)
{
  u.resize(n_variables_);

  for(int i=0; i<n_variables_; i++)
    u[i] = uniform_->Draw();

  std::vector<double> correlated_elastic_variables(n_variables_);

  correlated_elastic_variables = GetQuantileValues(u, s1, s2);

  return(correlated_elastic_variables);

}

void
Tabulated::CalculateSigmaSqrt(const NRLib::Grid2D<double> & Sigma)
{
  // Calculate square root of positive definite correlation matrix

  NRLib::Matrix corr_matrix(n_variables_,n_variables_);
  for(int i=0; i<n_variables_; i++) {
    for(int j=0; j<n_variables_; j++)
      corr_matrix(i,j) = Sigma(i,j);
  }

  NRLib::Vector eigen_values(n_variables_);
  NRLib::Matrix eigen_vectors(n_variables_,n_variables_);
  NRLib::ComputeEigenVectors(corr_matrix, eigen_values, eigen_vectors);

  NRLib::Matrix eigen_vectors_transpose(n_variables_,n_variables_);
  for(int i=0; i<n_variables_; i++) {
    for(int j=0; j<n_variables_; j++)
      eigen_vectors_transpose(j,i) = eigen_vectors(i,j);
  }

  NRLib::Matrix lambda_square(n_variables_,n_variables_);
  for(int i=0; i<n_variables_; i++) {
    for(int j=0; j<n_variables_; j++) {
      if(i == j)
        lambda_square(i,j) = std::sqrt(eigen_values(i));
      else
        lambda_square(i,j) = 0;
    }
  }

  NRLib::Matrix Sigma_sqrt1(n_variables_,n_variables_);
  NRLib::Matrix Sigma_sqrt(n_variables_,n_variables_);

  Sigma_sqrt1 = eigen_vectors * lambda_square;
  Sigma_sqrt  = Sigma_sqrt1   * eigen_vectors_transpose;

  for(int i=0; i<n_variables_; i++) {
    for(int j=0; j<n_variables_; j++)
      Sigma_sqrt_(i,j) = Sigma_sqrt(i,j);
  }
}
