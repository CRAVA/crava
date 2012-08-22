#include "rplib/distributionsrocktrinormal.h"
#include "rplib/multinormalrock.h"
#include "rplib/trinormalwith2dtrend.h"
#include "rplib/pdf3dgaussian.h"
#include "rplib/pdf3d.h"

#include "nrlib/grid/grid2d.hpp"


DistributionsRockTriNormal::DistributionsRockTriNormal(NRLib::Trend * mean_vp,
                                                       NRLib::Trend * mean_vs,
                                                       NRLib::Trend * mean_density,
                                                       NRLib::Trend * variance_vp,
                                                       NRLib::Trend * variance_vs,
                                                       NRLib::Trend * variance_density,
                                                       NRLib::Trend * correlation_vp_vs,
                                                       NRLib::Trend * correlation_vp_density,
                                                       NRLib::Trend * correlation_vs_density)
{
  mult_normal_distr_ = NULL;

  NRLib::Trend * covariance_vp_vs       = NULL;
  NRLib::Trend * covariance_vp_density  = NULL;
  NRLib::Trend * covariance_vs_density  = NULL;

  CalculateCovarianceFromCorrelation(correlation_vp_vs,
                                     variance_vp,
                                     variance_vs,
                                     covariance_vp_vs);
  CalculateCovarianceFromCorrelation(correlation_vp_density,
                                     variance_vp,
                                     variance_density,
                                     covariance_vp_density);
  CalculateCovarianceFromCorrelation(correlation_vs_density,
                                     variance_vs,
                                     variance_density,
                                     covariance_vs_density);

  NRLib::Grid2D<NRLib::Trend *> cov_mat(3,3,NULL);
  cov_mat(0,0) = variance_vp;
  cov_mat(1,1) = variance_vs;
  cov_mat(2,2) = variance_density;
  cov_mat(0,1) = covariance_vp_vs;
  cov_mat(0,2) = covariance_vp_density;
  cov_mat(1,2) = covariance_vs_density;

  std::vector<NRLib::Trend * > mean_vec(3);
  mean_vec[0] = mean_vp;
  mean_vec[1] = mean_vs;
  mean_vec[2] = mean_density;

  std::vector<NRLib::Trend *> log_mean(3);
  for(int i=0; i<3; i++)
    log_mean[i] = NULL;

  NRLib::Grid2D<NRLib::Trend *> log_cov(3,3,NULL);
  bool                          diagonal_element;

  for(int i=0; i<3; i++) {
    for(int j=i; j<3; j++) {
      if(i == j)
        diagonal_element = true;
      else
        diagonal_element = false;

      LogTransformExpectationAndCovariance(mean_vec[i],
                                           mean_vec[j],
                                           cov_mat(i,j),
                                           log_mean[i],
                                           log_cov(i,j),
                                           diagonal_element);
    }
  }

  for(int i=1; i<3; i++) {
    for(int j=0; j<i; j++)
      log_cov(i,j) = log_cov(j,i);
  }

  mult_normal_distr_ = new TriNormalWith2DTrend(log_mean[0],
                                                log_mean[1],
                                                log_mean[2],
                                                log_cov);

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(mean_vp               ->GetTrendDimension(), mean_vp               ->GetReference());
  FindUseTrendCube(mean_vs               ->GetTrendDimension(), mean_vs               ->GetReference());
  FindUseTrendCube(mean_density          ->GetTrendDimension(), mean_density          ->GetReference());
  FindUseTrendCube(variance_vp           ->GetTrendDimension(), variance_vp           ->GetReference());
  FindUseTrendCube(variance_vs           ->GetTrendDimension(), variance_vs           ->GetReference());
  FindUseTrendCube(variance_density      ->GetTrendDimension(), variance_density      ->GetReference());
  FindUseTrendCube(correlation_vp_vs     ->GetTrendDimension(), correlation_vp_vs     ->GetReference());
  FindUseTrendCube(correlation_vp_density->GetTrendDimension(), correlation_vp_density->GetReference());
  FindUseTrendCube(correlation_vs_density->GetTrendDimension(), correlation_vs_density->GetReference());

  delete covariance_vp_vs;
  delete covariance_vp_density;
  delete covariance_vs_density;

  for(int i=0; i<3; i++) {
    delete log_mean[i];
    for(int j=i; j<3; j++) {
      delete log_cov(i,j);
    }
  }
}

DistributionsRockTriNormal::~DistributionsRockTriNormal()
{
  delete mult_normal_distr_;
}

Rock*
DistributionsRockTriNormal::GenerateSample(const std::vector<double> & trend_params) const {
  std::vector<double>  param(3, 0);
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }

  NRLib::Normal vp0, vs0, density0;
  mult_normal_distr_->ReSample(vp0, vs0, density0, s1, s2, param[0], param[1], param[2]);
  return new MultiNormalRock(param);
}

std::vector<double>
DistributionsRockTriNormal::GetExpectation(const std::vector<double> & trend_params) const
{
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  std::vector<double> expectation(3);
  mult_normal_distr_->GetExpectation(s1, s2, expectation);
  return(expectation);
}

NRLib::Grid2D<double>
DistributionsRockTriNormal::GetCovariance(const std::vector<double> & trend_params) const
{
  double s1 = 0, s2 = 0;
  if (trend_params.size() >= 2) {
    s1 = trend_params[0];
    s2 = trend_params[1];
  }
  NRLib::Grid2D<double> covariance(3,3,0);
  mult_normal_distr_->GetCovariance(s1, s2, covariance);
  return(covariance);
}

Pdf3D*
DistributionsRockTriNormal::GeneratePdf(void) const
{

  return new Pdf3DGaussian(mult_normal_distr_);

}

bool
DistributionsRockTriNormal::HasDistribution() const
{
  bool has_distribution = false; //Always false as no variable is allowed to have a distribution
  return(has_distribution);
}

std::vector<bool>
DistributionsRockTriNormal::HasTrend() const
{
  return(use_trend_cube_);
}

void
DistributionsRockTriNormal::FindUseTrendCube(int dim, int reference)
{
  if(dim == 1) {
    use_trend_cube_[reference-1] = true;
  }
  else if(dim == 2) {
    use_trend_cube_[0] = true;
    use_trend_cube_[1] = true;
  }
}

void
DistributionsRockTriNormal::LogTransformExpectationAndCovariance(NRLib::Trend *  mean1,
                                                                 NRLib::Trend *  mean2,
                                                                 NRLib::Trend *  cov,
                                                                 NRLib::Trend *& log_mean,
                                                                 NRLib::Trend *& log_cov,
                                                                 bool          & diagonal_element) const
{
  std::vector<NRLib::Trend *> trender(3);
  trender[0] = mean1;
  trender[1] = mean2;
  trender[2] = cov;

  int new_dim = FindNewGridDimension(trender);

  std::vector<int>    size(2);
  std::vector<double> increment(2);

  FindNewGridSizeAndIncrement(size, increment, trender, new_dim);

  if(new_dim == 0) {
    double log_mu;
    double log_covariance;
    int    dummy = 0;

    LogTransformCovariance(trender[0]->GetTrendElement(dummy,dummy,dummy),
                           trender[1]->GetTrendElement(dummy,dummy,dummy),
                           trender[2]->GetTrendElement(dummy,dummy,dummy),
                           log_covariance);

    if(diagonal_element == true) {
      LogTransformExpectation(trender[0]->GetTrendElement(dummy,dummy,dummy),
                              log_covariance,
                              log_mu);

      log_mean = new NRLib::TrendConstant(log_mu);
    }
    log_cov = new NRLib::TrendConstant(log_covariance);
  }

  else if(new_dim == 1) {
    std::vector<int> dim(3);
    dim[0] = trender[0]->GetTrendDimension();
    dim[1] = trender[1]->GetTrendDimension();
    dim[2] = trender[2]->GetTrendDimension();

    int reference = 0;
    for(int i=0; i<3; i++) {
      if(dim[i] > 0) {
        reference = trender[i]->GetReference();
        break;
      }
    }

    std::vector<std::vector<double> > expanded_trend;
    std::vector<double>               log_mu(size[0]);
    std::vector<double>               log_covariance(size[0]);

    expanded_trend = ExpandGrids1D(trender,size);

    for(int i=0; i<size[0]; i++) {
      LogTransformCovariance(expanded_trend[0][i],
                             expanded_trend[1][i],
                             expanded_trend[2][i],
                             log_covariance[i]);

      if(diagonal_element == true)
        LogTransformExpectation(expanded_trend[0][i],
                                log_covariance[i],
                                log_mu[i]);
    }

    if(diagonal_element == true)
      log_mean = new NRLib::Trend1D(log_mu, reference, increment[0]);

    log_cov = new NRLib::Trend1D(log_covariance, reference, increment[0]);
  }

  else {
    std::vector<NRLib::Grid2D<double> > expanded_trend;
    NRLib::Grid2D<double>               log_mu(size[0],size[1],0);
    NRLib::Grid2D<double>               log_covariance(size[0],size[1],0);

    expanded_trend = ExpandGrids2D(trender,size);

    for(int i=0; i<size[0]; i++) {
      for(int j=0; j<size[1]; j++) {
        LogTransformCovariance(expanded_trend[0](i,j),
                               expanded_trend[1](i,j),
                               expanded_trend[2](i,j),
                               log_covariance(i,j));

        if(diagonal_element == true)
          LogTransformExpectation(expanded_trend[0](i,j),
                                  log_covariance(i,j),
                                  log_mu(i,j));
      }
    }

    if(diagonal_element == true)
      log_mean = new NRLib::Trend2D(log_mu, 1, 2, increment[0], increment[1]);

    log_cov = new NRLib::Trend2D(log_covariance, 1, 2, increment[0], increment[1]);
  }
}

void
DistributionsRockTriNormal::CalculateCovarianceFromCorrelation(NRLib::Trend *  corr,
                                                               NRLib::Trend *  var1,
                                                               NRLib::Trend *  var2,
                                                               NRLib::Trend *& cov) const
{
  std::vector<NRLib::Trend *> trender(3);
  trender[0] = corr;
  trender[1] = var1;
  trender[2] = var2;


  int new_dim = FindNewGridDimension(trender);

  std::vector<int>    size(2);
  std::vector<double> increment(2);

  FindNewGridSizeAndIncrement(size, increment, trender, new_dim);

  if(new_dim == 0) {
    double covariance;
    int    dummy = 0;

    CalculateCovariance(trender[0]->GetTrendElement(dummy,dummy,dummy),
                        trender[1]->GetTrendElement(dummy,dummy,dummy),
                        trender[2]->GetTrendElement(dummy,dummy,dummy),
                        covariance);

    cov = new NRLib::TrendConstant(covariance);
  }

  else if(new_dim == 1) {
    std::vector<double>               covariance(size[0]);
    std::vector<std::vector<double> > expanded_trend;

    expanded_trend = ExpandGrids1D(trender,size);

    for(int i=0; i<size[0]; i++)
      CalculateCovariance(expanded_trend[0][i],
                          expanded_trend[1][i],
                          expanded_trend[2][i],
                          covariance[i]);

    std::vector<int> dim(3);
    dim[0] = trender[0]->GetTrendDimension();
    dim[1] = trender[1]->GetTrendDimension();
    dim[2] = trender[2]->GetTrendDimension();

    int reference = 0;
    for(int i=0; i<3; i++)
      if(dim[i] > 0) {
        reference = trender[i]->GetReference();
        break;
      }

    cov = new NRLib::Trend1D(covariance, reference, increment[0]);
  }

  else {
    NRLib::Grid2D<double>               covariance(size[0],size[1],0);
    std::vector<NRLib::Grid2D<double> > expanded_trend;

    expanded_trend = ExpandGrids2D(trender,size);

    for(int i=0; i<size[0]; i++) {
      for(int j=0; j<size[1]; j++) {
        CalculateCovariance(expanded_trend[0](i,j),
                            expanded_trend[1](i,j),
                            expanded_trend[2](i,j),
                            covariance(i,j));
      }
    }

    cov = new NRLib::Trend2D(covariance, 1, 2, increment[0], increment[1]);
  }
}

int
DistributionsRockTriNormal::FindNewGridDimension(const std::vector<NRLib::Trend *> trender) const
{
  int new_dim;

  std::vector<int> dim;
  for(int i=0; i<static_cast<int>(trender.size()); i++)
    dim.push_back(trender[i]->GetTrendDimension());

  int number = 0;
  for(int i=0; i<static_cast<int>(trender.size()); i++) {
    if(dim[i] > number)
      number = dim[i];
  }

  if(number == 0)
    new_dim = 0;

  else if(number == 2)
    new_dim = 2;

  else {
    new_dim = 1;

    // new_dim = 2 if the 1D-trends have different reference-variable
    std::vector<int> reference;

    for(int i=0; i<static_cast<int>(trender.size()); i++) {
      if(dim[i] > 0)
        reference.push_back(trender[i]->GetReference());
    }

    if(reference.size() > 1) {
      int compare = reference[0];

      for(int i=1; i<static_cast<int>(reference.size()); i++) {
        if(reference[i] != compare)
          new_dim = 2;
      }
    }
  }
  return(new_dim);
}

void
DistributionsRockTriNormal::FindNewGridSizeAndIncrement(std::vector<int>                  & size,
                                                        std::vector<double>               & increment,
                                                        const std::vector<NRLib::Trend *> & trender,
                                                        const int                         & new_dim) const
{
  int n_trends = static_cast<int>(trender.size());

  // Find dimension of each of the trends in trender
  std::vector<int> dim;
  for(int i=0; i<n_trends; i++)
    dim.push_back(trender[i]->GetTrendDimension());

  for(int i=0; i<static_cast<int>(size.size()); i++)
    size[i] = 0;

  for(int i=0; i<static_cast<int>(increment.size()); i++)
    increment[i] = 0;

  if(new_dim == 1) {
    // Know that at least one trend is 1D, and that if more than one trend is 1D,
    // they all use the same trend cube reference

    for(int k=0; k<n_trends; k++) {

      if(dim[k] > 0) {
        int reference = trender[k]->GetReference();

        std::vector<int> trend_size = trender[k]->GetTrendSize();
        size[0] = trend_size[reference-1];

        std::vector<double> inc = trender[k]->GetIncrement();
        increment[0] = inc[reference-1];

        break;
      }
    }
  }

  else if(new_dim == 2) {
    std::vector<int>    trend_size(n_trends);
    std::vector<double> inc(n_trends);
    int                 reference;

    for(int k=0; k<n_trends; k++) {

      trend_size = trender[k]->GetTrendSize();
      inc        = trender[k]->GetIncrement();
      reference  = trender[k]->GetReference();

      if(dim[k] == 1) {
        // Build 2D grid from 1D grids with different reference parameter
        // All 1D trends with the same reference parameter have same sampling; hence the size is the same

        if(trend_size[reference-1] != 0) {
          size[reference-1]      = trend_size[reference-1];
          increment[reference-1] = inc[reference-1];
        }
      }

      else if(dim[k] == 2) {
        // Sufficient to find dimension of the first 2D trend,
        // as all 2D trends have same sampling

        for(int i=0; i<static_cast<int>(size.size()); i++)
          size[i] = trend_size[i];

        for(int i=0; i<static_cast<int>(increment.size()); i++)
          increment[i] = inc[i];

        break;
      }
    }
  }
}

std::vector<std::vector<double> >
DistributionsRockTriNormal::ExpandGrids1D(const std::vector<NRLib::Trend *> trender,
                                          const std::vector<int>         &   size) const
{
  std::vector<std::vector<double> > expanded_trend(3);

  std::vector<int> dim(3);
  dim[0] = trender[0]->GetTrendDimension();
  dim[1] = trender[1]->GetTrendDimension();
  dim[2] = trender[2]->GetTrendDimension();

  int dummy = 0;
  for(int k=0; k<3; k++) {
    std::vector<double> grid(size[0]);

    if(dim[k] == 0) {
      double value = trender[k]->GetTrendElement(dummy,dummy,dummy);

      for(int i=0; i<size[0]; i++)
        grid[i] = value;
    }

    else {
      for(int i=0; i<size[0]; i++)
        grid[i] = trender[k]->GetTrendElement(i,dummy,dummy);
    }

    expanded_trend[k] = grid;
  }

  return(expanded_trend);
}

std::vector<NRLib::Grid2D<double> >
DistributionsRockTriNormal::ExpandGrids2D(const std::vector<NRLib::Trend *> trender,
                                          const std::vector<int>         &  size) const
{
  std::vector<NRLib::Grid2D<double> > expanded_trend(3);

  std::vector<int> dim(3);
  dim[0] = trender[0]->GetTrendDimension();
  dim[1] = trender[1]->GetTrendDimension();
  dim[2] = trender[2]->GetTrendDimension();

  int dummy = 0;
  for(int k=0; k<3; k++) {
    NRLib::Grid2D<double> grid(size[0], size[1], 0);

    if(dim[k] == 0) {
      double value = trender[k]->GetTrendElement(dummy,dummy,dummy);

      for(int i=0; i<size[0]; i++) {
        for(int j=0; j<size[1]; j++)
          grid(i,j) = value;
      }
    }

    else if(dim[k] == 1) {
      std::vector<int> trend_size = trender[k]->GetTrendSize();
      int              reference  = trender[k]->GetReference();

      std::vector<double> values(trend_size[reference-1]);

      for(int i=0; i<trend_size[reference-1]; i++)
        values[i] = trender[k]->GetTrendElement(i,dummy,dummy);

      for(int i=0; i<size[0]; i++) {
        for(int j=0; j<size[1]; j++) {

          if(reference == 1)
            grid(i,j) = values[i];

          else
            grid(i,j) = values[j];
        }
      }
    }

    else {
      for(int i=0; i<size[0]; i++) {
        for(int j=0; j<size[1]; j++)
          grid(i,j) = trender[k]->GetTrendElement(i,j,0);
      }
    }

    expanded_trend[k] = grid;
  }

  return(expanded_trend);
}

void
DistributionsRockTriNormal::LogTransformExpectation(const double & expectation,
                                                    const double & log_variance,
                                                    double       & mu) const
{
  // mu = log(E(X))-0.5*log(1+Var(X)/E(X)^2) = log(E(X))-0.5*sigma^2

  mu = std::log(expectation)-0.5*log_variance;

}

void
DistributionsRockTriNormal::LogTransformCovariance(const double & expectation1,
                                                   const double & expectation2,
                                                   const double & covariance,
                                                   double       & s2) const
{
  // sigma^2 = log(1+Var(X)/E(X)^2)

  s2 = std::log(1+covariance/(expectation1*expectation2));

}

void
DistributionsRockTriNormal::CalculateCovariance(const double & corr,
                                                const double & var1,
                                                const double & var2,
                                                double       & cov) const
{
  // Cov(X,Y) = Corr(X,Y)*sqrt(Var(X)*Var(Y))

  cov = corr * std::sqrt(var1*var2);
}

