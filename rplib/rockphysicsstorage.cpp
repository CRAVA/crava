#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/trend/trend.hpp"
#include "nrlib/grid/grid2d.hpp"

#include "rplib/rockphysicsstorage.h"
#include "rplib/distributionsrockt0.h"
#include "rplib/trinormalwith2dtrend.h"
#include "rplib/multinormaldistributedrockt0.h"

RockPhysicsStorage::RockPhysicsStorage()
{
}

RockPhysicsStorage::~RockPhysicsStorage()
{
}

GaussianRockPhysicsStorage::GaussianRockPhysicsStorage(NRLib::TrendStorage *mean_vp,
                                                       NRLib::TrendStorage *mean_vs,
                                                       NRLib::TrendStorage *mean_density,
                                                       NRLib::TrendStorage *variance_vp,
                                                       NRLib::TrendStorage *variance_vs,
                                                       NRLib::TrendStorage *variance_density,
                                                       NRLib::TrendStorage *correlation_vp_vs,
                                                       NRLib::TrendStorage *correlation_vp_density,
                                                       NRLib::TrendStorage *correlation_vs_density)
: mean_vp_(mean_vp),
  mean_vs_(mean_vs),
  mean_density_(mean_density),
  variance_vp_(variance_vp),
  variance_vs_(variance_vs),
  variance_density_(variance_density),
  correlation_vp_vs_(correlation_vp_vs),
  correlation_vp_density_(correlation_vp_density),
  correlation_vs_density_(correlation_vs_density)
{
}

GaussianRockPhysicsStorage::~GaussianRockPhysicsStorage()
{
}

DistributionsRockT0 *
GaussianRockPhysicsStorage::GenerateRockPhysics(const std::string              & path,
                                                const std::vector<std::string> & trend_cube_names,
                                                std::string                    & errTxt) const
{
  NRLib::Trend * mean_vp_trend                = mean_vp_               ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * mean_vs_trend                = mean_vs_               ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * mean_density_trend           = mean_density_          ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * variance_vp_trend            = variance_vp_           ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * variance_vs_trend            = variance_vs_           ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * variance_density_trend       = variance_density_      ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * correlation_vp_vs_trend      = correlation_vp_vs_     ->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * correlation_vp_density_trend = correlation_vp_density_->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * correlation_vs_density_trend = correlation_vs_density_->GenerateTrend(path,trend_cube_names,errTxt);
  NRLib::Trend * covariance_vp_vs_trend       = NULL;
  NRLib::Trend * covariance_vp_density_trend  = NULL;
  NRLib::Trend * covariance_vs_density_trend  = NULL;

  //Marit: Resample trendverdiene til de som trengs, vha global min/max
  //Pass på å vri 2D-gridene i resamplingen, slik at alle får samme parametre langs samme akser
  //Ønsker samme oppløsning for alle trender etter resampling. Finner dette fra trendkubene, utfra min/max-verdi der

  CalculateCovarianceFromCorrelation(correlation_vp_vs_trend,
                                     variance_vp_trend,
                                     variance_vs_trend,
                                     covariance_vp_vs_trend);
  CalculateCovarianceFromCorrelation(correlation_vp_density_trend,
                                     variance_vp_trend,
                                     variance_density_trend,
                                     covariance_vp_density_trend);
  CalculateCovarianceFromCorrelation(correlation_vs_density_trend,
                                     variance_vs_trend,
                                     variance_density_trend,
                                     covariance_vs_density_trend);

  NRLib::Grid2D<NRLib::Trend *> cov_mat(3,3,NULL);
  cov_mat(0,0) = variance_vp_trend;
  cov_mat(1,1) = variance_vs_trend;
  cov_mat(2,2) = variance_density_trend;
  cov_mat(0,1) = covariance_vp_vs_trend;
  cov_mat(0,2) = covariance_vp_density_trend;
  cov_mat(1,2) = covariance_vs_density_trend;

  std::vector<NRLib::Trend * > mean_vec(3);
  mean_vec[0] = mean_vp_trend;
  mean_vec[1] = mean_vs_trend;
  mean_vec[2] = mean_density_trend;

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

  delete correlation_vp_vs_trend;
  delete correlation_vp_density_trend;
  delete correlation_vs_density_trend;

  for(int i=0; i<3; i++) {
    delete mean_vec[i];

    for(int j=0; j<3; j++) {
      delete cov_mat(i,j);
    }
  }

  TriNormalWith2DTrend multi(log_mean[0],
                             log_mean[1],
                             log_mean[2],
                             log_cov);

  DistributionsRockT0 * rock = new MultiNormalDistributedRockT0(multi);

  return(rock);
}

void
GaussianRockPhysicsStorage::LogTransformExpectationAndCovariance(NRLib::Trend *  mean1,
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

  int new_dim           = FindNewGridDimension(trender);
  std::vector<int> size = FindNewGridSize(trender, new_dim);

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

      log_mean = new NRLib::ConstantTrend(log_mu);
    }
    log_cov = new NRLib::ConstantTrend(log_covariance);
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
      log_mean = new NRLib::Trend1D(log_mu, reference);

    log_cov = new NRLib::Trend1D(log_covariance,reference);
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
      log_mean = new NRLib::Trend2D(log_mu,1,2);

    log_cov = new NRLib::Trend2D(log_covariance,1,2);
  }
}

void
GaussianRockPhysicsStorage::CalculateCovarianceFromCorrelation(NRLib::Trend *  corr,
                                                               NRLib::Trend *  var1,
                                                               NRLib::Trend *  var2,
                                                               NRLib::Trend *& cov) const
{
  std::vector<NRLib::Trend *> trender(3);
  trender[0] = corr;
  trender[1] = var1;
  trender[2] = var2;

  int new_dim           = FindNewGridDimension(trender);
  std::vector<int> size = FindNewGridSize(trender, new_dim);

  if(new_dim == 0) {
    double covariance;
    int    dummy = 0;

    CalculateCovariance(trender[0]->GetTrendElement(dummy,dummy,dummy),
                        trender[1]->GetTrendElement(dummy,dummy,dummy),
                        trender[2]->GetTrendElement(dummy,dummy,dummy),
                        covariance);

    cov = new NRLib::ConstantTrend(covariance);
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

    cov = new NRLib::Trend1D(covariance,reference);
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

    cov = new NRLib::Trend2D(covariance,1,2);
  }
}

int
GaussianRockPhysicsStorage::FindNewGridDimension(const std::vector<NRLib::Trend *> trender) const
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

std::vector<int>
GaussianRockPhysicsStorage::FindNewGridSize(const std::vector<NRLib::Trend *> trender,
                                            const int                         new_dim) const
{
  std::vector<int> dim;
  for(int i=0; i<static_cast<int>(trender.size()); i++)
    dim.push_back(trender[i]->GetTrendDimension());

  std::vector<int> size(2);
  for(int i=0; i<2; i++)
    size[i] = 0;

  if(new_dim == 1) {
    for(int k=0; k<static_cast<int>(trender.size()); k++) {

      if(dim[k] > 0) {
        std::vector<int> trend_size = trender[k]->GetTrendSize();
        int reference = trender[k]->GetReference();
        size[0] = trend_size[reference-1];
      }
    }
  }

  else if(new_dim == 2) {
    std::vector<int> trend_size(static_cast<int>(trender.size()));
    int              reference;

    for(int i=0; i<2; i++)
      size[i] = 0;

    for(int k=0; k<static_cast<int>(trender.size()); k++) {
      if(dim[k] == 1) { //Build 2D grid from two 1D grids with different reference parameter

        trend_size = trender[k]->GetTrendSize();
        reference  = trender[k]->GetReference();

        if(trend_size[reference-1] != 0)
          size[reference-1] = trend_size[reference-1];
      }

      else if(dim[k] == 2) {

        trend_size = trender[k]->GetTrendSize();

        for(int i=0; i<2; i++)
          size[i] = trend_size[i];
      }
    }
  }

  return(size);
}

std::vector<std::vector<double> >
GaussianRockPhysicsStorage::ExpandGrids1D(const std::vector<NRLib::Trend *> trender,
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
GaussianRockPhysicsStorage::ExpandGrids2D(const std::vector<NRLib::Trend *> trender,
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
GaussianRockPhysicsStorage::LogTransformExpectation(const double & expectation,
                                                    const double & log_variance,
                                                    double       & mu) const
{
  // mu = log(E(X))-0.5*log(1+Var(X)/E(X)^2) = log(E(X))-0.5*sigma^2

  mu = std::log(expectation)-0.5*log_variance;

}

void
GaussianRockPhysicsStorage::LogTransformCovariance(const double & expectation1,
                                                   const double & expectation2,
                                                   const double & covariance,
                                                   double       & s2) const
{
  // sigma^2 = log(1+Var(X)/E(X)^2)

  s2 = std::log(1+covariance/(expectation1*expectation2));

}

void
GaussianRockPhysicsStorage::CalculateCovariance(const double & corr,
                                                const double & var1,
                                                const double & var2,
                                                double       & cov) const
{
  // Cov(X,Y) = Corr(X,Y)*sqrt(Var(X)*Var(Y))

  cov = corr * std::sqrt(var1*var2);
}
