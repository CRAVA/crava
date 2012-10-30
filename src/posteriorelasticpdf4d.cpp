#include "src/posteriorelasticpdf4d.h"
#include "src/fftgrid.h"
#include "nrlib/grid/grid2d.hpp"
#include <math.h>
#include <algorithm>
#include <src/simbox.h>
#include "lib/lib_matr.h"
#include "nrlib/trend/trend.hpp"


PosteriorElasticPDF4D::PosteriorElasticPDF4D(const std::vector<double>  & d1,
                 const std::vector<double>                              & d2,
                 const std::vector<double>                              & d3,
                 const std::vector<double>                              & t1,
                 const std::vector<double>                              & t2,
                 double                                                ** sigma_prior,
                 double                                                ** sigma_posterior,
                 int                                                      nx,
                 int                                                      ny,
                 int                                                      nt1,
                 int                                                      nt2,
                 double                                                   t1_min,
                 double                                                   t1_max,
                 double                                                   t2_min,
                 double                                                   t2_max):
nx_(nx),
ny_(ny),
nt1_(nt1),
nt2_(nt2),
t1_min_(t1_min),
t1_max_(t1_max),
t2_min_(t2_min),
t2_max_(t2_max)
{

  // set size of 2D grid to ni*nj
  histogram_.Resize(nt1_,nt2_);

  // We assume that the input vectors are of the same length
  assert(d1.size()==d2.size());
  assert(d2.size()==d3.size());
  assert(d3.size()==t1.size());
  assert(t1.size() ==t2.size());

  int dim = static_cast<int>(d1.size());

  //Vectors holding the transformed variables
  std::vector<double> x(dim);
  std::vector<double> y(dim);
  //Linear transformation vectors
  std::vector<double> v1(3);
  std::vector<double> v2(3);

  // Create new 2D covariance matrix from the input smoothing variables
  double **sigma_2d = new double *[2];
  for(int i=0; i<2; i++)
    sigma_2d[i] = new double[2];

  // computes the optimal linear transforms v1 and v2
  this->SolveGEVProblem(sigma_prior, sigma_posterior, v1, v2);

  //computes x and y
  this->CalculateTransform2D(d1, d2, d3, x, y, v1, v2);

  // computes the resulting 2D variance
  this->CalculateVariance2D(sigma_posterior, sigma_2d, v1, v2);

  v1_ = std::vector<double>(3);
  v2_ = std::vector<double>(3);

  for(int i=0;i<3;i++){
    v1_[i] = v1[i];
    v2_[i] = v2[i];
  }

  // Establish minimum and maximum limits of the input data

  x_min_ = *min_element(x.begin(), x.end()) - 5.0f*sqrt(sigma_2d[0][0]);
  x_max_ = *max_element(x.begin(), x.end()) + 5.0f*sqrt(sigma_2d[0][0]);
  y_min_ = *min_element(y.begin(), y.end()) - 5.0f*sqrt(sigma_2d[1][1]);
  y_max_ = *max_element(y.begin(), y.end()) + 5.0f*sqrt(sigma_2d[1][1]);

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/nx_;
  dy_ = (y_max_ - y_min_)/ny_;
  dt1_ = (t1_max_ - t1_min_)/nt1_;
  dt2_ = (t2_max_ - t2_min_)/nt2_;


  // Matrix inversion of the covariance matrix sigma
  double **sigma_inv = new double *[2];
  for(int i=0; i<2; i++)
    sigma_inv[i] = new double [2];

  this->InvertSquareMatrix(sigma_2d, sigma_inv, 2);

  // For each combination of trend parameters, initialize a density grid as an FFTGrid

  for(int i=0; i<nt1_; i++){
    for (int j=0; j<nt2_; j++){
      histogram_(i,j) = new FFTGrid(nx_, ny_, 1, nx_, ny_, 1);
      histogram_(i,j)->fillInConstant(0); // initialize to zero
    }
  }

  // Loop over data points and place in bins in histogram_

  for (int ind = 0; ind < dim; ind++){

    int i = static_cast<int>(floor((t1[ind]-t1_min_)/dt1_));
    int j = static_cast<int>(floor((t2[ind]-t2_min_)/dt2_));
    int m = static_cast<int>(floor((x[ind]-x_min_)/dx_));
    int n = static_cast<int>(floor((y[ind]-y_min_)/dy_));

    // Counting data points in index (i,j,k)
    histogram_(i,j)->setRealValue(m,n,0, histogram_(i,j)->getRealValue(m,n,0) + 1);
  }

  // Setup of Gaussian smoother
  float *smooth = new float[nx_*ny_];

  int mm,nn,mmm,nnn;
  nnn=2;

  //sum is the normalizing constant
  float sum = 0.0;

  for(int n=0; n<ny_; n++) {
    mmm=2;
    if(n<=ny_/2)
      nn=n;
    else {
      nn = -(n-nnn);
      nnn+=2;
    }
    for(int m=0; m<nx_; m++) {
      if(m<=nx_/2)
        mm=m;
      else {
        mm = -(m-mmm);
        mmm+=2;
      }
      smooth[m+n*nx_] = float(exp(-0.5f*(mm*dx_ *mm*dx_*float(sigma_inv[0][0])
                        +nn*dy_  *nn*dy_*float(sigma_inv[1][1])
                        +2*mm*dx_*nn*dy_*float(sigma_inv[1][0]))));

      sum = sum+smooth[m+n*nx_];
    }
  }



  FFTGrid *smoother;
  smoother = new FFTGrid(nx_, ny_, 1, nx_, ny_, 1);
  smoother->fillInFromArray(smooth);
  //normalizing constant for the smoother
  float normalizer = float(1.0/(sum*dx_*dy_));
  smoother->multiplyByScalar(normalizer);
  smoother->fftInPlace();

  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  for (int i=0; i<nt1_; i++){
    for (int j=0; j<nt2_; j++){
      //multiply by normalizing constant for the pdf - dim is the total number of entries
      histogram_(i,j)->multiplyByScalar(float(1/(dim*dx_*dy_*dt1_*dt2_)));
      histogram_(i,j)->fftInPlace();
      histogram_(i,j)->multiply(smoother);
      histogram_(i,j)->invFFTInPlace();
      // normalizing constant for the fft
      histogram_(i,j)->multiplyByScalar(sqrt(float(nx_*ny_)));
    }
  }
  delete smoother;
  delete [] smooth;
  for(int i=0;i<2;i++)
    delete [] sigma_2d[i];
  delete [] sigma_2d;
  for(int i=0;i<2;i++)
    delete [] sigma_inv[i];
  delete [] sigma_inv;
}

PosteriorElasticPDF4D::PosteriorElasticPDF4D(int nx,
                                             int ny,
                                             int nt1,
                                             int nt2):
nx_(nx),
ny_(ny),
nt1_(nt1),
nt2_(nt2)
{
}

PosteriorElasticPDF4D::~PosteriorElasticPDF4D()
{
  for (size_t i=0; i<histogram_.GetNI(); i++){
    for (size_t j=0; j<histogram_.GetNJ(); j++){
      delete histogram_(i,j);
    }
  }
}

// Density function with trend parameters
/*double PosteriorElasticPDF4D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho,
                                      const double & s1) const
{
  (void) s1;
  (void) s2;
  return this->Density(vp,vs,rho);
}

double PosteriorElasticPDF4D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho) const
{
  (void) s1;
  return this->Density(vp,vs,rho);
}*/

// Density function with trend parameters
double PosteriorElasticPDF4D::Density(const double & vp,
                               const double & vs,
                               const double & rho,
                               const double & s1,
                               const double & s2
                               ) const
{
  double x = vp*v1_[0] + vs*v1_[1] + rho*v1_[2];
  double y = vp*v2_[0] + vp*v2_[1] + rho*v2_[2];

  double returnvalue;

  if (s1<t1_min_ || s1> t1_max_ || s2<t2_min_ || s2>t2_max_){
    // skal returnvalue være 0 eller missing hvis utenfor gridet?
    returnvalue = 0;
  }else{
    int i = static_cast<int>(floor((s1-t1_min_)/dt1_));
    int j = static_cast<int>(floor((s2-t2_min_)/dt2_));
    returnvalue = histogram_(i,j)->InterpolateTrilinear(x_min_, x_max_,
                                   y_min_,y_max_,0,0,x,y,0);
    // If the value is outside FFTGrid, return probability 0
    if (returnvalue == RMISSING)
      returnvalue = 0;

  }

  return returnvalue;
}

/*void PosteriorElasticPDF4D::WriteAsciiFile(std::string filename) const
{

}*/

void PosteriorElasticPDF4D::ResampleAndWriteDensity(const std::string & /*fileName*/,
                                                    const Simbox      * /*origVol*/,
                                                    Simbox            * /*volume*/,
                                                    int                 /*gridNo*/,
                                                    bool                /*writeSurface*/) const
{
}
