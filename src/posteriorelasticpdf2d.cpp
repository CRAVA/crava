#include "src/posteriorelasticpdf2d.h"
#include "src/fftgrid.h"
#include <src/simbox.h>
#include <math.h>
#include "lib/lib_matr.h"


PosteriorElasticPDF2D::PosteriorElasticPDF2D(const std::vector<double> & d1,
                               const std::vector<double> & d2,
                               const std::vector<double> & d3,
                               int n1,
                               int n2,
                               double **smoothvar,
                               const std::vector<double> & v1,
                               const std::vector<double> & v2)
                               :
n1_(n1),
n2_(n2)
{

  int dim = static_cast<int>(d1.size());

  // v1 and v2 should both be of size 3
  assert (v1.size()== 3);
  assert (v2.size()== 3);

  //Vectors holding the transformed variables
  std::vector<double> x(dim);
  std::vector<double> y(dim);

  // x and y are linear combinations of d1, d2 and d3 given by the 3d vectors v1 and v2.
  this->CalculateTransform2D(d1, d2, d3, x, y, v1, v2);

  for (int i = 0; i<3; i++){
    v1_[i] = v1[i];
    v2_[i] = v2[i];
  }

  // Create covariance matrix from the input smoothing variables
  double ** sigma_2d = new double *[2];
  for(int i=0; i<2; i++)
    sigma_2d[i] = new double[2];

  this->CalculateVariance2D(smoothvar, sigma_2d, v1, v2);

  // Establish minimum and maximum limits of the input data
  // (Should this be done before inverting sigma?)
  x_min_ = *min_element(x.begin(), x.end()) - 5.0f*sqrt(sigma_2d[0][0]);
  x_max_ = *max_element(x.begin(), x.end()) + 5.0f*sqrt(sigma_2d[0][0]);
  y_min_ = *min_element(y.begin(), y.end()) - 5.0f*sqrt(sigma_2d[1][1]);
  y_max_ = *max_element(y.begin(), y.end()) + 5.0f*sqrt(sigma_2d[1][1]);

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/n1_;
  dy_ = (y_max_ - y_min_)/n2_;

  // Matrix inversion of the covariance matrix sigma
  double **sigma_inv = new double *[2];
  for(int i=0; i<2; i++)
    sigma_inv[i] = new double [2];

  this->InvertSquareMatrix(sigma_2d,sigma_inv,2);

  for(int i=0; i<2; i++)
    delete [] sigma_2d[i];
  delete [] sigma_2d;

  // Create density grid as a FFTGrid
  histogram_ = new FFTGrid(n1_, n2_, 1, n1_, n2_, 1);
  histogram_->fillInConstant(0); // initialize to zero


  // Go through data points and place in bins in histogram
  for (int ind = 0; ind < static_cast<int>(x.size()); ind++){

    int i = static_cast<int>(floor((x[ind]-x_min_)/dx_));
    int j = static_cast<int>(floor((y[ind]-y_min_)/dy_));
    int k = 0;

    // Counting data points in index (i,j,k)
    histogram_->setRealValue(i,j,k, histogram_->getRealValue(i,j,k) + 1);
  }

  // Setup of Gaussian smoother

  FFTGrid *smoother = new FFTGrid(n1, n2_, 1, n1_, n2_, 1);
  this->SetupSmoothingGaussian2D(smoother,sigma_inv,n1,n2,dx_,dy_);
  smoother->fftInPlace();

  //multiply by normalizing constant for the PDF - dim is the total number of entries
  histogram_->multiplyByScalar(float(1/(dim*dx_*dy_)));
  histogram_->fftInPlace();
  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  histogram_->multiply(smoother);
  histogram_->invFFTInPlace();
  histogram_->multiplyByScalar(float(sqrt(double(n1_*n2_))));

  delete smoother;
  for(int i=0;i<2;i++)
    delete [] sigma_inv[i];
  delete [] sigma_inv;
}

PosteriorElasticPDF2D::PosteriorElasticPDF2D(const std::vector<double> & d1,
                               const std::vector<double> & d2,
                               const std::vector<double> & d3,
                               int n1,
                               int n2,
                               double ** sigma_prior,
                               double ** sigma_posterior)
                               :
n1_(n1),
n2_(n2)
{
  int dim = static_cast<int>(d1.size());

  //Linear transformation vectors
  std::vector<double> v1(3);
  std::vector<double> v2(3);
  //Vectors holding the transformed variables
  std::vector<double> x(dim);
  std::vector<double> y(dim);

  // find v1 and v2 automagically
  this->SolveGEVProblem(sigma_prior, sigma_posterior, v1, v2);

  for (int i = 0; i<3; i++){
    v1_[i] = v1[i];
    v2_[i] = v2[i];
  }

  // x and y are linear combinations of d1, d2 and d3 given by the 3d vectors v1 and v2.
  this->CalculateTransform2D(d1, d2, d3, x, y, v1, v2);

  // Create covariance matrix from the input smoothing variables
  double ** sigma_2d = new double *[2];
  for(int i=0; i<2; i++)
    sigma_2d[i] = new double[2];

  this->CalculateVariance2D(sigma_posterior, sigma_2d, v1, v2);

  // Establish minimum and maximum limits of the input data
  // (Should this be done before inverting sigma?)
  x_min_ = *min_element(x.begin(), x.end()) - 5.0f*sqrt(sigma_2d[0][0]);
  x_max_ = *max_element(x.begin(), x.end()) + 5.0f*sqrt(sigma_2d[0][0]);
  y_min_ = *min_element(y.begin(), y.end()) - 5.0f*sqrt(sigma_2d[1][1]);
  y_max_ = *max_element(y.begin(), y.end()) + 5.0f*sqrt(sigma_2d[1][1]);

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/n1_;
  dy_ = (y_max_ - y_min_)/n2_;

  // Matrix inversion of the covariance matrix sigma
  double **sigma_inv = new double *[2];
  for(int i=0; i<2; i++)
    sigma_inv[i] = new double [2];

  this->InvertSquareMatrix(sigma_2d,sigma_inv,2);

  for(int i=0; i<2; i++)
    delete [] sigma_2d[i];
  delete [] sigma_2d;

  // Create density grid as a FFTGrid
  histogram_ = new FFTGrid(n1_, n2_, 1, n1_, n2_, 1);
  histogram_->fillInConstant(0); // initialize to zero


  // Go through data points and place in bins in histogram
  for (int ind = 0; ind < static_cast<int>(x.size()); ind++){

    int i = static_cast<int>(floor((x[ind]-x_min_)/dx_));
    int j = static_cast<int>(floor((y[ind]-y_min_)/dy_));
    int k = 0;

    // Counting data points in index (i,j,k)
    histogram_->setRealValue(i,j,k, histogram_->getRealValue(i,j,k) + 1);
  }

  // Setup of Gaussian smoother
  FFTGrid *smoother = new FFTGrid(n1, n2, 1, n1, n2, 1);
  this->SetupSmoothingGaussian2D(smoother,sigma_inv,n1,n2,dx_,dy_);
  smoother->fftInPlace();

  //multiply by normalizing constant for the PDF - dim is the total number of entries
  histogram_->multiplyByScalar(float(1/(dim*dx_*dy_)));
  histogram_->fftInPlace();
  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  histogram_->multiply(smoother);
  histogram_->invFFTInPlace();
  histogram_->multiplyByScalar(float(sqrt(double(n1_*n2_))));

  delete smoother;
  for(int i=0;i<2;i++)
    delete [] sigma_inv[i];
  delete [] sigma_inv;
}


PosteriorElasticPDF2D::PosteriorElasticPDF2D(int n1, int n2)
                               :
n1_(n1),
n2_(n2)
{
}

PosteriorElasticPDF2D::~PosteriorElasticPDF2D()
{
  delete histogram_;
}

void PosteriorElasticPDF2D::SetupSmoothingGaussian2D(FFTGrid * smoother,
                                double ** smoothingVar,
                                int n1,
                                int n2,
                                double dx,
                                double dy)
{
  float *smooth = new float[n1*n2];
  int j,k,jj,kk,jjj,kkk;
  kkk=2;

  float sum = 0.0f;

  for(k=0; k<n2; k++) {
    jjj=2;
    if(k<=n2/2)
      kk=k;
    else {
      kk = -(k-kkk);
      kkk+=2;
    }
    for(j=0; j<n1; j++) {
      if(j<=n1/2)
        jj=j;
      else {
        jj = -(j-jjj);
        jjj+=2;
      }
      smooth[j+k*n1] = float(exp(-0.5f*(jj*dx*jj*dx*smoothingVar[0][0]
                        +kk*dy*kk*dy*smoothingVar[1][1]
                        +2*jj*dx*kk*dy*smoothingVar[1][0])));
      sum = sum+smooth[j+k*n1];
    }
  }

  smoother->fillInFromArray(smooth);
  //normalizing constant for the smoother
  smoother->multiplyByScalar(float(1/(sum*dx*dy)));
}


// Density function with trend parameters
double PosteriorElasticPDF2D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho,
                                      const double & s1,
                                      const double & s2) const
{
  (void) s1;
  (void) s2;
  return this->Density(vp,vs,rho);
}

double PosteriorElasticPDF2D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho,
                                      const double & s1) const
{
  (void) s1;
  return this->Density(vp,vs,rho);
}

/*void PosteriorElasticPDF2D::WriteAsciiFile(std::string filename) const
{

}*/

void PosteriorElasticPDF2D::ResampleAndWriteDensity(const std::string & /* fileName*/,
                                                    const Simbox      * /* origVol*/,
                                                    Simbox            * /* volume*/,
                                                    int                 /* gridNo*/,
                                                    bool                /* writeSurface*/) const
{
}

 // Density function without trend parameters
double PosteriorElasticPDF2D::Density(const double & vp,
                               const double & vs,
                               const double & rho) const
{

    // Transform the elastic variables to 2D
    double x = vp*v1_[0] + vs*v1_[1] + rho*v1_[2];
    double y = vp*v2_[0] + vp*v2_[1] + rho*v2_[2];

    // Trilinear interpolation with z = 0
    double returnvalue =  histogram_->InterpolateTrilinear(x_min_, x_max_,
                                   y_min_, y_max_, 0, 0, x, y, 0);


    return (returnvalue);
}

