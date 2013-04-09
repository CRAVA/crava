#include "src/posteriorelasticpdf4d.h"
#include "src/fftgrid.h"
#include "nrlib/grid/grid2d.hpp"
#include <math.h>
#include "src/modelsettings.h"
#include <src/simbox.h>
#include "lib/lib_matr.h"
#include "nrlib/trend/trend.hpp"


PosteriorElasticPDF4D::PosteriorElasticPDF4D(const std::vector<double>                   & d1, // first dimension of data points
                                             const std::vector<double>                   & d2, // second dimension of data points
                                             const std::vector<double>                   & d3, // third dimension of data points
                                             const std::vector<int>                      & t1, // fourth dimension (trend parameters)
                                             const std::vector<int>                      & t2, // fourth dimension (trend parameters)
                                             const std::vector<std::vector<double> >     & v,  // Transformation of elastic variables from 3D to 2D
                                             const double                     *const*const sigma, // Gaussian smoothing kernel in 2D
                                             int                                           n1,    // resolution of density grid in elastic dimension 1
                                             int                                           n2,    // resolution of density grid in elastic dimension 2
                                             int                                           nt1,    // resolution of density grid in trend dimension 1
                                             int                                           nt2,    // resolution of density grid in trend dimension 2
                                             double                                        d1_min,
                                             double                                        d1_max,
                                             double                                        d2_min,
                                             double                                        d2_max,
                                             double                                        t1_min,
                                             double                                        t1_max,
                                             double                                        t2_min,
                                             double                                        t2_max,
                                             int                                           ind):
nx_(n1),
ny_(n2),
nt1_(nt1),
nt2_(nt2),
x_min_(d1_min),
x_max_(d1_max),
y_min_(d2_min),
y_max_(d2_max),
t1_min_(t1_min),
t1_max_(t1_max),
t2_min_(t2_min),
t2_max_(t2_max)
{
  // We assume that the input vectors are of the same length
  if (d1.size()!=d2.size() || d2.size()!=d3.size() || d3.size()!=t1.size() || t1.size() != t2.size())
    throw NRLib::Exception("Facies probabilities: Size of input vectors do not match.");
  if (v.size()>2 || static_cast<int>(v[0].size()) != 3 ||  static_cast<int>(v[1].size()) != 3)
    throw NRLib::Exception("Facies probabilities: Transformation matrix v does not have the right dimensions");

  v1_.resize(3,0.0);
  v2_.resize(3,0.0);

  for(int i=0; i<3; i++){
    v1_[i] = v[0][i];
    v2_[i] = v[1][i];
  }

  // set size of 2D grid to ni*nj
  histogram_.Resize(nt1_,nt2_,NULL);
  int nData = 0;

  int dim = static_cast<int>(d1.size());

  std::vector<std::vector<double> > x(2);
  x[0].resize(dim,0.0);
  x[1].resize(dim,0.0);

  for(int i=0; i<nt1_; i++){
    for (int j=0; j<nt2_; j++){
      histogram_(i,j) = new FFTGrid(nx_, ny_, 1, nx_, ny_, 1);
      //histogram_(i,j)->createRealGrid(true);
      //int rnxp = histogram_(i,j)->getRNxp();
      histogram_(i,j)->setType(FFTGrid::PARAMETER);
      histogram_(i,j)->setAccessMode(FFTGrid::WRITE);

      histogram_(i,j)->fillInConstant(0.0);
      /*
      for(int k=0;k<ny_;k++){
        for(int j=0;j<rnxp;j++)
          histogram_(i,j)->setNextReal(0.0f);
      }
      */

      histogram_(i,j)->endAccess();
    }
  }

  //computes x and y from d1, d2 and d3
  CalculateTransform2D(d1, d2, d3, x, v);

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/nx_;
  dy_ = (y_max_ - y_min_)/ny_;
  dt1_ = (t1_max_ - t1_min_)/nt1_;
  dt2_ = (t2_max_ - t2_min_)/nt2_;

  // Loop over data points and place in bins in histogram_

  for (int l = 0; l < dim; l++){
    int i = t1[l];
    int j = t2[l];
    int m = static_cast<int>(floor((x[0][l]-x_min_)/dx_));
    int n = static_cast<int>(floor((x[1][l]-y_min_)/dy_));

    // Counting data points in index (i,j,k)
    histogram_(i,j)->setAccessMode(FFTGrid::RANDOMACCESS);
    histogram_(i,j)->setRealValue(m,n,0, histogram_(i,j)->getRealValue(m,n,0) + 1);
    histogram_(i,j)->endAccess();
    nData++;
  }

  //multiply by normalizing constant for the PDF - dim is the total number of entries
  for(int i=0; i<nt1_; i++){
    for (int j=0; j<nt2_; j++){
      histogram_(i,j)->setAccessMode(FFTGrid::READANDWRITE);
      histogram_(i,j)->multiplyByScalar(float(1.0f/nData));
      histogram_(i,j)->endAccess();
      if(ModelSettings::getDebugLevel() >= 1){
        std::string baseName = "Hist_" + NRLib::ToString(ind) + IO::SuffixAsciiFiles();
        std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
        //histogram_(i,j)->writeAsciiFile(fileName);
      }
      histogram_(i,j)->fftInPlace();

      double **sigma_tmp = new double *[2];
      for (int m=0;m<2;m++){
        sigma_tmp[m] = new double[2];
      }
      for(int m=0;m<2;m++){
        for(int n=0; n<2; n++)
          sigma_tmp[m][n] = sigma[m][n];
      }

      // Matrix inversion of the covariance matrix sigma
      double **sigma_inv = new double *[2];
      for(int m=0; m<2; m++)
        sigma_inv[m] = new double [2];

      InvertSquareMatrix(sigma_tmp,sigma_inv,2);

      FFTGrid * smoother = new FFTGrid(nx_, ny_, 1, nx_, ny_, 1);

      smoother->createRealGrid(false);
      smoother->setType(FFTGrid::PARAMETER);
      smoother->setAccessMode(FFTGrid::WRITE);
      for(int k=0;k<ny_;k++){
        for(int l=0;l<static_cast<int>(smoother->getRNxp());l++)
          smoother->setNextReal(0.0f);
      }
      smoother->endAccess();

      SetupSmoothingGaussian2D(smoother, sigma_inv, nx_, ny_, 1, dx_, dy_);

      // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
      smoother->fftInPlace();
      histogram_(i,j)->multiply(smoother);
      histogram_(i,j)->invFFTInPlace();
      histogram_(i,j)->multiplyByScalar(sqrt(float(nx_*ny_*1)));
      histogram_(i,j)->endAccess();

      delete smoother;
      for(int i=0;i<2;i++){
        delete [] sigma_inv[i];
        delete [] sigma_tmp[i];
      }
      delete [] sigma_tmp;
      delete [] sigma_inv;
    }
  }
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
double PosteriorElasticPDF4D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho,
                                      const double & s1) const
{
  (void) s1;
  return this->Density(vp,vs,rho);
}

double PosteriorElasticPDF4D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho) const
{
  (void) vp;
  (void) vs;
  (void) rho;
  return 0;
}

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

    /***/
    /*
    std::string baseName = "Smoothed_hist_" + NRLib::ToString(i) + NRLib::ToString(j) + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    histogram_(i,j)->writeAsciiFile(fileName);
    */
    /****/

    returnvalue = histogram_(i,j)->InterpolateTrilinear(x_min_, x_max_,
                                   y_min_, y_max_, 0, 0, x, y, 0);
    // If the value is outside FFTGrid, return probability 0
    if (returnvalue == RMISSING)
      returnvalue = 0;

  }

  return returnvalue;
}

void PosteriorElasticPDF4D::WriteAsciiFile(std::string filename,
                                           int         i,
                                           int         j) const
{
  if(i<nx_ && j<ny_)
    histogram_(i,j)->writeAsciiFile(filename);
}

double PosteriorElasticPDF4D::FindDensity(const double & vp,
                                          const double & vs,
                                          const double & rho,
                                          const double & s1,
                                          const double & s2,
                                          const Simbox * const /*volume*/) const{
  return Density(vp, vs, rho, s1, s2);
}

void PosteriorElasticPDF4D::ResampleAndWriteDensity(const std::string & /*fileName*/,
                                                    const Simbox      * /*origVol*/,
                                                    Simbox            * /*volume*/,
                                                    int                 /*gridNo*/,
                                                    bool                /*writeSurface*/) const{
  LogKit::LogFormatted(LogKit::Low,"\nWARNING: Not possible to resample and write a file for facies probabilities 4D density objects\n");
}
