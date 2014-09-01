#include <math.h>
#include <iostream>

#include "src/posteriorelasticpdf.h"
#include "src/posteriorelasticpdf3d.h"
#include "src/fftgrid.h"
#include "src/modelsettings.h"
#include <src/simbox.h>


PosteriorElasticPDF3D::PosteriorElasticPDF3D(const std::vector<double> & d1,    // first dimension of data points
                 const std::vector<double>                             & d2,    // second dimension of data points
                 const std::vector<double>                             & d3,    // third dimension of data points
                 const double                               *const*const sigma, // Gaussian smoothing kernel in 3D
                 int                                                     n1,    // resolution of density grid in dimension 1
                 int                                                     n2,    // resolution of density grid in dimension 2
                 int                                                     n3,    // resolution of density grid in dimension 3
                 double                                                  d1_min,
                 double                                                  d1_max,
                 double                                                  d2_min,
                 double                                                  d2_max,
                 double                                                  d3_min,
                 double                                                  d3_max,
                 int                                                     ind)
                 :
n1_(n1),
n2_(n2),
n3_(n3),
x_min_(d1_min),
x_max_(d1_max),
y_min_(d2_min),
y_max_(d2_max),
z_min_(d3_min),
z_max_(d3_max)
{
  // We assume that the input vectors are of the same length
  if (d1.size()!=d2.size() || d2.size()!=d3.size())
    throw NRLib::Exception("Facies probabilities: Size of input vectors do not match.");

  int dim = static_cast<int>(d1.size());

  histogram_ = new FFTGrid(n1_, n2_, n3_, n1_, n2_, n3_);
  histogram_->createRealGrid(false);
  int rnxp = histogram_->getRNxp();
  histogram_->setType(FFTGrid::PARAMETER);
  histogram_->setAccessMode(FFTGrid::WRITE);

  for(int l=0;l<n3_;l++){
    for(int k=0;k<n2_;k++){
      for(int j=0;j<rnxp;j++)
        histogram_->setNextReal(0.0f);
    }
  }
  histogram_->endAccess();

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/n1_;
  dy_ = (y_max_ - y_min_)/n2_;
  dz_ = (z_max_ - z_min_)/n3_;

  // Go through data points and place in bins in histogram

  for (int i = 0; i < dim; i++){
    //volume->getIndexes(d1[i], d2[i], d3[i], i_tmp, j_tmp, k_tmp);
    int i_tmp = static_cast<int>(floor((d1[i]-x_min_)/dx_));
    int j_tmp = static_cast<int>(floor((d2[i]-y_min_)/dy_));
    int k_tmp = static_cast<int>(floor((d3[i]-z_min_)/dz_));
    // Counting data points in index (i,j,k)
    histogram_->setAccessMode(FFTGrid::RANDOMACCESS);
    histogram_->setRealValue(i_tmp,j_tmp,k_tmp, histogram_->getRealValue(i_tmp,j_tmp,k_tmp) + 1.0f);
    histogram_->endAccess();
  }

  //multiply by normalizing constant for the PDF - dim is the total number of entries
  histogram_->setAccessMode(FFTGrid::READANDWRITE);
  histogram_->multiplyByScalar(float(1.0f/dim));
  histogram_->endAccess();

  if(ModelSettings::getDebugLevel() >= 1){
    std::string baseName = "Hist_" + NRLib::ToString(ind) + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    histogram_->writeAsciiFile(fileName);
  }

  histogram_->fftInPlace();

  NRLib::Matrix sigma_tmp(3,3);

  /*
  double **sigma_tmp = new double *[3];
  for (int i=0;i<3;i++){
    sigma_tmp[i] = new double[3];
  }
  */
  for(int i=0;i<3;i++){
    for(int j=0; j<3; j++)
      sigma_tmp(i,j) = sigma[i][j];
    }
  // Matrix inversion of the covariance matrix sigma
  NRLib::Matrix sigma_inv;
  /*
  double **sigma_inv = new double *[3];
  for(int i=0; i<3; i++)
    sigma_inv[i] = new double [3];
  */

  InvertSquareMatrix(sigma_tmp,sigma_inv,3);

  FFTGrid *smoother = new FFTGrid(n1, n2, n3, n1, n2, n3);

  smoother->createRealGrid(false);
  smoother->setType(FFTGrid::PARAMETER);
  smoother->setAccessMode(FFTGrid::WRITE);
  for(int l=0;l<n3_;l++){
      for(int k=0;k<n2_;k++){
        for(int j=0;j<static_cast<int>(smoother->getRNxp());j++)
          smoother->setNextReal(0.0f);
      }
  }
  smoother->endAccess();

  SetupSmoothingGaussian3D(smoother, sigma_inv);

  if(ModelSettings::getDebugLevel() >= 1) {
    std::string baseName = "Smoother" + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    smoother->writeAsciiFile(fileName);
  }

  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  smoother->fftInPlace();
  histogram_->multiply(smoother);
  histogram_->invFFTInPlace();
  histogram_->multiplyByScalar(sqrt(float(n1_*n2_*n3_)));
  histogram_->endAccess();

  delete smoother;
}

PosteriorElasticPDF3D::PosteriorElasticPDF3D(const std::vector<double>               & d1,        // first dimension of data points
                                             const std::vector<double>               & d2,        // second dimension of data points
                                             const std::vector<double>               & d3,        // third dimension of data points
                                             const std::vector<int>                  & t1,        // trend data points
                                             const NRLib::Matrix                     & v,        // Transformation of elastic variables from 3D to 2D
                                             const double             *const*const   sigma,     // Gaussian smoothing kernel in 2D
                                             int                                     n1,        // resolution of density grid in elastic dimension 1
                                             int                                     n2,        // resolution of density grid in elastic dimension 2
                                             int                                     n3,        // resolution of density grid in the trend dimension
                                             double                                  d1_min,
                                             double                                  d1_max,
                                             double                                  d2_min,
                                             double                                  d2_max,
                                             double                                  t1_min,
                                             double                                  t1_max,
                                             int                                     ind)
: n1_(n1),
  n2_(n2),
  n3_(n3),
  x_min_(d1_min),
  x_max_(d1_max),
  y_min_(d2_min),
  y_max_(d2_max),
  z_min_(t1_min),
  z_max_(t1_max)
{
  // We assume that the input vectors are of the same length
  if (d1.size()!=d2.size() || d2.size()!=d3.size() || d3.size()!=t1.size())
    throw NRLib::Exception("Facies probabilities: Size of input vectors do not match.");
  if (v.numRows()>2 )
    throw NRLib::Exception("Facies probabilities: Transformation matrix v does not have the right dimensions");

  v1_.resize(3);
  v2_.resize(3);

  for(int i=0; i<3; i++){
    v1_[i] = v(0,i);
    v2_[i] = v(1,i);
  }

  int dim = static_cast<int>(d1.size());

  std::vector<std::vector<double> > x(2);
  x[0].resize(dim);
  x[1].resize(dim);

  //computes x and y from d1, d2 and d3
  CalculateTransform2D(d1, d2, d3, x, v);

  histogram_ = new FFTGrid(n1_, n2_, n3_, n1_, n2_, n3_);
  histogram_->createRealGrid(false);
  int rnxp = histogram_->getRNxp();
  histogram_->setType(FFTGrid::PARAMETER);
  histogram_->setAccessMode(FFTGrid::WRITE);

  for(int l=0;l<n3_;l++){
    for(int k=0;k<n2_;k++){
      for(int j=0;j<rnxp;j++)
        histogram_->setNextReal(0.0f);
    }
  }
  histogram_->endAccess();

  // Spacing variables in the density grid
  dx_ = (x_max_ - x_min_)/n1_;
  dy_ = (y_max_ - y_min_)/n2_;
  dz_ = (z_max_ - z_min_)/n3_;

  // Go through data points and place in bins in histogram

  for (int i = 0; i < dim; i++){
    //volume->getIndexes(d1[i], d2[i], d3[i], i_tmp, j_tmp, k_tmp);
    int i_tmp = static_cast<int>(floor((x[0][i]-x_min_)/dx_));
    int j_tmp = static_cast<int>(floor((x[1][i]-y_min_)/dy_));
    int k_tmp = t1[i];
    // Counting data points in index (i,j,k)
    histogram_->setAccessMode(FFTGrid::RANDOMACCESS);
    histogram_->setRealValue(i_tmp, j_tmp, k_tmp, histogram_->getRealValue(i_tmp,j_tmp,k_tmp) + 1.0f);
    histogram_->endAccess();
  }

  //multiply by normalizing constant for the PDF - dim is the total number of entries
  histogram_->setAccessMode(FFTGrid::READANDWRITE);
  histogram_->multiplyByScalar(float(1.0f/dim));
  histogram_->endAccess();

  if(ModelSettings::getDebugLevel() >= 1){
    std::string baseName = "Hist_" + NRLib::ToString(ind) + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    histogram_->writeAsciiFile(fileName);
  }

  histogram_->fftInPlace();

  NRLib::Matrix sigma_tmp(2,3,0);// = sigma;
  /*
  double **sigma_tmp = new double *[2];
  for (int i=0;i<2;i++){
    sigma_tmp[i] = new double[2];
  }
  */
  for(int i=0;i<2;i++){
    for(int j=0; j<2; j++)
      sigma_tmp(i,j) = sigma[i][j];
  }

  // Matrix inversion of the covariance matrix sigma
  NRLib::Matrix sigma_inv;

  InvertSquareMatrix(sigma_tmp,sigma_inv,2);

  FFTGrid *smoother = new FFTGrid(n1, n2, n3, n1, n2, n3);

  smoother->createRealGrid(false);
  smoother->setType(FFTGrid::PARAMETER);
  smoother->setAccessMode(FFTGrid::WRITE);
  for(int l=0;l<n3_;l++){
      for(int k=0;k<n2_;k++){
        for(int j=0;j<static_cast<int>(smoother->getRNxp());j++)
          smoother->setNextReal(0.0f);
      }
  }
  smoother->endAccess();

  //SetupSmoothingGaussian2D(smoother, sigma_inv, n1, n2, n3, dx_, dy_);

  if(ModelSettings::getDebugLevel() >= 1) {
    std::string baseName = "Smoother" + IO::SuffixAsciiFiles();
    std::string fileName = IO::makeFullFileName(IO::PathToDebug(), baseName);
    smoother->writeAsciiFile(fileName);
  }

  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  smoother->fftInPlace();
  histogram_->multiply(smoother);
  histogram_->invFFTInPlace();
  histogram_->multiplyByScalar(sqrt(float(n1_*n2_*n3_)));
  histogram_->endAccess();

  delete smoother;
}


PosteriorElasticPDF3D::PosteriorElasticPDF3D(int n1,
                                             int n2,
                                             int n3)
: n1_(n1),
  n2_(n2),
  n3_(n3)
{
}

PosteriorElasticPDF3D::~PosteriorElasticPDF3D()
{
  delete histogram_;
}

/*void PosteriorElasticPDF3D::WriteAsciiFile(std::string filename) const
{
  histogram_->writeAsciiFile(filename);
}*/

void PosteriorElasticPDF3D::ResampleAndWriteDensity(const std::string & fileName,
                                                    const Simbox      * origVol,
                                                    Simbox            * volume,
                                                    int                 gridNo,
                                                    bool                writeSurface) const

{
  if(writeSurface == true) {
    int format = IO::STORM;
    std::string topSurf  = IO::PrefixSurface() + IO::PrefixTop()  + IO::PrefixDensity() + NRLib::ToString(gridNo);
    std::string baseSurf = IO::PrefixSurface() + IO::PrefixBase() + IO::PrefixDensity() + NRLib::ToString(gridNo);
    volume->WriteTopBaseSurfaceGrids(topSurf, baseSurf, IO::PathToInversionResults(), format);
    volume->setTopBotName(topSurf, baseSurf, format);
  }


  FFTGrid expDens(n1_, n2_, n3_, n1_, n2_, n3_);
  expDens.createRealGrid();

  double sum = 0;
  for(int k=0; k<n3_; k++) {
    for(int j=0; j<n2_; j++) {
      for(int i=0; i<n1_; i++) {
        double alpha, beta, rho;
        volume->getCoord(i, j, k, alpha, beta, rho);
        double aInd, bInd, rInd;
        origVol->getInterpolationIndexes(log(alpha), log(beta), log(rho), aInd, bInd, rInd);
        double ti = aInd-floor(aInd);
        double tj = bInd-floor(bInd);
        double tk = rInd-floor(rInd);
        int li  = static_cast<int>(floor(aInd));
        int li2 = (li == n1_-1) ? li : li+1;
        int lj = static_cast<int>(floor(bInd));
        int lj2 = (lj == n2_-1) ? lj : lj+1;
        int lk = static_cast<int>(floor(rInd));
        int lk2 = (lk == n3_-1) ? lk : lk+1;
        double tmpFrontLeft  = histogram_->getRealValue(li, lj, lk)*(1-tk)+histogram_->getRealValue(li, lj, lk2)*tk;
        double tmpFrontRight = histogram_->getRealValue(li2, lj, lk)*(1-tk)+histogram_->getRealValue(li2, lj, lk2)*tk;
        double tmpBackLeft   = histogram_->getRealValue(li, lj2, lk)*(1-tk)+histogram_->getRealValue(li, lj2, lk2)*tk;
        double tmpBackRight  = histogram_->getRealValue(li2, lj2, lk)*(1-tk)+histogram_->getRealValue(li2, lj2, lk2)*tk;
        double tmpLeft  = tmpFrontLeft*(1-tj)+tmpBackLeft*tj;
        double tmpRight = tmpFrontRight*(1-tj)+tmpBackRight*tj;
        double value = (tmpLeft*(1-ti)+tmpRight*ti)/alpha/beta/rho;
        expDens.setRealValue(i, j, k, static_cast<float>(value));
        sum += value;
      }
    }
  }
  expDens.multiplyByScalar(static_cast<float>(1.0/sum));
  expDens.writeFile(fileName, "", volume);
}

void PosteriorElasticPDF3D::SetupSmoothingGaussian3D(FFTGrid              * smoother,
                                                     const NRLib::Matrix  & sigma_inv)
{
  float *smooth = new float[n1_*n2_*n3_];
  int j,k,l,jj,jjj,kk,kkk,ll,lll;
  lll=2;

  float sum = 0.0f;
  for(l=0; l<n3_; l++) {
    kkk=2;
    if(l<=n3_/2)
      ll = l;
    else {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0; k<n2_; k++) {
      jjj=2;
      if(k<=n2_/2)
        kk=k;
      else {
        kk = -(k-kkk);
        kkk+=2;
      }
      for(j=0; j<n1_; j++) {
        if(j<=n1_/2)
          jj=j;
        else {
          jj = -(j-jjj);
          jjj+=2;
        }
        smooth[j+k*n1_+l*n1_*n2_] = float(exp(-0.5f*(jj*dx_*jj*dx_*sigma_inv(0,0)
                                                   +kk*dy_*kk*dy_*sigma_inv(1,1)
                                                   +ll*dz_*ll*dz_*sigma_inv(2,2)
                                                   +2*jj*dx_*kk*dy_*sigma_inv(1,0)
                                                   +2*jj*dx_*ll*dz_*sigma_inv(2,0)
                                                   +2*kk*dy_*ll*dz_*sigma_inv(2,1))));
        sum += smooth[j+k*n1_+l*n1_*n2_];
      }
    }
  }

  // normalize smoother
  for(l=0;l<n3_;l++)
    for(k=0;k<n2_;k++)
      for(j=0;j<n1_;j++)
        smooth[j+k*n1_+l*n1_*n2_]/=sum;

  smoother->fillInFromArray(smooth); //No mode/randomaccess
  //normalizing constant for the smoother
  //smoother->multiplyByScalar(static_cast<float>(1.0/sum)); //No mode required

  delete [] smooth;
}

void PosteriorElasticPDF3D::SetupSmoothingGaussian2D(FFTGrid    * smoother,
                                                     double    ** sigma_inv,
                                                     int          n1,
                                                     int          n2,
                                                     int          n3,
                                                     double       dx,
                                                     double       dy,
                                                     double       dz)
{
  float *smooth = new float[n1*n2*n3];
  int j,k,l,jj,jjj,kk,kkk,ll,lll;
  lll=2;

  float sum = 0.0f;
  for(l=0; l<n3; l++) {
    kkk=2;
    if(l<=n3/2)
      ll = l;
    else {
      ll = -(l-lll);
      lll+=2;
    }
    for(k=0; k<n2; k++) {
      jjj=2;
      if(k<=n2_/2)
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
        smooth[j+k*n1+l*n1*n2] = float(exp(-0.5f*(jj*dx *jj*dx*sigma_inv[0][0]
                                  +kk*dy  *kk*dy*sigma_inv[1][1]
                                  +2*jj*dx*kk*dy*sigma_inv[1][0])));
        sum = sum+smooth[j+k*n1+l*n1*n2];
      }
    }
  }


  smoother->fillInFromArray(smooth);
  //normalizing constant for the smoother
  smoother->multiplyByScalar(float(1/(sum*dx*dy*dz)));
}


double PosteriorElasticPDF3D::Density (const double & vp,
                                       const double & vs,
                                       const double & rho,
                                       const double & s1,
                                       const double & s2) const
{
  (void)s2;
  return this->Density(vp,vs,rho,s1);
}

// Density function with trend parameter
double PosteriorElasticPDF3D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho,
                                      const double & s1) const
{
  // check that there is an active transformation for vp, vs and rho
  // is there a better way to do this?
  assert (v1_.size() == 3);

  // Transform the elastic variables to 2D
  double x = vp*v1_[0] + vs*v1_[1] + rho*v1_[2];
  double y = vp*v2_[0] + vp*v2_[1] + rho*v2_[2];

  double returnvalue = histogram_->InterpolateTrilinear(x_min_, x_max_,
                                   y_min_, y_max_, z_min_, z_max_, x, y, s1);

  return returnvalue;
}


 // Density function without trend parameters
double PosteriorElasticPDF3D::Density(const double & vp,
                                      const double & vs,
                                      const double & rho) const
{

    // Trilinear interpolation with z = 0
    double returnvalue = histogram_->InterpolateTrilinear(x_min_, x_max_,
                                   y_min_, y_max_, z_min_, z_max_, vp, vs, rho);

    return returnvalue;
}

double PosteriorElasticPDF3D::FindDensity(const double & vp,
                                          const double & vs,
                                          const double & rho,
                                          const double & s1,
                                          const double & s2,
                                          const Simbox * const volume) const
{
  double value = 0.0;
  // if the size of v1 and v2 > 0 the PDF has been constructed with dimension reduction and one trend value
  if (v1_.size()>0 && v2_.size()>0){
    value = Density(vp, vs, rho, s1, s2);
    if (value == RMISSING)
      value = 0.0;
  }
  else{
    double jFull, kFull, lFull;

    volume->getInterpolationIndexes(vp, vs, rho, jFull, kFull, lFull);
    int j1,k1,l1;
    int j2,k2,l2;
    float wj,wk,wl;
    j1=k1=l1=j2=k2=l2=0;
    j1 = static_cast<int>(floor(jFull));
    if(j1<0) {
      j1 = 0;
      j2 = 0;
      wj = 0;
    }
    else if(j1>=volume->getnx()-1) {
      j1 = volume->getnx()-1;
      j2 = j1;
      wj = 0;
    }
    else {
      j2 = j1 + 1;
      wj = static_cast<float>(jFull-j1);
    }

    k1 = static_cast<int>(floor(kFull));
    if(k1<0) {
      k1 = 0;
      k2 = 0;
      wk = 0;
    }
    else if(k1>=volume->getny()-1) {
      k1 = volume->getny()-1;
      k2 = k1;
      wk = 0;
    }
    else {
      k2 = k1 + 1;
      wk = static_cast<float>(kFull-k1);
    }

    l1 = static_cast<int>(floor(lFull));
    if(l1<0) {
      l1 = 0;
      l2 = 0;
      wl = 0;
    }
    else if(l1>=volume->getnz()-1) {
      l1 = volume->getnz()-1;
      l2 = l1;
      wl = 0;
    }
    else {
      l2 = l1 + 1;
      wl = static_cast<float>(lFull-l1);
    }


      histogram_->setAccessMode(FFTGrid::RANDOMACCESS);
      float value1 = std::max<float>(0,histogram_->getRealValue(j1,k1,l1));
      float value2 = std::max<float>(0,histogram_->getRealValue(j1,k1,l2));
      float value3 = std::max<float>(0,histogram_->getRealValue(j1,k2,l1));
      float value4 = std::max<float>(0,histogram_->getRealValue(j1,k2,l2));
      float value5 = std::max<float>(0,histogram_->getRealValue(j2,k1,l1));
      float value6 = std::max<float>(0,histogram_->getRealValue(j2,k1,l2));
      float value7 = std::max<float>(0,histogram_->getRealValue(j2,k2,l1));
      float value8 = std::max<float>(0,histogram_->getRealValue(j2,k2,l2));
      histogram_->endAccess();

      value += (1.0f-wj)*(1.0f-wk)*(1.0f-wl)*value1;
      value += (1.0f-wj)*(1.0f-wk)*(     wl)*value2;
      value += (1.0f-wj)*(     wk)*(1.0f-wl)*value3;
      value += (1.0f-wj)*(     wk)*(     wl)*value4;
      value += (     wj)*(1.0f-wk)*(1.0f-wl)*value5;
      value += (     wj)*(1.0f-wk)*(     wl)*value6;
      value += (     wj)*(     wk)*(1.0f-wl)*value7;
      value += (     wj)*(     wk)*(     wl)*value8;
  }

  return value;
}
