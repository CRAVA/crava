#include "rplib/pdf3dempirical.h"
#include "src/fftgrid.h"
#include <math.h>
#include "lib/lib_matr.h"


Pdf3DEmpirical::Pdf3DEmpirical(const std::vector<double> & d1,
                               const std::vector<double> & d2,
                               const std::vector<double> & d3,
                               int n1,
                               int n2,
                               int n3,
                               double smooth_var1,
                               double smooth_var2,
                               double smooth_var3,
                               double smooth_corr12,
                               double smooth_corr13,
                               double smooth_corr23)
                               :
n1_(n1),
n2_(n2),
n3_(n3)
{
// Create covariance matrix from the input smoothing variables
  double ** sigma = new double *[3];
  for(int i=0; i<3; i++)
    sigma[i] = new double[3];

  sigma[0][0] = smooth_var1;
  sigma[0][1] = smooth_corr12;
  sigma[0][2] = smooth_corr13;
  sigma[1][0] = smooth_corr12;
  sigma[1][1] = smooth_var2;
  sigma[1][2] = smooth_corr23;
  sigma[2][0] = smooth_corr13;
  sigma[2][1] = smooth_corr23;
  sigma[2][2] = smooth_var3;

  // Establish minimum and maximum limits of the input data
  // (Is it a point to do this before inverting sigma?)
  v1_min_ = *min_element(d1.begin(), d1.end()) - 5.0f*sqrt(sigma[0][0]);
  v1_max_ = *max_element(d1.begin(), d1.end()) + 5.0f*sqrt(sigma[0][0]);
  v2_min_ = *min_element(d2.begin(), d2.end()) - 5.0f*sqrt(sigma[1][1]);
  v2_max_ = *max_element(d2.begin(), d2.end()) + 5.0f*sqrt(sigma[1][1]);
  v3_min_ = *min_element(d3.begin(), d3.end()) - 5.0f*sqrt(sigma[2][2]);
  v3_max_ = *max_element(d3.begin(), d3.end()) + 5.0f*sqrt(sigma[2][2]);

  // Spacing variables in the density grid
  dv1_ = (v1_max_ - v1_min_)/n1_;
  dv2_ = (v2_max_ - v2_min_)/n2_;
  dv3_ = (v3_max_ - v3_min_)/n3_;

  // Matrix inversion of the covariance matrix sigma
  double **sigmainv = new double *[3];
  for(int i=0; i<3; i++)
    sigmainv[i] = new double [3];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      if(i==j)
        sigmainv[i][j] = 1.0;
      else
        sigmainv[i][j] = 0.0;
  lib_matrCholR(3, sigma);
  lib_matrAXeqBMatR(3, sigma, sigmainv, 3);

  for(int i=0; i<3; i++)
    delete [] sigma[i];
  delete [] sigma;

  // Create density grid as a FFTGrid
  histogram_ = new FFTGrid(n1_, n2_, n3_, n1_, n2_, n3_);
  histogram_->fillInConstant(0); // initialize to zero

  // We assume in the for loop below that the size of the input datapoint vectors are the same
  assert(d1.size()==d2.size());
  assert(d2.size()==d3.size());


  // Go through data points and place in bins in histogram
  for (int ind = 0; ind < static_cast<int>(d1.size()); ind++){

    int i = static_cast<int>(floor((d1[ind]-v1_min_)/dv1_));
    int j = static_cast<int>(floor((d2[ind]-v2_min_)/dv2_));
    int k = static_cast<int>(floor((d3[ind]-v3_min_)/dv3_));

    // Counting data points in index (i,k,j)
    histogram_->setRealValue(i,j,k, histogram_->getRealValue(i,j,k) + 1);
  }

  // Set upt of Gaussian smoother
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
        smooth[j+k*n1_+l*n1_*n2_] = float(exp(-0.5f*(jj*dv1_ *jj*dv1_*sigmainv[0][0]
                                                   +kk*dv2_  *kk*dv2_*sigmainv[1][1]
                                                   +ll*dv3_  *ll*dv3_*sigmainv[2][2]
                                                   +2*jj*dv1_*kk*dv2_*sigmainv[1][0]
                                                   +2*jj*dv1_*ll*dv3_*sigmainv[2][0]
                                                   +2*kk*dv2_*ll*dv3_*sigmainv[2][1])));
        sum = sum+smooth[j+k*n1_+l*n1_*n2_];
      }
    }
  }
  // normalize smoother
  for(l=0;l<n3_;l++)
    for(k=0;k<n2_;k++)
      for(j=0;j<n1_;j++)
        smooth[j+k*n1_+l*n1_*n2_]/=sum;


  FFTGrid *smoother;
  smoother = new FFTGrid(n1_, n2_, n3_, n1_, n2_, n3_);
  smoother->fillInFromArray(smooth);
  smoother->fftInPlace();

  // Carry out multiplication of the smoother with the density grid (histogram) in the Fourier domain
  histogram_->fftInPlace();
  histogram_->multiply(smoother);
  histogram_->invFFTInPlace();
  histogram_->multiplyByScalar(float(sqrt(double(n1_*n2_*n3_))));

  delete smoother;
  delete [] smooth;
  for(int i=0;i<3;i++)
    delete [] sigmainv[i];
  delete [] sigmainv;
}

Pdf3DEmpirical::~Pdf3DEmpirical()
{
}

double
Pdf3DEmpirical::density(const double & vp,
                        const double & vs,
                        const double & rho,
                        const double & s1,
                        const double & s2) const
{
  // Disregarding trend variables s1 and s2
  return density(vp, vs, rho);

}

 // Call to density function without dummy trend input variables
double Pdf3DEmpirical::density(const double & vp,
                               const double & vs,
                               const double & rho) const
  {
  // If values are outside definition area.
  if(vp < v1_min_ || vp > v1_max_ || vs < v2_min_ || vs > v2_max_ || rho < v3_min_ || rho > v3_max_)
    return RMISSING;

  double returnvalue;

  // Find indices of input values and perform look-up in the grid
  int i = static_cast<int>(floor((vp - v1_min_)/dv1_));
  int j = static_cast<int>(floor((vs - v2_min_)/dv2_));
  int k = static_cast<int>(floor((rho- v3_min_)/dv3_));

  returnvalue = histogram_->getRealValue(i,j,k);
  return (returnvalue);
}
