#ifndef POSTERIORELASTICPDF2D_H
#define POSTERIORELASTICPDF2D_H

#include <vector>
#include <src/posteriorelasticpdf.h>
#include "src/simbox.h"

class FFTGrid;

// Class for holding a 2D posterior elastic PDF.

// The 2D empirical PDF is found by transforming the variables according to the input vectors and sampling the data points in
// bins in a grid with a different resolution (the resolution is here given by the input parameters n1 and n2).
// Denote this grid a density grid. We count the number of data points that are placed in each cell in the
// density grid. Next, a Gaussian 2D kernel is created to smooth the density grid. The smoothing
// is carried out in the Fourier domain, meaning we transform both the Gaussian 2D kernel and the
// density grid to the Fourier domain and carry out a multiplication, corresponding to a
// convolution of the kernel with each point in the density grid.

class PosteriorElasticPDF2D : public PosteriorElasticPDF {
public:

  // Constructor with given dimension reduction vectors v1 and v2
  PosteriorElasticPDF2D(const std::vector<double> & d1,    // first dimension of data points
                        const std::vector<double> & d2,    // second dimension of data points
                        const std::vector<double> & d3,    // third dimension of data points
                        int                         n1,    // resolution of density grid in dimension 1
                        int                         n2,    // resolution of density grid in dimension 2
                        double                   ** smoothvar,  // Gaussian smoothing kernel in 3D
                        const std::vector<double> & v1,    // Linear combination 1 of d1, d2 and d3
                        const std::vector<double> & v2);     // Linear combination 2 of d1, d2 and d3

  // Constructor with automatic dimension reduction
  PosteriorElasticPDF2D(const std::vector<double> & d1,    // first dimension of data points
                 const std::vector<double>        & d2,    // second dimension of data points
                 const std::vector<double>        & d3,    // third dimension of data points
                 int                                n1,    // resolution of density grid in dimension 1
                 int                                n2,    // resolution of density grid in dimension 2
                 double                          ** sigma_prior, //Covariance matrix, prior model
                 double                          ** sigma_posterior);  // Covariance matrix, posterior model

  PosteriorElasticPDF2D(int n1, int n2);

  ~PosteriorElasticPDF2D();

  // Returns missing if values are outside definition area.
  virtual double Density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const;

  virtual double Density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1) const;

  // Call to density function without dummy trend input variables
  virtual double Density(const double & vp,
                 const double & vs,
                 const double & rho) const;

  //virtual void WriteAsciiFile(std::string filename) const;

  virtual void ResampleAndWriteDensity(const std::string & fileName,
                                       const Simbox      * origVol,
                                       Simbox            * volume,
                                       int                 gridNo,
                                       bool                writeSurface) const;

  // Get functions for the limits of the variables
  double GetXMin() const { return x_min_; }
  double GetXMax() const { return x_max_; }
  double GetdX()   const { return dx_; }

  double GetYMin() const { return y_min_; }
  double GetYMax() const { return y_max_; }
  double GetdY()   const { return dy_; }

private:

  void SetupSmoothingGaussian2D(FFTGrid * smoother,      //2D Gaussian
                                double ** smoothingVar,  //
                                int       n1,            //Grid size n1
                                int       n2,
                                double    dx,
                                double    dy);

  FFTGrid * histogram_;   // Density grid of size n1_ \times n2_ \times 1

  std::vector<double> v1_; //Transform 1
  std::vector<double> v2_; //Transform 2

  int n1_;                // Grid resolution for each variable.
  int n2_;

  double x_min_; //Limits for variable 1
  double x_max_;
  double dx_;    //dv1 = (v1_max-v1_min)/n1. To find index, take floor((v-v_min)/dv)).

  double y_min_;  // Limits for variable 2
  double y_max_;
  double dy_;

};

#endif
