#ifndef PDF3DEMPIRICAL_H
#define PDF3DEMPIRICAL_H

#include <vector>
#include <rplib/pdf3d.h>

class FFTGrid;

// Class for holding a 3D empirical pdf.

// The 3D empirical PDF is found by sampling the 3D data points in bins in a grid with a lower
// resolution (the resolution is here given by the input parameters n1, n2 and n3). Denote this
// grid a density grid. We count the number of data points that are placed in each cell in the
// density grid. Next, a Gaussian 3D kernel is created to smooth the density grid. The smoothing
// is carried out in the Fourier domain, meaning we transform both the Gaussian 3D kernel and the
// density grid to the Fourier domain and carry out a multiplication, corresponding to a
// convolution of the kernel with each point in the density grid.

class Pdf3DEmpirical : public Pdf3D {
public:

  Pdf3DEmpirical(const std::vector<double> & d1,    // first dimension of datapoints
                 const std::vector<double> & d2,    // second dimension of datapoints
                 const std::vector<double> & d3,    // third dimension of datapoints
                 int n1,                            // resolution of density grid
                 int n2,                            //
                 int n3,                            //
                 double smooth_var1,                // Variance for the Gaussian smoothing kernel
                 double smooth_var2,
                 double smooth_var3,
                 double smooth_corr12,              // Covariance for the Gaussian smoothing kernel
                 double smooth_corr13,
                 double smooth_corr23
                 );
  virtual ~Pdf3DEmpirical();

  // Returns missing if values are outside definition area.
  virtual double density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const;

  // Call to density function without dummy trend input variables
  double density(const double & vp,
                 const double & vs,
                 const double & rho) const;

private:

  FFTGrid * histogram_;   // Density grid of size n1_ \times n2_ \times n3_

  int n1_;                // Grid resolution for each varaiable.
  int n2_;
  int n3_;

  double v1_min_; //Limits for variable 1
  double v1_max_;
  double dv1_;    //dv1 = (v1_max-v1_min)/n1. To find index, take floor((v-v_min)/dv)).

  double v2_min_;  // Limits for variable 2
  double v2_max_;
  double dv2_;

  double v3_min_;  // Limits for variable 3
  double v3_max_;
  double dv3_;
};

#endif
