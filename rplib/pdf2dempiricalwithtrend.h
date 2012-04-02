#ifndef PDF2DEMPIRICALWITHTREND_H
#define PDF2DEMPIRICALWITHTREND_H

#include <vector>
#include <rplib/pdf3d.h>
//#include <rplib/pdf3dempirical.h>
#include "nrlib/trend/trend.hpp"

class Pdf3DEmpirical;

// Class for holding a 3D empirical pdf with two trends.


class Pdf2DEmpiricalWithTrend : public Pdf3D {
public:

  Pdf2DEmpiricalWithTrend(const std::vector<double> & d1,    // first dimension of datapoints
                 const std::vector<double> & d2,             // second dimension of datapoints
                 const std::vector<double> & d3,             // third dimension of datapoints
                 const std::vector<double> & t1,
                 const std::vector<double> & t2,
                 const NRLib::Trend2D & a1,
                 const NRLib::Trend2D & b1,
                 const NRLib::Trend2D & c1,
                 const NRLib::Trend2D & a2,
                 const NRLib::Trend2D & b2,
                 const NRLib::Trend2D & c2,
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

  virtual ~Pdf2DEmpiricalWithTrend();

  virtual double density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const;

private:

  // Three dimensional probability density functions based on a linear combination of the data points and the trend vectors
  Pdf3DEmpirical * density_pdf1_;
  Pdf3DEmpirical * density_pdf2_;

  // Local copies of dimension reduction coefficients
  const NRLib::Trend2D a1_;
  const NRLib::Trend2D b1_;
  const NRLib::Trend2D c1_;
  const NRLib::Trend2D a2_;
  const NRLib::Trend2D b2_;
  const NRLib::Trend2D c2_;

  // Limits of trend values
  double t1_min_;
  double t1_max_;
  double dt1_;
  double t2_min_;
  double t2_max_;
  double dt2_;

  // Scaling coefficients for normalization of probability density function
  NRLib::Grid2D<double> scaling_coefficients_;
};

#endif
