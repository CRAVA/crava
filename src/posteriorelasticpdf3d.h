#ifndef POSTERIORELASTICPDF2D_H
#define POSTERIORELASTICPDF2D_H


#include <src/posteriorelasticpdf.h>
#include "src/simbox.h"

class FFTGrid;

// Class for a 3D posterior PDF with either  (i) three elastic variables
// or (ii) a reduction from three elastic parameters and one trend variable to a 3D PDF.


class PosteriorElasticPDF3D : public PosteriorElasticPDF {
public:

  // (i) Constructor with no dimension reduction: input: three elastic parameters
  PosteriorElasticPDF3D(const std::vector<double>                   & d1,          // first dimension of data points
                        const std::vector<double>                   & d2,          // second dimension of data points
                        const std::vector<double>                   & d3,          // third dimension of data points
                        const double                   *const *const  smoothvar,   // Gaussian smoothing kernel in 3D
                        int                                           n1,          // resolution of density grid in dimension 1
                        int                                           n2,          // resolution of density grid in dimension 2
                        int                                           n3,          // resolution of density grid in dimension 3
                        double                                        d1_min,
                        double                                        d1_max,
                        double                                        d2_min,
                        double                                        d2_max,
                        double                                        d3_min,
                        double                                        d3_max,
                        int                                           ind = 0);

  // (ii) Constructor with dimension reduction: input: three elastic parameters and one trend variable
  PosteriorElasticPDF3D(const std::vector<double>                   & d1, // first dimension of data points
                        const std::vector<double>                   & d2, // second dimension of data points
                        const std::vector<double>                   & d3, // third dimension of data points
                        const std::vector<int>                      & t1, // fourth dimension (trend parameters)
                        const NRLib::Matrix                         & v,  // Transformation of elastic variables from 3D to 2D
                        const double                     *const*const sigma, // Gaussian smoothing kernel in 2D
                        int                                           n1,    // resolution of density grid in elastic dimension 1
                        int                                           n2,    // resolution of density grid in elastic dimension 2
                        int                                           nt,    // resolution of density grid in the trend dimension
                        double                                        d1_min,
                        double                                        d1_max,
                        double                                        d2_min,
                        double                                        d2_max,
                        double                                        t1_min,
                        double                                        t1_max,
                        int                                           ind);

  PosteriorElasticPDF3D(int                                    n1,
                        int                                    n2,
                        int                                    n3);


  ~PosteriorElasticPDF3D();

  // Returns missing if values are outside definition area.

  virtual double Density (const double  & vp,
                          const double  & vs,
                          const double  & rho,
                          const double  & s1,
                          const double  & s2) const;

  virtual double Density(const double   & vp,
                         const double   & vs,
                         const double   & rho,
                         const double   & s1) const;

  // Call to density function without dummy trend input variables
  virtual double Density(const double   & vp,
                         const double   & vs,
                         const double   & rho) const;

  virtual double FindDensity(const double & vp,
                             const double & vs,
                             const double & rho,
                             const double & s1 = 0,
                             const double & s2 = 0,
                             const Simbox * const volume = 0) const;

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

  double GetZMin() const { return z_min_; }
  double GetZMax() const { return z_max_; }
  double GetdZ()   const { return dz_; }

private:

  FFTGrid *histogram_;   // Density grid of size n1_ \times n2_ \times n3_
  std::vector<double> v1_; //Transform 1
  std::vector<double> v2_; //Transform 2
  int n1_;                // Grid resolution for each variable.
  int n2_;
  int n3_;
  double x_min_; //Limits for variable 1
  double x_max_;
  double dx_;
  double y_min_;  // Limits for variable 2
  double y_max_;
  double dy_;
  double z_min_;  // Limits for variable 3
  double z_max_;
  double dz_;

  void SetupSmoothingGaussian3D(FFTGrid                 * smoother,
                                const NRLib::Matrix     & sigmainv);

  void SetupSmoothingGaussian2D(FFTGrid    * smoother,
                                double    ** sigmainv,
                                int          n1,
                                int          n2,
                                int          n3,
                                double       dx,
                                double       dy,
                                double       dz);
};

#endif
