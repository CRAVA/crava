#ifndef POSTERIORELASTICPDF_H
#define POSTERIORELASTICPDF_H

#include <vector>
#include "src/simbox.h"
#include "src/fftgrid.h"

class Simbox;

// Abstract class for holding a posterior pdf for elastic parameters, either in 2D, 3D or 4D.
// Provides the density for the given set of values.
// Returns missing if values are outside definition area.

class PosteriorElasticPDF {

public:

  PosteriorElasticPDF() {}

  virtual ~PosteriorElasticPDF() {}

  virtual double Density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1,
                         const double & s2) const = 0;

  virtual double Density(const double & vp,
                         const double & vs,
                         const double & rho,
                         const double & s1) const = 0;

  virtual double Density(const double & vp,
                         const double & vs,
                         const double & rho) const = 0;

  virtual double FindDensity(const double & vp,
                             const double & vs,
                             const double & rho,
                             const double & s1 = 0,
                             const double & s2 = 0,
                             const Simbox * const volume = 0) const = 0;

  virtual void  ResampleAndWriteDensity(const std::string & fileName,
                                        const Simbox      * origVol,
                                        Simbox            * volume,
                                        int                 gridNo,
                                        bool                writeSurface) const = 0;

  //virtual void WriteAsciiFile(std::string filename) const;

protected:

  void SolveGEVProblem(double               ** sigmaPrior,       //Covariance matrix prior model
                       double               ** sigmaPosterior,   //Covariance matrix posterior model
                       std::vector<double>   & v1,              //Linear transformation 1
                       std::vector<double>   & v2);             //Linear transformation 2

  void CalculateVariance2D(double                    ** sigmaSmooth,  //the smoothing 3D covariance matrix
                           double                    ** sigma2D,      //the resulting 2D covariance matrix
                           const std::vector<double>    v1,      //linear transformation 1
                           const std::vector<double>    v2);     //linear transformation 2

  void CalculateTransform2D(const std::vector<double>                & d1,
                            const std::vector<double>                & d2,
                            const std::vector<double>                & d3,
                            std::vector<std::vector<double> >        & x, // either 2 x 3 or 3 x 3
                            const std::vector<std::vector<double> >  & v);

  void InvertSquareMatrix(double                     ** matrix,     //matrix to be inverted
                          double                     ** invMatrix,  //inverted matrix
                          int                           n);         // size

  void SetupSmoothingGaussian2D(FFTGrid                   * smoother,
                                const double   *const*const sigma_inv,
                                int                         n1,
                                int                         n2,
                                int                         n3,
                                double                      dx,
                                double                      dy);
};

#endif
