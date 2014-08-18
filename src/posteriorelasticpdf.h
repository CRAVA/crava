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

  void SolveGEVProblem(NRLib::Matrix           & sigma_prior,       //Covariance matrix prior model
                       NRLib::Matrix           & sigma_post,   //Covariance matrix posterior model
                       std::vector<double>     & v1,              //Linear transformation 1
                       std::vector<double>     & v2);             //Linear transformation 2

  void CalculateVariance2D(NRLib::Matrix              & sigma_smooth,  //the smoothing 3D covariance matrix
                           NRLib::Matrix              & sigma_2d,      //the resulting 2D covariance matrix
                           const std::vector<double>  & v1_std,      //linear transformation 1
                           const std::vector<double>  & v2_std);     //linear transformation 2

  void CalculateTransform2D(const std::vector<double>                & d1,
                            const std::vector<double>                & d2,
                            const std::vector<double>                & d3,
                            std::vector<std::vector<double> >        & x, // either 2 x 3 or 3 x 3
                            const NRLib::Matrix                      & v);

  void InvertSquareMatrix(NRLib::Matrix               & matrix,     //matrix to be inverted
                          NRLib::Matrix               & inv_matrix,  //inverted matrix
                          int                           n);         // size

  void SetupSmoothingGaussian2D(FFTGrid                   * smoother,
                                const NRLib::Matrix       & sigma_inv,
                                int                         n1,
                                int                         n2,
                                int                         n3,
                                double                      dx,
                                double                      dy);
};

#endif
