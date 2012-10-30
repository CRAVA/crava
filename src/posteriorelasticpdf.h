#ifndef POSTERIORELASTICPDF_H
#define POSTERIORELASTICPDF_H

#include <vector>
#include "src/simbox.h"

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
                           double                    ** sigma2D,                 //the resulting 2D covariance matrix
                           const std::vector<double>    v1,      //linear transformation 1
                           const std::vector<double>    v2);     //linear transformation 2

  void CalculateTransform2D(const std::vector<double>   d1,     //data vector 1
                            const std::vector<double>   d2,     //data vector 2
                            const std::vector<double>   d3,     //data vector 3
                            std::vector<double>         x,      //output dimension 1
                            std::vector<double>         y,      //output dimension 2
                            const std::vector<double>   v1,     //linear transformation 1
                            const std::vector<double>   v2);    //linear transformation 2

  void InvertSquareMatrix(double                     ** matrix,     //matrix to be inverted
                          double                     ** invMatrix,  //inverted matrix
                          int                           n);         // size
};

#endif
