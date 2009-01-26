#ifndef KRIGING2D_H
#define KRIGING2D_H

#include "src/covgrid2d.h"
#include "nrlib/grid/grid2d.hpp"


class Kriging2D
{
public:
  static void      krigSurface(NRLib2::Grid2D<double> *trend, float *data, int *indexi, int *indexj, int md,       
                              CovGrid2D *cov);

private:
 // static void      locateValidData(float **data, 
 //                                  int *indexi, 
 //                                  int nx, 
 //                                  int *indexj, 
  //                                 int ny, 
  //                                 int &md);
    
  static void      subtractTrend(float * data,
                                 NRLib2::Grid2D<double> * trend,
                                 int   * indexi,
                                 int   * indexj,
                                 int     md);
 // static void      addTrend(float ** data,
 //                           float ** trend,
  //                          int     nx, 
  //                          int ny);
  static void      allocateSpaceForMatrixEq(double ** & K, 
                                            double ** & C,
                                            double  * & k,
                                            int         md);
  static void      deAllocateSpaceForMatrixEq(double ** & K, 
                                              double ** & C,
                                              double  * & k,
                                              int         md);
  static void      fillKrigingMatrix(double **K, 
                                     int     *indexi, 
                                     int     *indexj, 
                                     int     md, 
                                     CovGrid2D *cov);
  static void      cholesky(double ** K,
                            double ** C,
                            int       md);
  static double ** copyMatrix(double ** in,
                              double ** out,
                              int       md);
  static void      fillKrigingVector(double *k, 
                                     int *indexi,
                                     int *indexj,
                                     int md,
                                     int i, 
                                     int j, 
                                     CovGrid2D *cov);                             
                                    
};

#endif
