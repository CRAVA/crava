/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/


#ifndef KRIGING1D_H
#define KRIGING1D_H

class Kriging1D
{
public:
  static void      krigVector(float * data,
                              float * trend,
                              int     nd,
                              float   dz);

  static void      krigVector(double * data,
                              float  * trend,
                              int      nd,
                              float    dz);

private:
  static void      locateValidData(float * data,
                                   int   * index,
                                   int     nd,
                                   int   & md);

  static void      locateValidData(double * data,
                                   int    * index,
                                   int      nd,
                                   int    & md);

  static void      subtractTrend(float * data,
                                 float * trend,
                                 int   * index,
                                 int     md);

  static void      subtractTrend(double * data,
                                 float  * trend,
                                 int    * index,
                                 int      md);

  static void      addTrend(float * data,
                            float * trend,
                            int     nd);

  static void      addTrend(double * data,
                            float  * trend,
                            int      nd);

  static void      allocateSpaceForMatrixEq(double ** & K,
                                            double ** & C,
                                            double  * & k,
                                            int         md);
  static void      deAllocateSpaceForMatrixEq(double ** & K,
                                              double ** & C,
                                              double  * & k,
                                              int         md);
  static void      fillKrigingMatrix(double ** K,
                                     int     * index,
                                     int       md,
                                     double    range,
                                     double    power,
                                     double    dz);
  static void      cholesky(double ** K,
                            double ** C,
                            int       md);
  static double ** copyMatrix(double ** in,
                              double ** out,
                              int       md);
  static void      fillKrigingVector(double * k,
                                     int    * index,
                                     int      md,
                                     double   range,
                                     double   power,
                                     double   dz,
                                     int      krigK);
};

#endif
