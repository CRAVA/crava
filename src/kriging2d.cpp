/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>

#include "src/definitions.h"
#include "src/kriging2d.h"
#include "lib/lib_matr.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"


void Kriging2D::krigSurface(Grid2D              & trend,
                            const KrigingData2D & krigingData,
                            const CovGrid2D     & cov,
                            bool                  getResiduals)
{
  //
  // This routine by default returns z(x) = m(x) + k(x)K^{-1}(d - m). If only
  // residuals are wanted a copy of the inpuit trend
  //
  int md = krigingData.getNumberOfData();
  const std::vector<int> & indexi = krigingData.getIndexI();
  const std::vector<int> & indexj = krigingData.getIndexJ();
  std::vector<float>       data   = krigingData.getData();   // Take an editable copy

  int nx = static_cast<int>(trend.GetNI());
  int ny = static_cast<int>(trend.GetNJ());
  if (md < nx*ny) {
    subtractTrend(data, trend, indexi, indexj, md);

    double ** K;  // Kriging matrix
    double ** C;  // Kriging matrix cholesky decomposed
    double *  k;  // Kriging vector

    allocateSpaceForMatrixEq(K, C, k, md);
    fillKrigingMatrix(K, cov, indexi, indexj, md);
    cholesky(K, C, md);

    for (int i = 0 ; i < nx ; i++)
      for (int j = 0 ; j < ny ; j++)
      {
        fillKrigingVector(k, cov, indexi, indexj, md, i, j);
        lib_matrAxeqbR(md, C, k); // solve kriging equation
        if (getResiduals)
          for (int ii = 0 ; ii < md ; ii++)
            trend(i,j)  = k[ii] * data[ii];
        else
          for (int ii = 0 ; ii < md ; ii++)
            trend(i,j) += k[ii] * data[ii];
      }
    deAllocateSpaceForMatrixEq(K, C, k, md);
  }
}

void
Kriging2D::subtractTrend(std::vector<float>     & data,
                         const Grid2D           & trend,
                         const std::vector<int> & indexi,
                         const std::vector<int> & indexj,
                         int                      md)
{
  for (int i = 0 ; i < md ; i++)
    data[i] -= static_cast<float>(trend(indexi[i],indexj[i]));

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"\nData vector after trend subtraction:\n");
    for (int i = 0 ; i < md ; i++) {
      LogKit::LogFormatted(LogKit::Low," i indexi[i] indexj[i] data : %3d %3d %3d  %.5f\n",i,indexi[i],indexj[i],data[i]);
    }
  }
}

void
Kriging2D::fillKrigingMatrix(double                 ** K,
                             const CovGrid2D         & cov,
                             const std::vector<int>  & indexi,
                             const std::vector<int>  & indexj,
                             int                       md)
{
  for(int i=0;i<md;i++)
    for(int j=0;j<md;j++)
    {
      int deltai = indexi[i] - indexi[j];
      int deltaj = indexj[i] - indexj[j];
      K[i][j] = static_cast<double>(cov.getCov(deltai,deltaj));
    }
}

void
Kriging2D::fillKrigingVector(double                 * k,
                             const CovGrid2D        & cov,
                             const std::vector<int> & indexi,
                             const std::vector<int> & indexj,
                             int md, int i, int j)
{
  for(int ii=0;ii<md;ii++)
  {
    int deltai = indexi[ii]-i;
    int deltaj = indexj[ii]-j;
    k[ii] = static_cast<double>(cov.getCov(deltai,deltaj));
  }
}

void
Kriging2D::allocateSpaceForMatrixEq(double ** & K,
                                    double ** & C,
                                    double  * & k,
                                    int         md)
{
  K = new double * [md];
  C = new double * [md];
  k = new double[md];
  for (int i = 0 ; i < md ; i++) {
    K[i] = new double[md];
    C[i] = new double[md];
  }
}

void
Kriging2D::deAllocateSpaceForMatrixEq(double ** K,
                                      double ** C,
                                      double  * k,
                                      int       md)
{
  for (int i = 0 ; i < md ; i++) {
    delete [] K[i];
    delete [] C[i];
  }
  delete [] K;
  delete [] C;
  delete [] k;
  K = NULL;
  C = NULL;
  k = NULL;
}

void
Kriging2D::cholesky(double ** K,
                    double ** C,
                    int       md)
{
  static const double choleskyRepairFactor = 1.001;
  int count = 0;
  while ( lib_matrCholR(md, copyMatrix(K, C, md)) ) {
    for (int i = 0 ; i < md ; i++)
      K[i][i] *= choleskyRepairFactor;
    count++;
    if (count > 5)
      LogKit::LogFormatted(LogKit::Low,"\nERROR in Kriging1D::Cholesky(): Could not find cholesky factor\n");
  }
}

double **
Kriging2D::copyMatrix(double ** in,
                      double ** out,
                      int       md)
{
  for (int i = 0 ; i < md ; i++) {
    for (int j = 0 ; j < md ; j++)
      out[i][j] = in[i][j];
  }
  return out;
}

