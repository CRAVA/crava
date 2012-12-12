/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>

#include "src/definitions.h"
#include "src/kriging2d.h"
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

  if (md > 0 && md < nx*ny) {

    NRLib::SymmetricMatrix K(md);
    NRLib::Vector residual(md);
    NRLib::Vector k(md);
    NRLib::Vector x(md);

    subtractTrend(residual, data, trend, indexi, indexj);

    fillKrigingMatrix(K, cov, indexi, indexj);

    for (int i = 0 ; i < nx ; i++) {
      for (int j = 0 ; j < ny ; j++) {
        fillKrigingVector(k, cov, indexi, indexj, i, j);

        NRLib::CholeskySolve(K, k, x);

        if (getResiduals) {
          trend(i,j) = x * residual;
        }
        else {
          trend(i,j) += x * residual;
        }
      }
    }
  }
}

void
Kriging2D::subtractTrend(NRLib::Vector            & residual,
                         const std::vector<float> & data,
                         const Grid2D             & trend,
                         const std::vector<int>   & indexi,
                         const std::vector<int>   & indexj)
{
  int md = residual.length();
  for (int i = 0 ; i < md ; i++)
    residual(i) = data[i] - static_cast<float>(trend(indexi[i],indexj[i]));

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"\nData vector after trend subtraction:\n");
    for (int i = 0 ; i < md ; i++) {
      LogKit::LogFormatted(LogKit::Low," i indexi[i] indexj[i] residual : %3d %3d %3d  %.5f\n",i,indexi[i],indexj[i],residual(i));
    }
  }
}

void
Kriging2D::fillKrigingMatrix(NRLib::SymmetricMatrix  & K,
                             const CovGrid2D         & cov,
                             const std::vector<int>  & indexi,
                             const std::vector<int>  & indexj)
{
  int n = K.dim();
  for(int i=0 ; i < n  ; i++) {
    for(int j=0 ; j <= i ; j++) {
      int deltai = indexi[i] - indexi[j];
      int deltaj = indexj[i] - indexj[j];
      K(j,i) = static_cast<double>(cov.getCov(deltai,deltaj));
    }
  }
}

void
Kriging2D::fillKrigingVector(NRLib::Vector          & k,
                             const CovGrid2D        & cov,
                             const std::vector<int> & indexi,
                             const std::vector<int> & indexj,
                             int i, int j)
{
  for(int ii=0 ; ii < k.length() ; ii++) {
    int deltai = indexi[ii] - i;
    int deltaj = indexj[ii] - j;
    k(ii) = static_cast<double>(cov.getCov(deltai,deltaj));
  }
}

