/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <math.h>

#include "src/definitions.h"
#include "src/kriging2d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"
#include "src/covgrid2d.h"
#include "src/simbox.h"

void Kriging2D::krigSurface(Grid2D              & trend,
                            const KrigingData2D & krigingData,
                            const CovGrid2D     & cov,
                            bool                  getResiduals)
{
  //
  // This routine by default returns z(x) = m(x) + k(x)K^{-1}(d - m). If only
  // residuals are wanted a copy of the input trend
  //
  int md = krigingData.getNumberOfData();
  const std::vector<int> & indexi = krigingData.getIndexI();
  const std::vector<int> & indexj = krigingData.getIndexJ();
  std::vector<float>       data   = krigingData.getData();   // Take an editable copy

  int nx = static_cast<int>(trend.GetNI());
  int ny = static_cast<int>(trend.GetNJ());

  if (md > 0 && md <= nx*ny) {

    NRLib::SymmetricMatrix K(md);
    NRLib::Vector residual(md);
    NRLib::Vector k(md);
    NRLib::Vector x(md);

    subtractTrend(residual, data, trend, indexi, indexj);


    Grid2D              filled(nx,ny,0);

    for(int i=0;i<md;i++){
      if (getResiduals) {  // Only get the residuals
        trend(indexi[i],indexj[i]) = residual(i);
        filled(indexi[i],indexj[i])=1.0;
      }
      else {
        trend(indexi[i],indexj[i]) += residual(i);
        filled(indexi[i],indexj[i])=1.0;
      }
    }
    bool first=true;
    for (int i = 0 ; i < nx ; i++) {
      for (int j = 0 ; j < ny ; j++) {
        if(!(filled(i,j) > 0.0)) // if this is not a datapoint
        {
          if(first)
          {
             fillKrigingMatrix(K, cov, indexi, indexj);
             NRLib::CholeskySolve(K, residual, x);
             first = false;
          }
          fillKrigingVector(k, cov, indexi, indexj, i, j);

          if (getResiduals) {  // Only get the residuals
            trend(i,j) = k * x;
          }
          else {
            trend(i,j) += k * x;
          }
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

CovGrid2D &
Kriging2D::makeCovGrid2D(const Simbox * simbox,
                         Vario        * vario,
                         int            debugFlag)
{
  //
  // Pretabulate all needed covariances
  //
  const int    nx = simbox->getnx();
  const int    ny = simbox->getny();

  const float  dx = static_cast<float>(simbox->getdx());
  const float  dy = static_cast<float>(simbox->getdy());

  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);

  if(debugFlag == 1) {
    std::string baseName = IO::PrefixBackground() + "covGrid2D" + IO::SuffixAsciiIrapClassic();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    cov->writeToFile(fileName);
  }
  return (*cov);
}


