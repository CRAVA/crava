/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef KRIGING2D_H
#define KRIGING2D_H

#include "src/definitions.h"
#include "src/covgrid2d.h"
#include "src/krigingdata2d.h"
#include "nrlib/grid/grid2d.hpp"
#include "nrlib/flens/nrlib_flens.hpp"

class Kriging2D
{
public:
  static void  krigSurface(Grid2D              & trend,
                           const KrigingData2D & krigingData,
                           const CovGrid2D     & cov,
                           bool                  getResiduals = false);

private:
  static void  subtractTrend(NRLib::Vector            & d,
                             const std::vector<float> & data,
                             const Grid2D             & trend,
                             const std::vector<int>   & indexi,
                             const std::vector<int>   & indexj);

  static void  fillKrigingMatrix(NRLib::SymmetricMatrix  & K,
                                 const CovGrid2D         & cov,
                                 const std::vector<int>  & indexi,
                                 const std::vector<int>  & indexj);

  static void  fillKrigingVector(NRLib::Vector          & k,
                                 const CovGrid2D        & cov,
                                 const std::vector<int> & indexi,
                                 const std::vector<int> & indexj,
                                 int i,
                                 int j);
};
#endif
