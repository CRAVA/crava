#ifndef KRIGING2D_H
#define KRIGING2D_H

#include "src/definitions.h"
#include "src/covgrid2d.h"
#include "src/krigingdata2d.h"
#include "nrlib/grid/grid2d.hpp"

class Kriging2D
{
public:
  static void      krigSurface(Grid2D              & trend,
                               const KrigingData2D & krigingData,
                               const CovGrid2D     & cov,
                               bool                  getResiduals = false);

private:
  static void      subtractTrend(std::vector<float>     & data,
                                 const Grid2D           & trend,
                                 const std::vector<int> & indexi,
                                 const std::vector<int> & indexj,
                                 int                      md);
  static void      fillKrigingMatrix(double                 ** K,
                                     const CovGrid2D         & cov,
                                     const std::vector<int>  & indexi,
                                     const std::vector<int>  & indexj,
                                     int                       md);
  static void      fillKrigingVector(double                 * k,
                                     const CovGrid2D        & cov,
                                     const std::vector<int> & indexi,
                                     const std::vector<int> & indexj,
                                     int md,
                                     int i,
                                     int j);
  static void      allocateSpaceForMatrixEq(double ** & K,
                                            double ** & C,
                                            double  * & k,
                                            int         md);
  static void      deAllocateSpaceForMatrixEq(double ** K,
                                              double ** C,
                                              double  * k,
                                              int       md);
  static void      cholesky(double ** K,
                            double ** C,
                            int       md);
  static double ** copyMatrix(double ** in,
                              double ** out,
                              int       md);

};
#endif
