#include <math.h>

#include "src/definitions.h"
#include "src/kriging2d.h"
#include "lib/global_def.h"
#include "lib/lib_matr.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"


void Kriging2D::krigSurface(NRLib2::Grid2D<double> * trend, 
                            const KrigingData2D    & krigingData2D,
                            CovGrid2D              * cov)
{
  int md = krigingData2D.getNumberOfData();
  int   * indexi = new int[md];
  int   * indexj = new int[md];
  float * data   = new float[md];
  
  for (int i = 0 ; i < md ; i++) 
  {
    indexi[i] = krigingData2D.getIndexI(i);
    indexj[i] = krigingData2D.getIndexJ(i);
    data[i]   = krigingData2D.getData(i);
  }

  int nx = trend->GetNI();
  int ny = trend->GetNJ();
  if (md < nx*ny) {
    subtractTrend(data, trend, indexi, indexj, md);

    double ** K;  // Kriging matrix                     
    double ** C;  // Kriging matrix cholesky decomposed 
    double *  k;  // Kriging vector                     

    allocateSpaceForMatrixEq(K, C, k, md);
    fillKrigingMatrix(K, indexi, indexj, md, cov);
    cholesky(K, C, md);

    for (int i = 0 ; i < nx ; i++) 
      for (int j = 0 ; j < ny ; j++)
      {
        fillKrigingVector(k, indexi,indexj,md,i,j, cov);
        lib_matrAxeqbR(md, C, k); // solve kriging equation
        for (int ii = 0 ; ii < md ; ii++) {
          (*trend)(i,j) += k[ii] * data[ii];
        }
      }

    deAllocateSpaceForMatrixEq(K, C, k, md);
  }
  delete [] indexi;
  delete [] indexj;
  delete [] data;
}

void 
Kriging2D::subtractTrend(float                  * data,
                         NRLib2::Grid2D<double> * trend,
                         int                    * indexi,
                         int                    * indexj,
                         int                      md)
{  
  for (int i = 0 ; i < md ; i++) 
    data[i] -= float((*trend)(indexi[i],indexj[i]));
  
  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::LOW,"\nData vector after trend subtraction:\n");
    for (int i = 0 ; i < md ; i++) {
      LogKit::LogFormatted(LogKit::LOW," i indexi[i] indexj[i] data : %3d %3d %3d  %.5f\n",i,indexi[i],indexj[i],data[i]);
    }
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
Kriging2D::deAllocateSpaceForMatrixEq(double ** & K, 
                                      double ** & C,
                                      double  * & k,
                                      int         md)
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
Kriging2D::fillKrigingMatrix(double **K, int *indexi, int *indexj, int md, CovGrid2D *cov)
{
  for(int i=0;i<md;i++)
    for(int j=0;j<md;j++)
    {
      int deltai = indexi[i] - indexi[j];
      int deltaj = indexj[i] - indexj[j];
      K[i][j] = cov->getCov(deltai,deltaj);
    }
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
      LogKit::LogFormatted(LogKit::LOW,"\nERROR in Kriging1D::Cholesky(): Could not find cholesky factor\n");
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

void 
Kriging2D::fillKrigingVector(double *k, int *indexi,int *indexj,int md,int i, int j, CovGrid2D *cov)
{
  int ii, deltai, deltaj;
  for(ii=0;ii<md;ii++)
  {
    deltai = indexi[ii]-i;
    deltaj = indexj[ii]-j;
    k[ii] = cov->getCov(deltai,deltaj);
  }
}
