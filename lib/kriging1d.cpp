#include <math.h>

#include "src/definitions.h"
#include "lib/kriging1d.h"
#include "lib/lib_matr.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"

//-----------------------------------------------------------------------
void
Kriging1D::krigVector(float * data,
                      float * trend,
                      int     nd,
                      float   dz) 
{
  //
  // This class takes the incomplete vector 'data' and fills it using kriging.
  // 
  int * index = new int[nd];

  double range = 500.0;
  double power =   1.5;

  int md;
  locateValidData(data, index, nd, md);

  if (md < nd) {
    subtractTrend(data, trend, index, md);

    double ** K;  // Kriging matrix                     
    double ** C;  // Kriging matrix cholesky decomposed 
    double *  k;  // Kriging vector                     

    allocateSpaceForMatrixEq(K, C, k, md);
    fillKrigingMatrix(K, index, md, range, power, static_cast<double>(dz));
    cholesky(K, C, md);

    for (int krigK = 0 ; krigK < nd ; krigK++) {
      if (data[krigK] == RMISSING) {
        fillKrigingVector(k, index, md, range, power, static_cast<double>(dz), krigK);
        lib_matrAxeqbR(md, C, k); // solve kriging equation

        data[krigK] = 0.0f;
        for (int i = 0 ; i < md ; i++) {
          data[krigK] += static_cast<float>(k[i]) * data[index[i]];
        }
      }
    }
    deAllocateSpaceForMatrixEq(K, C, k, md);

    addTrend(data, trend, nd);

    bool debug = false;
    if (debug) {
      for (int i = 0 ; i < nd ; i++) {
        LogKit::LogFormatted(LogKit::LOW," i krigedData[i] : %3d  %.5f\n",i,data[i]);
      }
    }
  }
  delete [] index;
}

//-----------------------------------------------------------------------
void
Kriging1D::locateValidData(float * data,
                           int   * index,
                           int     nd,
                           int   & md) 
{
  int count = 0;
  for (int i = 0 ; i < nd ; i++) {
    if (data[i] != RMISSING) {
      index[count] = i;
      count++;
    }
  }
  if (count == 0) {
    LogKit::LogFormatted(LogKit::LOW,"\nWARNING in Kriging1D::locateDataIndices() : ");
    LogKit::LogFormatted(LogKit::LOW,"Only missing values found in data vector.\n");
  }
  md = count;

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::LOW,"\nData vector after trend subtraction:\n");
    for (int i = 0 ; i < md ; i++) {
      LogKit::LogFormatted(LogKit::LOW," i index[i]  data : %3d %3d  %.5f\n",i,index[i],data[index[i]]);
    }
  }
}

//-----------------------------------------------------------------------
void 
Kriging1D::subtractTrend(float * data,
                         float * trend,
                         int   * index,
                         int     md)
{  
  for (int i = 0 ; i < md ; i++) {
    data[index[i]] -= trend[index[i]];
  }
  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::LOW,"\nData vector after trend subtraction:\n");
    for (int i = 0 ; i < md ; i++) {
      LogKit::LogFormatted(LogKit::LOW," i index[i]  data : %3d %3d  %.5f\n",i,index[i],data[index[i]]);
    }
  }
}

//-----------------------------------------------------------------------
void 
Kriging1D::addTrend(float * data,
                    float * trend,
                    int     nd)
{  
  for (int i = 0 ; i < nd ; i++) {
    data[i] += trend[i];
  }
}

//-----------------------------------------------------------------------
void 
Kriging1D::allocateSpaceForMatrixEq(double ** & K, 
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

//-----------------------------------------------------------------------
void 
Kriging1D::deAllocateSpaceForMatrixEq(double ** & K, 
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

//-----------------------------------------------------------------------
void 
Kriging1D::fillKrigingMatrix(double ** K,
                             int     * index,
                             int       md,
                             double    range, 
                             double    power,
                             double    dz)
{
  for (int i = 0 ; i < md ; i++) {
    for (int j = 0 ; j < md ; j++) {
        double dist = abs(index[i] - index[j])*dz;
        double corr = exp(-3.0*pow(dist/range, power));
        K[i][j] = corr;
    }
  }
  bool debug = false;
  if (debug) {
    Utils::writeTitler("Kriging matrix");
    Utils::writeMatrix(K,md,md);
  }
}

//-----------------------------------------------------------------------
void
Kriging1D::cholesky(double ** K,
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

//-----------------------------------------------------------------------
double ** 
Kriging1D::copyMatrix(double ** in,
                      double ** out,
                      int       md)
{
  for (int i = 0 ; i < md ; i++) {
    for (int j = 0 ; j < md ; j++)
      out[i][j] = in[i][j];
  }
  return out;
}

//-----------------------------------------------------------------------
void 
Kriging1D::fillKrigingVector(double * k,                             
                             int    * index,
                             int      md,
                             double   range, 
                             double   power,
                             double   dz,
                             int      krigK)
{
  for (int i = 0 ; i < md ; i++) {
    double dist = abs(index[i] - krigK)*dz;
    double corr = exp(-3.0*pow(dist/range, power));
    k[i] = corr;
  }
  bool debug = false;
  if (debug) {
    Utils::writeTitler("Kriging vector");
    Utils::writeVector(k,md);
  }

}
