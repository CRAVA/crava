/*#include <math.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
*/
//#include <iostream>
#include <fstream>
using std::ifstream;
#include <string>

#include "lib/global_def.h"
#include "lib/sgri.h"
#include "lib/lib_misc.h"
#include "lib/log.h"
#include "lib/fileio.hpp"
#include "lib/exception.hpp"
#include "lib/grid.hpp"

using namespace std;
using namespace NRLib2;

Sgri::Sgri(char * fileName, char *errText, int &errCode) :
  Grid<float>(0,0,0),
  dX_(0.0),
  dY_(0.0),
  dZ_(0.0),
  x0_(RMISSING),
  y0_(RMISSING),
  z0_(RMISSING),
  scaleX_(1.0),
  scaleY_(1.0),
  scaleZ_(1.0),
  hasComplex_(0),
  rotAngle_(0.0),
  dipAngle_(0.0)
{
  if (!readHeaderFile(fileName, errText, errCode))
    return;

  ifstream binFile(binFileName_, ios::in | ios::binary);
  if(!binFile) {
    if (errCode == 0) 
      sprintf(errText, "Error: Could not open file %s for reading.\n", binFileName_);
    else 
      sprintf(errText, "%sError: Could not open file %s for reading.\n", errText, binFileName_);
    errCode = 1;
    return;
  }
  binFile.close();

  readBinaryFile(GetN(), errText, errCode);

  return;
}

Sgri::~Sgri() 
{
  delete [] axisLabels_;
  delete [] gridLabels_;
}

bool
Sgri::readHeaderFile(char * fileName, char *errText, int &errCode)
{
  ifstream headerFile(fileName, ios::in); //checkFileOpen performed previously
  int i;
  bool retValue = true;
  string tmpStr;

  //Reading record 1: Version header
  getline(headerFile, tmpStr);
  //Reading record 2: Grid dimension
  headerFile >> dim_;
  if (dim_ < 1 || dim_ > 3) {
    if (errCode == 0)
      sprintf(errText,"Error: Dimension of wavelet in file %s must be 1,2 or 3.\n", fileName);
    else 
      sprintf(errText,"%sError: Dimension of wavelet in file %s must be 1,2 or 3.\n", errText,fileName);
    errCode = 2;
    return (false);
  }
  else if (dim_ == 1) {
    LogKit::writeLog("\nWarning  !!!\n");
    LogKit::writeLog("Dimension for wavelet set to 1. Z-direction for 1D wavelet grid chosen.\n");
    LogKit::writeLog("\n");
  }
  getline(headerFile, tmpStr);
  //Reading record 3 ... 3+dim: Axis labels + grid value label
  axisLabels_ = new string[dim_];
  for (i=0; i<dim_; i++)
    getline(headerFile, axisLabels_[i]);
  getline(headerFile, gridValLabel_);
  int config;
  enum Planes {XZ, ZX, YZ, ZY};
  if (dim_ == 2) { //Must find the right axis configuration
    bool map[2][3];
    string dir[6] = {"X", "Y", "Z", "x", "y", "z"};
    for (i=0; i<2; i++) {
      int nHits = 0;
      for (int j=0; j<3; j++) {
        if ((axisLabels_[i].find(dir[j]) != string::npos) || (axisLabels_[i].find(dir[j+3]) != string::npos)) {
          map[i][j] = true;
          nHits++;
        }
        else
          map[i][j] = false;
      }
      if (nHits == 0) {
        if (errCode == 0)  
          sprintf(errText,"Error: Axis label %d does not contain either X, Y or Z. Must have exactly one direction specified.\n", i);
        else 
          sprintf(errText,"%sError: Axis label %d does not contain either X, Y or Z. Must have exactly one direction specified.\n", errText, i);
        errCode=2;
        return(false);
      }
      else if (nHits > 1) {
        if (errCode == 0)  
          sprintf(errText,"Error: Axis label %d contains at least to of the letters X, Y and Z. Must have exactly one direction specified.\n", i);
        else 
          sprintf(errText,"%sError: Axis label %d contains at least to of the letters X, Y and Z. Must have exactly one direction specified.\n", errText, i);
        errCode=2;
        return(false);
      }
    }
    if ((map[0][0] && map[1][0]) || (map[0][0] && map[1][1]) || (map[0][1] && map[1][0]) || (map[0][1] && map[1][1]) || (map[0][2] && map[1][2])) {
      if (errCode == 0)  
        sprintf(errText,"Error: Axis configuration can not be XX, YY, ZZ, XY or YX.\n");
      else 
        sprintf(errText,"%sError: Axis configuration can not be XX, YY, ZZ, XY or YX.\n", errText, i);
      errCode=2;
      return(false);
    }
    else {
      if (map[0][0] && map[1][2])
        config = XZ;
      else if (map[0][2] && map[1][0])
        config = ZX;
      else if (map[0][1] && map[1][0])
        config = YZ;
      else
        config = ZY;
    }
  }
  //Reading record 4+dim: Number of grids
  headerFile >> nGrid_;
  if (nGrid_ < 1) {
    if (errCode == 0)  
      sprintf(errText,"Error: Number of grids read from waveletfile %s must be > 0.\n", fileName);
    else 
      sprintf(errText,"%sError: Number of grids read from waveletfile %s must be > 0.\n", errText,fileName);
    errCode=2;
    return (false);
  }
  getline(headerFile, tmpStr);
  //Reading record 5+dim ... 5+dim+ngrid-1: Grid labels
  gridLabels_ = new string[nGrid_];
  for (i=0; i<nGrid_; i++)
    getline(headerFile, gridLabels_[i]);
  double *dValues1 = new double[dim_];
  double *dValues2 = new double[dim_];
  double *dValues3 = new double[dim_];
  int *iValues = new int[dim_];
  //Reading record 5+dim+ngrid: Scaling factor of grid values
  for (i=0; i<dim_; i++)
    headerFile >> dValues1[i];
  getline(headerFile,tmpStr);
  //Reading record 6+dim+ngrid: Number of samples in each dir.
  for (i=0; i<dim_; i++)
    headerFile >> iValues[i];
  getline(headerFile,tmpStr);
  //Reading record 7+dim+ngrid: Grid sampling in each dir.
  for (i=0; i<dim_; i++)
    headerFile >> dValues2[i];
  getline(headerFile,tmpStr);
  //Reading record 8+dim+ngrid: First point coord.
  for (i=0; i<dim_; i++)
    headerFile >> dValues3[i];
  getline(headerFile, tmpStr);

  int nX=1;
  int nY=1; 
  int nZ=1;
  bool xActive = false;
  bool yActive = false;
  if (dim_ == 3) {
    scaleX_ = dValues1[0];
    nX		  = iValues[0];
    dX_		  = dValues2[0];
    x0_		  = dValues3[0];
    scaleY_ = dValues1[1];
    nY		  = iValues[1];
    dY_		  = dValues2[1];
    y0_		  = dValues3[1];
    scaleZ_ = dValues1[2];
    nZ		  = iValues[2];
    dZ_		  = dValues2[2];
    z0_		  = dValues3[2];
  }
  else if (dim_ == 2) { 
    if (config == XZ) {
      xActive = true;
      scaleX_ = dValues1[0];
      nX		  = iValues[0];
      dX_		  = dValues2[0];
      x0_		  = dValues3[0];
      scaleZ_ = dValues1[1];
      nZ      = iValues[1];
      dZ_     = dValues2[1];
      z0_     = dValues3[1];
    }
    if (config == ZX) {
      xActive = true;
      scaleZ_ = dValues1[0];
      nZ		  = iValues[0];
      dZ_		  = dValues2[0];
      z0_		  = dValues3[0];
      scaleX_ = dValues1[1];
      nX      = iValues[1];
      dX_     = dValues2[1];
      x0_     = dValues3[1];
    }
    if (config == YZ) {
      yActive = true;
      scaleY_ = dValues1[0];
      nY		  = iValues[0];
      dY_		  = dValues2[0];
      y0_		  = dValues3[0];
      scaleZ_ = dValues1[1];
      nZ      = iValues[1];
      dZ_     = dValues2[1];
      z0_     = dValues3[1];
    }
    if (config == ZY) {
      yActive = true;
      scaleZ_ = dValues1[0];
      nZ		  = iValues[0];
      dZ_		  = dValues2[0];
      z0_		  = dValues3[0];
      scaleY_ = dValues1[1];
      nY      = iValues[1];
      dY_     = dValues2[1];
      y0_     = dValues3[1];
    }
  }
  else { //dim_ = 1 - currently not an active option
    scaleZ_ = dValues1[0];
    nZ		= iValues[0];
    dZ_		= dValues2[0];
    z0_		= dValues3[0];
  }
  if (nX < 1) {
    if (errCode == 0)  
      sprintf(errText,"Error: Number of samples in X-dir in waveletfile %s must be >= 1.\n", fileName);
    else 
      sprintf(errText,"%sError: Number of samples in X-dir in waveletfile %s must be >= 1.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  if (nY < 1) {
    if (errCode == 0)  
      sprintf(errText,"Error: Number of samples in Y-dir in waveletfile %s must be >= 1.\n", fileName);
    else 
      sprintf(errText,"%sError: Number of samples in Y-dir in waveletfile %s must be >= 1.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  if (nZ < 1) {
    if (errCode == 0)  
      sprintf(errText,"Error: Number of samples in Z-dir in waveletfile %s must be >= 1.\n", fileName);
    else 
      sprintf(errText,"%sError: Number of samples in Z-dir in waveletfile %s must be >= 1.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  if ((xActive) && (dX_ <= 0.0)) {
    if (errCode == 0)  
      sprintf(errText,"Error: Grid sampling in X-dir in waveletfile %s must be > 0.0.\n", fileName);
    else 
      sprintf(errText,"%sError: Grid sampling in X-dir in waveletfile %s must be > 0.0.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  if ((yActive) && (dY_ <= 0.0)) {
    if (errCode == 0)  
      sprintf(errText,"Error: Grid sampling in Y-dir in waveletfile %s must be > 0.0.\n", fileName);
    else 
      sprintf(errText,"%sError: Grid sampling in Y-dir in waveletfile %s must be > 0.0.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }  
  if (dZ_ <= 0.0) {
    if (errCode == 0)  
      sprintf(errText,"Error: Grid sampling in Z-dir in waveletfile %s must be > 0.0.\n", fileName);
    else 
      sprintf(errText,"%sError: Grid sampling in Z-dir in waveletfile %s must be > 0.0.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  delete [] dValues1;
  delete [] dValues2;
  delete [] dValues3;
  delete [] iValues;
  Resize(nX,nY,nZ);
  //Reading record 9+dim+ngrid: Angle of rotation
  if (dim_ > 1) {
    headerFile >> dipAngle_;
    if (dim_ > 2)
      headerFile >> rotAngle_;
  }
  getline(headerFile, tmpStr);
  //Reading record 10+dim+ngrid: Undef value
  headerFile >> unDef_;
  getline(headerFile, tmpStr);
  //Reading record 11+dim+ngrid: Filename of binary file
  headerFile >> binFileName_;
  getline(headerFile, tmpStr);
  //Reading record 12+dim+ngrid: Complex values
  headerFile >> hasComplex_;
  if (hasComplex_ != 0 && hasComplex_ != 1) {
    if (errCode == 0)  
      sprintf(errText,"Error: Code for complex grid in waveletfile %s must be 0 or 1.\n", fileName);
    else 
      sprintf(errText,"%sError: Code for complex grid in waveletfile %s must be 0 or 1.\n", errText,fileName);
    errCode=2;
    retValue = false;
  }
  //The remaining records are not relevant, stop reading

  return (retValue);
}

void
Sgri::readBinaryFile(int n, char *errText, int &errCode)
{
  ifstream binFile(binFileName_,ios::in | ios::binary); //Check opening of file before calling this function
 
  if (nGrid_ > 1) {
    LogKit::writeLog("\nWarning  !!!\n");
    LogKit::writeLog("%d %s", nGrid_, "grids read from Sgri-file. Only the first is used.\n");
    LogKit::writeLog("\n");
  }
  if (hasComplex_) {
    if (errCode == 0)  
      sprintf(errText,"Error: Complex grid given in waveletfile %s. Only handling of real grids implemented.\n", binFileName_);
    else 
      sprintf(errText,"%sError: Complex grid given in waveletfile %s. Only handling of real grids implemented.\n", errText, binFileName_);
    errCode=2;
    return;
  }
  
  try {
    ReadBinaryFloatArray(binFile, begin(), n);
  }
  catch (NRLib2::Exception& e) {
    if (errCode == 0)  
      sprintf(errText, "Error: Reading from binary waveletfile %s. %s\n", binFileName_, e.what().c_str());
    else 
      sprintf(errText, "%sError: Reading from binary waveletfile %s. %s\n", errText, binFileName_, e.what().c_str());
    errCode=2;
  }
  
  return;
}
