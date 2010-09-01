#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "nrlib/iotools/logkit.hpp"

#include "src/covgridseparated.h"
#include "src/fftgrid.h"
#include "src/model.h"
#include "src/definitions.h"
#include "src/io.h"

CovGridSeparated::CovGridSeparated(const FFTGrid & grid)
{
  nxp_ = grid.getNxp(); 
  nyp_ = grid.getNyp(); 
  nzp_ = grid.getNzp();
  gammaXY_ = new float[nxp_*nyp_];
  gammaZ_  = new float[nzp_];

  tabulateCorr_ = true;

  FFTGrid& gridTmp = const_cast<FFTGrid&>(grid);
  gridTmp.setAccessMode(FFTGrid::RANDOMACCESS);
  bool isTrans = gridTmp.getIsTransformed();
  if (isTrans) 
    gridTmp.invFFTInPlace();

  for (int j = 0; j < nyp_; j++) {
    for (int i = 0; i < nxp_; i++) {
      int index = Get2DIndex(i, j);
      gammaXY_[index] = gridTmp.getRealValue(i, j, 0, true);
    } 
  } 
  float maxval = fabs(gridTmp.getRealValue(0, 0, 0, true));
  int i = 1;
  float value = fabs(gridTmp.getRealValue(0, 0, i, true));
  while(value>maxval)
  {
    maxval = fabs(gridTmp.getRealValue(0, 0, i, true));
    i++;
    value = fabs(gridTmp.getRealValue(0, 0, i, true));
  }
  maxval = gridTmp.getRealValue(0, 0, i-1, true);

  for (int k = 0; k < nzp_; k++) {
    gammaZ_[k] = gridTmp.getRealValue(0, 0, k, true) / maxval;
  }

  if (isTrans)
    gridTmp.fftInPlace();
  gridTmp.endAccess();
}

CovGridSeparated::CovGridSeparated(int nxp, int nyp, int nzp, 
                                   float dx, float dy, float dz,
                                   float rangeX, float rangeY, float rangeZ, 
                                   float power, float rotAngle, 
                                   bool tabulateCor) : 
  nxp_(nxp), nyp_(nyp), nzp_(nzp), tabulateCorr_(tabulateCor),
  dx_(dx), dy_(dy), dz_(dz),
  rangeX_(rangeX), rangeY_(rangeY), rangeZ_(rangeZ), power_(power), rotAngle_(rotAngle) {

  InitRotMatrix();
  
  if (!tabulateCor) {
    gammaXY_ = gammaZ_ = NULL;
    return;
  }
  gammaXY_ = new float[nxp_*nyp_];
  gammaZ_ = new float[nzp_];

  const int nzp2 = nzp/2; const int nyp2 = nyp/2; const int nxp2 = nxp/2;
  const int uppernxp2 = (nxp % 2 == 0 ? nxp2 : nxp2 + 1);
  const int uppernyp2 = (nyp % 2 == 0 ? nyp2 : nyp2 + 1);
  const int uppernzp2 = (nzp % 2 == 0 ? nzp2 : nzp2 + 1);  
  //LogKit::LogFormatted(LogKit::DebugLow,": %d %d %d %d\n",nxp,nxp2,nyp,nyp2);
  int i, j;
  int countxy = 0;
  for (j = -nyp2; j < uppernyp2; j++) {
      int j1 = (j < 0 ? nyp + j : j);
      for (i = -nxp2; i < uppernxp2; i++) {
        int i1 = (i < 0 ? nxp + i : i);
        float deltaX = i*dx;
        float deltaY = j*dy;
        float deltaZ = 0.f; 
        if (rotAngle != 0.f)
          RotateVec(deltaX, deltaY, deltaZ, rotMatrix_);
        float h = float(sqrt(deltaX*deltaX/(rangeX*rangeX) + deltaY*deltaY/(rangeY*rangeY)));
        gammaXY_[i1 + j1*nxp_] = float(exp(-3.0*pow(h,power)));
        countxy++;
      } // end i
    } //end j

  int k;
  int countz = 0;
  for (k = -nzp2; k < uppernzp2; k++) {
    int k1 = (k < 0 ? nzp + k : k);
    float deltaZ = k*dz;
    float h = float(sqrt(deltaZ*deltaZ/(rangeZ*rangeZ)));
    gammaZ_[k1] = float(exp(-3.0*pow(h,power)));
    countz++;
  }
  
  if (countxy != nxp_*nyp_ || countz != nzp_) {
    LogKit::LogFormatted(LogKit::Low,"ERROR in CovGridSeparated constructor.");
    exit(1);
  }
}

CovGridSeparated::CovGridSeparated(int nxp, int nyp, int nzp) : 
  nxp_(nxp), nyp_(nyp), nzp_(nzp), tabulateCorr_(true) {
  gammaXY_ = new float[nxp_*nyp_];
  gammaZ_ = new float[nzp_];

  int index;
  for (index = 0; index < nxp_*nyp_; index++)
    gammaXY_[index] = 0.f;
    
  for (index = 0; index < nzp_; index++) {
    gammaZ_[index] = 0.f;
  }
}

CovGridSeparated::~CovGridSeparated(void)
{
  if (gammaXY_)
    delete [] gammaXY_;

  if (gammaZ_)
    delete [] gammaZ_;
}

void CovGridSeparated::InitRotMatrix() {
  float cosA = float(cos(rotAngle_));
  float sinA = float(sin(rotAngle_));
  rotMatrix_[0][0] =  cosA; rotMatrix_[0][1] = sinA; rotMatrix_[0][2] = 0.0f;
  rotMatrix_[1][0] = -sinA; rotMatrix_[1][1] = cosA; rotMatrix_[1][2] = 0.0f;
  rotMatrix_[2][0] =  0.0f; rotMatrix_[2][1] = 0.0f; rotMatrix_[2][2] = 1.0f;
}

void CovGridSeparated::EstimateRanges(int& rangeX, int& rangeY, int& rangeZ) const {
  EstimateRangeX(rangeX);
  EstimateRangeY(rangeY);
  EstimateRangeZ(rangeZ);
}

void CovGridSeparated::EstimateRangeX(int& rangeX) const {
  const float f = 0.05f;
  int i,j;
  const float sillXY = f * gammaXY_[0];
  //const float sillXY = f;
  const int nxp2 = nxp_/2;
  const int nyp2 = nyp_/2;
  const int nyp = nyp_;

  for (i = nxp2 - 1; i >= 0; i--) {
    for (j = nyp2 - 1; j >= 0; j--) {
      int index = Get2DIndex(i, j);
      //      LogKit::LogFormatted(LogKit::DebugHigh,"%d ", index);
      float absGamma = float(fabs(gammaXY_[index]));
      if(absGamma >= sillXY) {
        rangeX = i;
        LogKit::LogFormatted(LogKit::DebugHigh,": %d\n",i);
        return;
      }
    } // end backward j

    for (j = nyp2; j < nyp; j++) {
      int index = Get2DIndex(i, j);
      LogKit::LogFormatted(LogKit::DebugHigh,"%d ", index);
      float absGamma = float(fabs(gammaXY_[index]));
      if(absGamma >= sillXY) {
        rangeX = i;
        LogKit::LogFormatted(LogKit::DebugHigh,": %d\n", i);
        return;
      }
    } // end forward j
  } // end i  
}

void CovGridSeparated::EstimateRangeY(int& rangeY) const {
  const float f = 0.05f;
  const float sillXY = f * gammaXY_[0];
  const int nxp2 = nxp_/2;
  const int nyp2 = nyp_/2;
  const int nxp = nxp_;

  for (int j = nyp2 - 1; j >= 0; j--) {
    for (int i = nxp2 - 1; i >= 0; i--) {
      int index = Get2DIndex(i, j);
      float absGamma = float(fabs(gammaXY_[index]));
      if (absGamma >= sillXY) {
        rangeY = j; return;
      }
    } // end backward i

    for (int i = nxp2; i < nxp; i++) {
      int index = Get2DIndex(i, j);
      float absGamma = float(fabs(gammaXY_[index]));
      if (absGamma >= sillXY) {
        rangeY = j; return;
      }
    } // end forward i
  } // end j
}

void CovGridSeparated::EstimateRangeZ(int& rangeZ) const {
  const float sillZ = 0.05f;

//  for (int k = nzp_/2; k >= 0; k--) {
  for(int k = 0;k<=nzp_/2;k++){
    float absGamma = float(fabs(gammaZ_[k]));
    //if (absGamma >= sillZ) {
    if(absGamma<=sillZ){
      rangeZ = k-1; return;
    }
  }
}


void CovGridSeparated::findTaperRanges(float & rangeX, float & rangeY, float & rangeZ) const
{
  float sum = 0;
  //Handle z-direction
  for (int k = 0; k < nzp_/2; k++)
    sum += float(fabs(gammaZ_[k]));
  rangeZ = 2.5f*sum;

  //X and Y ranges estimated simultaneously, by sweeping a 180 degree angle.
  float angle = 0;
  float dAngle = float(M_PI)/20.0f;
  rangeX = 0;
  rangeY = 0;
  while(angle < M_PI)
  {
    sum = 2.5f*angleSum(angle);
    if(sum*cos(angle) > rangeX)
      rangeX = float(sum*cos(angle));
    if(sum*sin(angle) > rangeY)
      rangeY = float(sum*sin(angle));
    angle += dAngle;
  }
}

float CovGridSeparated::angleSum(float angle) const
{
  float fi, fj, di, dj;
  float cosang = float(cos(angle));
  float sinang = float(sin(angle));
  if(fabs(cosang) > fabs(sinang))
  {
    di = 1;
    dj = sinang/cosang;
  }
  else
  {
    dj = cosang/sinang;
    di = 1;
  }

  double maxi = nxp_/2 + 1;
  double maxj = nyp_/2 + 1;
  int i, j, nUsed = 0;
  fi = 0;
  fj = 0;
  float sum = 0;
  while(fabs(fi) < maxi && fabs(fj) < maxj)
  {
    nUsed++;
    i = (fi >= 0 ? int(floor(fi)) : nxp_ + int(floor(fi)));
    j = (fj >= 0 ? int(floor(fj)) : nyp_ + int(floor(fj)));
    int index = Get2DIndex(i, j);
    sum += gammaXY_[index];
    fi += di;
    fj += dj;
  }
  sum /= gammaXY_[0];
  return(sum);
}


void CovGridSeparated::performTapering(float rangeX, float rangeY, float rangeZ)
{
  int i, j, k;
  float sq3 = float(sqrt(3.0f));

  //Taper z with Gausskernel with given range.
  int nzp = nzp_;
  int nzp2 = int(ceil(float(nzp)/2.0f)); //Round up if odd
  for(k= 0;k < nzp2; k++)
  {
    gammaZ_[k]       *= float(exp(-pow(float(k  )/rangeZ*sq3,2)));
    gammaZ_[nzp-1-k] *= float(exp(-pow(float(k+1)/rangeZ*sq3,2)));
  }

  int nxp  = nxp_;
  int nxp2 = int(ceil(float(nxp)/2.0f)); //Round up if odd
  int nyp  = nyp_;
  int nyp2 = int(ceil(float(nyp)/2.0f)); //Round up if odd

  //Taper x
  for(j=0;j<nyp;j++)
  {
    for(i=0;i<nxp2;i++)
    {
      int index1 = Get2DIndex(i, j);
      int index2 = Get2DIndex(nxp-i-1,j);
      gammaXY_[index1] *= float(exp(-pow(float(i  )/rangeX*sq3,2)));
      gammaXY_[index2] *= float(exp(-pow(float(i+1)/rangeX*sq3,2)));
    }
  }

  //Taper y
  for(i=0;i<nxp;i++)
  {
    for(j=0;j<nyp2;j++)
    {
      int index1 = Get2DIndex(i, j);
      int index2 = Get2DIndex(i,nyp-j-1);
      gammaXY_[index1] *= float(exp(-pow(float(j  )/rangeY*sq3,2)));
      gammaXY_[index2] *= float(exp(-pow(float(j+1)/rangeY*sq3,2)));
    }
  }
}

float CovGridSeparated::GetGamma2(int i1, int j1, int k1, int i2, int j2, int k2) const {

  if (i1 == IMISSING || j1 == IMISSING || k1 == IMISSING || 
      i2 == IMISSING || j2 == IMISSING || k2 == IMISSING) 
    return RMISSING;

  int deltai = i2 - i1;
  int deltaj = j2 - j1;
  int deltak = k2 - k1;

  if (tabulateCorr_) {
    //NBNB fjellvoll should be commented in when NOT kriging BGB model
    //
    //LogKit::LogFormatted(LogKit::DebugHigh,"PAL: abs(deltai),abs(deltaj),abs(deltak) : nxp_/2 nyp_/2 nzp_/2  =  %d %d %d : %d %d %d\n",
    //                      abs(deltai),abs(deltaj),abs(deltak),nxp_/2,nyp_/2,nzp_/2);
    //
    //NBNB-PAL: The variogram covers only half the grid in x- and y-directions. This explains 
    //          the if below. This, however, is problematic if the grid is small and the ranges 
    //          large. Probably the variogram should be defined in a grid twice as large (area) 
    //          as the modelling grid, so that we can estimate the correct covariance between
    //          any pair of points in the grid.
    //
    if (abs(deltai) >= nxp_/2 || abs(deltaj) >= nyp_/2 || abs(deltak) >= nzp_/2)
      return 0.0f;

    deltai = (deltai >= 0 ? deltai : nxp_ + deltai);
    deltaj = (deltaj >= 0 ? deltaj : nyp_ + deltaj);
    deltak = (deltak >= 0 ? deltak : nzp_ + deltak); 
  }
  return GetGamma(deltai, deltaj, deltak);
}

float 
CovGridSeparated::GetGamma(int i, int j, int k) const {
  if (tabulateCorr_) {
    if (!IsIndexValid(i, j, k))
      return RMISSING;
    return gammaXY_[Get2DIndex(i, j)] * gammaZ_[k];
  }
 
  float deltaX = i*dx_;
  float deltaY = j*dy_;
  float deltaZ = k*dz_; 
  if (rotAngle_ != 0.f)
    RotateVec(deltaX, deltaY, deltaZ, rotMatrix_);
  float h = float(sqrt(deltaX*deltaX/(rangeX_*rangeX_) + deltaY*deltaY/(rangeY_*rangeY_) +
                       deltaZ*deltaZ/(rangeZ_*rangeZ_)));
  return float(exp(-3.0*pow(h,power_)));
  
}

void CovGridSeparated::RotateVec(float& rx, float& ry, float& rz, const float mat[][3]) const {
  float res[3] = {0.0f}; 

  float input[3]; 
  input[0] = rx; input[1] = ry; input[2] = rz;
  int i,j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      res[i] += mat[i][j]* input[j];
    }
  } // end i
  rx = res[0]; ry = res[1]; rz = res[2];
}

void CovGridSeparated::writeXYGrid(const std::string fName) const {
  if (!tabulateCorr_) //NBNB fjellvoll add support for writing in this case also
    return;

  int nx = nxp_;
  int ny = nyp_;

  std::string baseName = fName + IO::SuffixAsciiIrapClassic();
  std::string fileName = IO::makeFullFileName(IO::PathToCorrelations(), baseName);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  file << std::fixed 
       << nx << " " << ny << " 1.0 1.0\n"
       << std::setprecision(2)
       << "0.0 " << static_cast<float>(nx-1) << " " 
       << "0.0 " << static_cast<float>(ny-1) << "\n";
  for(int j=0 ; j < ny ; j++)
    for(int i=0 ; i < nx ; i++)
      file << gammaXY_[i+j*nx] << "\n";
  file.close();
}
