/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef CCOVGRIDSEPARATED_H
#define COVGRIDSEPARATED_H

class FFTGrid;

class CovGridSeparated
{
public:
  CovGridSeparated(const FFTGrid& grid);
  CovGridSeparated(int nxp, int nyp, int nzp,
                   float dx, float dy, float dz,
                   float rangeX, float rangeY,
                   float rangeZ, float power,
                   float rotAngle = 0.f,
                   bool tabulateCorr = true); // for BG model
  CovGridSeparated(int nxp, int nyp, int nzp); // for CrCorr = 0
  ~CovGridSeparated(void);

  float   GetGamma2(int i1, int j1, int k1, int i2, int j2, int k2) const; // returns RMISSING if it fails
  float   GetGamma(int i, int j, int k) const; // returns RMISSING if outside of simbox
  bool    IsIndexValid(int i, int j, int k) const;
  void    EstimateRanges(int& rangeX, int& rangeY, int& rangeZ) const;
  void    findTaperRanges(float & rangeX, float & rangeY, float & rangeZ) const;
  void    performTapering(float rangeX, float rangeY, float rangeZ);
  void    writeXYGrid(const std::string fName) const;

private:
  int     Get2DIndex(int i, int j)        const { return i + j*nxp_;}
  int     Get3DIndex(int i, int j, int k) const { return i + j*nxp_ + k*nxp_*nyp_; }
  void    EstimateRangeX(int& rangeX) const;
  void    EstimateRangeY(int& rangeY) const;
  void    EstimateRangeZ(int& rangeZ) const;
  float   angleSum(float angle) const;
  void    InitRotMatrix();
  void    RotateVec(float& rx, float& ry, float& rz, const float mat[][3]) const;

private:
  int     nxp_;
  int     nyp_;
  int     nzp_;
  float * gammaXY_;
  float * gammaZ_;
  float   rotMatrix_[3][3];
  bool    tabulateCorr_;
  // these next variables are only init if tabulateCorr_ is false
  float   dx_;
  float   dy_;
  float   dz_;
  float   rangeX_;
  float   rangeY_;
  float   rangeZ_;
  float   power_;
  float   rotAngle_;
};

inline bool
CovGridSeparated::IsIndexValid(int i, int j, int k) const // returns true if inside grid
{
  if (i < 0 || i >= nxp_ || j < 0 || j >= nyp_ || k < 0 || k >= nzp_)
    return false;
  return true;
}

#endif
