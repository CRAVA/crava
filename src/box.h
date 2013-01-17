/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef CBOX_H
#define CBOX_H

class Simbox;

// a simple helping class: a box is used to define a neighbourhood in cokriging.
// The indexes are relative to the simbox

class CBox
{
public:
  CBox(int iMin = 0, int jMin = 0, int kMin = 0, int iMax = 0, int jMax = 0, int kMax = 0,
    const Simbox* pSBox = 0);
  ~CBox(void);
  bool IsInside(int i, int j, int k) const; // returns true if inside box, ijk is rotated indexes relative to simbox
  int GetNx() const { return iMax_ - iMin_ + 1;}
  int GetNy() const { return jMax_ - jMin_ + 1;}
  int GetNz() const { return kMax_ - kMin_ + 1;}
  void GetMin(int& iMin, int& jMin, int& kMin) const { iMin = iMin_; jMin = jMin_; kMin = kMin_;}
  void GetMax(int& iMax, int& jMax, int& kMax) const { iMax = iMax_; jMax = jMax_; kMax = kMax_;}
  void ModifyBox(const CBox& boxS, const Simbox* pSBox = 0,
    float rangeX = 0.0f, float rangeY = 0.0f, float rangeZ = 0.0f); // a special function
  bool operator==(const CBox& boxS) const
  {
    return (iMin_ == boxS.iMin_ && jMin_ == boxS.jMin_ && kMin_ == boxS.kMin_ &&
      iMax_ == boxS.iMax_ && jMax_ == boxS.jMax_ && kMax_ == boxS.kMax_);
  }

private:
  void TruncateIndex(int& i, int& j, int& k, const Simbox* pSBox);
  void CheckMinMax(int iMin, int jMin, int kMin, int iMax, int jMax, int kMax);

protected:
  int iMin_, jMin_, kMin_, iMax_, jMax_, kMax_;
};

#endif
