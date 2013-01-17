/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include "src/box.h"
#include "src/simbox.h"

#include <math.h>

CBox::CBox(int iMin, int jMin, int kMin, int iMax, int jMax, int kMax, const Simbox* pSBox) {
  CheckMinMax(iMin, jMin, kMin, iMax, jMax, kMax);

  if (!pSBox)
    return;

  TruncateIndex(iMin_, jMin_, kMin_, pSBox);
  TruncateIndex(iMax_, jMax_, kMax_, pSBox);

}

void CBox::CheckMinMax(int iMin, int jMin, int kMin, int iMax, int jMax, int kMax) {
  if (iMax >= iMin) {
    iMax_ = iMax; iMin_ = iMin;
  }
  else {
    iMin_ = iMax; iMax_ = iMin;
  }

  if (jMax >= jMin) {
    jMax_ = jMax; jMin_ = jMin;
  }
  else {
    jMin_ = jMax; jMax_ = jMin;
  }

  if (kMax >= kMin) {
    kMax_ = kMax; kMin_ = kMin;
  }
  else {
    kMin_ = kMax; kMax_ = kMin;
  }
}

CBox::~CBox() {

}

bool
CBox::IsInside(int i, int j, int k) const {
  return (iMin_ <= i && iMax_ >= i && jMin_ <= j && jMax_ >= j && kMin_ <= k && kMax_ >= k);
}

void
CBox::TruncateIndex(int& i, int& j, int& k, const Simbox* pSBox) {
  if (i < 0) i = 0; if (j < 0) j = 0; if (k < 0) k = 0;
  if (i >= pSBox->getnx()) i = pSBox->getnx() - 1;
  if (j >= pSBox->getny()) j = pSBox->getny() - 1;
  if (k >= pSBox->getnz()) k = pSBox->getnz() - 1;

}

void CBox::ModifyBox(const CBox& boxS, const Simbox* pSBox, float rangeX, float rangeY, float rangeZ) {
  iMax_ = int((iMax_ + boxS.iMax_ + rangeX)/2.0f);
  jMax_ = int((jMax_ + boxS.jMax_ + rangeY)/2.0f);
  kMax_ = int((kMax_ + boxS.kMax_ + rangeZ)/2.0f);

  iMin_ = int((iMin_ + boxS.iMin_ - rangeX)/2.0f);
  jMin_ = int((jMin_ + boxS.jMin_ - rangeY)/2.0f);
  kMin_ = int((kMin_ + boxS.kMin_ - rangeZ)/2.0f);

  if (pSBox) {
    TruncateIndex(iMin_, jMin_, kMin_, pSBox);
    TruncateIndex(iMax_, jMax_, kMax_, pSBox);
  }
}
