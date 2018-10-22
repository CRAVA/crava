/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "nrlib/iotools/logkit.hpp"

#include "src/krigingadmin.h"
#include "src/fftgrid.h"
#include "src/modelgeneral.h"
#include "src/bwellpt.h"
#include "src/simbox.h"
#include "src/covgridseparated.h"
#include "src/definitions.h"

CKrigingAdmin::CKrigingAdmin(const Simbox      & simbox,
                             CBWellPt         ** pBWellPt,
                             int                 noData,
                             CovGridSeparated  & covAlpha,
                             CovGridSeparated  & covBeta,
                             CovGridSeparated  & covRho,
                             CovGridSeparated  & covCrAlphaBeta,
                             CovGridSeparated  & covCrAlphaRho,
                             CovGridSeparated  & covCrBetaRho,
                             int                 dataTarget,
                             bool                backgroundModel) :
  simbox_(simbox),
  trendAlpha_(0),
  trendBeta_(0),
  trendRho_(0),
  covAlpha_(covAlpha),
  covBeta_(covBeta),
  covRho_(covRho),
  covCrAlphaBeta_(covCrAlphaBeta),
  covCrAlphaRho_(covCrAlphaRho),
  covCrBetaRho_(covCrBetaRho),
  pBWellPt_(pBWellPt),
  noData_(noData),
  dataTarget_(dataTarget),
  backgroundModel_(backgroundModel)
{
  Init(); // Common init
}

CKrigingAdmin::~CKrigingAdmin(void)
{
  delete pBWellGrid_;

  delete [] pIndexAlpha_;
  delete [] pIndexBeta_;
  delete [] pIndexRho_;
  int i;
  for (i = 0; i < GetSmoothBlockNx() - 2; i++) {
    delete [] ppKrigSmoothWeightsX_[i];
  }

  for (i = 0; i < GetSmoothBlockNy() - 2; i++) {
    delete [] ppKrigSmoothWeightsY_[i];
  }

  for (i = 0; i < GetSmoothBlockNz() - 2; i++) {
    delete [] ppKrigSmoothWeightsZ_[i];
  }
  delete [] ppKrigSmoothWeightsX_;
  delete [] ppKrigSmoothWeightsY_;
  delete [] ppKrigSmoothWeightsZ_;
}

void CKrigingAdmin::Init() {
  i_ = j_ = k_ = 0;

  //
  // Create indicator grid having 1.0f if data in cell and -1.0f if no data in cell
  // I wonder why Bjørn didn't choose and int grid with 1s and 0s instead?
  //
  pBWellGrid_ = CreateValidGrid();

  //
  // Find number of valid data
  //
  noValidAlpha_ = 0;
  noValidBeta_ = 0;
  noValidRho_ = 0;
  for (int m = 0 ; m < noData_ ; m++) {
    bool validA, validB, validR;
    pBWellPt_[m]->IsValidObs(validA, validB, validR);
    if (validA) noValidAlpha_++;
    if (validB) noValidBeta_++;
    if (validR) noValidRho_++;
  }
  noValid_ = noValidAlpha_ + noValidBeta_ + noValidRho_;

  noKrigedCells_ = noKrigedVariables_ = totalNoDataInCurrKrigBlock_= noEmptyDataBlocks_ = 0;
  sizeAlpha_ = sizeBeta_ = sizeRho_ = 0;
  rangeAlphaX_ = rangeAlphaY_ = rangeAlphaZ_ = 0;
  rangeBetaX_ = rangeBetaY_ = rangeBetaZ_ = 0;
  rangeRhoX_ = rangeRhoY_ = rangeRhoZ_ = 0;
  noSolvedMatrixEq_ = noCholeskyDecomp_ = noRMissing_ = 0;
  dxSmoothBlock_ = dySmoothBlock_ = dzSmoothBlock_ = 0;
  ppKrigSmoothWeightsX_ = ppKrigSmoothWeightsY_ = ppKrigSmoothWeightsZ_ = 0;
  failed2EstimateRange_ = failed2EstimateDefaultDataBoxAndBlock_ = false;

  if (backgroundModel_) {
    rangeAlphaX_ = rangeBetaX_ = rangeRhoX_  = simbox_.getnx();
    rangeAlphaY_ = rangeBetaY_ = rangeRhoY_  = simbox_.getny();
    rangeAlphaZ_ = rangeBetaZ_ = rangeRhoZ_  = 0;

    rangeX_ = static_cast<float>(std::max(rangeAlphaX_, rangeBetaX_));
    rangeX_ = static_cast<float>(std::max(static_cast<float>(rangeRhoX_), rangeX_));
    rangeY_ = static_cast<float>(std::max(rangeAlphaY_, rangeBetaY_));
    rangeY_ = static_cast<float>(std::max(static_cast<float>(rangeRhoY_), rangeY_));
    rangeZ_ = static_cast<float>(std::max(rangeAlphaZ_, rangeBetaZ_));
    rangeZ_ = static_cast<float>(std::max(static_cast<float>(rangeRhoZ_), rangeZ_));

    dxBlock_ = simbox_.getnx();
    dyBlock_ = simbox_.getny();
    dzBlock_ = 1;
    dxBlockExt_ = simbox_.getnx();
    dyBlockExt_ = simbox_.getny();
    dzBlockExt_ = 0;
  }
  else {
    EstimateSizeOfBlock();
    if (dxBlock_ < 1) dxBlock_ = 1;
    if (dyBlock_ < 1) dyBlock_ = 1;
    if (dzBlock_ < 1) dzBlock_ = 1;
  }

  dxSmoothBlock_ = (dxBlock_ <= 4 ? 1 : 2);
  dySmoothBlock_ = (dyBlock_ <= 4 ? 1 : 2);
  dzSmoothBlock_ = (dzBlock_ <= 4 ? 1 : 2);

  int i;
  ppKrigSmoothWeightsX_ = new double*[GetSmoothBlockNx() - 2];
  for (i = 0; i < GetSmoothBlockNx() - 2; i++) {
    ppKrigSmoothWeightsX_[i] = new double[GetSmoothBlockNx()];
  }

  ppKrigSmoothWeightsY_ = new double*[GetSmoothBlockNy() - 2];
  for (i = 0; i < GetSmoothBlockNy() - 2 ; i++) {
    ppKrigSmoothWeightsY_[i] = new double[GetSmoothBlockNy()];
  }

  ppKrigSmoothWeightsZ_ = new double*[GetSmoothBlockNz()- 2];
  for (i = 0; i < GetSmoothBlockNz()- 2; i++) {
    ppKrigSmoothWeightsZ_[i] = new double[GetSmoothBlockNz()];
  }

  const int sizeMaxBlock =
    (dxBlock_ + 2*static_cast<int>(ceil(rangeX_))) *
    (dyBlock_ + 2*static_cast<int>(ceil(rangeY_))) *
    (dzBlock_ + 2*static_cast<int>(ceil(rangeZ_)));

  int maxAlphaData2 = std::min(noValidAlpha_, sizeMaxBlock);
  int maxBetaData2  = std::min(noValidBeta_, sizeMaxBlock);
  int maxRhoData2   = std::min(noValidRho_, sizeMaxBlock);

  pIndexAlpha_ = new int[maxAlphaData2];
  pIndexBeta_ = new int[maxBetaData2];
  pIndexRho_  = new int[maxRhoData2];

  Require(dxBlockExt_ <= rangeX_ && dyBlockExt_ <= rangeY_ && dzBlockExt_ <= rangeZ_,
    "dxBlockExt_ <= rangeX_ && dyBlockExt_ <= rangeY_ && dzBlockExt_ <= rangeZ_");

  WriteDebugOutput();
}

void CKrigingAdmin::KrigAll(Gamma gamma, bool doSmoothing) {
  // basic set of neighbourhoods
  noCholeskyDecomp_ = noSolvedMatrixEq_ = 0;
  noRMissing_ = 0;
  const int nxBlock = NBlocks(dxBlock_, simbox_.getnx());
  const int nyBlock = NBlocks(dyBlock_, simbox_.getny());
  const int nzBlock = NBlocks(dzBlock_, simbox_.getnz());

  monitorSize_ = int(3*simbox_.getnx()*simbox_.getny()*simbox_.getnz()*0.02);
  monitorSize_ = std::max(1,monitorSize_);

  // loop over all kriging blocks
  int i,j,k;
  for (k = 0; k < nzBlock; k++) {
    int k1 = k*dzBlock_;
    for (j = 0; j < nyBlock; j++) {
      int j1 = j*dyBlock_;
      for (i = 0; i < nxBlock; i++) {
   //     LogKit::LogFormatted(LogKit::DebugLow,"%d %d %d : %d %d %d\n",i, j, k, nxBlock, nyBlock, nzBlock);
        int i1 = i*dxBlock_;
        currBlock_ = CBox(i1, j1, k1, i1 + dxBlock_ - 1, j1 + dyBlock_ - 1, k1 + dzBlock_ - 1, &simbox_);
        currDataBox_ = CBox(i1 - dxBlockExt_, j1 - dyBlockExt_, k1 - dzBlockExt_,
          i1 + dxBlock_ + dxBlockExt_ - 1, j1 + dyBlock_ + dyBlockExt_ - 1, k1 + dzBlock_ + dzBlockExt_ - 1,
          &simbox_);
        KrigBlock(gamma);
      } // end i
    } // end j
  } // end k
  noKrigedVariables_++;
  if (!backgroundModel_ && doSmoothing==true) {
    //LogKit::LogFormatted(LogKit::Low,"SmoothKrigedResult start\n");
    SmoothKrigedResult(gamma);
    //LogKit::LogFormatted(LogKit::Low,"SmoothKrigedResult end\n");
  }


}

void CKrigingAdmin::KrigAll(FFTGrid& trendAlpha, FFTGrid& trendBeta, FFTGrid& trendRho, SeismicParametersHolder & seismicParameters,
                            bool trendsAlreadySubtracted, int debugflag, bool doSmoothing) {
  Require(!trendAlpha.getIsTransformed()
          && !trendBeta.getIsTransformed()
          && !trendRho.getIsTransformed(),
          "!trendAlpha.getIsTransformed() && !trendBeta.getIsTransformed() && !trendRho.getIsTransformed()");

  if (!trendsAlreadySubtracted)
    SubtractTrends(trendAlpha, trendBeta, trendRho);
  if(debugflag>0)
  {
    FFTGrid *blockGrid = new FFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(), simbox_.getnx(), simbox_.getny(), simbox_.getnz());
    LogKit::LogFormatted(LogKit::DebugLow,"Grid for reporting block size created. Memory requirements may be wrongly calculated.\n");
    blockGrid->fillInConstant(0.0, false);
    int i,j,k;
    for(i=0;i<simbox_.getnx();i++)
    {
      int i1 =  div(i,dxBlock_).quot+1;
      for(j=0;j<simbox_.getny();j++)
      {
        int j1 = div(j,dyBlock_).quot +1;
        for(k=0;k<simbox_.getnz();k++)
        {
          int k1 = div(k,dzBlock_).quot+1;
          if(div(i1+j1+k1,2).rem==0)
            blockGrid->setRealValue(i,j,k,0.0);
          else
            blockGrid->setRealValue(i,j,k,1.0);
        }
      }
    }

    seismicParameters.SetBlockGrid(blockGrid);
  }

  trendAlpha_ = &trendAlpha;
  trendBeta_  = &trendBeta;
  trendRho_   = &trendRho;

  printf("\n  0%%       20%%       40%%       60%%       80%%      100%%");
  printf("\n  |    |    |    |    |    |    |    |    |    |    |  ");
  printf("\n  ^");

  trendAlpha_->setAccessMode(FFTGrid::RANDOMACCESS);
  LogKit::LogFormatted(LogKit::DebugHigh,"Start CKrigingAdminKrigAll: Alpha\n");
  KrigAll(ALPHA_KRIG, doSmoothing);
  WriteDebugOutput2();
  LogKit::LogFormatted(LogKit::DebugHigh,"End CKrigingAdminKrigAll: Alpha\n");
  trendAlpha_->endAccess();

  trendBeta_->setAccessMode(FFTGrid::RANDOMACCESS);
  LogKit::LogFormatted(LogKit::DebugHigh,"Start CKrigingAdminKrigAll: Beta\n");
  KrigAll(BETA_KRIG, doSmoothing);
  WriteDebugOutput2();
  LogKit::LogFormatted(LogKit::DebugHigh,"End CKrigingAdminKrigAll: Beta\n");
  trendBeta_->endAccess();

  trendRho_->setAccessMode(FFTGrid::RANDOMACCESS);
  LogKit::LogFormatted(LogKit::DebugHigh,"Start CKrigingAdminKrigAll: Rho\n");
  KrigAll(RHO_KRIG, doSmoothing);
  WriteDebugOutput2();
  LogKit::LogFormatted(LogKit::DebugHigh,"End CKrigingAdminKrigAll: Rho\n");
  trendRho_->endAccess();

  printf("\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"KrigAll finished\n");
}

void CKrigingAdmin::KrigBlock(Gamma gamma)
{
  // search for neighbours
  LogKit::LogFormatted(LogKit::DebugHigh,"FindDataInDataBlockLoop(gamma) called next\n");
  FindDataInDataBlockLoop(gamma);
  LogKit::LogFormatted(LogKit::DebugHigh,"sizeAlpha_: %d\n", sizeAlpha_);
  LogKit::LogFormatted(LogKit::DebugHigh,"sizeBeta_: %d\n", sizeBeta_);
  LogKit::LogFormatted(LogKit::DebugHigh,"sizeRho_: %d\n", sizeRho_);
  LogKit::LogFormatted(LogKit::DebugHigh,"totalNoDataInCurrKrigBlock_: %d\n", totalNoDataInCurrKrigBlock_);
  //FindDataInDataBlock(gamma); //DEBUG

  int iMin, jMin, kMin, iMax, jMax, kMax;
  currBlock_.GetMin(iMin, jMin, kMin); currBlock_.GetMax(iMax, jMax, kMax);
  if (!totalNoDataInCurrKrigBlock_) {
    int cellsInBlock = (kMax - kMin + 1)*(jMax - jMin + 1)*(iMax - iMin + 1);
    int i;
    for (i=noKrigedCells_ ; i<noKrigedCells_+cellsInBlock ; i++) {
      if (noKrigedCells_ > 0 && i%monitorSize_ == 0) {
        printf("^");
        fflush(stdout);
      }
    }
    noKrigedCells_ += cellsInBlock;
    noEmptyDataBlocks_++;
    return;
  }

  int n = sizeAlpha_ + sizeBeta_ + sizeRho_;

  NRLib::Matrix krigMatrix(n, n);
  NRLib::Vector residual(n);

  // Set kriging matrix based on data finds

  SetMatrix(krigMatrix,
            residual,
            gamma);

  NRLib::SymmetricMatrix K(n);

  for (int j = 0 ; j < n ; j++) {
    for (int i = 0 ; i <= j ; i++) {
      K(i,j) = krigMatrix(i,j);
    }
  }

  FFTGrid * pGrid = 0;

  switch(gamma) {
  case ALPHA_KRIG :
    pGrid = trendAlpha_;
    break;
  case BETA_KRIG :
    pGrid = trendBeta_;
    break;
  case RHO_KRIG :
    pGrid = trendRho_;
    break;
  default :
    Require(false, "switch failed");
  } // end switch

  NRLib::Vector x(n);
  // NBNB-PAL: Add try/catch loop around CholeskySolve call with a regularization term.
  NRLib::CholeskySolve(K, residual, x);

  NRLib::Vector kVec(n);

  for (k_ = kMin; k_ <= kMax; k_++) {
    for (j_ = jMin; j_ <= jMax; j_++) {
      for (i_ = iMin; i_ <= iMax; i_++) {

        // set kriging vector
        SetKrigVector(kVec, gamma);

        // kriging;
        float result = pGrid->getRealValue(i_, j_, k_);
        if (result == RMISSING) {
          noRMissing_++;
        }
        else {
          result += static_cast<float>(kVec * x);

          if(pGrid->setRealValue(i_, j_, k_, result))
            Require(false, "pGrid->setRealValue failed"); // something is serious wrong...

          noSolvedMatrixEq_++;
        }

        noKrigedCells_++;
        // return result
        if (noKrigedCells_%monitorSize_ == 0) {
          printf("^");
          fflush(stdout);
        }
      } // end for i
    } // end for j
  } // end for k

}

FFTGrid* CKrigingAdmin::CreateValidGrid() const
{
  //FFTGrid* pGrid = new FFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(),
  //                             simbox_.getnx(), simbox_.getny(), simbox_.getnz());
  FFTGrid* pGrid = ModelGeneral::CreateFFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(),
                                               simbox_.getnx(), simbox_.getny(), simbox_.getnz(), false);
  pGrid->fillInConstant(-1.0f, false);
  pGrid->setAccessMode(FFTGrid::RANDOMACCESS);
  int i;
  for (i = 0; i < noData_; i++) {
    bool validA, validB, validR;
    pBWellPt_[i]->IsValidObs(validA, validB, validR);
    if (validA || validB || validR) {
      int i1, j1, k1;
      pBWellPt_[i]->GetIJK(i1, j1, k1);
      pGrid->setRealValue(i1, j1, k1, 1.0f);
    }
  }
  pGrid->endAccess();
  return pGrid;
}

void
CKrigingAdmin::SubtractTrends(FFTGrid& trend_alpha, FFTGrid& trend_beta, FFTGrid& trend_rho)
{
  FFTGrid* pGrid = 0;
  CBWellPt::Gamma type = CBWellPt::ERROR;
  for (int l = 0; l < 3; l++) {
    switch (l) {
      case 0 :
        pGrid = &trend_alpha;
        type  = CBWellPt::ALPHA;
        break;
      case 1 :
        pGrid = &trend_beta;
        type  = CBWellPt::BETA;
        break;
      case 2 :
        pGrid = &trend_rho;
        type  = CBWellPt::RHO;
        break;
      default :
        assert(1 == 0);
        break;
    }
    pGrid->setAccessMode(FFTGrid::RANDOMACCESS);
    bool isTrans = pGrid->getIsTransformed();
    if (isTrans)
      pGrid->invFFTInPlace();
    for (int m = 0 ; m < noData_ ; m++) {
      int i,j,k;
      pBWellPt_[m]->GetIJK(i, j, k);
      float value = pGrid->getRealValue(i, j, k);
      pBWellPt_[m]->SubtractOnly(type, value);

    } // end loop over well obs
    if (isTrans)
      pGrid->fftInPlace();
    pGrid->endAccess();
  } // end l
}

CKrigingAdmin::DataBoxSize
CKrigingAdmin::FindDataInDataBlock(Gamma gamma, const CBox & dataBox) {
  sizeAlpha_ = sizeBeta_ = sizeRho_ = totalNoDataInCurrKrigBlock_ = 0;
  const int countTotalMin = int(dataTarget_*(1.0f - maxDataTolerance_/100.0f));
  const int countTotalMax = int(dataTarget_*(1.0f + maxDataTolerance_/100.0f));

  for (int i = 0; i < noData_ ; i++) {
    int i1, j1, k1;
    pBWellPt_[i]->GetIJK(i1, j1, k1);
    if (dataBox.IsInside(i1, j1, k1)) {
      bool validA, validB, validR;
      pBWellPt_[i]->IsValidObs(validA, validB, validR);
      switch (gamma) {
      case ALPHA_KRIG :
        if (validA && ++totalNoDataInCurrKrigBlock_)
          pIndexAlpha_[sizeAlpha_++] = i;
        else {
          if (validB && ++totalNoDataInCurrKrigBlock_)
            pIndexBeta_[sizeBeta_++] = i;

          if (validR && ++totalNoDataInCurrKrigBlock_)
            pIndexRho_[sizeRho_++] = i;
        }
        break;

      case BETA_KRIG :
        if (validB && ++totalNoDataInCurrKrigBlock_)
          pIndexBeta_[sizeBeta_++] = i;
        else {
          if (validA && ++totalNoDataInCurrKrigBlock_)
            pIndexAlpha_[sizeAlpha_++] = i;

          if (validR && ++totalNoDataInCurrKrigBlock_)
            pIndexRho_[sizeRho_++] = i;
        }
        break;
      case RHO_KRIG :
        if (validR && ++totalNoDataInCurrKrigBlock_)
          pIndexRho_[sizeRho_++] = i;
        else {
          if (validA && ++totalNoDataInCurrKrigBlock_)
            pIndexAlpha_[sizeAlpha_++] = i;

          if (validB && ++totalNoDataInCurrKrigBlock_)
            pIndexBeta_[sizeBeta_++] = i;
        }
        break;

      default:
        // should never happen
        Require(false, "switch failed");
      } // end switch
      // early exit
    } // end if
  } // end i
  LogKit::LogFormatted(LogKit::DebugHigh,"Found %d data. (%d, %d)\n", totalNoDataInCurrKrigBlock_,
    countTotalMin, countTotalMax);
  if (totalNoDataInCurrKrigBlock_ <= countTotalMax && totalNoDataInCurrKrigBlock_ >= countTotalMin)
    return DBS_RIGHT;
  if (totalNoDataInCurrKrigBlock_ < countTotalMin)
    return DBS_TOO_SMALL;
  else {//(totalNoDataInCurrKrigBlock_ > countTotalMax)
    return DBS_TOO_BIG;
  }
}


void CKrigingAdmin::FindDataInDataBlockLoop(Gamma gamma) {
  int counter = 0;
  DataBoxSize currDataBoxSize, startDataboxSize, testDataBoxSize;
  currDataBoxSize = FindDataInDataBlock(gamma, currDataBox_);
  startDataboxSize = currDataBoxSize;
  //CBox minDataBox = currBlock_;
  CBox minDataBox = currDataBox_;
  int iMin,iMax,jMin,jMax,kMin,kMax;
  currBlock_.GetMin(iMin,jMin,kMin);
  currBlock_.GetMax(iMax,jMax,kMax);
  CBox maxDataBox(iMin-int(rangeX_),jMin-int(rangeY_),kMin-int(rangeZ_),
    iMax+int(rangeX_),jMax+int(rangeY_),kMax+int(rangeZ_));

  switch (currDataBoxSize) {
  case DBS_RIGHT:
    // NBNB-PAL: Nothing to do here? I put in this switch option to avoid a crash (CRA-75)
    break;
  case DBS_TOO_SMALL:
    testDataBoxSize = FindDataInDataBlock(gamma, maxDataBox);
    if(testDataBoxSize != DBS_TOO_BIG)
    {
      currDataBox_ = maxDataBox;
      currDataBoxSize = DBS_RIGHT;
    }
    break;
  case DBS_TOO_BIG:
    testDataBoxSize = FindDataInDataBlock(gamma, minDataBox);
    if(testDataBoxSize != DBS_TOO_SMALL)
    {
      currDataBox_ = minDataBox;
      currDataBoxSize = DBS_RIGHT;
    }
    break;
  default :
    Require(false, "switch failed");
    break;

  }

  while (currDataBoxSize != DBS_RIGHT) {
    switch (currDataBoxSize) {
    case DBS_TOO_SMALL :
      minDataBox = currDataBox_;
      currDataBox_.ModifyBox(maxDataBox);
      break;
    case DBS_TOO_BIG :
      maxDataBox = currDataBox_;
      currDataBox_.ModifyBox(minDataBox);
      break;
    default :
      Require(false, "switch failed");
      break;

    } // end switch
    counter++;
    //if (currDataBoxSize != startDataboxSize || counter++ >= maxDataBlockLoopCounter_ || prevDataBox == currDataBox_)
    //if (currDataBoxSize != startDataboxSize || prevDataBox == currDataBox_)
    if(currDataBox_ == maxDataBox || currDataBox_ == minDataBox)
      break;

    currDataBoxSize = FindDataInDataBlock(gamma, currDataBox_);

  } // end while
  currDataBox_.ModifyBox(currDataBox_, &simbox_); //Does not modify, only truncates.

  LogKit::LogFormatted(LogKit::DebugHigh,"FindDataInDataBlock iterations: %d\n", counter);
}


int CKrigingAdmin::NBlocks(int dBlocks, int lSBox) const {
  if (lSBox % dBlocks == 0)
    return lSBox/dBlocks;

  return lSBox/dBlocks + 1;
}

void CKrigingAdmin::SetMatrix(NRLib::Matrix & krigMatrix,
                              NRLib::Vector & residual,
                              Gamma           gamma) {
  assert(gamma >= 0);
  if (!totalNoDataInCurrKrigBlock_)
    return;
  int a, b, r;

  // set Kriging Matrix
  // for alpha kriging
  int a2, b2, r2;
  // first row
  for (a = 0; a < sizeAlpha_; a++) {
    int krigRowIndex = a;
    int indexA = pIndexAlpha_[a];
    int i,j,k;
    pBWellPt_[indexA]->GetIJK(i, j, k);
    // K_aa
    for (a2 = 0; a2 < sizeAlpha_; a2++) {
      int indexA2 = pIndexAlpha_[a2];
      int i2, j2, k2;
      pBWellPt_[indexA2]->GetIJK(i2, j2, k2);

      krigMatrix(krigRowIndex,a2) = covAlpha_.GetGamma2(i, j, k, i2, j2, k2);
    } // end a2

    // K_ab
    for (b2 = 0; b2 < sizeBeta_; b2++) {
      int indexB2 = pIndexBeta_[b2];
      int i2, j2, k2;
      pBWellPt_[indexB2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, b2 + sizeAlpha_) = covCrAlphaBeta_.GetGamma2(i, j, k, i2, j2, k2);
    } // end b2

    // K_ar
    for (r2 = 0; r2 < sizeRho_; r2++) {
      int indexR2 = pIndexRho_[r2];
      int i2, j2, k2;
      pBWellPt_[indexR2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, r2 + sizeAlpha_ + sizeBeta_) = covCrAlphaRho_.GetGamma2(i, j, k, i2, j2, k2);
    } // end r2
  }// end a

  // second row
  for (b = 0; b < sizeBeta_; b++) {
    int krigRowIndex = b + sizeAlpha_;
    int indexB = pIndexBeta_[b];
    int i,j,k;
    pBWellPt_[indexB]->GetIJK(i, j, k);
    // K_ba
    for (a2 = 0; a2 < sizeAlpha_; a2++) {
      int indexA2 = pIndexAlpha_[a2];
      int i2, j2, k2;
      pBWellPt_[indexA2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex,a2) = covCrAlphaBeta_.GetGamma2(i2, j2, k2, i, j, k); // flip
    } // end a2

    // K_bb
    for (b2 = 0; b2 < sizeBeta_; b2++) {
      int indexB2 = pIndexBeta_[b2];
      int i2, j2, k2;
      pBWellPt_[indexB2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, b2 + sizeAlpha_) = covBeta_.GetGamma2(i, j, k, i2, j2, k2);
    } // end b2

    // K_br
    for (r2 = 0; r2 < sizeRho_; r2++) {
      int indexR2 = pIndexRho_[r2];
      int i2, j2, k2;
      pBWellPt_[indexR2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex,r2 + sizeAlpha_ + sizeBeta_) = covCrBetaRho_.GetGamma2(i, j, k, i2, j2, k2);
    } // end r2
  }// end b
  // third row
  for (r = 0; r < sizeRho_; r++) {
    int krigRowIndex = r + sizeAlpha_ + sizeBeta_;
    int indexR = pIndexRho_[r];
    int i,j,k;
    pBWellPt_[indexR]->GetIJK(i, j, k);
    // K_ra
    for (a2 = 0; a2 < sizeAlpha_; a2++) {
      int indexA2 = pIndexAlpha_[a2];
      int i2, j2, k2;
      pBWellPt_[indexA2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, a2) = covCrAlphaRho_.GetGamma2(i2, j2, k2, i, j, k); // flip
    } // end a2

    // K_rb
    for (b2 = 0; b2 < sizeBeta_; b2++) {
      int indexB2 = pIndexBeta_[b2];
      int i2, j2, k2;
      pBWellPt_[indexB2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, b2  + sizeAlpha_) = covCrBetaRho_.GetGamma2(i2, j2, k2, i, j, k); // flip
    } // end b2

    // K_rr
    for (r2 = 0; r2 < sizeRho_; r2++) {
      int indexR2 = pIndexRho_[r2];
      int i2, j2, k2;
      pBWellPt_[indexR2]->GetIJK(i2, j2, k2);
      krigMatrix(krigRowIndex, r2 + sizeAlpha_ + sizeBeta_) = covRho_.GetGamma2(i, j, k, i2, j2, k2);
    } // end r2
  }// end r

  // Also calulates the kriging data vector
  for (a = 0; a < sizeAlpha_; a++) {
    int indexA = pIndexAlpha_[a];
    residual(a) = pBWellPt_[indexA]->GetAlpha();
  } // end a

  for (b = 0; b < sizeBeta_; b++) {
    int indexB = pIndexBeta_[b];
    residual(sizeAlpha_ + b) = pBWellPt_[indexB]->GetBeta();
  } // end b

  for (r = 0; r < sizeRho_; r++) {
    int indexR = pIndexRho_[r];
    residual(sizeAlpha_ + sizeBeta_ + r) = pBWellPt_[indexR]->GetRho();
  } // end r

}

void CKrigingAdmin::SetKrigVector(NRLib::Vector & k,
                                  Gamma           gamma)
{
  int offsetB1, offsetR1;
  offsetB1 = sizeAlpha_; offsetR1 = sizeAlpha_ + sizeBeta_;
  const CovGridSeparated *pA = NULL, *pB = NULL, *pR = NULL;
  bool flipA = false, flipB = false, flipR = false;
  switch(gamma) {
  case ALPHA_KRIG :
    pA = &covAlpha_; pB = &covCrAlphaBeta_; pR = &covCrAlphaRho_;
    flipA = flipB = flipR = false;
    break;
  case BETA_KRIG :
    pA = &covCrAlphaBeta_; pB = &covBeta_; pR = &covCrBetaRho_;
    flipA = flipB = true; flipR = false;
    break;
  case RHO_KRIG :
    pA = &covCrAlphaRho_; pB = &covCrBetaRho_; pR = &covRho_;
    flipA = flipB = flipR = true;
    break;
  default :
    // should never arrive here
    Require(false, "switch failed");

  } // end switch

  // set krig vector
  // for alpha kriging

  //LogKit::LogFormatted(LogKit::DebugHigh,"\n");

  // k_a
  int a2;
  for (a2 = 0; a2 < sizeAlpha_; a2++) {
    int indexA2 = pIndexAlpha_[a2];
    int i2, j2, k2;
    pBWellPt_[indexA2]->GetIJK(i2, j2, k2);
    k(a2) = (!flipA ? pA->GetGamma2(i_, j_, k_, i2, j2, k2) : pA->GetGamma2(i2, j2, k2, i_, j_, k_));
  } // end a2

  // k_b
  int b2;
  for (b2 = 0; b2 < sizeBeta_; b2++) {
    int indexB2 = pIndexBeta_[b2];
    int i2, j2, k2;
    pBWellPt_[indexB2]->GetIJK(i2, j2, k2);
    k(b2 + offsetB1) = (!flipB ? pB->GetGamma2(i_, j_, k_, i2, j2, k2) : pB->GetGamma2(i2, j2, k2, i_, j_, k_));
  } // end b2

  // k_r
  int r2;
  for (r2 = 0; r2 < sizeRho_; r2++) {
    int indexR2 = pIndexRho_[r2];
    int i2, j2, k2;
    pBWellPt_[indexR2]->GetIJK(i2, j2, k2);
    k(r2 + offsetR1) = (!flipR ? pR->GetGamma2(i_, j_, k_, i2, j2, k2) : pR->GetGamma2(i2, j2, k2, i_, j_, k_));
  } // end r2
}

void CKrigingAdmin::EstimateSizeOfBlock() {
  float rAx, rAy, rAz, rBx, rBy, rBz, rRx, rRy, rRz, rx, ry, rz;
  covAlpha_.findTaperRanges(rAx, rAy, rAz);
  covBeta_.findTaperRanges(rBx, rBy, rBz);
  covRho_.findTaperRanges(rRx, rRy, rRz);
  rx = std::max(rAx, rBx);
  rx = std::max(rx, rRx);
  ry = std::max(rAy, rBy);
  ry = std::max(ry, rRy);
  rz = std::max(rAz, rBz);
  rz = std::max(rz, rRz);
  covAlpha_.performTapering(rx, ry, rz);
  covBeta_.performTapering(rx, ry, rz);
  covRho_.performTapering(rx, ry, rz);
  covCrAlphaBeta_.performTapering(rx, ry, rz);
  covCrAlphaRho_.performTapering(rx, ry, rz);
  covCrBetaRho_.performTapering(rx, ry, rz);
  covAlpha_.EstimateRanges(rangeAlphaX_, rangeAlphaY_, rangeAlphaZ_);
  covBeta_.EstimateRanges(rangeBetaX_, rangeBetaY_, rangeBetaZ_);
  covRho_.EstimateRanges(rangeRhoX_, rangeRhoY_, rangeRhoZ_);

  // Should check ranges

  if (!(rangeAlphaX_ > 0 && rangeAlphaY_ > 0 && rangeAlphaZ_ > 0 &&
        rangeBetaX_ > 0 && rangeBetaY_ > 0 && rangeBetaZ_ > 0 &&
        rangeRhoX_ > 0 && rangeRhoY_ > 0 && rangeRhoZ_ > 0)) {
    failed2EstimateRange_ = true;
  }

  rangeX_ = static_cast<float>(std::max(rangeAlphaX_, rangeBetaX_));
  rangeX_ = static_cast<float>(std::max(static_cast<float>(rangeRhoX_), rangeX_));
  rangeY_ = static_cast<float>(std::max(rangeAlphaY_, rangeBetaY_));
  rangeY_ = static_cast<float>(std::max(static_cast<float>(rangeRhoY_), rangeY_));
  rangeZ_ = static_cast<float>(std::max(rangeAlphaZ_, rangeBetaZ_));
  rangeZ_ = static_cast<float>(std::max(static_cast<float>(rangeRhoZ_), rangeZ_));

  LogKit::LogFormatted(LogKit::Low,"Estimated ranges(grid cells) from covariance cubes:\n");
  LogKit::LogFormatted(LogKit::Low,"             rangeX     rangeY     rangeZ\n");
  LogKit::LogFormatted(LogKit::Low,"-----------------------------------------\n");
  LogKit::LogFormatted(LogKit::Low,"Vp   :     %8d   %8d   %8d\n", rangeAlphaX_, rangeAlphaY_, rangeAlphaZ_);
  LogKit::LogFormatted(LogKit::Low,"Vs   :     %8d   %8d   %8d\n", rangeBetaX_, rangeBetaY_, rangeBetaZ_);
  LogKit::LogFormatted(LogKit::Low,"Rho  :     %8d   %8d   %8d\n", rangeRhoX_, rangeRhoY_, rangeRhoZ_);
  LogKit::LogFormatted(LogKit::Low,"Used :     %8.0f   %8.0f   %8.0f\n", rangeX_, rangeY_, rangeZ_);

  if (noData_ <= dataTarget_) {
     dxBlock_ = simbox_.getnx();
     dyBlock_ = simbox_.getny();
     dzBlock_ = simbox_.getnz();
     dxBlockExt_ = static_cast<int>(rangeX_);
     dyBlockExt_ = static_cast<int>(rangeY_);
     dzBlockExt_ = static_cast<int>(rangeZ_);
     return;
  }

  if(rangeX_==0.0 && rangeY_ ==0.0)
  {
    dxBlock_ = 1;
    dyBlock_ = 1;
    dxBlockExt_ = 0;
    dyBlockExt_ = 0;
    dzBlock_ = static_cast<int>(ceil(rangeZ_));
    dzBlockExt_ = dzBlock_;
    return;
  }
  else if(rangeX_==0.0)
  {
    dxBlock_ = 1;
    dxBlockExt_ = 0;
    dyBlock_ = static_cast<int>(ceil(rangeY_));
    dyBlockExt_ = dyBlock_;
    dzBlock_ = static_cast<int>(ceil(rangeZ_));
    dzBlockExt_ = dzBlock_;
    return;
  }
  else if(rangeY_ == 0.0)
  {
    dyBlock_ = 1;
    dyBlockExt_ = 0;
    dxBlock_ = static_cast<int>(ceil(rangeX_));
    dxBlockExt_ = dxBlock_;
    dzBlock_ = static_cast<int>(ceil(rangeZ_));
    dzBlockExt_ = dzBlock_;
    return;

  }
  int dxBlock;
  float nd;

  int dxBlockT = 0;
  bool rapidInc = false;

  // find minimum time given that dyBlockExt = rangeY_
  int dyBlockExt = static_cast<int>(rangeY_);
  int dyBlockExtT = dyBlockExt;
  float tMin = CalcCPUTime(1, static_cast<float>(dyBlockExt), nd, rapidInc);
  float ndT = nd;
  for (dxBlock = 2; dxBlock <= simbox_.getnx(); dxBlock += 1) {
    float t1 = CalcCPUTime(static_cast<float>(dxBlock),static_cast<float>(dyBlockExt), nd, rapidInc);
    if (t1 < tMin) {
      tMin = t1;
      dxBlockT = dxBlock;
      ndT = nd;
    }
  } // end loop dxBlock

  if (ndT < dataTarget_) {
    dxBlock_ = dxBlockT;
    dyBlockExt_ = dyBlockExtT;

    EstimateSizeOfBlock2();
    return;
  }
  // first estimate failed
  bool test = false;
  float sumdyBExt = 0.0f, sumdxBlock = 0.0f;
  int count = 0;
  CalcCPUTime(1.0f,1.0f, nd, rapidInc);
  CalcCPUTime(static_cast<float>(simbox_.getnx()),rangeY_, nd, rapidInc);
  //int maxNd = (int)nd; // LATER better to use dataTarget_ here ?
  int maxNd = dataTarget_;
  int nd2;
  const float V = static_cast<float>(simbox_.getnx() * simbox_.getny() * simbox_.getnz());
  const float R = rangeX_*rangeY_*rangeZ_;

  // for (nd2 = minNd; nd2 <= maxNd; nd2++) {
  for (nd2 = maxNd; nd2 <= maxNd; nd2++) {
    const float a = static_cast<float>(pow(nd2*V/(noValid_*R),0.333333333f));
    int maxdyBlockExt = static_cast<int>(0.5f*(a - 1.0f/rangeX_)*rangeY_);
    for (dyBlockExt = 1; dyBlockExt <= maxdyBlockExt; dyBlockExt++) {
      float dxBlock2 = rangeX_ * (a - 2*dyBlockExt/rangeY_);
      if(dxBlock2 < 1)
        continue; // should never happen, but...

      //float t1 = CalcCPUTime(dxBlock2,(float)dyBlockExt, nd, rapidInc); //NBNB Bjorn: Used before, keep if change of mind?
      if (rapidInc) {
        test = true;
        sumdyBExt += dyBlockExt;
        sumdxBlock += dxBlock2;
        count++;
        break;
      }

    } // end loop dyBlockExt
  } // end loop nd2
  if (test) {
    dxBlock_    = static_cast<int>(ceil(sumdxBlock/count));
    dyBlockExt_ = static_cast<int>(ceil(sumdyBExt/count));

    EstimateSizeOfBlock2();
    return;
  }

  // if we are here we failed to estimate, should never happen, but we put in worst case values
  failed2EstimateDefaultDataBoxAndBlock_ = true;
  // speed will slow down
  dxBlock_ = simbox_.getnx();
  dyBlock_ = simbox_.getny();
  dzBlock_ = simbox_.getnz();
  dxBlockExt_ = static_cast<int>(rangeX_);
  dyBlockExt_ = static_cast<int>(rangeY_);
  dzBlockExt_ = static_cast<int>(rangeZ_);
}


void CKrigingAdmin::EstimateSizeOfBlock2() {
  dyBlock_ = static_cast<int>(ceil((dxBlock_*rangeY_) / rangeX_));
  dzBlock_ = static_cast<int>(ceil((dxBlock_*rangeZ_) / rangeX_));

  if (dyBlockExt_ > static_cast<int>(rangeY_)) dyBlockExt_ = static_cast<int>(rangeY_);
  dxBlockExt_ = static_cast<int>(ceil((dyBlockExt_*rangeX_) / rangeY_));
  dzBlockExt_ = static_cast<int>(ceil((dyBlockExt_*rangeZ_) / rangeY_));

  if (dxBlockExt_ > static_cast<int>(rangeX_)) dxBlockExt_ = static_cast<int>(rangeX_);
  if (dzBlockExt_ > static_cast<int>(rangeZ_)) dzBlockExt_ = static_cast<int>(rangeZ_);
}


float CKrigingAdmin::CalcCPUTime(float dxBlock, float dyBlockExt, float& nd, bool& rapidInc) {

  //static const float t_bigK = 150.0f; static const float t_chol = 1.90f;
  //static const float t_solve = 70.0f; static const float t_smallk = 164.0f;

  Require(dxBlock >= 1.0f && dyBlockExt >= 0.0f, "dxBlock >= 1.0f && dyBlockExt >= 0.0f");
  static const float t_bigK   = static_cast<float>(1.50E-6);
  static const float t_chol   = static_cast<float>(1.90E-8);
  static const float t_solve  = static_cast<float>(7.00E-7);
  static const float t_smallk = static_cast<float>(1.64E-6);
  float dxBlockExt = (dyBlockExt*rangeX_) / rangeY_;
  float dzBlockExt = (dyBlockExt*rangeZ_) / rangeY_;
  float dyBlock = (dxBlock*rangeY_) / rangeX_;
  float dzBlock = (dxBlock*rangeZ_) / rangeX_;

  const float V = static_cast<float>(simbox_.getnx() * simbox_.getny() * simbox_.getnz());
  const float v = dxBlock * dyBlock * dzBlock;
  const float Nss = V/v;
  nd = (dxBlock + 2*dxBlockExt)*(dyBlock + 2*dyBlockExt)*(dzBlock + 2*dzBlockExt)* noValid_ / V;
  const float T_BigK = t_bigK*Nss*nd*nd;
  const float T_chol = t_chol*Nss*nd*nd*nd;
  const float T_solve = t_solve*Nss*nd*nd;
  const float T_smallk = t_smallk*V*nd;
  rapidInc = (T_chol >= 0.5*T_smallk);
  return T_BigK + T_chol + T_solve + T_smallk;
  //return t_bigK*Nss*nd*nd + t_chol*Nss*nd*nd*nd + t_solve*Nss*nd*nd + t_smallk*V*nd;
  //return nd * (Nss*nd*(t_bigK + t_chol*nd + t_solve) + t_smallk*V);
}


const FFTGrid&
CKrigingAdmin::CreateAndFillFFTGridWithCov(int i1) {
  FFTGrid* pGrid = new FFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(),
    simbox_.getnx(), simbox_.getny(), simbox_.getnz());

  pGrid->createRealGrid();
  pGrid->setAccessMode(FFTGrid::RANDOMACCESS);
  float rangeXf = 0.0f, rangeYf = 0.0f, rangeZf = 0.0f;
  float power = 1.0f;
  switch (i1) {
  case 1:
    rangeXf = 2.0f; rangeYf = 4.0f; rangeZf = 2.0f;//2.5f;
    power = 1.0f;
    break;
  case 2:
    rangeXf = 5.0f; rangeYf = 8.0f; rangeZf = 3.5f;
    power = 1.0f;
    break;
  case 3:
    rangeXf = 3.0f; rangeYf = 16.0f; rangeZf = 2.5f;
    power = 1.0f;
    break;
  default:
    Require(false, "switch failed");
  }// end switch
  //rangeXf = 2.0f; rangeYf = 4.0f; rangeZf = 2.0f;//2.5f;
  int i, j, k;
  float rangeX = static_cast<float>(simbox_.getlx())/rangeXf;
  float rangeY = static_cast<float>(simbox_.getly())/rangeYf;
  float rangeZ = static_cast<float>(simbox_.getlz())/rangeZf;

  // need to take into acccount practical range
  rangeX /= 3.0f;
  rangeY /= 3.0f;
  rangeZ /= 3.0f;

  const int nzp = pGrid->getNzp(); const int nzp2 = nzp/2;
  const int nyp = pGrid->getNyp(); const int nyp2 = nyp/2;
  const int nxp = pGrid->getNxp(); const int nxp2 = nxp/2;

  for (k = -nzp2; k < nzp2; k++) {
    int k1 = (k < 0 ? nzp + k : k);
    float deltaZ = static_cast<float>(k*simbox_.getdz());

    for (j = -nyp2; j < nyp2; j++) {
      int j1 = (j < 0 ? nyp + j : j);
      float deltaY = static_cast<float>(j*simbox_.getdy());

      for (i = -nxp2; i < nxp2; i++) {
        int ii1 = (i < 0 ? nxp + i : i);
        float deltaX = static_cast<float>(i*simbox_.getdx());
        const float h = float(sqrt(deltaX*deltaX/(rangeX*rangeX) + deltaY*deltaY/(rangeY*rangeY) +
          deltaZ*deltaZ/(rangeZ*rangeZ)));
        const float gamma = float(exp(-pow(h,power)));

        pGrid->setRealValue(ii1, j1, k1, gamma);

      } // end i
    } //end j
  } // endk;
  pGrid->endAccess();
  return *pGrid;
}


const FFTGrid&
CKrigingAdmin::CreateAndFillFFTGridWithCovRot(int i1) {
  FFTGrid* pGrid = new FFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(),
    simbox_.getnx(), simbox_.getny(), simbox_.getnz());

  pGrid->createRealGrid();
  pGrid->setAccessMode(FFTGrid::RANDOMACCESS);
  float rangeXf = 0.0f, rangeYf = 0.0f, rangeZf = 0.0f;
  float power = 1.0f;
  switch (i1) {
  case 1:
    rangeXf = 2.0f; rangeYf = 4.0f; rangeZf = 2.0f;//2.5f;
    power = 1.0f;
    break;
  case 2:
    rangeXf = 5.0f; rangeYf = 8.0f; rangeZf = 3.5f;
    power = 1.0f;
    break;
  case 3:
    rangeXf = 3.0f; rangeYf = 16.0f; rangeZf = 2.5f;
    power = 1.0f;
    break;
  default:
    Require(false, "switch failed");
  }// end switch
  //rangeXf = 2.0f; rangeYf = 4.0f; rangeZf = 2.0f;//2.5f;
  int i, j, k;
  float rangeX = static_cast<float>(simbox_.getlx())/rangeXf;
  float rangeY = static_cast<float>(simbox_.getly())/rangeYf;
  float rangeZ = static_cast<float>(simbox_.getlz())/rangeZf;

  // need to take into acccount practical range
  rangeX /= 3.0f;
  rangeY /= 3.0f;
  rangeZ /= 3.0f;

  const int nzp = pGrid->getNzp(); const int nzp2 = nzp/2;
  const int nyp = pGrid->getNyp(); const int nyp2 = nyp/2;
  const int nxp = pGrid->getNxp(); const int nxp2 = nxp/2;

  float rotAngle = 45.0f;
  rotAngle *= 3.14159f/180.0f;
  float cosA = float(cos(rotAngle));
  float sinA = float(sin(rotAngle));
  float rotMatrix[3][3];
  rotMatrix[0][0] = cosA; rotMatrix[0][1] = sinA; rotMatrix[0][2] = 0.0f;
  rotMatrix[1][0] = -sinA; rotMatrix[1][1] = cosA; rotMatrix[1][2] = 0.0f;
  rotMatrix[2][0] = 0.0f; rotMatrix[2][1] = 0.0f; rotMatrix[2][2] = 1.0f;

  for (k = -nzp2; k < nzp2; k++) {
    int k1 = (k < 0 ? nzp + k : k);

    for (j = -nyp2; j < nyp2; j++) {
      int j1 = (j < 0 ? nyp + j : j);
      for (i = -nxp2; i < nxp2; i++) {
        int ii1 = (i < 0 ? nxp + i : i);
        float deltaX = static_cast<float>(i*simbox_.getdx());
        float deltaY = static_cast<float>(j*simbox_.getdy());
        float deltaZ = static_cast<float>(k*simbox_.getdz());
        //deltaX -= nxp2; deltaY -= nyp2; deltaZ -= nzp2;
        RotateVec(deltaX, deltaY, deltaZ, rotMatrix);
        //deltaX += nxp2; deltaY += nyp2; deltaZ += nzp2;
        const float h = float(sqrt(deltaX*deltaX/(rangeX*rangeX) + deltaY*deltaY/(rangeY*rangeY) +
          deltaZ*deltaZ/(rangeZ*rangeZ)));
        const float gamma = float(exp(-pow(h,power)));

        pGrid->setRealValue(ii1, j1, k1, gamma);

      } // end i
    } //end j
  } // endk;
  pGrid->endAccess();
  return *pGrid;
}

const FFTGrid& CKrigingAdmin::CreateAndFillFFTGridWithCrCov() {

  FFTGrid* pGrid = new FFTGrid(simbox_.getnx(), simbox_.getny(), simbox_.getnz(),
    simbox_.getnx(), simbox_.getny(), simbox_.getnz());

  pGrid->createRealGrid();
  pGrid->setAccessMode(FFTGrid::RANDOMACCESS);

  int i, j, k;

  const int nzp = pGrid->getNzp(); const int nzp2 = nzp/2;
  const int nyp = pGrid->getNyp(); const int nyp2 = nyp/2;
  const int nxp = pGrid->getNxp(); const int nxp2 = nxp/2;

  for (k = -nzp2; k < nzp2; k++) {
    int k1 = (k < 0 ? nzp + k : k);
    //float deltaZ = k*simbox_.getdz();

    for (j = -nyp2; j < nyp2; j++) {
      int j1 = (j < 0 ? nyp + j : j);
      //float deltaY = j*simbox_.getdy();

      for (i = -nxp2; i < nxp2; i++) {
        int i1 = (i < 0 ? nxp + i : i);
        //float deltaX = i*simbox_.getdx();
        const float gamma2 = 0.0f;//exp(-pow(h,power));

        pGrid->setRealValue(i1, j1, k1, gamma2);

      } // end i
    } //end j
  } // endk;
  pGrid->endAccess();
  return *pGrid;
}

void CKrigingAdmin::RotateVec(float& rx, float& ry, float& rz, const float mat[][3]) {
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

void CKrigingAdmin::CalcSmoothWeights(Gamma gamma, int direction) {
  //FFTGrid* pGrid = 0;
  const CovGridSeparated* pCov = 0;
  switch (gamma) {
  case ALPHA_KRIG :
    pCov = &covAlpha_;
    break;
  case BETA_KRIG :
    pCov = &covBeta_;
    break;
  case RHO_KRIG :
    pCov = &covRho_;
    break;
  default :
    // should never arrive here
    Require(false, "switch failed");
  } // end switch gamma

 int size = 0;
  double** ppKrigSmoothWeights = 0;
  switch (direction) {
  case 1 :
    size = GetSmoothBlockNx();
    ppKrigSmoothWeights = ppKrigSmoothWeightsX_;
    break;
  case 2 :
    size = GetSmoothBlockNy();
    ppKrigSmoothWeights = ppKrigSmoothWeightsY_;
    break;
  case 3 :
    size = GetSmoothBlockNz();
    ppKrigSmoothWeights = ppKrigSmoothWeightsZ_;
    break;
  default :
    // should never arrive here
    Require(false, "switch failed");
  } // end switch direction

  // allocate
  int i;

  NRLib::Vector pVec(size);

  for (i = 0; i < size; i++) {
    //if (i > 0) {
    switch (direction) {
    case 1 :
      pVec(i) = pCov->GetGamma2(0, 0, 0, i, 0, 0);
      break;
    case 2 :
      pVec(i) = pCov->GetGamma2(0, 0, 0, 0, i, 0);
      break;
    case 3 :
      pVec(i) = pCov->GetGamma2(0, 0, 0, 0, 0, i);
      break;
    default :
      // should never arrive here
      Require(false, "switch failed");
    } // end switch direction
    //} // end if
  } // end i
  //pVec[0] = 1.0;
  //pVec[0] = pCov


  NRLib::SymmetricMatrix ppMatrix(size);

  for (int a = 0 ; a < size ; a++) {
    for (int b = 0 ; b <= a ; b++) {
      ppMatrix(b, a) = pVec(std::abs(a - b));
      // noise
      if (a == b && a > 0 && a < size - 1)
        ppMatrix(a,a) *= 1.2;
    }
  }

  int counter = 0;
  static const double cholesky_repair_factor = 1.01;
  bool robustify = false;

  NRLib::SymmetricMatrix ppMatrix2 = ppMatrix;
  NRLib::Matrix          ppInv     = NRLib::IdentityMatrix(size);

  try {
    NRLib::CholeskySolve(ppMatrix2, ppInv);
  }
  catch (NRLib::IOError e) {
    robustify = true;
  }

  while (robustify && counter < maxCholeskyLoopCounter_) {
    for (i = 0 ; i < size ; i++) {
      ppMatrix(i,i) *= cholesky_repair_factor;
    }
    ppMatrix2 = ppMatrix;
    robustify = false;
    try {
      NRLib::CholeskySolve(ppMatrix2, ppInv);
    }
    catch (NRLib::Exception & e) {
      (void) e;
      robustify = true;
    }
    counter++;
  }

  Require(counter <= maxCholeskyLoopCounter_, "counter <= maxCholeskyLoopCounter_");

  NRLib::Vector p(size);

  for (int krigI = 0 ; krigI < size - 2 ; krigI++) {
    for (int b = 0 ; b < size ; b++)
      p(b) = pVec(std::abs(b - krigI - 1));

    NRLib::Vector x = ppInv * p;

    double * const pWeights = ppKrigSmoothWeights[krigI];
    for (int b = 0; b < size; b++)
      pWeights[b] = x(b);
  }
}

void CKrigingAdmin::SmoothKrigedResult(Gamma gamma) {
  CalcSmoothWeights(gamma, 1);
  CalcSmoothWeights(gamma, 2);
  CalcSmoothWeights(gamma, 3);

  const int nxBlock = NBlocks(dxBlock_, simbox_.getnx());
  const int nyBlock = NBlocks(dyBlock_, simbox_.getny());
  const int nzBlock = NBlocks(dzBlock_, simbox_.getnz());

  FFTGrid* pGrid = 0;
  switch(gamma) {
  case ALPHA_KRIG :
    pGrid = trendAlpha_;
    break;
  case BETA_KRIG :
    pGrid = trendBeta_;
    break;
  case RHO_KRIG :
    pGrid = trendRho_;
    break;
  default :
    // should never arrive here
    Require(false, "switch failed");
  } // end switch

  // loop over all kriging blocks
  int i,j,k;
  for (k = 0; k < nzBlock; k++) {
    int k1 = k*dzBlock_;
    const int k2Max = std::min(k1 + dzBlock_, simbox_.getnz());
    for (j = 0; j < nyBlock; j++) {
      int j1 = j*dyBlock_;
      const int j2Max = std::min(j1 + dyBlock_, simbox_.getny());
      for (i = 0; i < nxBlock; i++) {
        int i1 = i*dxBlock_;
        const int i2Max = std::min(i1 + dxBlock_, simbox_.getnx());
        int j2, k2,i2;
        // smooth in X direction
        if(dxBlock_>1)
        {
        const int i2Start = i1 + dxBlock_ - dxSmoothBlock_;
        const int i2End = std::min(i1 + dxBlock_ + dxSmoothBlock_, simbox_.getnx());
        const int c2EndX = i2End;

        for (k2 = k1; k2 < k2Max; k2++) {
          for (j2 = j1; j2 < j2Max; j2++) {
            // check for well obs
        //    bool foundObs = false;
       //     for (i2 = i2Start; i2 < i2End; i2++) {
       //       if (pBWellGrid_->getRealValue(i2, j2, k2) == 1.0f) {
       //         foundObs = true; break;
       //       }
       //     } // end i2

            // actually smoothing
          //  if (!foundObs) {
              for (i2 = i2Start; i2 < i2End; i2++) {
                if (pBWellGrid_->getRealValue(i2, j2, k2) != 1.0f) {
                float result = pGrid->getRealValue(i2, j2, k2);
                if (result != RMISSING) {
                  const float trend = result;
                  int c;
                  for (c = i2Start - 1; c <= c2EndX; c++) {
                 // for(c = i2 - dxSmoothBlock_; c<= i2+dxSmoothBlock_;c++){
                    int c2 = (c >= simbox_.getnx() ? 2*simbox_.getnx() - c - 1 : c);
                    c2 = (c2<0 ? 0 : c2);
                    result += static_cast<float>(ppKrigSmoothWeightsX_[i2-i2Start][c - i2Start + 1])
                  //  result += static_cast<float>(ppKrigSmoothWeightsX_[1][c - i2 + dxSmoothBlock_])
                      * (pGrid->getRealValue(c2, j2, k2) - trend);
                  } // end c
                  if (pGrid->setRealValue(i2, j2, k2, result))
                    //if (pGrid->setRealValue(i2, j2, k2, 1.0f))
                    Require(false, "pGrid->setRealValue failed"); // something is serious wrong...
                }
              } // end if
            } // end i2
          } // end j2
        } // end k2
        }
        if(dyBlock_>1)
        {
        // smooth in y direction
        int j2Start = j1 + dyBlock_ - dySmoothBlock_;
        const int j2End = std::min(j1 + dyBlock_ + dySmoothBlock_, simbox_.getny());
        const int c2EndY = j2End;

        for (k2 = k1; k2 < k2Max; k2++) {
          for (i2 = i1; i2 < i2Max; i2++) {
            // check for well obs
        //    bool foundObs = false;
       //     for (j2 = j2Start; j2 < j2End; j2++) {
        //      if (pBWellGrid_->getRealValue(i2, j2, k2) == 1.0f) {
       //         foundObs = true; break;
        //      }
        //    } // end j2

            // actually smoothing
          //  if (!foundObs) {
              for (j2 = j2Start; j2 < j2End; j2++) {
              if (pBWellGrid_->getRealValue(i2, j2, k2) != 1.0f) {
                float result = pGrid->getRealValue(i2, j2, k2);
                if (result != RMISSING) {
                  const float trend = result;
                  int c;
                  for (c = j2Start - 1; c <= c2EndY; c++) {
                    //for(c = j2 - dySmoothBlock_; c<= j2+dySmoothBlock_;c++){
                    int c2 = (c >= simbox_.getny() ? 2*simbox_.getny() - c - 1 : c);
                    c2 = (c2<0 ? 0 : c2);
                    result += static_cast<float>(ppKrigSmoothWeightsY_[j2 - j2Start][c - j2Start + 1])
                   // result += static_cast<float>(ppKrigSmoothWeightsY_[1][c - j2 + dySmoothBlock_])
                      * (pGrid->getRealValue(i2, c2, k2) - trend);
                  } // end c
                  if (pGrid->setRealValue(i2, j2, k2, result))
                    //if (pGrid->setRealValue(i2, j2, k2, 1.0f))
                    Require(false, "pGrid->setRealValue failed"); // something is serious wrong...
                }
              } // end j2
            } // end if
          } // end i2
        } // end k2

        }
        if(dzBlock_>1)
        {
        // smooth in z direction
        const int k2Start = k1 + dzBlock_ - dzSmoothBlock_;
        const int k2End = std::min(k1 + dzBlock_ + dzSmoothBlock_, simbox_.getnz());
        const int c2EndZ = k2End;
        for (j2 = j1; j2 < j2Max; j2++) {
          for (i2 = i1; i2 < i2Max; i2++) {
            // check for well obs
         //   bool foundObs = false;
         //   for (k2 = k2Start; k2 <= k2End; k2++) {
          //    if (pBWellGrid_->getRealValue(i2, j2, k2) == 1.0f) {
         //       foundObs = true; break;
         //     }
         //   } // end k2

            // actually smoothing
          //  if (!foundObs) {
              for (k2 = k2Start; k2 < k2End; k2++) {
              if (pBWellGrid_->getRealValue(i2, j2, k2) != 1.0f) {
                float result = pGrid->getRealValue(i2, j2, k2);
                if (result != RMISSING) {
                  const float trend = result;
                  int c;
                  for (c = k2Start - 1; c <= c2EndZ; c++) {
                  //for(c = k2 - dzSmoothBlock_; c<= k2+dzSmoothBlock_;c++){
                    int c2 = (c >= simbox_.getnz() ? 2*simbox_.getnz() - c - 1 : c);
                    c2 = (c2<0 ? 0 : c2);
                    result += static_cast<float>(ppKrigSmoothWeightsZ_[k2 - k2Start][c - k2Start + 1])
                   // result += static_cast<float>(ppKrigSmoothWeightsZ_[1][c - k2 + dzSmoothBlock_])
                      * (pGrid->getRealValue(i2, j2, c2) - trend);
                  } // end c
                  if (pGrid->setRealValue(i2, j2, k2, result))
                    //if (pGrid->setRealValue(i2, j2, k2, 1.0f))
                    Require(false, "pGrid->setRealValue failed"); // something is serious wrong...
                }
              } // end k2
            } // end if
          } // end i2
        } // end j2
        }
      } // end i
    } // end j
  } // end k
}

void CKrigingAdmin::WriteDebugOutput() const {
  LogKit::LogFormatted(LogKit::DebugHigh,"Inside CKrigingAdmin::WriteDebugOutput, before KrigAll is called\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"number of cells to define a kriging block\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"dxBlock_: %d\n", dxBlock_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dyBlock_: %d\n", dyBlock_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dzBlock_: %d\n", dzBlock_);

  LogKit::LogFormatted(LogKit::DebugHigh,"number of additional cells to reach data neighbourhood\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"dxBlockExt_: %d\n", dxBlockExt_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dyBlockExt_: %d\n", dyBlockExt_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dzBlockExt_: %d\n", dzBlockExt_);

  LogKit::LogFormatted(LogKit::DebugHigh,"dataTarget_: %d\n", dataTarget_);

  LogKit::LogFormatted(LogKit::DebugHigh,"Estimated ranges from covariance cubes\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeAlphaX_: %d\n", rangeAlphaX_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeAlphaY_: %d\n", rangeAlphaY_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeAlphaZ_: %d\n", rangeAlphaZ_);

  LogKit::LogFormatted(LogKit::DebugHigh,"rangeBetaX_: %d\n", rangeBetaX_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeBetaY_: %d\n", rangeBetaY_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeBetaZ_: %d\n", rangeBetaZ_);

  LogKit::LogFormatted(LogKit::DebugHigh,"rangeRhoX_: %d\n", rangeRhoX_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeRhoY_: %d\n", rangeRhoY_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeRhoZ_: %d\n", rangeRhoZ_);

  LogKit::LogFormatted(LogKit::DebugHigh,"rangeX_: %f\n", rangeX_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeY_: %f\n", rangeY_);
  LogKit::LogFormatted(LogKit::DebugHigh,"rangeZ_: %f\n", rangeZ_);

  LogKit::LogFormatted(LogKit::DebugHigh,"Bool flags\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"failed2EstimateRange_: %d\n", static_cast<int>(failed2EstimateRange_));
  LogKit::LogFormatted(LogKit::DebugHigh,"failed2EstimateDefaultDataBoxAndBlock_: %d\n", static_cast<int>(failed2EstimateDefaultDataBoxAndBlock_));

  LogKit::LogFormatted(LogKit::DebugHigh,"Smoothing\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"dxSmoothBlock_: %d\n", dxSmoothBlock_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dySmoothBlock_: %d\n", dySmoothBlock_);
  LogKit::LogFormatted(LogKit::DebugHigh,"dzSmoothBlock_: %d\n", dzSmoothBlock_);

  LogKit::LogFormatted(LogKit::DebugHigh,"Kriging Data (from blocked wells)\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"NumberOfObs: %d\n", noData_);
  LogKit::LogFormatted(LogKit::DebugHigh,"NoValidData: %d\n", noValid_);
  LogKit::LogFormatted(LogKit::DebugHigh,"NoValidDataAlpha: %d\n", noValidAlpha_);
  LogKit::LogFormatted(LogKit::DebugHigh,"NoValidDataBeta: %d\n", noValidBeta_);
  LogKit::LogFormatted(LogKit::DebugHigh,"NoValidDataRho: %d\n", noValidRho_);
}

void CKrigingAdmin::WriteDebugOutput2() const {
  LogKit::LogFormatted(LogKit::DebugHigh,"Inside CKrigingAdmin::WriteDebugOutput2, Kriging for a variable is finished\n");
  LogKit::LogFormatted(LogKit::DebugHigh,"noKrigedCells_: %d\n", noKrigedCells_);
  LogKit::LogFormatted(LogKit::DebugHigh,"noEmptyDataBlocks_: %d\n", noEmptyDataBlocks_);
  LogKit::LogFormatted(LogKit::DebugHigh,"noKrigedVariables_: %d\n", noKrigedVariables_);
  LogKit::LogFormatted(LogKit::DebugHigh,"noCholeskyDecomp_: %d\n", noCholeskyDecomp_);
  LogKit::LogFormatted(LogKit::DebugHigh,"noSolvedMatrixEq_: %d\n", noSolvedMatrixEq_);
  LogKit::LogFormatted(LogKit::DebugHigh,"noRMissing_: %d\n", noRMissing_);
}

void CKrigingAdmin::Require(bool test, const std::string msg) const {
  if (!test) {
    LogKit::LogFormatted(LogKit::Low,"Requirement failed: %s\n", msg.c_str());
    exit(1);
  }
}
