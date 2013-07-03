/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <iomanip>
#include <math.h>

#include "nrlib/iotools/logkit.hpp"
#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "src/definitions.h"
#include "src/blockedlogsforzone.h"
#include "src/welldata.h"
#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/io.h"

BlockedLogsForZone::BlockedLogsForZone(WellData            * well,
                                       const StormContGrid & background_grid,
                                       const StormContGrid & eroded_grid)
: firstM_(IMISSING),
  lastM_(IMISSING)
{
  nLayers_ = static_cast<int>(background_grid.GetNK());
  dz_      = static_cast<double>(background_grid.GetLZ()/background_grid.GetNK());

  std::vector<int> bInd(well->getNd()); // Gives which block each well log entry contributes to

  findSizeAndBlockPointers(well, background_grid, eroded_grid, bInd);
  findBlockIJK(well, background_grid, bInd);

  int dummy;

  blockContinuousLog(bInd, well->getAlpha(dummy), alpha_);
  blockContinuousLog(bInd, well->getBeta(dummy), beta_);
  blockContinuousLog(bInd, well->getRho(dummy), rho_);

  blockContinuousLog(bInd, well->getAlphaBackgroundResolution(dummy), alpha_highcut_background_);
  blockContinuousLog(bInd, well->getBetaBackgroundResolution(dummy), beta_highcut_background_);
  blockContinuousLog(bInd, well->getRhoBackgroundResolution(dummy), rho_highcut_background_);

}

//------------------------------------------------------------------------------

BlockedLogsForZone::~BlockedLogsForZone()
{
  delete [] ipos_;
  delete [] jpos_;
  delete [] kpos_;
}

//------------------------------------------------------------------------------

void
BlockedLogsForZone::findBlockIJK(WellData               * well,
                                 const StormContGrid    & stormgrid,
                                 const std::vector<int>   bInd)
{
  ipos_ = new int[nBlocks_];
  jpos_ = new int[nBlocks_];
  kpos_ = new int[nBlocks_];

  int   dummy;
  const double * x = well->getXpos(dummy);
  const double * y = well->getYpos(dummy);
  const double * z = well->getZpos(dummy);

  //
  // Set IJK for virtual part of well in upper part of stormgrid
  //
  int b = -1; // block counter;
  size_t firstI;
  size_t firstJ;
  size_t firstK;
  stormgrid.FindIndex(x[firstM_], y[firstM_], z[firstM_], firstI, firstJ, firstK);

  for (size_t k = 0 ; k < firstK ; k++) {
    b++;
    ipos_[b] = static_cast<int>(firstI);
    jpos_[b] = static_cast<int>(firstJ);
    kpos_[b] = static_cast<int>(k);
  }

  //
  // Set IJK for the defined part of the well
  //
  b = static_cast<int>(firstK);
  ipos_[b] = static_cast<int>(firstI);
  jpos_[b] = static_cast<int>(firstJ);
  kpos_[b] = static_cast<int>(firstK);
  size_t i, j, k;
  for (int m = firstM_ + 1 ; m < lastM_ + 1 ; m++) {
    if (bInd[m] != bInd[m - 1]) {
      b++;
      stormgrid.FindIndex(x[m], y[m], z[m], i, j, k);
      ipos_[b] = static_cast<int>(i);
      jpos_[b] = static_cast<int>(j);
      kpos_[b] = static_cast<int>(k);
    }
  }
  firstB_ = static_cast<int>(firstK);
  lastB_  = b;

  //
  // Set IJK for the virtual part of well in lower part of simbox
  //
  size_t lastI,  lastJ,  lastK;
  stormgrid.FindIndex(x[lastM_], y[lastM_], z[lastM_], lastI, lastJ, lastK);

  for (int k = static_cast<int>(lastK) + 1 ; k < nLayers_ ; k++) {
    b++;
    ipos_[b] = static_cast<int>(lastI);
    jpos_[b] = static_cast<int>(lastJ);
    kpos_[b] = k;
  }
}

//------------------------------------------------------------------------------

void
BlockedLogsForZone::blockContinuousLog(const std::vector<int>   bInd,
                                       const float            * wellLog,
                                       std::vector<float>     & blockedLog)
{
  if (wellLog != NULL) {
    blockedLog.resize(nBlocks_);
    std::vector<int> count(nBlocks_);
    //
    // Initialise arrays
    //
    for (int l = 0 ; l < nBlocks_ ; l++) {
      blockedLog[l] = 0.0;
      count[l] = 0;
    }
    //
    // Block log
    //
    for (int m = firstM_ ; m < lastM_ + 1 ; m++) {
      if(m >= first_eroded_M_ && m <= last_eroded_M_) {
        if (wellLog[m] != RMISSING) {
          blockedLog[bInd[m]] += log(wellLog[m]);
          count[bInd[m]]++;
        }
      }
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      if (count[l] > 0) {
        blockedLog[l] /= count[l];
      }
      else
        blockedLog[l] = RMISSING;
    }
  }
}

//------------------------------------------------------------------------------

void
BlockedLogsForZone::findSizeAndBlockPointers(WellData            * well,
                                             const StormContGrid & background_grid,
                                             const StormContGrid & eroded_grid,
                                             std::vector<int>    & bInd)
{
  int            dummy;
  int            missing = 99999;
  const int      nd = well->getNd();
  const double * x  = well->getXpos(dummy);
  const double * y  = well->getYpos(dummy);
  const double * z  = well->getZpos(dummy);

  //
  // Find first cell in eroded_grid that the well hits
  //
  bool   inside = false;
  size_t firstI = missing;
  size_t firstJ = missing;
  size_t firstK = missing;

  for(int m=0; m<nd; m++) {
    inside = eroded_grid.IsInside(x[m], y[m], z[m]);
    if(inside == true) {
      eroded_grid.FindIndex(x[m], y[m], z[m], firstI, firstJ, firstK);
      first_eroded_M_ = m;
      break;
    }
  }

  //
  // Find last cell in eroded_grid that the well hits
  //
  inside       = false;
  size_t lastI = missing;
  size_t lastJ = missing;
  size_t lastK = missing;

  for (int m=nd-1; m>0; m--) {
    inside = eroded_grid.IsInside(x[m], y[m], z[m]);
    if(inside == true) {
      eroded_grid.FindIndex(x[m], y[m], z[m], lastI, lastJ, lastK);
      last_eroded_M_ = m;
      break;
    }
  }

  //
  // Find first cell in background_grid that the well hits
  //
  inside = false;
  firstI = missing;
  firstJ = missing;
  firstK = missing;

  for(int m=0; m<nd; m++) {
    inside = background_grid.IsInside(x[m], y[m], z[m]);
    if(inside == true) {
      background_grid.FindIndex(x[m], y[m], z[m], firstI, firstJ, firstK);
      firstM_ = m;
      break;
    }
  }
  size_t oldI = firstI;
  size_t oldJ = firstJ;
  size_t oldK = firstK;

  //
  // Find last cell in background_grid that the well hits
  //
  inside = false;
  lastI  = missing;
  lastJ  = missing;
  lastK  = missing;

  for (int m=nd-1; m>0; m--) {
    inside = background_grid.IsInside(x[m], y[m], z[m]);
    if(inside == true) {
      background_grid.FindIndex(x[m], y[m], z[m], lastI, lastJ, lastK);
      lastM_ = m;
      break;
    }
  }

  //
  // Count number of blocks needed for the defined part of well.
  //
  for (int m = 0 ; m < nd ; m++) {
    bInd[m] = IMISSING;
  }
  int nDefinedBlocks = 0;
  bInd[firstM_] = static_cast<int>(firstK); // The first defined well log entry contributes to this block.

  std::vector<int> stormInd(nd);
  const int nx    = static_cast<int>(background_grid.GetNI());
  const int ny    = static_cast<int>(background_grid.GetNJ());
  stormInd[0] = nx*ny*static_cast<int>(oldK) + nx*static_cast<int>(oldJ)+static_cast<int>(oldI);

  size_t newI = missing;
  size_t newJ = missing;
  size_t newK = missing;

  for (int m = firstM_ + 1 ; m < lastM_ + 1 ; m++) {
    background_grid.FindIndex(x[m], y[m], z[m], newI, newJ, newK);

    if (newI != oldI || newJ != oldJ || newK != oldK) {

      int  thisInd = nx*ny*static_cast<int>(newK) + nx*static_cast<int>(newJ)+static_cast<int>(newI);
      bool blockNotListed = true;

      for (int l = 0 ; l < nDefinedBlocks ; l++) {
        if (thisInd == stormInd[l]) {
          blockNotListed = false;
          break;
        }
      }

      if (blockNotListed) {
        stormInd[nDefinedBlocks+1] = thisInd;
        oldI = newI;
        oldJ = newJ;
        oldK = newK;
        nDefinedBlocks++;
      }
    }
    bInd[m] = static_cast<int>(firstK) + nDefinedBlocks;
  }
  nDefinedBlocks++;

  //
  // Why we cannot use nBlocks_ = nDefinedBlocks:
  //
  // When we calculate the background model for each parameter we first
  // estimate a vertical trend in the total volume, and then we interpolate
  // the blocked log in this trend volume. To avoid sharp contrast we
  // ensure that the blocked log is defined from top to base of the volume.
  // In regions where the log is undefined we generate it by kriging from
  // the rest of the log. Likewise, in regions where there is no blocked
  // log at all because the well was too short, we have to make a virtual
  // well.
  //
  nBlocks_ = static_cast<int>(firstK) + nDefinedBlocks + (nLayers_ - static_cast<int>(lastK) - 1);
}

//------------------------------------------------------------------------------

void
BlockedLogsForZone::getVerticalTrend(const std::vector<float> & blockedLog,
                                     float                    * trend) const
{
  if (blockedLog.size() != 0 && trend != NULL) {
    std::vector<int> count(nLayers_);

    for (int k = 0 ; k < nLayers_ ; k++) {
      trend[k] = 0.0;
      count[k] = 0;
    }

    for (int m = 0 ; m < nBlocks_ ; m++) {
      if (blockedLog[m] != RMISSING) {
        trend[kpos_[m]] += blockedLog[m];
        count[kpos_[m]]++;
      }
    }

    for (int k = 0 ; k < nLayers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
  }
  else {
    if (blockedLog.size() == 0)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsForZone::getVerticalTrend(): Trying to use an undefined blocked log\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogsForZone::getVerticalTrend(): Trying to use an undefined trend\n");
    exit(1);
  }
}

