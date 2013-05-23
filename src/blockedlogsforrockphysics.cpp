/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#include <iostream>
#include <iomanip>
#include <math.h>

#include "nrlib/flens/nrlib_flens.hpp"

#include "nrlib/iotools/logkit.hpp"
#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "src/definitions.h"
#include "src/blockedlogsforrockphysics.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/io.h"


BlockedLogsForRockPhysics::BlockedLogsForRockPhysics(WellData  * well,
                                                     Simbox    * simbox)

: firstM_(IMISSING),
  lastM_(IMISSING),
  nLayers_(simbox->getnz())
{
  int * bInd = new int[well->getNd()]; // Gives which block each well log entry contributes to

  BlockedLogs::findSizeAndBlockPointers(well, simbox, bInd, nLayers_, firstM_, lastM_, nBlocks_);

  BlockedLogs::findBlockIJK(well, simbox, bInd, firstM_, lastM_, nLayers_, nBlocks_, ipos_, jpos_, kpos_, dz_, firstB_, lastB_);

  int dummy;

  int nFacies = well->getNFacies();

  alpha_.resize(nFacies, NULL);
  beta_.resize(nFacies, NULL);
  rho_.resize(nFacies, NULL);

  alpha_highcut_background_.resize(nFacies, NULL);
  beta_highcut_background_.resize(nFacies, NULL);
  rho_highcut_background_.resize(nFacies, NULL);

  blockContinuousLog(bInd, well->getAlpha(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, alpha_);
  blockContinuousLog(bInd, well->getBeta(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, beta_);
  blockContinuousLog(bInd, well->getRho(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, rho_);

  blockContinuousLog(bInd, well->getAlphaBackgroundResolution(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, alpha_highcut_background_);
  blockContinuousLog(bInd, well->getBetaBackgroundResolution(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, beta_highcut_background_);
  blockContinuousLog(bInd, well->getRhoBackgroundResolution(dummy), well->getFacies(dummy), well->getFaciesNr(), nFacies, rho_highcut_background_);

  delete [] bInd;
}

//------------------------------------------------------------------------------
BlockedLogsForRockPhysics::~BlockedLogsForRockPhysics(void)
{
  if (ipos_ != NULL)
    delete [] ipos_;
  if (jpos_ != NULL)
    delete [] jpos_;
  if (kpos_ != NULL)
    delete [] kpos_;

  for(size_t i = 0; i< alpha_.size(); i++) {
    if (alpha_[i] != NULL)
      delete [] alpha_[i];
    if (beta_[i] != NULL)
      delete [] beta_[i];
    if (rho_[i] != NULL)
      delete [] rho_[i];

    if (alpha_highcut_background_[i] != NULL)
      delete [] alpha_highcut_background_[i];
    if (beta_highcut_background_[i] != NULL)
      delete [] beta_highcut_background_[i];
    if (rho_highcut_background_[i] != NULL)
      delete [] rho_highcut_background_[i];
  }

}
//------------------------------------------------------------------------------
void
BlockedLogsForRockPhysics::blockContinuousLog(const int            * bInd,
                                              const float          * wellLog,
                                              const int            * faciesLog,
                                              const int            * faciesNumbers,
                                              const int            & nFacies,
                                              std::vector<float *> & blockedLog)
{
  if (wellLog != NULL) {
    std::vector<std::vector<int> > count(nFacies, std::vector<int>(nBlocks_,0));
    for(int i=0; i<nFacies; i++) {
      float * blockedFaciesLog = new float[nBlocks_];
      //
      // Initialise arrays
      //
      for (int l = 0 ; l < nBlocks_ ; l++)
        blockedFaciesLog[l] = 0.0f;

      blockedLog[i] = blockedFaciesLog;
    }

    //
    // Block log
    //
    for (int m = firstM_ ; m < lastM_ + 1 ; m++) {
      if (wellLog[m] != RMISSING && faciesLog[m] != IMISSING) {
        for(int j=0; j<nFacies; j++) {
          if(faciesNumbers[j] == faciesLog[m]) {
            blockedLog[j][bInd[m]] += log(wellLog[m]);
            count[j][bInd[m]]++;
          }
        }
      }
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      for(int j=0; j<nFacies; j++) {
        if(faciesNumbers[j] == faciesLog[l]) {
          if (count[j][l] > 0)
            blockedLog[j][l] /= count[j][l];
          else
            blockedLog[j][l] = RMISSING;
        }
      }
    }
  }
}
