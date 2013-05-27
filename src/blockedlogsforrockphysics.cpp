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
{
  int nFacies = well->getNFacies();
  facies_names_.resize(nFacies);
  for(int i=0; i<nFacies; i++)
    facies_names_[i] = well->getFaciesName(i);

  int * bInd    = new int[well->getNd()]; // Gives which block each well log entry contributes to
  int   nLayers = simbox->getnz();
  int   nBlocks;
  int   firstM;
  int   lastM;
  int   firstB;
  int   lastB;
  float dz;
  int * ipos;                     ///<
  int * jpos;                     ///< Simbox IJK value for block
  int * kpos;                     ///<


  BlockedLogs::findSizeAndBlockPointers(well, simbox, bInd, nLayers, firstM, lastM, nBlocks);

  BlockedLogs::findBlockIJK(well, simbox, bInd, firstM, lastM, nLayers, nBlocks, ipos, jpos, kpos, dz, firstB, lastB);

  delete [] ipos;
  delete [] jpos;
  delete [] kpos;

  int dummy;

  alpha_.resize(nFacies, std::vector<float>(nBlocks,0));
  beta_.resize(nFacies, std::vector<float>(nBlocks,0));
  rho_.resize(nFacies, std::vector<float>(nBlocks,0));

  blockContinuousLog(bInd, well->getAlpha(dummy), well->getFacies(dummy), well->getFaciesNr(), firstM, lastM, alpha_);
  blockContinuousLog(bInd, well->getBeta(dummy), well->getFacies(dummy), well->getFaciesNr(), firstM, lastM, beta_);
  blockContinuousLog(bInd, well->getRho(dummy), well->getFacies(dummy), well->getFaciesNr(), firstM, lastM, rho_);


  delete [] bInd;
}

//------------------------------------------------------------------------------
BlockedLogsForRockPhysics::~BlockedLogsForRockPhysics(void)
{
}
//------------------------------------------------------------------------------
std::vector<float>
BlockedLogsForRockPhysics::getAlphaForFacies(const std::string & facies_name)
{
  std::vector<float> alpha_given_facies;
  for(size_t i=0; i<facies_names_.size(); i++) {
    if(facies_name == facies_names_[i])
      alpha_given_facies = alpha_[i];
  }

  return alpha_given_facies;
}
//------------------------------------------------------------------------------
void
BlockedLogsForRockPhysics::blockContinuousLog(const int                        * bInd,
                                              const float                      * wellLog,
                                              const int                        * faciesLog,
                                              const int                        * faciesNumbers,
                                              const int                        & firstM,
                                              const int                        & lastM,
                                              std::vector<std::vector<float> > & blockedLog)
{
  if (wellLog != NULL) {
    int nFacies = blockedLog.size();
    int nBlocks = blockedLog[0].size();

    std::vector<std::vector<int> > count(nFacies, std::vector<int>(nBlocks,0));
    //
    // Block log
    //
    for (int m = firstM; m < lastM + 1; m++) {
      if (wellLog[m] != RMISSING && faciesLog[m] != IMISSING) {
        for(int j=0; j<nFacies; j++) {
          if(faciesNumbers[j] == faciesLog[m]) {
            blockedLog[j][bInd[m]] += log(wellLog[m]);
            count[j][bInd[m]]++;
          }
        }
      }
    }
    for (int l = 0 ; l < nBlocks; l++) {
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
