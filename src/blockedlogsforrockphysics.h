/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKED_LOGS_FOR_ROCK_PHYSICS_H
#define BLOCKED_LOGS_FOR_ROCK_PHYSICS_H

#include "nrlib/iotools/logkit.hpp"
#include <stdlib.h>
#include <string.h>
#include "fftw.h"
#include "lib/utils.h"

class ModelSettings;
class FFTGrid;
class WellData;
class Simbox;
class Wavelet;

class BlockedLogsForRockPhysics
{
public:
  BlockedLogsForRockPhysics(WellData  * well,
                            Simbox    * simbox);

  ~BlockedLogsForRockPhysics(void);

private:
  void                     blockContinuousLog(const int            * bInd,
                                              const float          * wellLog,
                                              const int            * faciesLog,
                                              const int            * faciesNumbers,
                                              const int            & nFacies,
                                              std::vector<float *> & blockedLog);


  int          * ipos_;                     ///<
  int          * jpos_;                     ///< Simbox IJK value for block
  int          * kpos_;                     ///<

  float          dz_;                       ///< Simbox dz value for block

  std::vector<float *> alpha_;                    ///<
  std::vector<float *> beta_;                     ///< Raw logs (log-domain)
  std::vector<float *> rho_;                      ///<

  std::vector<float *> alpha_highcut_background_; ///<
  std::vector<float *> beta_highcut_background_;  ///< Logs high-cut filtered to background resolution (log-domain)
  std::vector<float *> rho_highcut_background_;   ///<


  int            firstM_;                   ///< First well log entry contributing to blocked well
  int            lastM_;                    ///< Last well log entry contributing to blocked well

  int            firstB_;                   ///< First block with contribution from well log
  int            lastB_;                    ///< Last block with contribution from well log

  int            nBlocks_;                  ///< Number of blocks
  int            nLayers_;                  ///< Number of vertical layers (nz)
};

#endif
