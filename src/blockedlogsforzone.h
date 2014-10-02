/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKEDLOGSFORZONE_H
#define BLOCKEDLOGSFORZONE_H

#include "nrlib/iotools/logkit.hpp"
#include "nrlib/stormgrid/stormcontgrid.hpp"
#include <stdlib.h>
#include <string.h>
#include "fftw.h"
#include "lib/utils.h"

//#include "nrlib/well/well.hpp"
//#include "src/blockedlogscommon.h"

class WellData;
class BlockedLogs;

class BlockedLogsForZone
{
public:
  BlockedLogsForZone(WellData            * well,
                     const StormContGrid & background_grid,
                     const NRLib::Volume & eroded_grid);

  //BlockedLogsForZone(const NRLib::Well   & well,
  //                   const StormContGrid & stormgrid,
  //                   BlockedLogsCommon   * mapped_blocked_logs);

  ~BlockedLogsForZone();

  int                getNumberOfBlocks()          { return nBlocks_                  ;}

  const int *        getIpos()                    { return ipos_                     ;}
  const int *        getJpos()                    { return jpos_                     ;}
  const int *        getKpos()                    { return kpos_                     ;}

  std::vector<float> getAlpha()                   { return alpha_                    ;}
  std::vector<float> getBeta()                    { return beta_                     ;}
  std::vector<float> getRho()                     { return rho_                      ;}

  std::vector<float> getAlphaHighCutBackground()  { return alpha_highcut_background_ ;}
  std::vector<float> getBetaHighCutBackground()   { return beta_highcut_background_  ;}
  std::vector<float> getRhoHighCutBackground()    { return rho_highcut_background_   ;}

  void               getVerticalTrend(const std::vector<float> & blockedLog,
                                      float                    * trend) const;

private:
  void       blockContinuousLog(const std::vector<int>   bInd,
                                const float            * wellLog,
                                std::vector<float>     & blockedLog);

  void       blockContinuousLog(const std::vector<int>    bInd,
                                const double *  wellLog,
                                std::vector<float>     & blockedLog);

  void       findBlockIJK(WellData               * well,
                          const StormContGrid    & stormgrid,
                          const std::vector<int>   bInd);

  //void       findBlockIJK(const StormContGrid    & stormgrid,
  //                        const std::vector<int>   bInd,
  //                        BlockedLogsCommon      * mapped_blocked_log);

  void       findSizeAndBlockPointers(WellData            * well,
                                      const StormContGrid & background_grid,
                                      const NRLib::Volume & eroded_grid,
                                      std::vector<int>    & bInd);

  //void       findSizeAndBlockPointers(const StormContGrid  & stormgrid,
  //                                    std::vector<int>     & bInd,
  //                                    BlockedLogsCommon    * mapped_blocked_log,
  //                                    int                    number_of_nonmissing_welldata);

  std::vector<float>        alpha_;                    ///<
  std::vector<float>        beta_;                     ///< Logs (log-domain)
  std::vector<float>        rho_;                      ///<

  std::vector<float>        alpha_highcut_background_; ///<
  std::vector<float>        beta_highcut_background_;  ///< Logs high-cut filtered to background resolution (log-domain)
  std::vector<float>        rho_highcut_background_;   ///<

  int                     * ipos_;                     ///<
  int                     * jpos_;                     ///< Stormgrid IJK value for block
  int                     * kpos_;                     ///<

  int                       firstM_;                   ///< First well log entry contributing to blocked well
  int                       lastM_;                    ///< Last well log entry contributing to blocked well

  int                       first_eroded_M_;           ///< First well log entry to be used from the eroded grid
  int                       last_eroded_M_;            ///< Last well log entry to be used from the eroded grid

  int                       firstB_;                   ///< First block with contribution from well log
  int                       lastB_;                    ///< Last block with contribution from well log

  int                       nBlocks_;                  ///< Number of blocks
  int                       nLayers_;                  ///< Number of vertical layers (nz)

  double                    dz_;                       ///< Stormgrid dz value for block

};

#endif
