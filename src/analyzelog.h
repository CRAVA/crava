/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef ANALYZELOG_H
#define ANALYZELOG_H

#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/welldata.h"
#include "nrlib/grid/grid.hpp"

class BlockedLogsCommon;
class Background;

class Analyzelog {
public:
  Analyzelog(const std::vector<WellData *>                      & wells,
             const std::map<std::string, BlockedLogsCommon *>   & mapped_blocked_logs_for_correlation,
             const std::vector<std::vector<NRLib::Grid<double> > >        & background,
             const Simbox                                       * simbox,
             const ModelSettings                                * model_settings,
             std::string                                        & errTxt);

  Analyzelog(std::vector<WellData *>  & wells,
             Background               * background,
             const Simbox             * simbox,
             const ModelSettings      * modelSettings,
             std::string              & errTxt);

  ~Analyzelog(void);

  float        ** getVar0(void)         const { return Var0_         ;}
  float        ** getPointVar0(void)    const { return pointVar0_    ;}
  float         * getCorrT(void)        const { return CorrT_        ;}
  int             getNumberOfLags(void) const { return numberOfLags_ ;}

private:
  void            estimate(const ModelSettings * modelSettings,
                           Background    * background,
                           std::string   & errTxt);

  void            estimateLnData(float      **& lnData,
                                 FFTGrid      * background,
                                 int            logNr,
                                 std::string  & errTxt);

  void            estimatePointVar0(float       ** Var0,
                                    float       ** lnDataAlpha,
                                    float       ** lnDataBeta,
                                    float       ** lnDataRho,
                                    std::string & errTxt);

  void            estimateCorrTAndVar0(float       * CorrT,
                                       float      ** Var0,
                                       float      ** lnDataAlpha,
                                       float      ** lnDataBeta,
                                       float      ** lnDataRho,
                                       bool          allVsLogsAreSynthetic,
                                       float         dt,
                                       int           n,
                                       int           maxnd,
                                       std::string & errTxt);

  void            readMeanData(FFTGrid *cube, int nd, const double *xpos, const double *ypos,
                               const double *zpos, float *meanValue);

  void            calculateNumberOfLags(int         & numberOfLags,
                                        int         & maxnd,
                                        std::string & errTxt);

  void            findConstructedVsLogs(void);

  void            checkVariances(const ModelSettings  * modelSettings,
                                 const float  * const * pointVar0,
                                 const float  * const * Var0,
                                 float            dt,
                                 std::string    & errTxt);

  const Simbox          * simbox_;
  std::vector<WellData *> wells_;
  int                     nwells_;       // Number of wells

  float        ** Var0_;
  float        ** pointVar0_;
  float         * CorrT_;
  int             numberOfLags_;
};
#endif
