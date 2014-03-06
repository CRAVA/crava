/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef ANALYZELOG_H
#define ANALYZELOG_H

#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/welldata.h"
#include "nrlib/grid/grid.hpp"
#include "nrlib/well/well.hpp"
#include "src/multiintervalgrid.h"

class BlockedLogsCommon;
class Background;

class Analyzelog {

public:
  Analyzelog(const std::vector<NRLib::Well>                           & wells,
             const std::map<std::string, BlockedLogsCommon *>         & mapped_blocked_logs_for_correlation,
             const std::vector<NRLib::Grid<double> >                  & background,
             const std::vector<NRLib::Grid<double> >                  & background_max_Hz,
             const Simbox                                             * simbox,
             const ModelSettings                                      * model_settings,
             std::string                                              & err_txt);

  Analyzelog(const std::vector<NRLib::Well>                           & wells,
             const std::map<std::string, BlockedLogsCommon *>         & mapped_blocked_logs_for_correlation,
             const std::vector<std::vector<NRLib::Grid<double> > >    & background,
             const std::vector<std::vector<NRLib::Grid<double> > >    & background_max_Hz,
             const std::vector<Simbox>                                & interval_simboxes,
             const ModelSettings                                      * model_settings,
             std::string                                              & err_txt);

  /*
  Analyzelog(std::vector<WellData *>  & wells,
             Background               * background,
             const Simbox             * simbox,
             const ModelSettings      * modelSettings,
             std::string              & errTxt);
  */

  ~Analyzelog(void);

  //NRLib::Grid2D<double> & GetVar0()     const { return Var0_         ;}
  const NRLib::Matrix       & GetVar0(void)         const { return var_0_                           ;}
  const NRLib::Matrix       & GetPointVar0(void)    const { return point_var_0_                     ;}
  const std::vector<float>  & GetCorrT(void)        const { return corr_T_                          ;}
  int                         GetNumberOfLags(void) const { return n_lags_                  ;}
  bool                        GetEnoughData()       const { return enough_data_for_corr_estimation_ ;}

private:

  void            EstimateLnData(std::vector<std::vector<float> >                      & log_data,
                                 const NRLib::Grid<double>                             * background,
                                 const std::vector<std::string>                        & well_names,
                                 const std::map<std::string, BlockedLogsCommon *>      & mapped_blocked_logs_for_correlation,
                                 const std::string                                     & interval_name,
                                 const std::string                                     & log_name,
                                 std::string                                           & err_txt);

  void            EstimateCorrelation(const ModelSettings                                      * model_settings,
                                      const std::vector<NRLib::Well>                           & wells,
                                      const std::vector<std::string>                           & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>         & mapped_blocked_logs_for_correlation,
                                      std::string                                              & interval_name,
                                      const Simbox                                             * interval_simbox,
                                      bool                                                     & enough_data_for_corr_estimation,
                                      const std::vector<NRLib::Grid<double> >                  & background,
                                      std::string                                              & errTxt);

    void          EstimateCorrelation(const ModelSettings                                      * model_settings,
                                      const std::vector<NRLib::Well>                           & wells,
                                      const std::vector<std::string>                           & well_names,
                                      const std::map<std::string, BlockedLogsCommon *>         & mapped_blocked_logs_for_correlation,
                                      const std::vector<Simbox>                                & interval_simboxes,
                                      bool                                                     & enough_data_for_corr_estimation,
                                      const std::vector<std::vector<NRLib::Grid<double> > >    & background,
                                      std::string                                              & errTxt);

  void            estimate(const ModelSettings  * modelSettings,
                           Background           * background,
                           std::string          & errTxt);

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

  void            CalculateNumberOfLags(int                               & numberOfLags,
                                        const std::vector<NRLib::Well>    & wells,
                                        int                               & maxnd,
                                        const Simbox                      * simbox,
                                        std::string                       & errTxt);

  void            findConstructedVsLogs(void);

  void            checkVariances(const ModelSettings  * modelSettings,
                                 const float  * const * pointVar0,
                                 const float  * const * Var0,
                                 float            dt,
                                 std::string    & errTxt);

  int                                                   min_blocks_with_data_for_corr_estim_;
  bool                                                  enough_data_for_corr_estimation_;       ///< Vector with size == n_intervals
  //std::map<std::string, BlockedLogsCommon *>            mapped_blocked_logs_for_correlation_;
  //const Simbox                                        * simbox_;
  //std::vector<WellData *>                               wells_;
  std::vector<std::string>                              well_names_;
  int                                                   n_wells_;       // Number of wells

  std::string            interval_name_;
  NRLib::Matrix          var_0_;
  NRLib::Matrix          point_var_0_;
  std::vector<float>     corr_T_;
  int                    n_lags_;

  //float        ** Var0_;
  //float        ** pointVar0_;
  //float         * CorrT_;
  //int             numberOfLags_;
};
#endif
