/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef SPATIALWELLFILTER_H
#define SPATIALWELLFILTER_H

#include "src/definitions.h"
#include "rplib/syntwelldata.h"

class FFTGrid;
class AVOInversion;
//class WellData;
//class BlockedLogs;
class BlockedLogsCommon;
class SeismicParametersHolder;

class SpatialWellFilter
{
public:
  SpatialWellFilter(int nwells);

  ~SpatialWellFilter();

  void                     setPriorSpatialCorrSyntWell(FFTGrid             * parSpatialCorr,
                                                       SyntWellData        * well,
                                                       int                   wellnr);

  //void                     setPriorSpatialCorr(FFTGrid    * parSpatialCorr,
  //                                             WellData   * well,
  //                                             int          wellnr);

  void                     setPriorSpatialCorr(FFTGrid           * parSpatialCorr,
                                               BlockedLogsCommon * blocked_log,
                                               int                 wellnr);

  //void                     doFiltering(std::vector<WellData *>         wells,
  //                                     int                             nWells,
  //                                     bool                            useVpRhoFilter,
  //                                     int                             nAngles,
  //                                     const Crava                   * cravaResult,
  //                                     const std::vector<Grid2D *>   & noiseScale,
  //                                     SeismicParametersHolder       & seismicParameters);

  void                    doFiltering(std::map<std::string, BlockedLogsCommon *> blocked_logs,
                                      //int                                        nWells,
                                      bool                                       useVpRhoFilter,
                                      int                                        nAngles,
                                      const AVOInversion                       * avoInversionResult,
                                      const std::vector<Grid2D *>              & noiseScale,
                                      SeismicParametersHolder                  & seismicParameters);

  void                     doFilteringSyntWells(std::vector<SyntWellData *>              & syntWellData,
                                                const std::vector<std::vector<double> >  & v,
                                                SeismicParametersHolder                  & seismicParameters,
                                                int                                        nWells,
                                                const NRLib::Matrix                      & priorVar0);

  //void                         doFiltering(WellData                   ** wells,
  //                                         int                           nWells,
  //                                         bool                          useVpRhoFilter,
  //                                         int                           nAngles,
  //                                         const AVOInversion          * avoInversionResult,
  //                                         const std::vector<Grid2D *> & noiseScale);

  std::vector<NRLib::Matrix> & getSigmae(void)     { return sigmae_ ;}
  std::vector<double **>     & getSigmaeSynt()     {return sigmaeSynt_;}


private:

  //void doVpRhoFiltering(std::vector<NRLib::Matrix> &  sigmaeVpRho,
  //                      double                     ** sigmapri,
  //                      double                     ** sigmapost,
  //                      const int                     n,
  //                      BlockedLogs                 * blockedLogs);

  void doVpRhoFiltering(std::vector<NRLib::Matrix> &  sigmaeVpRho,
                        double                     ** sigmapri,
                        double                     ** sigmapost,
                        const int                     n,
                        BlockedLogsCommon           * blockedLogs);

  void updateSigmaE(NRLib::Matrix       & sigmae,
                    const NRLib::Matrix & filter,
                    const NRLib::Matrix & Spost,
                    int                   n);

  void updateSigmaeSynt(double ** filter, double ** postCov,  int n);

  void completeSigmaE(std::vector<NRLib::Matrix>  & sigmae,
                      int                           lastn,
                      const AVOInversion          * avoInversionResult,
                      const std::vector<Grid2D *> & noiseScale);

  void updateSigmaEVpRho(std::vector<NRLib::Matrix> & sigmaeVpRho,
                         const NRLib::Matrix        & Aw,
                         const NRLib::Matrix        & Spost,
                         int                          nDim,
                         int                          n);

  void completeSigmaEVpRho(std::vector<NRLib::Matrix>  & sigmaeVpRho,
                           int                           lastn,
                           const AVOInversion          * avoInversionResult,
                           const std::vector<Grid2D *> & noiseScale);

  void fillValuesInSigmapostSyntWell(double     ** sigmapost,
                                     const int  *  ipos,
                                     const int  *  jpos,
                                     const int  *  kpos,
                                     FFTGrid    *  covgrid,
                                     int           n,
                                     int           ni,
                                     int           nj);

  void computeSigmaEAdjusted(NRLib::Matrix & sigmae,
                             NRLib::Matrix & sigmaE0,
                             NRLib::Matrix & sigmaETmp,
                             NRLib::Matrix & sigmaEAdj);

  void adjustDiagSigma(NRLib::Matrix & sigmae);

  //void calculateFilteredLogs(const NRLib::Matrix & Aw,
  //                           BlockedLogs         * blockedlogs,
  //                           int                   n,
  //                           bool                  useVs);

  void calculateFilteredLogs(const NRLib::Matrix & Aw,
                             BlockedLogsCommon   * blockedlogs,
                             int                   n,
                             bool                  useVs);

  //void MakeInterpolatedResiduals(const float   * bwLog,
  //                               const float   * bwLogBG,
  //                               const int       n,
  //                               const int       offset,
  //                               NRLib::Vector & residuals);

  void MakeInterpolatedResiduals(const std::vector<double> & bwLog,
                                 const std::vector<double> & bwLogBG,
                                 const int                   n,
                                 const int                   offset,
                                 NRLib::Vector &             residuals);

  void fillValuesInSigmapost(double    ** sigmapost,
                             const int *  ipos,
                             const int *  jpos,
                             const int *  kpos,
                             FFTGrid   *  covgrid,
                             int          n,
                             int          ni,
                             int          nj);

  std::vector<NRLib::Matrix> sigmae_;
  std::vector<double **> sigmaeSynt_;

  int                        nData_;   ///< sum no blocks in all wells
  int                        nWells_;
  int                    *   n_;
  double                 *** priorSpatialCorr_;

  std::vector<double **>     sigmaeVpRho_;

};
#endif

