/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef ANALYZELOG_H
#define ANALYZELOG_H

class  FFTGrid;
class  Simbox;
class  WellData;
class  Background;
class  ModelSettings;
struct irapgrid;

class Analyzelog {
public:
  Analyzelog(WellData      ** wells,
             Background     * background,
             Simbox         * simbox,
             ModelSettings  * modelSettings,
             std::string    & errTxt);
  ~Analyzelog(void);

  float        ** getVar0(void)         const { return Var0_         ;}
  float        ** getPointVar0(void)    const { return pointVar0_    ;}
  float         * getCorrT(void)        const { return CorrT_        ;}
  int             getNumberOfLags(void) const { return numberOfLags_ ;}

private:
  void            estimate(ModelSettings * modelSettings,
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

  void            checkVariances(ModelSettings  * modelSettings,
                                 float         ** pointVar0,
                                 float         ** Var0,
                                 float            dt,
                                 std::string    & errTxt);

  const Simbox  * simbox_;
  WellData     ** wells_;
  int             nwells_;       // Number of wells

  float        ** Var0_;
  float        ** pointVar0_;
  float         * CorrT_;
  int             numberOfLags_;
};
#endif
