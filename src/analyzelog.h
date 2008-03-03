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
             ModelSettings  * modelSettings); //beregn logdataene her
  ~Analyzelog(void);

  float        ** getVar0(void)         const { return Var0_         ;}
  float        ** getPointVar0(void)    const { return pointVar0_    ;}
  float         * getCorrT(void)        const { return CorrT_        ;}
  int             getNumberOfLags(void) const { return numberOfLags_ ;}

private:
  void            estimate(ModelSettings * modelSettings,
                           Background    * background);
  void            estimateLnData(float   **& lnData, 
                                 FFTGrid   * background, 
                                 int         logNr);
  void            estimatePointVar0(float ** Var0,
                                    float ** lnDataAlpha,
                                    float ** lnDataBeta,
                                    float ** lnDataRho); //Returns 1 if error. NBNB prints own error message
  void            estimateCorrTAndVar0(float  * CorrT, 
                                       float ** Var0, 
                                       float ** lnDataAlpha,
                                       float ** lnDataBeta,
                                       float ** lnDataRho, 
                                       bool     allVsLogsAreSynthetic,
                                       float    dt, 
                                       int      n, 
                                       int      maxnd);
  void            readMeanData(FFTGrid *cube, int nd, const double *xpos, const double *ypos, 
                               const double *zpos, float *meanValue);
  void            calculateNumberOfLags(int & numberOfLags,
                                        int & maxnd);
  void            findConstructedVsLogs(void);
  void            checkVariances(ModelSettings  * modelSettings,
                                 float         ** pointVar0,
                                 float         ** Var0,
                                 float            dt);
  
  const Simbox  * simbox_;
  WellData     ** wells_;
  int             nwells_;       // Number of wells

  float        ** Var0_;
  float        ** pointVar0_;
  float         * CorrT_;
  int             numberOfLags_;
  bool            failed_;
};
#endif
