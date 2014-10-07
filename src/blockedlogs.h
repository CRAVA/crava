/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef BLOCKEDLOGS_H
#define BLOCKEDLOGS_H

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

class BlockedLogs
{
public:
  BlockedLogs(WellData  * well,
              Simbox    * simbox,
              bool interpolate = false);

  ~BlockedLogs(void);

  const std::string         getWellname(void)                  const { return wellname_ ;}
  const int                 getNumberOfBlocks(void)            const { return nBlocks_  ;}
  const double            * getXpos(void)                      const { return xpos_     ;}
  const double            * getYpos(void)                      const { return ypos_     ;}
  const double            * getZpos(void)                      const { return zpos_     ;}
  const int               * getIpos(void)                      const { return ipos_     ;}
  const int               * getJpos(void)                      const { return jpos_     ;}
  const int               * getKpos(void)                      const { return kpos_     ;}
  const std::vector<int>    getIposVector(void)                const;
  const std::vector<int>    getJposVector(void)                const;
  const std::vector<double> getXposVector(void)                const;
  const std::vector<double> getYposVector(void)                const;
  const std::vector<double> getZposVector(void)                const;
  const float             * getAlpha(void)                     const { return alpha_    ;}
  const float             * getBeta(void)                      const { return beta_     ;}
  const float             * getRho(void)                       const { return rho_      ;}
  const int               * getFacies(void)                    const { return facies_   ;}
  const float               getDz(void)                        const { return dz_       ;}
  const float             * getAlphaHighCutBackground(void)    const { return alpha_highcut_background_ ;}
  const float             * getBetaHighCutBackground(void)     const { return beta_highcut_background_  ;}
  const float             * getRhoHighCutBackground(void)      const { return rho_highcut_background_   ;}
  const float             * getAlphaHighCutSeismic(void)       const { return alpha_highcut_seismic_    ;}
  const float             * getBetaHighCutSeismic(void)        const { return beta_highcut_seismic_     ;}
  const float             * getRhoHighCutSeismic(void)         const { return rho_highcut_seismic_      ;}
  const float             * getAlphaSeismicResolution(void)    const { return alpha_seismic_resolution_ ;}
  const float             * getBetaSeismicResolution(void)     const { return beta_seismic_resolution_  ;}
  const float             * getRhoSeismicResolution(void)      const { return rho_seismic_resolution_   ;}
  const float             * getAlphaForFacies(void)            const { return alpha_for_facies_         ;}
  const float             * getRhoForFacies(void)              const { return rho_for_facies_           ;}
  const float             * getAlphaPredicted(void)            const { return alpha_predicted_          ;}
  const float             * getBetaPredicted(void)             const { return beta_predicted_           ;}
  const float             * getRhoPredicted(void)              const { return rho_predicted_            ;}
  float                  ** getRealSeismicData(void)           const { return real_seismic_data_        ;}
  float                  ** getSyntSeismicData(void)           const { return actual_synt_seismic_data_        ;}
  float                  ** getCpp(void)                       const { return cpp_ ;}
  void                      getVerticalTrend(const float * blockedLog, float * trend);
  void                      getVerticalTrendLimited(const float * blockedLog, float * trend, const std::vector<Surface *> & limits);
  void                      getVerticalTrend(const int * blockedLog,int * trend);
  void                      getBlockedGrid(const FFTGrid * grid, float * blockedLog, int iOffset = 0, int jOffset = 0);
  void                      setLogFromVerticalTrend(float * vertical_trend, double z0, double dz,
                                         int nz, std::string type, int iAngle = IMISSING);
  void                      setLogFromGrid(FFTGrid * grid, int iAngle, int nAngles, std::string type);

  void                      writeWell(ModelSettings * modelSettings, std::vector<std::string> facies_name, std::vector<int> facies_label);       //Use this, will handle formats automatically.
  void                      writeRMSWell(ModelSettings * modelSettings, std::vector<std::string> facies_name, std::vector<int> facies_label);
  void                      writeNorsarWell(ModelSettings * modelSettings);

  void                      setSpatialFilteredLogs(float * filteredlog, int nData, std::string type, const float *bg);

  void                      fillInCpp(const float * coeff,int start,int length,fftw_real* cpp_r,int nzp);
  void                      fillInSeismic(float* seismicData,int start,int length,fftw_real* seis_r,int nzp) const;
  void                      estimateCor(fftw_complex* var1_c,fftw_complex* var2_c,fftw_complex* ccor_1_2_c,int cnzp) const;
  void                      findContinuousPartOfData(const std::vector<bool> & hasData,int nz,int &start,int &length) const;
  void                      findOptimalWellLocation(FFTGrid                   ** seisCube,
                                                    Simbox                     * timeSimbox,
                                                    float                     ** reflCoef,
                                                    int                          nAngles,
                                                    const std::vector<float>   & angleWeight,
                                                    float                        maxShift,
                                                    int                          iMaxOffset,
                                                    int                          jMaxOffset,
                                                    const std::vector<Surface *> limits,
                                                    int                        & iMove,
                                                    int                        & jMove,
                                                    float                      & kMove);
  void                      generateSyntheticSeismic(const float * const * reflCoef,
                                                     int nAngles,
                                                     Wavelet ** wavelet,
                                                     int nz,
                                                     int nzp,
                                                     const Simbox * timeSimbox);
  void                      setSeismicGradient(  double v0,
                                                 const NRLib::Grid2D<float>   &    structureDepthGradX,
                                                 const NRLib::Grid2D<float>   &    structureDepthGradY,
                                                 const NRLib::Grid2D<float>   &    refTimeGradX ,
                                                 const NRLib::Grid2D<float>   &    refTimeGradY,
                                                 std::vector<double>        & xGradient,
                                                 std::vector<double>        & yGradient);

  void                      setTimeGradientSettings(float distance, float sigma_m);
  void                      findSeismicGradient(const FFTGrid     * const * seisCube,
                                                const Simbox              * timeSimbox,
                                                int                        nAngles,
                                                std::vector<double>       & xGradient,
                                                std::vector<double>       & yGradient,
                                                std::vector<std::vector<double> > &Sigma_gradient);

  void           deleteDynamicBlockedLogs();

  static void               findSizeAndBlockPointers(WellData  *  well,
                                                     Simbox    *  simbox,
                                                     int       *& bInd,
                                                     const int &  nLayers,
                                                     int       &  firstM,
                                                     int       &  lastM,
                                                     int       &  nBlocks);

  static void               findBlockIJK(WellData   * well,
                                         Simbox     * simbox,
                                         const int  * bInd,
                                         const int  & firstM,
                                         const int  & lastM,
                                         const int  & nLayers,
                                         const int  & nBlocks,
                                         int       *& ipos,
                                         int       *& jpos,
                                         int       *& kpos,
                                         float      & dz,
                                         int        & firstB,
                                         int        & lastB);

  static void               blockContinuousLog(const int    *  bInd,
                                               const float  *  wellLog,
                                               const int    & firstM,
                                               const int    & lastM,
                                               const int    & nBlocks,
                                               float       *& blockedLog);

  static void               blockDiscreteLog(const int * bInd,
                                             const int *  wellLog,
                                             const int *  faciesNumbers,
                                             int          nFacies,
                                             const int  & firstM,
                                             const int  & lastM,
                                             const int  & nBlocks,
                                             int       *& blockedLog);

  static int                findMostProbable(const int * count,
                                             int         nFacies,
                                             int         blockIndex);
private:
  void                      setLogFromVerticalTrend(float *& log, double * zpos, int nBlocks,
                                                    float * vertical_trend, double z0, double dz, int nz);

  void                      blockWell(WellData * well,
                                      Simbox   * simbox,
                                      bool       interpolate = false);

  void                      blockCoordinateLog(const int    *  bInd,
                                               const double *  coord,
                                               double       *& blockedCoord);

  void                      interpolateContinuousLog(double * blockedLog, int start, int end, int index, float rel);
  void                      interpolateContinuousLog(float * blockedLog, int start, int end, int index, float rel);
  void                      interpolateTrend(const float * blockedLog, float * trend);
  void                      interpolateTrend(const float * blockedLog, float * trend, const std::vector<Surface *> & limits);
  void                      interpolateTrend(const int * blockedLog, int * trend);

  void                      findBlockXYZ(Simbox * simbox);

  void                      findXYZforVirtualPart(Simbox * simbox);

  float                     computeElasticImpedance(float         vp,
                                                    float         vs,
                                                    float         rho,
                                                    const float * coeff) const;

  void                      smoothTrace(std::vector<float> & trace);

  void                      findPeakTrace(std::vector<float>  & trace,
                                          std::vector<double> & zPeak,
                                          std::vector<double> & peak,
                                          std::vector<double> & b,
                                          double dz,
                                          double ztop);

  void                      peakMatch(std::vector<double> & zPeak,
                                      std::vector<double> & peak,
                                      std::vector<double> & b,
                                      std::vector<double> & zPeakW,
                                      std::vector<double> & peakW,
                                      std::vector<double> & bW);

  double                    computeShift(std::vector<double> & zPeak,
                                         std::vector<double> & zPeakW,
                                         double z0);

  void                      computeGradient(std::vector<double> & Qepsilon,
                                            std::vector<double> & Qepsilon_data,
                                            std::vector<double> & zShift,
                                            int nx,
                                            int ny,
                                            double dx,
                                            double dy);

  void                      smoothGradient(std::vector<double>               & xGradient,
                                           std::vector<double>               & yGradient,
                                           std::vector<double>               & Qepsilon,
                                           std::vector<double>               & Qepsilon_data,
                                           std::vector<std::vector<double> > & Sigma_gradient);
  void                      computePrecisionMatrix(double & a,
                                                   double & b,
                                                   double & c);

  std::string    wellname_;                 ///< Name of well

  double       * xpos_;                     ///<
  double       * ypos_;                     ///< Simbox XYZ value for block
  double       * zpos_;                     ///<

  double       * md_;                       ///< Preserve this when handling NORSAR-wells.

  int          * ipos_;                     ///<
  int          * jpos_;                     ///< Simbox IJK value for block
  int          * kpos_;                     ///<

  float          dz_;                       ///< Simbox dz value for block

  float        * alpha_;                    ///<
  float        * beta_;                     ///< Raw logs (log-domain)
  float        * rho_;                      ///<
  int          * facies_;                   ///< Facies numbers using *internal* numbering
  float       ** facies_prob_;              ///< Facies probabilities calculated in wells

  float        * alpha_highcut_background_; ///<
  float        * beta_highcut_background_;  ///< Logs high-cut filtered to background resolution (log-domain)
  float        * rho_highcut_background_;   ///<

  float        * alpha_highcut_seismic_;    ///<
  float        * beta_highcut_seismic_;     ///< Logs high-cut filtered to approx. seismic resolution (log-domain)
  float        * rho_highcut_seismic_;      ///<

  float        * alpha_seismic_resolution_; ///<
  float        * beta_seismic_resolution_;  ///< Logs filtered to resolution of inversion result
  float        * rho_seismic_resolution_;   ///<

  float        * alpha_predicted_;          ///< Predicted P-wave
  float        * beta_predicted_;           ///< Predicted S-wave
  float        * rho_predicted_;            ///< Predicted density

  float        * alpha_for_facies_;         ///< As above, but omit Vs from filter. Used for facies probabilities.
  float        * rho_for_facies_;           ///<

  float       ** real_seismic_data_;        ///< Seismic data
  float       ** actual_synt_seismic_data_; ///< Forward modelled seismic data using local wavelet
  float       ** well_synt_seismic_data_;   ///< Forward modelled seismic data using wavelet estimated in well
  float       ** cpp_;                      ///< Reflection coefficients
  int            nAngles_;                  ///< Number of AVA stacks

  int            firstM_;                   ///< First well log entry contributing to blocked well
  int            lastM_;                    ///< Last well log entry contributing to blocked well

  int            firstB_;                   ///< First block with contribution from well log
  int            lastB_;                    ///< Last block with contribution from well log

  int            nBlocks_;                  ///< Number of blocks
  int            nLayers_;                  ///< Number of vertical layers (nz)

  const int    * faciesNumbers_;
  int            nFacies_;

  float          lateralThresholdGradient_;  //Minimum lateral distance where gradient lines must not cross
  float          sigma_m_;             //Smoothing factor for the gradients

  bool           interpolate_;         //Should vertical interpolation be done? (for Roxar)
};

#endif
