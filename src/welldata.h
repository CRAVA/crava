#ifndef WELLDATA_H
#define WELLDATA_H

#include "lib/global_def.h"
#include "src/blockedlogs.h"
#include "src/simbox.h"
#include <string.h>
#include <stdio.h>

class ModelSettings;

class WellData 
{
public:
  WellData(const char     * wellFileName, 
           ModelSettings  * modelSettings,
           char          ** headerList,
           bool             faciesLogGiven);
  ~WellData(void);

  const float   * getAlpha(int &nData) const;
  const float   * getBeta(int &nData) const;
  const float   * getRho(int &nData) const;
  const int     * getFacies(int &nData) const;
  const float   * getAlphaBackgroundResolution(int &nData) const; 
  const float   * getBetaBackgroundResolution(int &nData) const;
  const float   * getRhoBackgroundResolution(int &nData) const;
  const float   * getAlphaSeismicResolution(int &nData) const;
  const float   * getBetaSeismicResolution(int &nData) const; 
  const float   * getRhoSeismicResolution(int &nData) const;  
  const double  * getXpos(int &nData)      const;
  const double  * getYpos(int &nData)      const;
  const double  * getZpos(int &nData)      const;
  const char    * getWellname(void)        const { return wellname_       ;} 
  const bool      hasSyntheticVsLog(void)  const { return syntheticVsLog_ ;}
  const bool      isDeviated(void)         const { return isDeviated_     ;}
  bool            isFaciesLogDefined(void) const;
  int             getNFacies(void)         const { return nFacies_        ;}
  int           * getFaciesNr(void)        const { return faciesNr_       ;}
  int             getFaciesNr(int i)       const { return faciesNr_[i]    ;}
  char          * getFaciesName(int i)     const { return faciesNames_[i] ;}
  void            getMinMaxFnr(int &min, int &max) const;

  BlockedLogs   * getBlockedLogsPropThick(void)  const { return blockedLogsPropThick_ ;}
  BlockedLogs   * getBlockedLogsConstThick(void) const { return blockedLogsConstThick_ ;} 

  void            setBlockedLogsPropThick(BlockedLogs * bl)  { blockedLogsPropThick_  = bl ;}
  void            setBlockedLogsConstThick(BlockedLogs * bl) { blockedLogsConstThick_ = bl ;}

  int             getNd() const;
  int             checkError(char * errText);
  int             checkSimbox(Simbox *simbox);  
  void            removeDuplicateLogEntries(int & nMerges);
  void            setWrongLogEntriesUndefined(int & count_alpha, int & count_beta, int & count_rho);
  void            filterLogs(void);
  void            lookForSyntheticVsLog(float & rank_correlation);
  void            calculateDeviation(float  & devAngle,
                                     Simbox * timeSimbox);
  void            countFacies(Simbox *simbox, int * faciesCount);
  void            writeRMSWell(void);

  static void     applyFilter(float *log_filtered, float *log_interpolated, int nt, double dt_milliseconds, float maxHz);
  static void     interpolateLog(float *log_interpolated, const float *log_raw, int nd);
  int             isFaciesOk(){return faciesok_;};
private:
  void            readRMSWell(const char * wellFileName, char **headerList, bool faciesLogGiven);
  void            mergeCells(const char * name, double * log_resampled, double * log, 
                             int ii, int istart, int iend, bool debug);
  void            mergeCells(const char * name, float * log_resampled, float * log, 
                            int ii, int istart, int iend, bool debug);
  void            mergeCells(const char * name, int * log_resampled, int * log, 
                             int ii, int istart, int iend, bool debug);
  void            mergeCellsDiscrete(const char * name, int * log_resampled, int * log, int ii, 
                                     int istart, int iend, bool printToScreen);
  void            resampleTime(double * time_resampled, int nd, double & dt);
  void            resampleLog(float * log_resampled, const float * log_interpolated, 
                              const double * time, const double * time_resampled, 
                              int nd, double dt);

  ModelSettings * modelSettings_;

  char          * wellname_;                    
  char          * wellfilename_;                // wellname given in RMS well file

  double          xpos0_;                       // x-coordinate from well file header
  double          ypos0_;                       // y-coordinate from well file header
  double        * xpos_;                        // x-coord. in well
  double        * ypos_;                        // y-coord in well
  double        * zpos_;                        // time step

  float         * alpha_; 
  float         * beta_;
  float         * rho_;
  int 	        * facies_;

  float         * alpha_background_resolution_; // Vp  - filtered to background resolution 
  float         * beta_background_resolution_;  // 
  float         * rho_background_resolution_;   // 

  float         * alpha_seismic_resolution_;    // Vp  - filtered to seismic resolution
  float         * beta_seismic_resolution_;     // 
  float         * rho_seismic_resolution_;      // 

  char            errTxt_[MAX_STRING];
  int             error_;
  int             timemissing_;
  bool            syntheticVsLog_;
  bool            isDeviated_;
  int             nFacies_;
  char          * faciesLogName_;
  char         ** faciesNames_;
  int           * faciesNr_;
  int             nd_;                          // number of obs. in well
  int             faciesok_;                    // all faciesnumbers read are present in header

  BlockedLogs   * blockedLogsPropThick_;
  BlockedLogs   * blockedLogsConstThick_;
};

inline const float* WellData::getAlpha(int &nData) const
{
  nData = nd_;
  return alpha_;
}

inline const float* WellData::getBeta(int &nData) const
{
  nData = nd_;
  return beta_;
}

inline const float* WellData::getRho(int &nData) const
{
  nData = nd_;
  return rho_;
}
inline const int* WellData::getFacies(int &nData) const
{
  nData = nd_;
  return facies_;
}
inline const float* WellData::getAlphaBackgroundResolution(int &nData) const
{
  nData = nd_;
  return alpha_background_resolution_;
}

inline const float* WellData::getBetaBackgroundResolution(int &nData) const
{
  nData = nd_;
  return beta_background_resolution_;
}

inline const float* WellData::getRhoBackgroundResolution(int &nData) const
{
  nData = nd_;
  return rho_background_resolution_;
}

inline const float* WellData::getAlphaSeismicResolution(int &nData) const
{
  nData = nd_;
  return alpha_seismic_resolution_;
}

inline const float* WellData::getBetaSeismicResolution(int &nData) const
{
  nData = nd_;
  return beta_seismic_resolution_;
}

inline const float* WellData::getRhoSeismicResolution(int &nData) const
{
  nData = nd_;
  return rho_seismic_resolution_;
}


inline const double * WellData::getXpos(int &nData) const
{
  nData = nd_;
  return xpos_;
}
inline const double * WellData::getYpos(int &nData) const
{
  nData = nd_;
  return ypos_;
}
inline const double * WellData::getZpos(int &nData) const
{
  nData = nd_;
  return zpos_;
}
inline int WellData::getNd() const
{
  return nd_;
}
inline bool WellData::isFaciesLogDefined() const
{
  if(nFacies_==0)
    return 0;
  else
    return 1;
}
#endif
