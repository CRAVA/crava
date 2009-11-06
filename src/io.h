#ifndef CRAVA_SRC_IO_H
#define CRAVA_SRC_IO_H

#include <string>

namespace IO
{
  // Directory names

  inline static  std::string    TopDirectory(void)                    { return std::string("")                           ;}
  inline static  std::string    InputDirectory(void)                  { return std::string("")                           ;}
  inline static  std::string    OutputDirectory(void)                 { return std::string("output")                     ;}
                                                                    
  inline static  std::string    PathToWells(void)                     { return std::string("wells/")                     ;}
  inline static  std::string    PathToBackground(void)                { return std::string("background/")                ;}
  inline static  std::string    PathToWavelets(void)                  { return std::string("wavelets/")                  ;}
  inline static  std::string    PathToSeismicData(void)               { return std::string("seismic/")                   ;}
  inline static  std::string    PathToSurfaces(void)                  { return std::string("surfaces/")                  ;}
  inline static  std::string    PathToVelocity(void)                  { return std::string("velocity/")                  ;}
  inline static  std::string    PathToCorrelations(void)              { return std::string("correlations/")              ;}
  inline static  std::string    PathToInversionResult(void)           { return std::string("inversionresult/")           ;}

  // File names
                                                                                                  
  inline static  std::string    FileLog(void)                         { return std::string("logFile.txt")                ;}
  inline static  std::string    FileDebug(void)                       { return std::string("debug.txt")                  ;}
  inline static  std::string    FilePriorParameterCov(void)           { return std::string("Prior_Var0")                 ;}
  inline static  std::string    FilePriorLateralCorr(void)            { return std::string("Prior_CorrXY")               ;}
  inline static  std::string    FilePriorTemporalCorr(void)           { return std::string("Prior_CorrT")                ;}
  inline static  std::string    FilePriorTemporalCorrUnfiltered(void) { return std::string("Prior_CorrT_Unfiltered")     ;}
  inline static  std::string    FilePosteriorParameterCov(void)       { return std::string("Posterior_Var0")             ;}

  // Prefixes

  inline static  std::string    PrefixDirect(void)                    { return std::string("Direct_")                    ;}
  inline static  std::string    PrefixSeismic(void)                   { return std::string("Seis_")                      ;}
  inline static  std::string    PrefixWells(void)                     { return std::string("Well_")                      ;}
  inline static  std::string    PrefixVelocity(void)                  { return std::string("Velocity_")                  ;}
  inline static  std::string    PrefixBlockedWells(void)              { return std::string("BW_")                        ;}
  inline static  std::string    PrefixBackground(void)                { return std::string("BG_")                        ;}
  inline static  std::string    PrefixBackgroundTrends(void)          { return std::string("BG_verticalTrend_")          ;}
  inline static  std::string    PrefixPosteriorTemporalCorr(void)     { return std::string("Posterior_CorrT_")           ;}
  inline static  std::string    PrefixPosteriorCov(void)              { return std::string("Posterior_Cov_")             ;}
  inline static  std::string    PrefixPosteriorCrCov(void)            { return std::string("Posterior_CrCov_")           ;}
  inline static  std::string    PrefixKrigingData(void)               { return std::string("KrigingData")                ;}

  // Suffixes

  inline static  std::string    SuffixDirectData(void)                { return std::string(".bin")                       ;}
  inline static  std::string    SuffixGeneralData(void)               { return std::string(".dat")                       ;}
  inline static  std::string    SuffixAsciiIrapClassic(void)          { return std::string(".irap")                      ;}
  inline static  std::string    SuffixStormBinary(void)               { return std::string(".storm")                     ;}
  inline static  std::string    SuffixRmsWells(void)                  { return std::string(".rmswell")                   ;}
  inline static  std::string    SuffixNorsarWells(void)               { return std::string(".nwh")                       ;}
  inline static  std::string    SuffixNorsarLog(void)                 { return std::string(".n00")                       ;}
  inline static  std::string    SuffixNorsarTrack(void)               { return std::string(".nwt")                       ;}
  inline static  std::string    SuffixNorsarWavelet(void)             { return std::string(".Swav")                      ;}
  inline static  std::string    SuffixJasonWavelet(void)              { return std::string(".wlt")                       ;}
}
#endif
