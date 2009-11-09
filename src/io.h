#ifndef CRAVA_SRC_IO_H
#define CRAVA_SRC_IO_H

#include <string>

class IO
{
public:
  IO(void);
  ~IO(void);

  // Directory names

  inline static  std::string    TopDirectory(void)                 { return std::string("")                    ;}
  inline static  std::string    InputDirectory(void)               { return std::string("")                    ;}
  inline static  std::string    OutputDirectory(void)              { return std::string("output")              ;}
                                                                   
  inline static  std::string    PathToWells(void)                  { return std::string("wells/")              ;}
  inline static  std::string    PathToBackground(void)             { return std::string("background/")         ;}
  inline static  std::string    PathToWavelets(void)               { return std::string("wavelets/")           ;}
  inline static  std::string    PathToSeismicData(void)            { return std::string("seismic/")            ;}
  inline static  std::string    PathToVelocity(void)               { return std::string("velocity/")           ;}
  inline static  std::string    PathToCorrelations(void)           { return std::string("correlations/")       ;}
  inline static  std::string    PathToInversionResults(void)       { return std::string("inversionresults/")   ;}
  inline static  std::string    PathToTmpFiles(void)               { return std::string("")                    ;}
  inline static  std::string    PathToDebug(void)                  { return std::string("")                    ;}

  // File names
                                                                                               
  inline static  std::string    FileLog(void)                      { return std::string("logFile")             ;}
  inline static  std::string    FileDebug(void)                    { return std::string("debug")               ;}
  inline static  std::string    FileParameterCov(void)             { return std::string("Var0")                ;}
  inline static  std::string    FileLateralCorr(void)              { return std::string("CorrXY")              ;}
  inline static  std::string    FileTemporalCorr(void)             { return std::string("CorrT")               ;}
  inline static  std::string    FileTemporalCorrUnfiltered(void)   { return std::string("CorrT_Unfiltered")    ;}

  // Prefixes

  inline static  std::string    PrefixDirect(void)                 { return std::string("Direct_")             ;}
  inline static  std::string    PrefixSeismic(void)                { return std::string("Seis_")               ;}
  inline static  std::string    PrefixWells(void)                  { return std::string("Well_")               ;}
  inline static  std::string    PrefixVelocity(void)               { return std::string("Velocity_")           ;}
  inline static  std::string    PrefixBlockedWells(void)           { return std::string("BW_")                 ;}
  inline static  std::string    PrefixBackground(void)             { return std::string("BG_")                 ;}
  inline static  std::string    PrefixTrend(void)                  { return std::string("Trend_")              ;}
  inline static  std::string    PrefixPrior(void)                  { return std::string("Prior_")              ;}
  inline static  std::string    PrefixPosterior(void)              { return std::string("Posterior_")          ;}
  inline static  std::string    PrefixTemporalCorr(void)           { return std::string("CorrT_")              ;}
  inline static  std::string    PrefixCovariance(void)             { return std::string("Cov_")                ;}
  inline static  std::string    PrefixCrossCovariance(void)        { return std::string("CrCov_")              ;}
  inline static  std::string    PrefixKrigingData(void)            { return std::string("KrigingData")         ;}
  inline static  std::string    PrefixReflectionCoefficients(void) { return std::string("Refl_")               ;}
  inline static  std::string    PrefixResiduals(void)              { return std::string("Residuals_")          ;}
  inline static  std::string    PrefixOriginalSeismic(void)        { return std::string("Original_Seismic_")   ;}
  inline static  std::string    PrefixSyntheticSeismic(void)       { return std::string("Synthetic_Seismic_")  ;}
  inline static  std::string    PrefixWavelet(void)                { return std::string("Wavelet_")            ;}
  inline static  std::string    PrefixLocalWaveletGain(void)       { return std::string("Local_Wavelet_Gain_") ;}
  inline static  std::string    PrefixLocalWaveletShift(void)      { return std::string("Local_Wavelet_Shift_");}
  inline static  std::string    PrefixLocalNoise(void)             { return std::string("Local_Noise_")        ;}
  inline static  std::string    PrefixSurface(void)                { return std::string("Surface_")            ;}
  inline static  std::string    PrefixTop(void)                    { return std::string("Top_")                ;}
  inline static  std::string    PrefixBase(void)                   { return std::string("Base_")               ;}
  inline static  std::string    PrefixTime(void)                   { return std::string("Time")                ;}
  inline static  std::string    PrefixDepth(void)                  { return std::string("Depth")               ;}
  inline static  std::string    PrefixTmpGrids(void)               { return std::string("tmpGrid_")            ;}

  // Suffixes

  inline static  std::string    SuffixDirectData(void)             { return std::string(".bin")                ;}
  inline static  std::string    SuffixGeneralData(void)            { return std::string(".dat")                ;}
  inline static  std::string    SuffixTextFiles(void)              { return std::string(".txt")                ;}
  inline static  std::string    SuffixAsciiFiles(void)             { return std::string(".ascii")              ;}
  inline static  std::string    SuffixAsciiIrapClassic(void)       { return std::string(".irap")               ;}
  inline static  std::string    SuffixStormBinary(void)            { return std::string(".storm")              ;}
  inline static  std::string    SuffixRmsWells(void)               { return std::string(".rmswell")            ;}
  inline static  std::string    SuffixNorsarWells(void)            { return std::string(".nwh")                ;}
  inline static  std::string    SuffixNorsarLog(void)              { return std::string(".n00")                ;}
  inline static  std::string    SuffixNorsarTrack(void)            { return std::string(".nwt")                ;}
  inline static  std::string    SuffixNorsarWavelet(void)          { return std::string(".Swav")               ;}
  inline static  std::string    SuffixJasonWavelet(void)           { return std::string(".wlt")                ;}
  inline static  std::string    SuffixSegy(void)                   { return std::string(".segy")               ;}
  inline static  std::string    SuffixSgriHeader(void)             { return std::string(".Sgrh")               ;}
  inline static  std::string    SuffixSgri(void)                   { return std::string(".Sgri")               ;}

  //Note: By convention, input path is added to input file names at end of parsing.
  //      Output path and prefix is added to output file name by call to makeFullFileName
  //      just before writing.
  static         void           setFilePrefix(const std::string & filePrefix);
  static         void           setOutputPath(const std::string & outputPath);
  static         std::string    getOutputPath(void) { return outputPath_ ;}

  static         std::string    makeFullFileName(const std::string & subDir, 
                                                 const std::string & fileName);

  static         std::string    outputPath_;   // Path for all output files.
  static         std::string    filePrefix_;   // Prefix for all output files
};
#endif

