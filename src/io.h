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
  inline static  std::string    PathToSeismicData(void)            { return std::string("seismicdata/")        ;}
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
  inline static  std::string    FileTimeToDepthVelocity(void)      { return std::string("Time-To-Depth_Velocity");}

  // Prefixes

  inline static  std::string    PrefixPredicted(void)              { return std::string("Predicted_")          ;}
  inline static  std::string    PrefixSimulated(void)              { return std::string("Simulated_")          ;}
  inline static  std::string    PrefixWells(void)                  { return std::string("Well_")               ;}
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
  inline static  std::string    PrefixLocalWaveletGain(void)       { return std::string("Local_Wavelet_Scale_");}
  inline static  std::string    PrefixLocalWaveletShift(void)      { return std::string("Local_Wavelet_Shift_");}
  inline static  std::string    PrefixLocalNoise(void)             { return std::string("Local_Noise_Scaled_")      ;}
  inline static  std::string    PrefixSurface(void)                { return std::string("Surface_")            ;}
  inline static  std::string    PrefixTop(void)                    { return std::string("Top_")                ;}
  inline static  std::string    PrefixBase(void)                   { return std::string("Base_")               ;}
  inline static  std::string    PrefixTime(void)                   { return std::string("Time")                ;}
  inline static  std::string    PrefixDepth(void)                  { return std::string("Depth")               ;}
  inline static  std::string    PrefixTmpGrids(void)               { return std::string("tmpGrid_")            ;}

  // Suffixes

  inline static  std::string    SuffixGeneralData(void)            { return std::string(".dat")                ;}
  inline static  std::string    SuffixTextFiles(void)              { return std::string(".txt")                ;}
  inline static  std::string    SuffixCrava(void)                  { return std::string(".crava")              ;}
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

  static         int            findGridType(const std::string & fileName);

  //Note: By convention, input path is added to input file names at end of parsing.
  //      Output path and prefix is added to output file name by call to makeFullFileName
  //      just before writing.
  static         void           setFilePrefix(const std::string & filePrefix);
  static         void           setOutputPath(const std::string & outputPath);
  static         std::string    getOutputPath(void) { return outputPath_ ;}

  static         std::string    makeFullFileName(const std::string & subDir, 
                                                 const std::string & fileName);

  static         void           writeSurfaceToFile(const Surface     & surface,
                                                   const std::string & baseName,
                                                   const std::string & path,
                                                   int                 format);

  enum           domains{TIMEDOMAIN  = 1, 
                         DEPTHDOMAIN = 2};

  enum           gridFormats{UNKNOWN =  0,
                             SEGY    =  1, 
                             STORM   =  2, 
                             ASCII   =  4, 
                             SGRI    =  8,
                             CRAVA   = 16};

  enum           wellFormats{RMSWELL    = 1, 
                             NORSARWELL = 2};

  enum           outputGrids{CORRELATION            = 1, 
                             RESIDUAL               = 2, 
                             VP                     = 4, 
                             VS                     = 8, 
                             RHO                    = 16,
                             LAMELAMBDA             = 32, 
                             LAMEMU                 = 64, 
                             POISSONRATIO           = 128, 
                             AI                     = 256,
                             SI                     = 512, 
                             VPVSRATIO              = 1024, 
                             MURHO                  = 2048, 
                             LAMBDARHO              = 4096, 
                             BACKGROUND             = 8192, 
                             BACKGROUND_TREND       = 16384, 
                             FACIESPROB             = 32768, 
                             FACIESPROBRELATIVE     = 65536, 
                             EXTRA_GRIDS            = 131072,
                             SEISMIC_DATA           = 262144,
                             TIME_TO_DEPTH_VELOCITY = 524288};

  enum           outputWells{WELLS              = 1,
                             BLOCKED_WELLS      = 2,
                             BLOCKED_LOGS       = 4};

  enum           outputOther{WAVELETS            = 1,
                             EXTRA_SURFACES      = 2,
                             PRIORCORRELATIONS   = 4,
                             BACKGROUND_TREND_1D = 8};

private:
  static         bool           IsCravaBinaryFile(const std::string & fileName);
  static         bool           IsStormBinaryFile(const std::string & fileName);

  static         std::string    outputPath_;   // Path for all output files.
  static         std::string    filePrefix_;   // Prefix for all output files
};
#endif

