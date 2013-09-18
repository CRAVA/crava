/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef CRAVA_SRC_IO_H
#define CRAVA_SRC_IO_H

#include <string>
#include "src/definitions.h"

class IO
{
public:
  IO(void);
  ~IO(void);

  // Directory names

  inline static  std::string    TopDirectory(void)                 { return std::string("")                         ;}
  inline static  std::string    InputDirectory(void)               { return std::string("")                         ;}
  inline static  std::string    OutputDirectory(void)              { return std::string("output")                   ;}

  inline static  std::string    PathToNoise(void)                  { return std::string("noise/")                   ;}
  inline static  std::string    PathToWells(void)                  { return std::string("wells/")                   ;}
  inline static  std::string    PathToBackground(void)             { return std::string("background/")              ;}
  inline static  std::string    PathToWavelets(void)               { return std::string("wavelets/")                ;}
  inline static  std::string    PathToSeismicData(void)            { return std::string("seismicdata/")             ;}
  inline static  std::string    PathToTravelTimeData()             { return std::string("traveltimedata/")          ;}
  inline static  std::string    PathToVelocity(void)               { return std::string("velocity/")                ;}
  inline static  std::string    PathToCorrelations(void)           { return std::string("correlations/")            ;}
  inline static  std::string    PathToInversionResults(void)       { return std::string("inversionresults/")        ;}
  inline static  std::string    PathToTmpFiles(void)               { return std::string("")                         ;}
  inline static  std::string    PathToDebug(void)                  { return std::string("")                         ;}

  // File names

  inline static  std::string    FileLog(void)                      { return std::string("logFile")                  ;}
  inline static  std::string    FileDebug(void)                    { return std::string("debug")                    ;}
  inline static  std::string    FileError(void)                    { return std::string("error")                    ;}
  inline static  std::string    FileTasks(void)                    { return std::string("tasks")                    ;}
  inline static  std::string    FileParameterCov(void)             { return std::string("Parameter_Covariance")     ;}
  inline static  std::string    FileLateralCorr(void)              { return std::string("Lateral_Correlation")      ;}
  inline static  std::string    FileTemporalCorr(void)             { return std::string("Temporal_Correlation")     ;}
  inline static  std::string    FileTimeToDepthVelocity(void)      { return std::string("Time-To-Depth_Velocity")   ;}
  inline static  std::string    FileTemporarySeismic(void)         { return std::string("Temp_seis")                ;}

  // Prefixes

  inline static  std::string    PrefixPredictions(void)            { return std::string("Predicted_")               ;}
  inline static  std::string    PrefixSimulations(void)            { return std::string("Simulated_")               ;}
  inline static  std::string    PrefixWells(void)                  { return std::string("Well_")                    ;}
  inline static  std::string    PrefixBlockedWells(void)           { return std::string("Blocked_Well_")            ;}
  inline static  std::string    PrefixBackground(void)             { return std::string("Background_")              ;}
  inline static  std::string    PrefixTrend(void)                  { return std::string("Trend_")                   ;}
  inline static  std::string    PrefixPrior(void)                  { return std::string("Prior_")                   ;}
  inline static  std::string    PrefixPosterior(void)              { return std::string("Posterior_")               ;}
  inline static  std::string    PrefixTemporalCorr(void)           { return std::string("Temporal_Correlation_")    ;}
  inline static  std::string    PrefixCovariance(void)             { return std::string("Covariance_")              ;}
  inline static  std::string    PrefixCrossCovariance(void)        { return std::string("Cross_Covariance_")        ;}
  inline static  std::string    PrefixKrigingData(void)            { return std::string("KrigingData")              ;}
  inline static  std::string    PrefixReflectionCoefficients(void) { return std::string("Reflection_Coefficients_") ;}
  inline static  std::string    PrefixResiduals(void)              { return std::string("Residuals_")               ;}
  inline static  std::string    PrefixSyntheticResiduals(void)     { return std::string("Synthetic_Residuals_")     ;}
  inline static  std::string    PrefixOriginalSeismicData(void)    { return std::string("Original_Seismic_Data_")   ;}
  inline static  std::string    PrefixSyntheticSeismicData(void)   { return std::string("Synthetic_Seismic_Data_")  ;}
  inline static  std::string    PrefixTravelTimeData()             { return std::string("Travel_time_data_")        ;}
  inline static  std::string    PrefixTrendCubes()                 { return std::string("Trend_cube_")              ;}
  inline static  std::string    PrefixFaciesProbability(void)      { return std::string("Facies_Probability_")      ;}
  inline static  std::string    PrefixLikelihood(void)             { return std::string("Seismic_Likelihood_")      ;}
  inline static  std::string    PrefixWavelet(void)                { return std::string("Wavelet_")                 ;}
  inline static  std::string    PrefixWellWavelet(void)            { return std::string("Wavelet_Well_")            ;}
  inline static  std::string    PrefixLocalWaveletGain(void)       { return std::string("Local_Wavelet_Scale_")     ;}
  inline static  std::string    PrefixLocalWaveletShift(void)      { return std::string("Local_Wavelet_Shift_")     ;}
  inline static  std::string    PrefixLocalNoise(void)             { return std::string("Local_Noise_Scaled_")      ;}
  inline static  std::string    PrefixSurface(void)                { return std::string("Surface_")                 ;}
  inline static  std::string    PrefixTop(void)                    { return std::string("Top_")                     ;}
  inline static  std::string    PrefixBase(void)                   { return std::string("Base_")                    ;}
  inline static  std::string    PrefixTime(void)                   { return std::string("Time")                     ;}
  inline static  std::string    PrefixDepth(void)                  { return std::string("Depth")                    ;}
  inline static  std::string    PrefixTmpGrids(void)               { return std::string("tmpGrid_")                 ;}
  inline static  std::string    PrefixDensity(void)                { return std::string("Density_")                 ;}

  // Suffixes

  inline static  std::string    SuffixGeneralData(void)            { return std::string(".dat")                     ;}
  inline static  std::string    SuffixTextFiles(void)              { return std::string(".txt")                     ;}
  inline static  std::string    SuffixCrava(void)                  { return std::string(".crava")                   ;}
  inline static  std::string    SuffixAsciiFiles(void)             { return std::string(".ascii")                   ;}
  inline static  std::string    SuffixAsciiIrapClassic(void)       { return std::string(".irap")                    ;}
  inline static  std::string    SuffixStormBinary(void)            { return std::string(".storm")                   ;}
  inline static  std::string    SuffixRmsWells(void)               { return std::string(".rmswell")                 ;}
  inline static  std::string    SuffixNorsarWells(void)            { return std::string(".nwh")                     ;}
  inline static  std::string    SuffixNorsarLog(void)              { return std::string(".n00")                     ;}
  inline static  std::string    SuffixNorsarTrack(void)            { return std::string(".nwt")                     ;}
  inline static  std::string    SuffixNorsarWavelet(void)          { return std::string(".Swav")                    ;}
  inline static  std::string    SuffixJasonWavelet(void)           { return std::string(".wlt")                     ;}
  inline static  std::string    SuffixSegy(void)                   { return std::string(".segy")                    ;}
  inline static  std::string    SuffixSgriHeader(void)             { return std::string(".Sgrh")                    ;}
  inline static  std::string    SuffixSgri(void)                   { return std::string(".Sgri")                    ;}

  static         int            findGridType(const std::string & fileName);

  //Note: By convention, input path is added to input file names at end of parsing.
  //      Output path and prefix is added to output file name by call to makeFullFileName
  //      just before writing.
  static         void           setFilePrefix(const std::string & filePrefix);
  static         std::string    getFilePrefix(void) { return filePrefix_ ;}
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

  enum           waveletFormats{JASONWAVELET     = 1,
                                NORSARWAVELET    = 2};

  enum           outputGridsElastic{VP                     = 1,
                                    VS                     = 2,
                                    RHO                    = 4,
                                    LAMELAMBDA             = 8,
                                    LAMEMU                 = 16,
                                    POISSONRATIO           = 32,
                                    AI                     = 64,
                                    SI                     = 128,
                                    VPVSRATIO              = 256,
                                    MURHO                  = 512,
                                    LAMBDARHO              = 1024,
                                    BACKGROUND             = 2048,
                                    BACKGROUND_TREND       = 4096};

  enum           outputGridsOther{FACIESPROB             = 1,
                                  FACIESPROB_WITH_UNDEF  = 2,
                                  TIME_TO_DEPTH_VELOCITY = 4,
                                  EXTRA_GRIDS            = 8,
                                  CORRELATION            = 16,
                                  SEISMIC_QUALITY_GRID   = 32,
                                  FACIES_LIKELIHOOD      = 64,
                                  RMS_VELOCITIES         = 128,
                                  TREND_CUBES            = 256};

  enum           outputGridsSeismic{ORIGINAL_SEISMIC_DATA  = 1,
                                    SYNTHETIC_SEISMIC_DATA = 2,
                                    RESIDUAL               = 4,
                                    SYNTHETIC_RESIDUAL     = 8};

  enum           outputWells{WELLS              = 1,
                             BLOCKED_WELLS      = 2,
                             BLOCKED_LOGS       = 4};

  enum           outputOther{EXTRA_SURFACES      =  1,
                             PRIORCORRELATIONS   =  2,
                             BACKGROUND_TREND_1D =  4,
                             LOCAL_NOISE         =  8,
                             ROCK_PHYSICS        = 16,
                             ERROR_FILE          = 32,
                             TASK_FILE           = 64,
                             ROCK_PHYSICS_TRENDS = 128};

  enum           outputWavelets{WELL_WAVELETS    = 1,
                                GLOBAL_WAVELETS  = 2,
                                LOCAL_WAVELETS   = 4};

private:
  static         bool           IsCravaBinaryFile(const std::string & fileName);
  static         bool           IsStormBinaryFile(const std::string & fileName);

  static         std::string    outputPath_;   // Path for all output files.
  static         std::string    filePrefix_;   // Prefix for all output files
};
#endif

