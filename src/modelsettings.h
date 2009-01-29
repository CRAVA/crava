
#ifndef MODELSETTINGS_H
#define MODELSETTINGS_H

#include <stdio.h>
#include <string.h>

#include "lib/global_def.h"
#include "src/definitions.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/traceheader.hpp"
#include "nrlib/segy/segy.hpp"

class Vario;

class ModelSettings
{
public:
  ModelSettings(void);
  ~ModelSettings(void); 

  Vario             * getAngularCorr(void)          const { return angularCorr_         ;} 
  Vario             * getLateralCorr(void)          const { return lateralCorr_         ;}
  Vario             * getBackgroundVario(void)      const { return backgroundVario_     ;} 
  Vario             * getLocalWaveletVario(void)    const { return localWaveletVario_   ;} 
  SegyGeometry      * getAreaParameters(void)       const { return geometry_            ;}    
  TraceHeaderFormat * getTraceHeaderFormat(void)    const { return traceHeaderFormat_   ;}
  float             * getKrigingParameters(void)    const { return krigingParams_       ;}
  float               getConstBackValue(int i)      const { return constBackValue_[i]   ;}
  int                 getNumberOfAngles(void)       const { return angle_.size()        ;} 
  int                 getSeismicType(int i)         const { return seismicType_[i]      ;}
  float               getAngle(int i)               const { return angle_[i]            ;}
  float               getWaveletScale(int i)        const { return waveletScale_[i]     ;} 
  float               getSNRatio(int i)             const { return SNRatio_[i]          ;} 
  bool                getMatchEnergies(int i)       const { return matchEnergies_[i]    ;} 
  bool                getEstimateWavelet(int i)     const { return estimateWavelet_[i]  ;}
  bool                getEstimateSNRatio(int i)     const { return estimateSNRatio_[i]  ;}
  char              * getFaciesName(int i)          const { return faciesNames_[i]      ;}
  int                 getFaciesLabel(int i)         const { return faciesLabels_[i]     ;}
  int                 getIndicatorBGTrend(int i)    const { return indBGTrend_[i]       ;}
  int                 getIndicatorWavelet(int i)    const { return indWavelet_[i]       ;}
  int                 getIndicatorFacies(int i)     const { return indFacies_[i]        ;}
  int                 getNumberOfFacies(void)       const { return nFacies_             ;}
  int                 getNumberOfWells(void)        const { return nWells_              ;}
  int                 getNumberOfSimulations(void)  const { return nSimulations_        ;}
  float               getAlphaMin(void)             const { return alpha_min_           ;}
  float               getAlphaMax(void)             const { return alpha_max_           ;}
  float               getBetaMin(void)              const { return beta_min_            ;}
  float               getBetaMax(void)              const { return beta_max_            ;}
  float               getRhoMin(void)               const { return rho_min_             ;}
  float               getRhoMax(void)               const { return rho_max_             ;}
  float               getVarAlphaMin(void)          const { return var_alpha_min_       ;}
  float               getVarAlphaMax(void)          const { return var_alpha_max_       ;}
  float               getVarBetaMin(void)           const { return var_beta_min_        ;}
  float               getVarBetaMax(void)           const { return var_beta_max_        ;}
  float               getVarRhoMin(void)            const { return var_rho_min_         ;}
  float               getVarRhoMax(void)            const { return var_rho_max_         ;}
  float               getMaxHzBackground(void)      const { return maxHz_background_    ;}
  float               getMaxHzSeismic(void)         const { return maxHz_seismic_       ;}
  float               getMaxRankCorr(void)          const { return maxRankCorr_         ;}
  float               getMaxMergeDist(void)         const { return maxMergeDist_        ;}
  float               getMaxDevAngle(void)          const { return maxDevAngle_         ;}
  float               getLowCut(void)               const { return lowCut_              ;}
  float               getHighCut(void)              const { return highCut_             ;}
  float               getWNC(void)                  const { return wnc_                 ;}
  float               getEnergyThreshold(void)      const { return energyThreshold_     ;}
  float               getMinRelWaveletAmp(void)     const { return minRelWaveletAmp_    ;}
  float               getMaxWaveletShift(void)      const { return maxWaveletShift_     ;}
  float               getWaveletTaperingL(void)     const { return waveletTaperingL_    ;}
  float               getXpad(void)                 const { return xPad_                ;}
  float               getYpad(void)                 const { return yPad_                ;}
  float               getZpad(void)                 const { return zPad_                ;}
  int                 getNXpad(void)                const { return nxPad_               ;}
  int                 getNYpad(void)                const { return nyPad_               ;}
  int                 getNZpad(void)                const { return nzPad_               ;}
  float               getSegyOffset(void)           const { return segyOffset_          ;}
  float               getPundef(void)               const { return p_undef_             ;}
  double              getLzLimit(void)              const { return lzLimit_             ;}
  double              getTimeDTop(void)             const { return time_dTop_           ;}
  double              getTimeLz(void)               const { return time_lz_             ;}
  double              getTimeDz(void)               const { return time_dz_             ;}
  int                 getTimeNz(void)               const { return time_nz_             ;}
  int                 getOutputFlag(void)           const { return outputFlag_          ;}
  int                 getFormatFlag(void)           const { return formatFlag_          ;}
  int                 getDebugFlag(void)            const { return debugFlag_           ;}
  static int          getDebugLevel(void)                 { return debugFlag_           ;}
  int                 getFileGrid(void)             const { return fileGrid_            ;}
  bool                getGenerateSeismic(void)      const { return generateSeismic_     ;}
  bool                getGenerateBackground(void)   const { return generateBackground_  ;}
  bool                getFaciesLogGiven(void)       const { return faciesLogGiven_      ;}
  bool                getDoDepthConversion(void)    const { return doDepthConversion_   ;}
  bool                getParallelTimeSurfaces(void) const { return parallelTimeSurfaces_;}
  bool                getUseLocalWavelet(void)      const { return useLocalWavelet_     ;}
  int                 getLogLevel(void)             const { return logLevel_            ;}
  int                 getSeed(void)                 const { return seed_                ;}
  bool                getDoInversion(void);                                            
  //NBNB Ragnar: Not active yet.                                                       
  Surface           * getCorrelationSurface()       const { return NULL                 ;}
  void                rotateVariograms(float angle);

  void                setAngularCorr(Vario * vario);    
  void                setLateralCorr(Vario * vario);    
  void                setBackgroundVario(Vario * vario);
  void                setLocalWaveletVario(Vario * vario);
  void                copyBackgroundVarioToLocalWaveletVario(void);
  void                setAreaParameters(const SegyGeometry * geometry);
  void                setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat);
  void                setKrigingParameters(float * krigingParams, int nParams);
  void                setConstBackValue(int i, float constBackValue){ constBackValue_[i]    = constBackValue      ;}
  void                addSeismicType(int seismicType)               { seismicType_.push_back(seismicType)         ;}
  void                addAngle(float angle)                         { angle_.push_back(angle)                     ;}
  void                addWaveletScale(float waveletScale)           { waveletScale_.push_back(waveletScale)       ;}
  void                setSNRatio(int i, float SNRatio)              { SNRatio_[i]           = SNRatio             ;}
  void                addSNRatio(float SNRatio)                     { SNRatio_.push_back(SNRatio)                 ;}
  void                addMatchEnergies(bool matchEnergies)          { matchEnergies_.push_back(matchEnergies)     ;}
  void                addEstimateWavelet(bool estimateWavelet)      { estimateWavelet_.push_back(estimateWavelet) ;}
  void                addEstimateSNRatio(bool estimateSNRatio)      { estimateSNRatio_.push_back(estimateSNRatio) ;}
  void                setAllIndicatorsTrue(int nWells);
  void                setIndicatorBGTrend(int * indBGTrend, int nWells);
  void                setIndicatorWavelet(int * indWavelet, int nWells);
  void                setIndicatorFacies(int * indFacies, int nWells);
  void                setFaciesLabels(int * faciesLabels, int nFacies);
  void                setFaciesNames(char ** faciesNames, int nFacies);
  void                setNumberOfFacies(int nFacies)                { nFacies_              = nFacies            ;}
  void                setNumberOfWells(int nWells)                  { nWells_               = nWells             ;} 
  void                setNumberOfSimulations(int nSimulations)      { nSimulations_         = nSimulations       ;} 
  void                setAlphaMin(float alpha_min)                  { alpha_min_            = alpha_min          ;}
  void                setAlphaMax(float alpha_max)                  { alpha_max_            = alpha_max          ;}
  void                setBetaMin(float beta_min)                    { beta_min_             = beta_min           ;}
  void                setBetaMax(float beta_max)                    { beta_max_             = beta_max           ;}
  void                setRhoMin(float rho_min)                      { rho_min_              = rho_min            ;}
  void                setRhoMax(float rho_max)                      { rho_max_              = rho_max            ;}
  void                setVarAlphaMin(float var_alpha_min)           { var_alpha_min_        = var_alpha_min      ;}
  void                setVarAlphaMax(float var_alpha_max)           { var_alpha_max_        = var_alpha_max      ;}
  void                setVarBetaMin(float var_beta_min)             { var_beta_min_         = var_beta_min       ;}
  void                setVarBetaMax(float var_beta_max)             { var_beta_max_         = var_beta_max       ;}
  void                setVarRhoMin(float var_rho_min)               { var_rho_min_          = var_rho_min        ;}
  void                setVarRhoMax(float var_rho_max)               { var_rho_max_          = var_rho_max        ;}
  void                setMaxHzBackground(float maxHz_background)    { maxHz_background_     = maxHz_background   ;}
  void                setMaxHzSeismic(float maxHz_seismic)          { maxHz_seismic_        = maxHz_seismic      ;}
  void                setMaxRankCorr(float maxRankCorr)             { maxRankCorr_          = maxRankCorr        ;}
  void                setMaxMergeDist(float maxMergeDist)           { maxMergeDist_         = maxMergeDist       ;}
  void                setMaxDevAngle(float maxDevAngle)             { maxDevAngle_          = maxDevAngle        ;}
  void                setLowCut(float lowCut)                       { lowCut_               = lowCut             ;}
  void                setHighCut(float highCut)                     { highCut_              = highCut            ;}
  void                setWNC(float wnc)                             { wnc_                  = wnc                ;}
  void                setEnergyThreshold(float energyThreshold)     { energyThreshold_      = energyThreshold    ;}
  void                setinRelWaveletAmp(float minRelWaveletAmp)    { minRelWaveletAmp_     = minRelWaveletAmp   ;}
  void                setMaxWaveletShift(float maxWaveletShift)     { maxWaveletShift_      = maxWaveletShift    ;}
  void                setWaveletTaperingL(float waveletTaperingL)   { waveletTaperingL_     = waveletTaperingL   ;}
  void                setXpad(float xPad)                           { xPad_                 = xPad               ;}
  void                setYpad(float yPad)                           { yPad_                 = yPad               ;}
  void                setZpad(float zPad)                           { zPad_                 = zPad               ;}
  void                setNXpad(int nxPad)                           { nxPad_                = nxPad              ;}
  void                setNYpad(int nyPad)                           { nyPad_                = nyPad              ;}
  void                setNZpad(int nzPad)                           { nzPad_                = nzPad              ;}
  void                setSegyOffset(float segyOffset)               { segyOffset_           = segyOffset         ;}
  void                setPundef(float p_undef)                      { p_undef_              = p_undef            ;}
  void                setLzLimit(double lzLimit)                    { lzLimit_              = lzLimit            ;}
  void                setTimeDTop(double time_dTop)                 { time_dTop_            = time_dTop          ;}
  void                setTimeLz(double time_lz)                     { time_lz_              = time_lz            ;}
  void                setTimeDz(double time_dz)                     { time_dz_              = time_dz            ;}
  void                setTimeNz(int time_nz)                        { time_nz_              = time_nz            ;}
  void                setOutputFlag(int outputFlag);                                        
  void                setFormatFlag(int formatFlag)                 { formatFlag_           = formatFlag         ;}
  void                setDebugFlag(int debugFlag)                   { debugFlag_            = debugFlag          ;}
  void                setFileGrid(int fileGrid)                     { fileGrid_             = fileGrid           ;}
  void                setGenerateSeismic(bool generateSeismic)      { generateSeismic_      = generateSeismic    ;}
  void                setGenerateBackground(bool generateBackgr)    { generateBackground_   = generateBackgr     ;}
  void                setFaciesLogGiven(bool faciesLogGiven)        { faciesLogGiven_       = faciesLogGiven     ;}
  void                setDoDepthConversion(bool doDepthConversion)  { doDepthConversion_    = doDepthConversion  ;}
  void                setParallelTimeSurfaces(bool pTimeSurfaces)   { parallelTimeSurfaces_ = pTimeSurfaces      ;}
  void                setUseLocalWavelet(bool useLocalWavelet)      { useLocalWavelet_      = useLocalWavelet    ;}
  void                setLogLevel(int logLevel)                     { logLevel_             = logLevel           ;}
  void                setSeed(int seed)                             { seed_                 = seed               ;}

  enum                outputGrids{PREDICTION         = 1, 
                                  CORRELATION        = 2, 
                                  RESIDUAL           = 4, 
                                  VP                 = 8, 
                                  VS                 = 16,
                                  RHO                = 32, 
                                  LAMELAMBDA         = 64, 
                                  LAMEMU             = 128, 
                                  POISSONRATIO       = 256,
                                  AI                 = 512, 
                                  SI                 = 1024, 
                                  VPVSRATIO          = 2048, 
                                  MURHO              = 4096, 
                                  LAMBDARHO          = 8192, 
                                  PRIORCORRELATIONS  = 16384, 
                                  BACKGROUND         = 32768, 
                                  WELLS              = 65536, 
                                  WAVELETS           = 131072, 
                                  NOTIME             = 262144,
                                  FACIESPROB         = 524288,
                                  FACIESPROBRELATIVE = 1048576,
                                  BLOCKED_WELLS      = 2097152,
                                  BLOCKED_LOGS       = 4194304,
                                  EXTRA_SURFACES     = 8388608,
                                  EXTRA_GRIDS        = 16777216,
                                  BACKGROUND_TREND   = 33554432};
                   
  enum                sseismicTypes{STANDARDSEIS = 0, PSSEIS = 1};
                   
  static void         setFilePrefix(char * filePrefix);
  static std::string  makeFullFileName(const std::string name, const std::string postfix = "");
                   
private:           
                   
  Vario             * angularCorr_;           // Variogram for lateral error correlation
  Vario             * lateralCorr_;           // Variogram for lateral parameter correlation 
  Vario             * backgroundVario_;       // Used for lateral background correlation.
  Vario             * localWaveletVario_;     // Used for local wavelet (gain and shift) and local noise.
                 
  SegyGeometry      * geometry_;              // area parameters
  TraceHeaderFormat * traceHeaderFormat_;     // 
  float             * krigingParams_;   
                   
  std::vector<int>    seismicType_;           // PP- or PS- seismic
  std::vector<float>  angle_;                 // Angles
  std::vector<float>  waveletScale_;          // Signal-to-noise ratio
  std::vector<float>  SNRatio_;               // Signal-to-noise ratio
  std::vector<bool>   matchEnergies_;         // Let dataVariance_ = signalVariance_
  std::vector<bool>   estimateWavelet_;       // 
  std::vector<bool>   estimateSNRatio_;       //

  std::vector<float>  constBackValue_;        // Values set for constant background model
                                              // Negative value ==> read from file (actual value gives format).
  int               * indBGTrend_;            // Use well to estimate background trend? (1=yes,0=no)
  int               * indWavelet_;            // Use well to estimate wavelet? (1=yes,0=no)
  int               * indFacies_;             // Use well to estimate facies? (1=yes,0=no)
                   
  char             ** faciesNames_;           // Facies names
  int               * faciesLabels_;          // Facies labels
  int                 nFacies_;
                   
  int                 nWells_;
  int                 nSimulations_;
                   
  float               alpha_min_;             // Vp - smallest allowed value
  float               alpha_max_;             // Vp - largest allowed value
  float               beta_min_;              // Vs - smallest allowed value
  float               beta_max_;              // Vs - largest allowed value
  float               rho_min_;               // Rho - smallest allowed value
  float               rho_max_;               // Rho - largest allowed value
                      
  float               var_alpha_min_;         //| These min and max values are used for consistency check. If  
  float               var_alpha_max_;         //| variances are outside these ranges there is probably a
  float               var_beta_min_;          //| problem with the logs.
  float               var_beta_max_;          //| 
  float               var_rho_min_;           //| The limits are for point variances. The minimum allowed variance 
  float               var_rho_max_;           //| for parameters will be scaled with 1/dt*dt
                      
  float               maxHz_background_;      // Background resolution (high cut frequency)
  float               maxHz_seismic_;         // Seismic resolution (high cut frequency)
                      
  float               maxRankCorr_;           // Vp-Vs correlation threshold for regarding Vs log synthetic
  float               maxMergeDist_;          // log entries closer than this will be merged
  float               maxDevAngle_;           // Wells with a local deviation larger than this is treated as deviated
                      
  float               lowCut_;                // lower limit for frequency to be inverted
  float               highCut_;               // upper limit for frecuency to be inverted
                      
  float               wnc_;                   // White noise component, see crava.h  
                      
  float               energyThreshold_;       // If energy in reflection trace divided by mean energy
                                              // in reflection traces is lower than this, the reflections
                                              // will be interpolated. Default 0.
                      
  float               minRelWaveletAmp_;      // Minimum relative wavelet amplitude. Smaller amplitudes are disregarded.
  float               maxWaveletShift_;       // Largest allowed shift when estimating wavelet
  float               waveletTaperingL_;      // Til Odds waveletestimering
  
  float               xPad_;                  // Padding factor in x direction
  float               yPad_;
  float               zPad_; 
                      
  int                 nxPad_;                 // Number of cells to pad in x direction
  int                 nyPad_;
  int                 nzPad_; 
                      
  float               segyOffset_;            // Starttime for SegY cubes.
                      
  float               p_undef_;               // NBNB-PAL: Hva gjør denne?
                      
  double              lzLimit_;               // Minimum allowed value for (min interval thickness)/(max interval thickness)
  double              time_dTop_;             // Used when top and base surfaces are parallel
  double              time_lz_;               // Used when top and base surfaces are parallel
  double              time_dz_;               // Used when top and base surfaces are parallel
  int                 time_nz_;               // Used when top and base surfaces are parallel

  int                 outputFlag_;            // Decides which grids to write (except simulation)
  int                 formatFlag_;            // Decides output format, see fftgird.h
  int                 fileGrid_;              // Indicator telling if grids are to be kept on file
                      
  bool                generateSeismic_;       // Forward modelling
  bool                generateBackground_;    // Make background model
  bool                faciesLogGiven_;
  bool                doDepthConversion_;     // 
  bool                parallelTimeSurfaces_;
  bool                useLocalWavelet_;       // Wavelets are mul;tiplied with gain and shift maps

  int                 logLevel_;      

  int                 seed_;                  // Random seed.
                      
  static int          debugFlag_;
  static std::string  filePrefix_;            // Prefix (including path) for all output files
                      
};

#endif
