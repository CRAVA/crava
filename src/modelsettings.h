
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

  Vario                          * getAngularCorr(void)                const { return angularCorr_           ;} 
  Vario                          * getLateralCorr(void)                const { return lateralCorr_           ;}
  Vario                          * getBackgroundVario(void)            const { return backgroundVario_       ;} 
  Vario                          * getLocalWaveletVario(void)          const { return localWaveletVario_     ;} 
  SegyGeometry                   * getAreaParameters(void)             const { return geometry_              ;}    
  TraceHeaderFormat              * getTraceHeaderFormat(void)          const { return traceHeaderFormat_     ;}
  TraceHeaderFormat              * getTraceHeaderFormat(int i)         const { return localTHF_[i]           ;}
  int                              getKrigingParameter(void)           const { return krigingParameter_      ;}
  float                            getConstBackValue(int i)            const { return constBackValue_[i]     ;}
  int                              getNumberOfAngles(void)             const { return angle_.size()          ;} 
  int                              getSeismicType(int i)               const { return seismicType_[i]        ;}
  float                            getAngle(int i)                     const { return angle_[i]              ;}
  float                            getWaveletScale(int i)              const { return waveletScale_[i]       ;} 
  float                            getSNRatio(int i)                   const { return SNRatio_[i]            ;} 
  bool                             getMatchEnergies(int i)             const { return bool(matchEnergies_[i])      ;} 
  bool                             getEstimateWavelet(int i)           const { return bool(estimateWavelet_[i])    ;}
  bool                             getEstimateSNRatio(int i)           const { return bool(estimateSNRatio_[i])    ;}
  bool                             getEstimateLocalShift(int i)        const { return bool(estimateLocalShift_[i]) ;}
  bool                             getEstimateLocalScale(int i)        const { return bool(estimateLocalScale_[i]) ;}
  bool                             getEstimateBackground()             const { return estimateBackground_    ;}
  bool                             getEstimateCorrelations()           const { return estimateCorrelations_  ;}
  bool                             getEstimateWaveletNoise()           const { return estimateWaveletNoise_  ;}
  const std::vector<std::string> & getLogNames(void)                   const { return logNames_              ;}
  const std::vector<bool>        & getInverseVelocity(void)            const { return inverseVelocity_       ;}
  int                              getNumberOfFacies(void)             const { return faciesNames_.size()    ;}
  const std::string              & getFaciesName(int i)                const { return faciesNames_[i]        ;}
  int                              getFaciesLabel(int i)               const { return faciesLabels_[i]       ;}
  int                              getIndicatorBGTrend(int i)          const { return indBGTrend_[i]         ;}
  int                              getIndicatorWavelet(int i)          const { return indWavelet_[i]         ;}
  int                              getIndicatorFacies(int i)           const { return indFacies_[i]          ;}
  int                              getNumberOfWells(void)              const { return nWells_                ;}
  int                              getNumberOfSimulations(void)        const { return nSimulations_          ;}
  float                            getAlphaMin(void)                   const { return alpha_min_             ;}
  float                            getAlphaMax(void)                   const { return alpha_max_             ;}
  float                            getBetaMin(void)                    const { return beta_min_              ;}
  float                            getBetaMax(void)                    const { return beta_max_              ;}
  float                            getRhoMin(void)                     const { return rho_min_               ;}
  float                            getRhoMax(void)                     const { return rho_max_               ;}
  float                            getVarAlphaMin(void)                const { return var_alpha_min_         ;}
  float                            getVarAlphaMax(void)                const { return var_alpha_max_         ;}
  float                            getVarBetaMin(void)                 const { return var_beta_min_          ;}
  float                            getVarBetaMax(void)                 const { return var_beta_max_          ;}
  float                            getVarRhoMin(void)                  const { return var_rho_min_           ;}
  float                            getVarRhoMax(void)                  const { return var_rho_max_           ;}
  float                            getMaxHzBackground(void)            const { return maxHz_background_      ;}
  float                            getMaxHzSeismic(void)               const { return maxHz_seismic_         ;}
  float                            getMaxRankCorr(void)                const { return maxRankCorr_           ;}
  float                            getMaxMergeDist(void)               const { return maxMergeDist_          ;}
  float                            getMaxDevAngle(void)                const { return maxDevAngle_           ;}
  float                            getLowCut(void)                     const { return lowCut_                ;}
  float                            getHighCut(void)                    const { return highCut_               ;}
  float                            getWNC(void)                        const { return wnc_                   ;}
  float                            getEnergyThreshold(void)            const { return energyThreshold_       ;}
  float                            getMinRelWaveletAmp(void)           const { return minRelWaveletAmp_      ;}
  float                            getMaxWaveletShift(void)            const { return maxWaveletShift_       ;}
  float                            getWaveletTaperingL(void)           const { return waveletTaperingL_      ;}
  double                           getXPadFac(void)                    const { return xPadFac_               ;}
  double                           getYPadFac(void)                    const { return yPadFac_               ;}
  double                           getZPadFac(void)                    const { return zPadFac_               ;}
  int                              getNXpad(void)                      const { return nxPad_                 ;}
  int                              getNYpad(void)                      const { return nyPad_                 ;}
  int                              getNZpad(void)                      const { return nzPad_                 ;}
  float                            getSegyOffset(void)                 const { return segyOffset_            ;}
  float                            getLocalSegyOffset(int i)           const { return localSegyOffset_[i]    ;}
  float                            getPundef(void)                     const { return p_undef_               ;}
  double                           getLzLimit(void)                    const { return lzLimit_               ;}
  double                           getTimeDTop(void)                   const { return time_dTop_             ;}
  double                           getTimeLz(void)                     const { return time_lz_               ;}
  double                           getTimeDz(void)                     const { return time_dz_               ;}
  int                              getTimeNz(void)                     const { return time_nz_               ;}
  bool                             getVelocityFromInversion()          const { return velocityFromInv_       ;}
  bool                             getWritePrediction(void)            const { return writePrediction_       ;}
  int                              getGridOutputFlag(void)             const { return gridFlag_              ;}
  bool                             getDefaultGridOutputInd(void)       const { return defaultGridOutput_     ;}
  int                              getOutputFormatFlag(void)           const { return formatFlag_            ;}
  int                              getOutputDomainFlag(void)           const { return domainFlag_            ;}
  int                              getWellOutputFlag(void)             const { return wellFlag_              ;}
  int                              getOtherOutputFlag(void)            const { return otherFlag_             ;}
  bool                             getDirectBGOutput(void)             const { return directBGOutput_        ;}
  bool                             getDirectSeisOutput(void)           const { return directSeisOutput_      ;}
  bool                             getDirectVelOutput(void)            const { return directVelOutput_       ;}
  bool                             getDirectBGInput(void)              const { return directBGInput_         ;}
  bool                             getDirectSeisInput(void)            const { return directSeisInput_       ;}
  bool                             getDirectVelInput(void)             const { return directVelInput_        ;}
  int                              getDebugFlag(void)                  const { return debugFlag_             ;}
  static int                       getDebugLevel(void)                       { return debugFlag_             ;}
  static std::string               getOutputPath(void)                       { return outputPath_            ;}
  int                              getFileGrid(void)                   const { return fileGrid_              ;}
  bool                             getEstimationMode(void)             const { return estimationMode_        ;}
  bool                             getGenerateSeismic(void)            const { return generateSeismic_       ;}
  bool                             getGenerateBackground(void)         const { return generateBackground_    ;}
  bool                             getEstimateFaciesProb(void)         const { return estimateFaciesProb_    ;}
  bool                             getFaciesLogGiven(void)             const { return faciesLogGiven_        ;}
  bool                             getDepthDataOK(void)                const { return depthDataOk_           ;}
  bool                             getParallelTimeSurfaces(void)       const { return parallelTimeSurfaces_  ;}
  bool                             getUseLocalWavelet(void)            const { return useLocalWavelet_       ;}
  int                              getLogLevel(void)                   const { return logLevel_              ;}
  int                              getSeed(void)                       const { return seed_                  ;}
  bool                             getDoInversion(void);                                            
  bool                             getDoDepthConversion(void)          const;
                                                               
  Surface                        * getCorrelationSurface()             const { return NULL                   ;}
  void                             rotateVariograms(float angle);
                                   
  void                             setAngularCorr(Vario * vario);    
  void                             setLateralCorr(Vario * vario);    
  void                             setBackgroundVario(Vario * vario);
  void                             setLocalWaveletVario(Vario * vario);
  void                             copyBackgroundVarioToLocalWaveletVario(void);
  void                             setAreaParameters(const SegyGeometry * geometry);
  void                             setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat);
  void                             addTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat);
  void                             setKrigingParameter(int krigingParameter)     { krigingParameter_     = krigingParameter     ;}
  void                             setConstBackValue(int i, float constBackValue){ constBackValue_[i]    = constBackValue       ;}
  void                             addSeismicType(int seismicType)               { seismicType_.push_back(seismicType)          ;}
  void                             addAngle(float angle)                         { angle_.push_back(angle)                      ;}
  void                             addWaveletScale(float waveletScale)           { waveletScale_.push_back(waveletScale)        ;}
  void                             setWaveletScale(int i, float waveletScale)    { waveletScale_[i] = waveletScale              ;} 
  void                             setSNRatio(int i, float SNRatio)              { SNRatio_[i]           = SNRatio              ;}
  void                             addSNRatio(float SNRatio)                     { SNRatio_.push_back(SNRatio)                  ;}
                                   
  void                             addMatchEnergies(int matchEnergies)           { matchEnergies_.push_back(matchEnergies)      ;}
  void                             addEstimateWavelet(int estimateWavelet)       { estimateWavelet_.push_back(estimateWavelet)  ;}
  void                             addEstimateSNRatio(int estimateSNRatio)       { estimateSNRatio_.push_back(estimateSNRatio)  ;}
  void                             addEstimateLocalShift(int estimateShift)      { estimateLocalShift_.push_back(estimateShift) ;}
  void                             addEstimateLocalScale(int estimateScale)      { estimateLocalScale_.push_back(estimateScale) ;}
                                   
  void                             setEstimateBackground(bool estimateBackground)     { estimateBackground_   = estimateBackground   ;}
  void                             setEstimateCorrelations(bool estimateCorrelations) { estimateCorrelations_ = estimateCorrelations ;}
  void                             setEstimateWaveletNoise(bool estimateWaveletNoise) { estimateWaveletNoise_ = estimateWaveletNoise ;}
                                   
  void                             setAllIndicatorsTrue(int nWells);                  
  void                             setIndicatorBGTrend(int * indBGTrend, int nWells); //NBNB kill when xml-model ok.
  void                             setIndicatorWavelet(int * indWavelet, int nWells); //NBNB kill when xml-model ok. 
  void                             setIndicatorFacies(int * indFacies, int nWells);   //NBNB kill when xml-model ok.
  void                             addIndicatorBGTrend(int indicator)            { indBGTrend_.push_back(indicator)           ;}
  void                             addIndicatorWavelet(int indicator)            { indWavelet_.push_back(indicator)           ;}
  void                             addIndicatorFacies(int indicator)             { indFacies_.push_back(indicator)            ;}
  void                             setLogName(int i, const std::string & logName){ logNames_[i]          = logName            ;}
  void                             setInverseVelocity(int i, bool inverse)       { inverseVelocity_[i]          = inverse     ;}
  void                             addFaciesLabel(int faciesLabel)               { faciesLabels_.push_back(faciesLabel)       ;}
  void                             addFaciesName(const std::string & faciesName) { faciesNames_.push_back(faciesName)         ;}
  void                             setNumberOfWells(int nWells)                  { nWells_               = nWells             ;} 
  void                             setNumberOfSimulations(int nSimulations)      { nSimulations_         = nSimulations       ;} 
  void                             setAlphaMin(float alpha_min)                  { alpha_min_            = alpha_min          ;}
  void                             setAlphaMax(float alpha_max)                  { alpha_max_            = alpha_max          ;}
  void                             setBetaMin(float beta_min)                    { beta_min_             = beta_min           ;}
  void                             setBetaMax(float beta_max)                    { beta_max_             = beta_max           ;}
  void                             setRhoMin(float rho_min)                      { rho_min_              = rho_min            ;}
  void                             setRhoMax(float rho_max)                      { rho_max_              = rho_max            ;}
  void                             setVarAlphaMin(float var_alpha_min)           { var_alpha_min_        = var_alpha_min      ;}
  void                             setVarAlphaMax(float var_alpha_max)           { var_alpha_max_        = var_alpha_max      ;}
  void                             setVarBetaMin(float var_beta_min)             { var_beta_min_         = var_beta_min       ;}
  void                             setVarBetaMax(float var_beta_max)             { var_beta_max_         = var_beta_max       ;}
  void                             setVarRhoMin(float var_rho_min)               { var_rho_min_          = var_rho_min        ;}
  void                             setVarRhoMax(float var_rho_max)               { var_rho_max_          = var_rho_max        ;}
  void                             setMaxHzBackground(float maxHz_background)    { maxHz_background_     = maxHz_background   ;}
  void                             setMaxHzSeismic(float maxHz_seismic)          { maxHz_seismic_        = maxHz_seismic      ;}
  void                             setMaxRankCorr(float maxRankCorr)             { maxRankCorr_          = maxRankCorr        ;}
  void                             setMaxMergeDist(float maxMergeDist)           { maxMergeDist_         = maxMergeDist       ;}
  void                             setMaxDevAngle(float maxDevAngle)             { maxDevAngle_          = maxDevAngle        ;}
  void                             setLowCut(float lowCut)                       { lowCut_               = lowCut             ;}
  void                             setHighCut(float highCut)                     { highCut_              = highCut            ;}
  void                             setWNC(float wnc)                             { wnc_                  = wnc                ;}
  void                             setEnergyThreshold(float energyThreshold)     { energyThreshold_      = energyThreshold    ;}
  void                             setMinRelWaveletAmp(float minRelWaveletAmp)   { minRelWaveletAmp_     = minRelWaveletAmp   ;}
  void                             setMaxWaveletShift(float maxWaveletShift)     { maxWaveletShift_      = maxWaveletShift    ;}
  void                             setWaveletTaperingL(float waveletTaperingL)   { waveletTaperingL_     = waveletTaperingL   ;}
  void                             setXPadFac(double xPadFac)                    { xPadFac_              = xPadFac            ;}
  void                             setYPadFac(double yPadFac)                    { yPadFac_              = yPadFac            ;}
  void                             setZPadFac(double zPadFac)                    { zPadFac_              = zPadFac            ;}
  void                             setNXpad(int nxPad)                           { nxPad_                = nxPad              ;}
  void                             setNYpad(int nyPad)                           { nyPad_                = nyPad              ;}
  void                             setNZpad(int nzPad)                           { nzPad_                = nzPad              ;}
  void                             setSegyOffset(float segyOffset)               { segyOffset_           = segyOffset         ;}
  void                             addLocalSegyOffset(float segyOffset)          { localSegyOffset_.push_back(segyOffset)     ;}
  void                             setPundef(float p_undef)                      { p_undef_              = p_undef            ;}
  void                             setLzLimit(double lzLimit)                    { lzLimit_              = lzLimit            ;}
  void                             setTimeDTop(double time_dTop)                 { time_dTop_            = time_dTop          ;}
  void                             setTimeLz(double time_lz)                     { time_lz_              = time_lz            ;}
  void                             setTimeDz(double time_dz)                     { time_dz_              = time_dz            ;}
  void                             setTimeNz(int time_nz)                        { time_nz_              = time_nz            ;}
  void                             setVelocityFromInversion(bool fromInversion)  { velocityFromInv_      = fromInversion      ;}
  void                             setWritePrediction(bool write)                { writePrediction_      = write              ;}
  void                             setGridOutputFlag(int gridFlag)               { gridFlag_             = gridFlag           ;}
  void                             setDefaultGridOutputInd(bool ind)             { defaultGridOutput_    = ind                ;}
  void                             setOutputFormatFlag(int formatFlag)           { formatFlag_           = formatFlag         ;}
  void                             setOutputDomainFlag(int domainFlag)           { domainFlag_           = domainFlag         ;}
  void                             setWellOutputFlag(int wellFlag)               { wellFlag_             = wellFlag           ;}
  void                             setOtherOutputFlag(int otherFlag)             { otherFlag_            = otherFlag          ;}
  void                             setDirectBGOutput(bool directBGOutput)        { directBGOutput_       = directBGOutput     ;}
  void                             setDirectSeisOutput(bool directSeisOutput)    { directSeisOutput_     = directSeisOutput   ;}
  void                             setDirectVelOutput(bool directVelOutput)      { directVelOutput_      = directVelOutput    ;}
  void                             setDirectBGInput(bool directBGInput)          { directBGInput_        = directBGInput      ;}
  void                             setDirectSeisInput(bool directSeisInput)      { directSeisInput_      = directSeisInput    ;}
  void                             setDirectVelInput(bool directVelInput)        { directVelInput_       = directVelInput    ;}
  void                             setDebugFlag(int debugFlag)                   { debugFlag_            = debugFlag          ;}
  void                             setFileGrid(int fileGrid)                     { fileGrid_             = fileGrid           ;}
  void                             setEstimationMode(bool estimationMode)        { estimationMode_       = estimationMode     ;}
  void                             setGenerateSeismic(bool generateSeismic)      { generateSeismic_      = generateSeismic    ;}
  void                             setGenerateBackground(bool generateBackgr)    { generateBackground_   = generateBackgr     ;}
  void                             setEstimateFaciesProb(bool estFaciesProb)     { estimateFaciesProb_   = estFaciesProb      ;}
  void                             setFaciesLogGiven(bool faciesLogGiven)        { faciesLogGiven_       = faciesLogGiven     ;}
  void                             setDepthDataOk(bool depthDataOk)              { depthDataOk_          = depthDataOk        ;}
  void                             setParallelTimeSurfaces(bool pTimeSurfaces)   { parallelTimeSurfaces_ = pTimeSurfaces      ;}
  void                             setUseLocalWavelet(bool useLocalWavelet)      { useLocalWavelet_      = useLocalWavelet    ;}
  void                             setLogLevel(int logLevel)                     { logLevel_             = logLevel           ;}
  void                             setSeed(int seed)                             { seed_                 = seed               ;}

  enum                             outputGrids{CORRELATION        = 1, 
                                               RESIDUAL           = 2, 
                                               VP                 = 4, 
                                               VS                 = 8, 
                                               RHO                = 16,
                                               LAMELAMBDA         = 32, 
                                               LAMEMU             = 64, 
                                               POISSONRATIO       = 128, 
                                               AI                 = 256,
                                               SI                 = 512, 
                                               VPVSRATIO          = 1024, 
                                               MURHO              = 2048, 
                                               LAMBDARHO          = 4096, 
                                               BACKGROUND         = 8192, 
                                               BACKGROUND_TREND   = 16384, 
                                               FACIESPROB         = 32768, 
                                               FACIESPROBRELATIVE = 65536, 
                                               EXTRA_GRIDS        = 131072};
                                   
  enum                             outputWells{WELLS              = 1,
                                               BLOCKED_WELLS      = 2,
                                               BLOCKED_LOGS       = 4};
                                   
  enum                             outputOther{WAVELETS            = 1,
                                               EXTRA_SURFACES      = 2,
                                               PRIORCORRELATIONS   = 4,
                                               BACKGROUND_TREND_1D = 8};
                                   
  enum                             domains{TIMEDOMAIN = 1, DEPTHDOMAIN = 2};
                                   
  enum                             gridFormats{SEGY = 1, STORM = 2, ASCII = 4, SGRI = 8};
                                   
  enum                             sseismicTypes{STANDARDSEIS = 0, PSSEIS = 1};
                                     
  //Note: By convention, input path is added to input file names at end of parsing.
  //      Output path and prefix is added to output file name by call to makeFullFileName
  //      just before writing.
  static void                      setFilePrefix(const std::string & filePrefix);
  static void                      setOutputPath(const std::string & outputPath);
  static std::string               makeFullFileName(const std::string name, const std::string postfix = "");
                                   
private:                           
                   
  Vario                          * angularCorr_;           ///< Variogram for lateral error correlation
  Vario                          * lateralCorr_;           ///< Variogram for lateral parameter correlation 
  Vario                          * backgroundVario_;       ///< Used for lateral background correlation.
  Vario                          * localWaveletVario_;     ///< Used for local wavelet (gain and shift) and local noise.
                          
  SegyGeometry                   * geometry_;              ///< area parameters
  float                            segyOffset_;            ///< Starttime for SegY cubes.
  std::vector<float>               localSegyOffset_;       ///< Starttime for SegY cubes per angle.
  TraceHeaderFormat              * traceHeaderFormat_;     ///< traceheader
  std::vector<TraceHeaderFormat*>  localTHF_;              ///< traceheader per angle

  int                              krigingParameter_;   
                          
  std::vector<int>                 seismicType_;           ///< PP- or PS- seismic
  std::vector<float>               angle_;                 ///< Angles
  std::vector<float>               waveletScale_;          ///< Signal-to-noise ratio
  std::vector<float>               SNRatio_;               ///< Signal-to-noise ratio

  // NBNB-PAL: I have temporarily converted the five arrays below from bool ==> int to avoid annoying UMRs in Purify
  std::vector<int>                 matchEnergies_;         ///< Let dataVariance_ = signalVariance_
  std::vector<int>                 estimateWavelet_;       ///< 
  std::vector<int>                 estimateSNRatio_;       ///<
  std::vector<int>                 estimateLocalShift_;    ///< Estimate local wavelet shift
  std::vector<int>                 estimateLocalScale_;    ///< Estimate local wavelet scale
                                 
  bool                             estimateBackground_;    ///< In estimation mode, skip estimation of background if false
  bool                             estimateCorrelations_;  ///< As above, but correlations.
  bool                             estimateWaveletNoise_;  ///< As above, but for wavelet and noise parameters.
                                 
  std::vector<float>               constBackValue_;        ///< Values set for constant background model
                                                           ///< Negative value ==> read from file (actual value gives format).
  std::vector<int>                 indBGTrend_;            ///< Use well to estimate background trend? (1=yes,0=no)
  std::vector<int>                 indWavelet_;            ///< Use well to estimate wavelet? (1=yes,0=no)
  std::vector<int>                 indFacies_;             ///< Use well to estimate facies? (1=yes,0=no)

  std::vector<std::string>         logNames_;              ///</< The keywords to look for for time, sonic, shear sonic and density
  std::vector<bool>                inverseVelocity_;       ///< If element 0 is true, vp comes from dt, if 1 is true, vs comes from dts in well.

  std::vector<int>                 faciesLabels_;          ///< Facies labels
  std::vector<std::string>         faciesNames_;           ///< Facies names   (nFacies = faciesNames.size())
                   
  int                              nWells_;
  int                              nSimulations_;
                           
  float                            alpha_min_;             ///< Vp - smallest allowed value
  float                            alpha_max_;             ///< Vp - largest allowed value
  float                            beta_min_;              ///< Vs - smallest allowed value
  float                            beta_max_;              ///< Vs - largest allowed value
  float                            rho_min_;               ///< Rho - smallest allowed value
  float                            rho_max_;               ///< Rho - largest allowed value
                                  
  float                            var_alpha_min_;         ///<| These min and max values are used for consistency check. If  
  float                            var_alpha_max_;         ///<| variances are outside these ranges there is probably a
  float                            var_beta_min_;          ///<| problem with the logs.
  float                            var_beta_max_;          ///<| 
  float                            var_rho_min_;           ///<| The limits are for point variances. The minimum allowed variance 
  float                            var_rho_max_;           ///<| for parameters will be scaled with 1/dt*dt
                                   
  float                            maxHz_background_;      ///< Background resolution (high cut frequency)
  float                            maxHz_seismic_;         ///< Seismic resolution (high cut frequency)
                                  
  float                            maxRankCorr_;           ///< Vp-Vs correlation threshold for regarding Vs log synthetic
  float                            maxMergeDist_;          ///< log entries closer than this will be merged
  float                            maxDevAngle_;           ///< Wells with a local deviation larger than this is treated as deviated
                                  
  float                            lowCut_;                ///< lower limit for frequency to be inverted
  float                            highCut_;               ///< upper limit for frecuency to be inverted
                                  
  float                            wnc_;                   ///< White noise component, see crava.h  
                                  
  float                            energyThreshold_;       ///< If energy in reflection trace divided by mean energy
                                                           ///< in reflection traces is lower than this, the reflections
                                                           ///< will be interpolated. Default 0.
                                  
  float                            minRelWaveletAmp_;      ///< Minimum relative wavelet amplitude. Smaller amplitudes are disregarded.
  float                            maxWaveletShift_;       ///< Largest allowed shift when estimating wavelet
  float                            waveletTaperingL_;      ///< Til Odds waveletestimering
                                  
  double                           xPadFac_;               ///< Padding factor/fraction in x direction
  double                           yPadFac_;               ///< Padding factor/fraction in y direction
  double                           zPadFac_;               ///< Padding factor/fraction in z direction
                                  
  int                              nxPad_;                 ///< Number of cells to pad in x direction
  int                              nyPad_;
  int                              nzPad_; 
                                  
  float                            p_undef_;               ///< NBNB-PAL: Hva gjør denne?
                                  
  double                           lzLimit_;               ///< Minimum allowed value for (min interval thickness)/(max interval thickness)
  double                           time_dTop_;             ///< Used when top and base surfaces are parallel
  double                           time_lz_;               ///< Used when top and base surfaces are parallel
  double                           time_dz_;               ///< Used when top and base surfaces are parallel
  int                              time_nz_;               ///< Used when top and base surfaces are parallel
  bool                             velocityFromInv_;       ///< Velocity for time depth from inverted Vs.
                                   
  bool                             writePrediction_;       ///< Determines whether prediction is written.
  int                              gridFlag_;              ///< Decides which grids to write (except simulation)
  int                              domainFlag_;            ///< Decides writing in time and/or depth.
  int                              formatFlag_;            ///< Decides output format, see above.
  int                              wellFlag_;              ///< Decides well output.
  int                              otherFlag_;             ///< Decides output beyond grids and wells.
  int                              fileGrid_;              ///< Indicator telling if grids are to be kept on file
  bool                             defaultGridOutput_;     ///< Indicator telling whether grid output has been actively controlled.
                                   
  bool                             directBGOutput_;        ///< Write raw background, can be read without resampling.
  bool                             directSeisOutput_;      ///< Write raw seismic, can be read without resampling.
  bool                             directVelOutput_;       ///< Write raw time-to-depth velocity, can be read without resampling.
  bool                             directBGInput_;         ///< Read raw background into grid without resampling.
  bool                             directSeisInput_;       ///< Read raw seismic into grid without resampling.
  bool                             directVelInput_;        ///< Read raw time-to-depth velocity into grid without resampling.
                                   
  bool                             generateSeismic_;       ///< Forward modelling
  bool                             estimationMode_;        ///< Estimation
                                   
  bool                             generateBackground_;    ///< Make background model
  bool                             estimateFaciesProb_;    ///< Shall facies probabilites be estimated?
  bool                             faciesLogGiven_;
  bool                             depthDataOk_;           ///< We have what we need to do depth conversion
  bool                             parallelTimeSurfaces_;
  bool                             useLocalWavelet_;       ///< Wavelets are multiplied with gain and shift maps

  int                              logLevel_;      
                           
  int                              seed_;                  ///< Random seed.
                           
  static int                       debugFlag_;

  static std::string               outputPath_;            ///< Path for all output files.
  static std::string               filePrefix_;            ///< Prefix for all output files
};

#endif
