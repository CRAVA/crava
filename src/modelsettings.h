#ifndef MODELSETTINGS_H
#define MODELSETTINGS_H

#include <stdio.h>
#include <string.h>
#include <map>

#include "src/definitions.h"
#include "src/io.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/traceheader.hpp"
#include "nrlib/segy/segy.hpp"

class Simbox;
class Vario;

class ModelSettings
{
public:
  ModelSettings(void);
  ~ModelSettings(void); 

  Vario                          * getAngularCorr(void)                 const { return angularCorr_                               ;} 
  Vario                          * getLateralCorr(void)                 const { return lateralCorr_                               ;}
  Vario                          * getBackgroundVario(void)             const { return backgroundVario_                           ;} 
  Vario                          * getLocalWaveletVario(void)           const { return localWaveletVario_                         ;} 
  SegyGeometry                   * getAreaParameters(void)              const { return geometry_                                  ;}    
  TraceHeaderFormat              * getTraceHeaderFormat(void)           const { return traceHeaderFormat_                         ;}
  TraceHeaderFormat              * getTraceHeaderFormatOutput(void)     const { return traceHeaderFormatOutput_                   ;}
  TraceHeaderFormat              * getTraceHeaderFormat(int i)          const { return localTHF_[i]                               ;}
  int                              getNumberOfTraceHeaderFormats(void)  const { return static_cast<int>(localTHF_.size())         ;}
  int                              getKrigingParameter(void)            const { return krigingParameter_                          ;}
  float                            getConstBackValue(int i)             const { return constBackValue_[i]                         ;}
  bool                             getUseAIBackground(void)             const { return useAIBackground_                           ;}
  bool                             getUseVpVsBackground(void)           const { return useVpVsBackground_                         ;}
  int                              getNumberOfAngles(void)              const { return static_cast<int>(angle_.size())            ;} 
  int                              getSeismicType(int i)                const { return seismicType_[i]                            ;}
  float                            getAngle(int i)                      const { return angle_[i]                                  ;}
  float                            getWaveletScale(int i)               const { return waveletScale_[i]                           ;} 
  float                            getSNRatio(int i)                    const { return SNRatio_[i]                                ;} 
  float                            getWellMoveAngle(int i,int j)        const { return wellMoveAngle_[i][j]                       ;}
  float                            getWellMoveWeight(int i,int j)       const { return wellMoveWeight_[i][j]                      ;}
  int                              getNumberOfWellAngles(int i)         const { return static_cast<int>(wellMoveAngle_[i].size()) ;} 
  bool                             getMatchEnergies(int i)              const { return matchEnergies_[i]      == 1                ;} 
  bool                             getEstimateWavelet(int i)            const { return estimateWavelet_[i]    == 1                ;}
  bool                             getUseRickerWavelet(int i)           const { return useRickerWavelet_[i]   == 1                ;}
  bool                             getEstimateSNRatio(int i)            const { return estimateSNRatio_[i]    == 1                ;}
  bool                             getEstimateLocalScale(int i)         const { return estimateLocalScale_[i] == 1                ;}
  bool                             getEstimateLocalShift(int i)         const { return estimateLocalShift_[i] == 1                ;}
  bool                             getEstimateLocalNoise(int i)         const { return estimateLocalNoise_[i] == 1                ;}
  bool                             getEstimateGlobalWaveletScale(int i) const { return estimateGlobalWaveletScale_[i]==1          ;}
  bool                             getEstimateBackground(void)          const { return estimateBackground_                        ;}
  bool                             getEstimateCorrelations(void)        const { return estimateCorrelations_                      ;}
  bool                             getEstimateWaveletNoise(void)        const { return estimateWaveletNoise_                      ;}
  bool                             getEstimate3DWavelet(void)           const { return estimate3DWavelet_                         ;}
  bool                             getHasTime3DMapping(void)            const { return hasTime3DMapping_                          ;}
  int                              getWaveletDim(int i)                 const { return waveletDim_[i]                             ;}
  float                            getStretchFactor(int i)              const { return stretchFactor_[i]                          ;}
  float                            getEstRangeX(int i)                  const { return estRangeX_[i]                              ;}
  float                            getEstRangeY(int i)                  const { return estRangeY_[i]                              ;}
  float                            getRickerPeakFrequency(int i)        const { return rickerPeakFrequency_[i]                    ;}
  const std::string                getBackgroundType(void)              const { return backgroundType_                            ;}            
  const std::vector<std::string> & getLogNames(void)                    const { return logNames_                                  ;}
  const std::vector<bool>        & getInverseVelocity(void)             const { return inverseVelocity_                           ;}
  int                              getNumberOfFacies(void)              const { return static_cast<int>(faciesNames_.size())      ;}
  const std::vector<std::string> & getFaciesNames(void)                 const { return faciesNames_                               ;}
  const std::string              & getFaciesName(int i)                 const { return faciesNames_[i]                            ;}
  int                              getFaciesLabel(int i)                const { return faciesLabels_[i]                           ;}
  int                              getIndicatorBGTrend(int i)           const { return indBGTrend_[i]                             ;}
  int                              getIndicatorWavelet(int i)           const { return indWavelet_[i]                             ;}
  int                              getIndicatorFacies(int i)            const { return indFacies_[i]                              ;}
  int                              getIndicatorRealVs(int i)            const { return indRealVs_[i]                              ;}
  int                              getIndicatorFilter(int i)            const { return indFilter_[i]                              ;}
  const std::vector<int>         & getIndicatorFilter(void)             const { return indFilter_                                 ;}
  int                              getNumberOfWells(void)               const { return nWells_                                    ;}
  int                              getNumberOfSimulations(void)         const { return nSimulations_                              ;}
  float                            getAlphaMin(void)                    const { return alpha_min_                                 ;}
  float                            getAlphaMax(void)                    const { return alpha_max_                                 ;}
  float                            getBetaMin(void)                     const { return beta_min_                                  ;}
  float                            getBetaMax(void)                     const { return beta_max_                                  ;}
  float                            getRhoMin(void)                      const { return rho_min_                                   ;}
  float                            getRhoMax(void)                      const { return rho_max_                                   ;}
  float                            getVarAlphaMin(void)                 const { return var_alpha_min_                             ;}
  float                            getVarAlphaMax(void)                 const { return var_alpha_max_                             ;}
  float                            getVarBetaMin(void)                  const { return var_beta_min_                              ;}
  float                            getVarBetaMax(void)                  const { return var_beta_max_                              ;}
  float                            getVarRhoMin(void)                   const { return var_rho_min_                               ;}
  float                            getVarRhoMax(void)                   const { return var_rho_max_                               ;}
  float                            getVpVsRatio(void)                   const { return vp_vs_ratio_                               ;}
  float                            getVpVsRatioMin(void)                const { return vp_vs_ratio_min_                           ;}
  float                            getVpVsRatioMax(void)                const { return vp_vs_ratio_max_                           ;}
  float                            getRefDepth(void)                    const { return ref_depth_                                 ;}
  float                            getAverageVelocity(void)             const { return average_velocity_                          ;}
  float                            getMaxHzBackground(void)             const { return maxHz_background_                          ;}
  float                            getMaxHzSeismic(void)                const { return maxHz_seismic_                             ;}
  float                            getMaxRankCorr(void)                 const { return maxRankCorr_                               ;}
  float                            getMaxMergeDist(void)                const { return maxMergeDist_                              ;}
  float                            getMaxDevAngle(void)                 const { return maxDevAngle_                               ;}
  float                            getLowCut(void)                      const { return lowCut_                                    ;}
  float                            getHighCut(void)                     const { return highCut_                                   ;}
  float                            getWNC(void)                         const { return wnc_                                       ;}
  float                            getEnergyThreshold(void)             const { return energyThreshold_                           ;}
  float                            getMaxWellOffset(void)               const { return maxWellOffset_                             ;}
  float                            getMaxWellShift(void)                const { return maxWellShift_                              ;}
  float                            getMinRelWaveletAmp(void)            const { return minRelWaveletAmp_                          ;}
  float                            getMaxWaveletShift(void)             const { return maxWaveletShift_                           ;}
  float                            getWaveletTaperingL(void)            const { return waveletTaperingL_                          ;}
  float                            getMinSamplingDensity(void)          const { return minSamplingDensity_                        ;}
  float                            getMinHorizontalRes(void)            const { return minHorizontalRes_                          ;}
  double                           getXPadFac(void)                     const { return xPadFac_                                   ;}
  double                           getYPadFac(void)                     const { return yPadFac_                                   ;}
  double                           getZPadFac(void)                     const { return zPadFac_                                   ;}
  int                              getNXpad(void)                       const { return nxPad_                                     ;}
  int                              getNYpad(void)                       const { return nyPad_                                     ;}
  int                              getNZpad(void)                       const { return nzPad_                                     ;}
  bool                             getEstimateXYPadding(void)           const { return estimateXYPadding_                         ;}
  bool                             getEstimateZPadding(void)            const { return estimateZPadding_                          ;}
  float                            getSegyOffset(void)                  const { return segyOffset_                                ;}
  float                            getLocalSegyOffset(int i)            const { return localSegyOffset_[i]                        ;}
  float                            getPundef(void)                      const { return p_undef_                                   ;}
  double                           getLzLimit(void)                     const { return lzLimit_                                   ;}
  double                           getTimeDTop(void)                    const { return time_dTop_                                 ;}
  double                           getTimeLz(void)                      const { return time_lz_                                   ;}
  double                           getTimeDz(void)                      const { return time_dz_                                   ;}
  int                              getTimeNz(void)                      const { return time_nz_                                   ;}
  std::vector<int>                 getAreaILXL(void)                    const { return areaILXL_                                  ;}
  int                              getAreaSpecification(void)           const { return areaSpecification_                         ;}     
  bool                             getVelocityFromInversion(void)       const { return velocityFromInv_                           ;}
  bool                             getWritePrediction(void)             const { return writePrediction_                           ;}
  int                              getOutputGridsElastic(void)          const { return outputGridsElastic_                        ;}
  int                              getOutputGridsOther(void)            const { return outputGridsOther_                          ;}
  int                              getOutputGridsSeismic(void)          const { return outputGridsSeismic_                        ;}
  int                              getOutputGridFormat(void)            const { return formatFlag_                                ;}
  int                              getOutputGridDomain(void)            const { return domainFlag_                                ;}
  bool                             getOutputGridsDefaultInd(void)       const { return outputGridsDefault_                        ;}
  int                              getWellOutputFlag(void)              const { return wellFlag_                                  ;}
  int                              getWellFormatFlag(void)              const { return wellFormatFlag_                            ;}
  int                              getWaveletOutputFlag(void)           const { return waveletFlag_                               ;}
  int                              getWaveletFormatFlag(void)           const { return waveletFormatFlag_                         ;}
  bool                             getWaveletFormatManual(void)         const { return waveletFormatManual_                       ;}
  int                              getOtherOutputFlag(void)             const { return otherFlag_                                 ;}
  int                              getDebugFlag(void)                   const { return debugFlag_                                 ;}
  static int                       getDebugLevel(void)                        { return debugFlag_                                 ;}
  bool                             getFileGrid(void)                    const { return fileGrid_                                  ;}
  bool                             getEstimationMode(void)              const { return estimationMode_                            ;}
  bool                             getForwardModeling(void)             const { return forwardModeling_                           ;}
  bool                             getGenerateSeismicAfterInv(void)     const { return generateSeismicAfterInv_                   ;}
  bool                             getGenerateBackground(void)          const { return generateBackground_                        ;}
  bool                             getEstimateFaciesProb(void)          const { return estimateFaciesProb_                        ;}
  bool                             getFaciesProbRelative(void)          const { return faciesProbRelative_                        ;}
  bool                             getNoVsFaciesProb(void)              const { return noVsFaciesProb_                            ;}
  bool                             getUseFilterForFaciesProb()          const { return useFilterForProb_                          ;}
  bool                             getFaciesLogGiven(void)              const { return faciesLogGiven_                            ;}
  std::map<std::string,float>      getPriorFaciesProb(void)             const { return priorFaciesProb_                           ;} 
  int                              getIsPriorFaciesProbGiven(void)      const { return priorFaciesProbGiven_                      ;}
  bool                             getDepthDataOK(void)                 const { return depthDataOk_                               ;}
  bool                             getParallelTimeSurfaces(void)        const { return parallelTimeSurfaces_                      ;}
  bool                             getUseLocalWavelet(void)             const { return useLocalWavelet_                           ;}
  bool                             getUseLocalNoise(void)               const { return useLocalNoise_                             ;}
  bool                             getOptimizeWellLocation(void)        const { return optimizeWellLocation_                      ;}
  bool                             getNoWellNedded(void)                const { return noWellNeeded_                              ;}
  bool                             getNoSeismicNeeded(void)             const { return noSeismicNeeded_                           ;}
  int                              getLogLevel(void)                    const { return logLevel_                                  ;}
  bool                             getErrorFileFlag()                   const { return ((otherFlag_ & IO::ERROR_FILE)>0)          ;}
  bool                             getTaskFileFlag()                    const { return ((otherFlag_ & IO::TASK_FILE)>0)           ;}
  int                              getSeed(void)                        const { return seed_                                      ;}
  bool                             getDoInversion(void);                                            
  bool                             getDoDepthConversion(void)           const;
  bool                             getDoSmoothKriging(void)             const { return smoothKrigedParameters_ ;}
  void                             getTimeGradientSettings(float &distance, float &alpha);      
  int                              getEstimateNumberOfWavelets(void)    const;

  void rotateVariograms(float angle);
  void setAngularCorr(Vario * vario);    
  void setLateralCorr(Vario * vario);    
  void setBackgroundVario(Vario * vario);
  void setLocalWaveletVario(Vario * vario);
  void copyBackgroundVarioToLocalWaveletVario(void);
  void setAreaParameters(const SegyGeometry * geometry);
  void setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat);
  void setTraceHeaderFormatOutput(TraceHeaderFormat * traceHeaderFormat);
  void addTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat);
  void setKrigingParameter(int krigingParameter)          { krigingParameter_         = krigingParameter         ;}
  void setConstBackValue(int i, float constBackValue)     { constBackValue_[i]        = constBackValue           ;}
  void setUseAIBackground(bool useAIBackground)           { useAIBackground_          = useAIBackground          ;}
  void setUseVpVsBackground(bool useVpVsBackground)       { useVpVsBackground_        = useVpVsBackground        ;}
  void addSeismicType(int seismicType)                    { seismicType_.push_back(seismicType)                  ;}
  void addAngle(float angle)                              { angle_.push_back(angle)                              ;}
  void addWaveletScale(float waveletScale)                { waveletScale_.push_back(waveletScale)                ;}
  void setWaveletScale(int i, float waveletScale)         { waveletScale_[i]          = waveletScale             ;} 
  void setSNRatio(int i, float SNRatio)                   { SNRatio_[i]               = SNRatio                  ;}
  void addSNRatio(float SNRatio)                          { SNRatio_.push_back(SNRatio)                          ;}
  void setBackgroundType(std::string type)                { backgroundType_           = type                     ;}
                                                                                                                
  void addMatchEnergies(int matchEnergies)                { matchEnergies_.push_back(matchEnergies)              ;}
  void addEstimateWavelet(int estimateWavelet)            { estimateWavelet_.push_back(estimateWavelet)          ;}
  void addEstimateSNRatio(int estimateSNRatio)            { estimateSNRatio_.push_back(estimateSNRatio)          ;}
  void setEstimateSNRatio(int i, int estimateSNRatio)     { estimateSNRatio_[i] = estimateSNRatio                ;}
  void addEstimateLocalShift(int estimateShift)           { estimateLocalShift_.push_back(estimateShift)         ;}
  void addEstimateLocalScale(int estimateScale)           { estimateLocalScale_.push_back(estimateScale)         ;}
  void addEstimateLocalNoise(int estimateNoise)           { estimateLocalNoise_.push_back(estimateNoise)         ;}
  void addEstimateGlobalWaveletScale(int estimateScale)   { estimateGlobalWaveletScale_.push_back(estimateScale) ;}
  void addUseRickerWavelet(int useRicker)                 { useRickerWavelet_.push_back(useRicker)               ;}
  void addRickerPeakFrequency(float pf)                   { rickerPeakFrequency_.push_back(pf)                   ;}

  void addWaveletDim(int waveletDim)                      { waveletDim_.push_back(waveletDim)                    ;}
  void addStretchFactor(float stretchFactor)              { stretchFactor_.push_back(stretchFactor)              ;}
  void addEstRangeX(float estRangeX)                      { estRangeX_.push_back(estRangeX)                      ;}
  void addEstRangeY(float estRangeY)                      { estRangeY_.push_back(estRangeY)                      ;}
                                                                                                                 
  void setEstimateBackground(bool estimateBackground)     { estimateBackground_       = estimateBackground       ;}
  void setEstimateCorrelations(bool estimateCorrelations) { estimateCorrelations_     = estimateCorrelations     ;}
  void setEstimateWaveletNoise(bool estimateWaveletNoise) { estimateWaveletNoise_     = estimateWaveletNoise     ;}
  void setEstimate3DWavelet(bool estimate3DWavelet)       { estimate3DWavelet_        = estimate3DWavelet        ;}
  void setHasTime3DMapping(bool hasTime3DMapping)         { hasTime3DMapping_         = hasTime3DMapping         ;}
                                                                                                                 
  void addMoveAngle(float moveAngle)                      { moveAngle_.push_back(moveAngle)                      ;} //Local temporary variable
  void addMoveWeight(float moveWeight)                    { moveWeight_.push_back(moveWeight)                    ;} //Local temporary variable
  void addMoveWell(void)                                  { wellMoveAngle_.push_back(moveAngle_);                
                                                            wellMoveWeight_.push_back(moveWeight_)               ;}                                         
  void clearMoveWell(void)                                { moveAngle_.clear();                                  
                                                            moveWeight_.clear()                                  ;}
  void addIndicatorBGTrend(int indicator)                 { indBGTrend_.push_back(indicator)                     ;}
  void addIndicatorWavelet(int indicator)                 { indWavelet_.push_back(indicator)                     ;}
  void addIndicatorFacies(int indicator)                  { indFacies_.push_back(indicator)                      ;}
  void addIndicatorRealVs(int indicator)                  { indRealVs_.push_back(indicator)                      ;}
  void addIndicatorFilter(int indicator)                  { indFilter_.push_back(indicator)                      ;}
  void setIndicatorFilter(int i ,int indicator)           { indFilter_[i]             = indicator                ;} 
  void setLogName(int i, const std::string & logName)     { logNames_[i]              = NRLib::Uppercase(logName);}
  void setInverseVelocity(int i, bool inverse)            { inverseVelocity_[i]       = inverse                  ;}
  void addFaciesLabel(int faciesLabel)                    { faciesLabels_.push_back(faciesLabel)                 ;}
  void addFaciesName(const std::string & faciesName)      { faciesNames_.push_back(faciesName)                   ;}
  void setNumberOfWells(int nWells)                       { nWells_                   = nWells                   ;} 
  void setNumberOfSimulations(int nSimulations)           { nSimulations_             = nSimulations             ;} 
  void setAlphaMin(float alpha_min)                       { alpha_min_                = alpha_min                ;}
  void setAlphaMax(float alpha_max)                       { alpha_max_                = alpha_max                ;}
  void setBetaMin(float beta_min)                         { beta_min_                 = beta_min                 ;}
  void setBetaMax(float beta_max)                         { beta_max_                 = beta_max                 ;}
  void setRhoMin(float rho_min)                           { rho_min_                  = rho_min                  ;}
  void setRhoMax(float rho_max)                           { rho_max_                  = rho_max                  ;}
  void setVarAlphaMin(float var_alpha_min)                { var_alpha_min_            = var_alpha_min            ;}
  void setVarAlphaMax(float var_alpha_max)                { var_alpha_max_            = var_alpha_max            ;}
  void setVarBetaMin(float var_beta_min)                  { var_beta_min_             = var_beta_min             ;}
  void setVarBetaMax(float var_beta_max)                  { var_beta_max_             = var_beta_max             ;}
  void setVarRhoMin(float var_rho_min)                    { var_rho_min_              = var_rho_min              ;}
  void setVarRhoMax(float var_rho_max)                    { var_rho_max_              = var_rho_max              ;}
  void setVpVsRatio(float vp_vs_ratio)                    { vp_vs_ratio_              = vp_vs_ratio              ;}
  void setVpVsRatioMin(float vp_vs_ratio_min)             { vp_vs_ratio_min_          = vp_vs_ratio_min          ;}
  void setVpVsRatioMax(float vp_vs_ratio_max)             { vp_vs_ratio_max_          = vp_vs_ratio_max          ;}
  void setRefDepth(float ref_depth)                       { ref_depth_                = ref_depth                ;}
  void setAverageVelocity(float average_velocity)         { average_velocity_         = average_velocity         ;}
  void setMaxHzBackground(float maxHz_background)         { maxHz_background_         = maxHz_background         ;}
  void setMaxHzSeismic(float maxHz_seismic)               { maxHz_seismic_            = maxHz_seismic            ;}
  void setMaxRankCorr(float maxRankCorr)                  { maxRankCorr_              = maxRankCorr              ;}
  void setMaxMergeDist(float maxMergeDist)                { maxMergeDist_             = maxMergeDist             ;}
  void setMaxDevAngle(float maxDevAngle)                  { maxDevAngle_              = maxDevAngle              ;}
  void setLowCut(float lowCut)                            { lowCut_                   = lowCut                   ;}
  void setHighCut(float highCut)                          { highCut_                  = highCut                  ;}
  void setWNC(float wnc)                                  { wnc_                      = wnc                      ;}
  void setEnergyThreshold(float energyThreshold)          { energyThreshold_          = energyThreshold          ;}
  void setMinRelWaveletAmp(float minRelWaveletAmp)        { minRelWaveletAmp_         = minRelWaveletAmp         ;}
  void setMaxWaveletShift(float maxWaveletShift)          { maxWaveletShift_          = maxWaveletShift          ;}
  void setMinSamplingDensity(float minSamplingDens)       { minSamplingDensity_       = minSamplingDens          ;}
  void setMinHorizontalRes(float minHorizontalRes)        { minHorizontalRes_         = minHorizontalRes         ;}
  void setMaxWellOffset(float maxWellOffset)              { maxWellOffset_            = maxWellOffset            ;}
  void setMaxWellShift(float maxWellShift)                { maxWellShift_             = maxWellShift             ;}
  void setWaveletTaperingL(float waveletTaperingL)        { waveletTaperingL_         = waveletTaperingL         ;}
  void setXPadFac(double xPadFac)                         { xPadFac_                  = xPadFac                  ;}
  void setYPadFac(double yPadFac)                         { yPadFac_                  = yPadFac                  ;}
  void setZPadFac(double zPadFac)                         { zPadFac_                  = zPadFac                  ;}
  void setNXpad(int nxPad)                                { nxPad_                    = nxPad                    ;}
  void setNYpad(int nyPad)                                { nyPad_                    = nyPad                    ;}
  void setNZpad(int nzPad)                                { nzPad_                    = nzPad                    ;}
  void setEstimateXYPadding(bool estimateXYPadding)       { estimateXYPadding_        = estimateXYPadding        ;}
  void setEstimateZPadding(bool estimateZPadding)         { estimateZPadding_         = estimateZPadding         ;}
  void setSegyOffset(float segyOffset)                    { segyOffset_               = segyOffset               ;}
  void addLocalSegyOffset(float segyOffset)               { localSegyOffset_.push_back(segyOffset)               ;}
  void setPundef(float p_undef)                           { p_undef_                  = p_undef                  ;}
  void setLzLimit(double lzLimit)                         { lzLimit_                  = lzLimit                  ;}
  void setTimeDTop(double time_dTop)                      { time_dTop_                = time_dTop                ;}
  void setTimeLz(double time_lz)                          { time_lz_                  = time_lz                  ;}
  void setTimeDz(double time_dz)                          { time_dz_                  = time_dz                  ;}
  void setTimeNz(int time_nz)                             { time_nz_                  = time_nz                  ;}
  void setVelocityFromInversion(bool fromInversion)       { velocityFromInv_          = fromInversion            ;}
  void setAreaILXLParameters(std::vector<int> ilxl)       { areaILXL_                 = ilxl                     ;}
  void setAreaSpecification(int areaSpecification)        { areaSpecification_        = areaSpecification        ;}
  void setWritePrediction(bool write)                     { writePrediction_          = write                    ;}
  void setOutputGridFormat(int formatFlag)                { formatFlag_               = formatFlag               ;}
  void setOutputGridDomain(int domainFlag)                { domainFlag_               = domainFlag               ;}
  void setOutputGridsElastic(int outputGridsElastic)      { outputGridsElastic_       = outputGridsElastic       ;}
  void setOutputGridsOther(int outputGridsOther)          { outputGridsOther_         = outputGridsOther         ;}
  void setOutputGridsSeismic(int outputGridsSeismic)      { outputGridsSeismic_       = outputGridsSeismic       ;}
  void setOutputGridsDefaultInd(bool ind)                 { outputGridsDefault_       = ind                      ;}
  void setWellOutputFlag(int wellFlag)                    { wellFlag_                 = wellFlag                 ;}
  void setWellFormatFlag(int formatFlag)                  { wellFormatFlag_           = formatFlag               ;}
  void setWaveletOutputFlag(int waveletFlag)              { waveletFlag_              = waveletFlag              ;}
  void setWaveletFormatFlag(int formatFlag)               { waveletFormatFlag_        = formatFlag               ;}
  void setWaveletFormatManual(bool waveletManual)         { waveletFormatManual_      = waveletManual            ;}
  void setOtherOutputFlag(int otherFlag)                  { otherFlag_                = otherFlag                ;}
  void setDebugFlag(int debugFlag)                        { debugFlag_                = debugFlag                ;}
  void setFileGrid(bool fileGrid)                         { fileGrid_                 = fileGrid                 ;}
  void setEstimationMode(bool estimationMode)             { estimationMode_           = estimationMode           ;}
  void setForwardModeling(bool forwardModeling)           { forwardModeling_          = forwardModeling          ;}
  void setGenerateSeismicAfterInv( bool generateSeismic)  { generateSeismicAfterInv_  = generateSeismic          ;}
  void setGenerateBackground(bool generateBackgr)         { generateBackground_       = generateBackgr           ;}
  void setEstimateFaciesProb(bool estFaciesProb)          { estimateFaciesProb_       = estFaciesProb            ;}
  void setFaciesProbRelative(bool faciesProbRel)          { faciesProbRelative_       = faciesProbRel            ;}
  void setNoVsFaciesProb(bool noVsFaciesProb)             { noVsFaciesProb_           = noVsFaciesProb           ;}
  void setUseFilterForFaciesProb(bool useFilterForProb)   { useFilterForProb_         = useFilterForProb         ;}
  void setFaciesLogGiven(bool faciesLogGiven)             { faciesLogGiven_           = faciesLogGiven           ;}
  void addPriorFaciesProb(std::string name, float value)  { priorFaciesProb_[name]    = value                    ;}
  void setPriorFaciesProbGiven(int fpg)                   { priorFaciesProbGiven_     = fpg                      ;}
  void setDepthDataOk(bool depthDataOk)                   { depthDataOk_              = depthDataOk              ;}
  void setParallelTimeSurfaces(bool pTimeSurfaces)        { parallelTimeSurfaces_     = pTimeSurfaces            ;}
  void setUseLocalWavelet(bool useLocalWavelet)           { useLocalWavelet_          = useLocalWavelet          ;}
  void setUseLocalNoise(bool useLocalNoise)               { useLocalNoise_            = useLocalNoise            ;}
  void setOptimizeWellLocation(bool optimizeWellLoc)      { optimizeWellLocation_     = optimizeWellLoc          ;}
  void setLogLevel(int logLevel)                          { logLevel_                 = logLevel                 ;}
  void setSeed(int seed)                                  { seed_                     = seed                     ;}
  void setDoSmoothKriging(bool smooth)                    { smoothKrigedParameters_  = smooth                   ;}
  void setNoWellNeeded(bool wellNeeded)                   { noWellNeeded_             = wellNeeded               ;}
  void setTimeGradientSettings(float distance, float alpha);
  void setNoSeismicNeeded(bool seismicNeeded)             { noSeismicNeeded_          = seismicNeeded            ;}
                        
  enum          priorFacies{FACIES_FROM_WELLS,
                            FACIES_FROM_MODEL_FILE,
                            FACIES_FROM_CUBES};

  enum          sseismicTypes{STANDARDSEIS = 0, PSSEIS = 1};                          

  enum          indicators{NO, YES, NOTSET};

  enum          areaSpecification{AREA_FROM_GRID_DATA, 
                                  AREA_FROM_UTM, 
                                  AREA_FROM_SURFACE};
private:           
                   
  Vario                           * angularCorr_;                // Variogram for lateral error correlation
  Vario                           * lateralCorr_;                // Variogram for lateral parameter correlation 
  Vario                           * backgroundVario_;            // Used for lateral background correlation.
  Vario                           * localWaveletVario_;          // Used for local wavelet (gain and shift) and local noise.
                                                                 
  SegyGeometry                    * geometry_;                   // area parameters
  float                             segyOffset_;                 // Starttime for SegY cubes.
  std::vector<float>                localSegyOffset_;            // Starttime for SegY cubes per angle.
  TraceHeaderFormat               * traceHeaderFormat_;          // traceheader of input
  std::vector<TraceHeaderFormat*>   localTHF_;                   // traceheader per angle
  TraceHeaderFormat               * traceHeaderFormatOutput_;    // traceheader for output files
  int                               krigingParameter_;           
                                                                 
  std::vector<int>                  seismicType_;                // PP- or PS- seismic
  std::vector<float>                angle_;                      // Angles
  std::vector<float>                waveletScale_;               // Signal-to-noise ratio
  std::vector<float>                SNRatio_;                    // Signal-to-noise ratio
                                                                 
  std::vector<float>                moveAngle_;                  // Angles for moving wells, local temporary variable
  std::vector<float>                moveWeight_;                 // Angle weights for moving wells, local temporary variable
  std::vector<std::vector<float> >  wellMoveAngle_;              // moveAngle_ collected for all wells
  std::vector<std::vector<float> >  wellMoveWeight_;             // moveWeight_ collected for all wells

  // --------- Start Purify-motivated use of 'int' instead of 'bool' -------------
  // NBNB-PAL: I have temporarily converted the arrays below from bool ==> int to avoid annoying UMRs in Purify
  std::vector<int>                  matchEnergies_;              // Let dataVariance_ = signalVariance_
  std::vector<int>                  estimateWavelet_;            // 
  std::vector<int>                  estimateSNRatio_;            //
  std::vector<int>                  estimateLocalShift_;         // Estimate local wavelet shift
  std::vector<int>                  estimateLocalScale_;         // Estimate local wavelet scale
  std::vector<int>                  estimateLocalNoise_;         // Estimate local noise 
  std::vector<int>                  estimateGlobalWaveletScale_;  
  std::vector<int>                  useRickerWavelet_;
  // --------- End Purify-motivated use of 'int' instead of 'bool' -------------

  std::vector<float>                rickerPeakFrequency_;

  std::vector<int>                  waveletDim_;                 // Holds if 1D-wavelet (=0) or 3D-wavelet (=1)
  std::vector<float>                stretchFactor_;              // Stretch factor for pulse in 3D-wavelet
  std::vector<float>                estRangeX_;                  // Estimation range in x-direction for 3D-wavelet 
  std::vector<float>                estRangeY_;                  // Estimation range in y-direction for 3D-wavelet 

  bool                              estimateBackground_;         // In estimation mode, skip estimation of background if false
  bool                              estimateCorrelations_;       // As above, but correlations.
  bool                              estimateWaveletNoise_;       // As above, but for wavelet and noise parameters.
  bool                              estimate3DWavelet_;          // True if a 3D wavelet is estimated for at least one angle.
  bool                              hasTime3DMapping_;           // True if command time-3D-mapping is used

  std::vector<float>                constBackValue_;             // Values set for constant background model
                                                                 // Negative value ==> read from file (actual value gives format).
  bool                              useAIBackground_;            // Read in file for AI background instead of Vp background
  bool                              useVpVsBackground_;          // Read in file for VpVs background instead of Vs background
  std::string                       backgroundType_;             // background or earth model

  //The following indicators use the indicators enum above. (2 = yes, but may override in QC, 1=yes, 0=no)
  std::vector<int>                  indBGTrend_;                 // Use well to estimate background trend?
  std::vector<int>                  indWavelet_;                 // Use well to estimate wavelet? 
  std::vector<int>                  indFacies_;                  // Use well to estimate facies? 
  std::vector<int>                  indRealVs_;                  // Treat Vs log as real? 
  std::vector<int>                  indFilter_;                  // Filter elastic logs using spatial multi-parameter filter?

  std::vector<std::string>          logNames_;                   // The keywords to look for for time, sonic, shear sonic and density
  std::vector<bool>                 inverseVelocity_;            // If element 0 is true, vp comes from dt, if 1 is true, vs comes from dts in well.

  std::vector<int>                  faciesLabels_;               // Facies labels
  std::vector<std::string>          faciesNames_;                // Facies names   (nFacies = faciesNames.size())
  int                               priorFaciesProbGiven_;
  std::map<std::string, float>      priorFaciesProb_;
                   
  int                               nWells_;
  int                               nSimulations_;
                           
  float                             alpha_min_;                  // Vp - smallest allowed value
  float                             alpha_max_;                  // Vp - largest allowed value
  float                             beta_min_;                   // Vs - smallest allowed value
  float                             beta_max_;                   // Vs - largest allowed value
  float                             rho_min_;                    // Rho - smallest allowed value
  float                             rho_max_;                    // Rho - largest allowed value
                                                                 
  float                             var_alpha_min_;              //| These min and max values are used for consistency check. If  
  float                             var_alpha_max_;              //| variances are outside these ranges there is probably a
  float                             var_beta_min_;               //| problem with the logs.
  float                             var_beta_max_;               //| 
  float                             var_rho_min_;                //| The limits are for point variances. The minimum allowed variance 
  float                             var_rho_max_;                //| for parameters will be scaled with 1/dt*dt 
 
  float                             vp_vs_ratio_min_;            // Smallest Vp/Vs-ratio regarded as likely
  float                             vp_vs_ratio_max_;            // Largest Vp/Vs-ratio regarded as likely
  float                             vp_vs_ratio_;                // Vp/Vs-ratio from input (for wavelet estimation if requested)

  float                             ref_depth_;                  // z0 - reference depth for target area
  float                             average_velocity_;           // v0 - average velocity in target area
                                                                 
  float                             maxHz_background_;           // Background resolution (high cut frequency)
  float                             maxHz_seismic_;              // Seismic resolution (high cut frequency)
                                                                 
  float                             maxRankCorr_;                // Vp-Vs correlation threshold for regarding Vs log synthetic
  float                             maxMergeDist_;               // log entries closer than this will be merged
  float                             maxDevAngle_;                // Wells with a local deviation larger than this is treated as deviated
                                                                 
  float                             lowCut_;                     // lower limit for frequency to be inverted
  float                             highCut_;                    // upper limit for frecuency to be inverted
                                                                 
  float                             wnc_;                        // White noise component, see crava.h  
                                                                 
  float                             energyThreshold_;            // If energy in reflection trace divided by mean energy
                                                                 // in reflection traces is lower than this, the reflections
                                                                 // will be interpolated. Default 0.
  float                             maxWellOffset_;              // Maximum offset for moving of wells
  float                             maxWellShift_;               // Maximum vertical shift for moving of wells
                                                                 
  float                             minRelWaveletAmp_;           // Minimum relative wavelet amplitude. Smaller amplitudes are disregarded.
  float                             maxWaveletShift_;            // Largest allowed shift when estimating wavelet
  float                             waveletTaperingL_;           // Til Odds waveletestimering
  
  float                             minSamplingDensity_;         // Threshold value for minimum sampling density in dz
  float                             minHorizontalRes_;           // Threshold value for minimum horizontal resolution in dx and dy

  double                            xPadFac_;                    // Padding factor/fraction in x direction
  double                            yPadFac_;                    // Padding factor/fraction in y direction
  double                            zPadFac_;                    // Padding factor/fraction in z direction
                                                                 
  int                               nxPad_;                      // Number of cells to pad in x direction
  int                               nyPad_;                    
  int                               nzPad_;                    
                                                               
  bool                              estimateXYPadding_;          // Estimate the z-padding from ranges
  bool                              estimateZPadding_;           // Estimate the z-padding from wavelet length

  float                             p_undef_;                    // Level for undefined facies
                                                                 
  double                            lzLimit_;                    // Minimum allowed value for (min interval thickness)/(max interval thickness)
  double                            time_dTop_;                  // Used when top and base surfaces are parallel
  double                            time_lz_;                    // Used when top and base surfaces are parallel
  double                            time_dz_;                    // Used when top and base surfaces are parallel
  int                               time_nz_;                    // Used when top and base surfaces are parallel
  bool                              velocityFromInv_;            // Velocity for time depth from inverted Vp.

  int                               areaSpecification_;          // Specifying whether are is taken from UTM-coord, seismic or surface
  std::vector<int>                  areaILXL_;                   // Vector with 6 elements (if used), in this order:
                                                                 // [0] = IL start
                                                                 // [1] = IL end
                                                                 // [2] = XL start
                                                                 // [3] = XL end
                                                                 // [4] = IL step
                                                                 // [5] = XL step
                                                                 
  bool                              writePrediction_;            // Determines whether prediction is written.
  int                               outputGridsElastic_;         // Decides which elastic grids to be written to file.
  int                               outputGridsOther_;           // Decides other grid output to be written to file
  int                               outputGridsSeismic_;         // Decides seismic grid output to be written to file.
  int                               domainFlag_;                 // Decides writing in time and/or depth.
  int                               formatFlag_;                 // Decides output format, see above.
  int                               wellFlag_;                   // Decides well output.
  int                               wellFormatFlag_;             // Decides well output format.
  int                               waveletFlag_;                // Decides wavelet output
  int                               waveletFormatFlag_;          // Decides wavelet output format
  int                               otherFlag_;                  // Decides output beyond grids and wells.
  bool                              fileGrid_;                   // Indicator telling if grids are to be kept on file
  bool                              outputGridsDefault_;         // Indicator telling if grid output has been actively controlled 
  bool                              waveletFormatManual_;        // True if wavelet format is decided in the model file

  bool                              forwardModeling_;            // Forward modelling
  bool                              estimationMode_;             // Estimation
  bool                              generateSeismicAfterInv_;    // Synthetic seismic from inversion result
  bool                              generateBackground_;         // Make background model
  bool                              estimateFaciesProb_;         // Shall facies probabilites be estimated?
  bool                              faciesProbRelative_;         // Use relative elastic parameters for facies prob estimation?
  bool                              noVsFaciesProb_;             // Do not use Vs for faciesprob.
  bool                              useFilterForProb_;           // Use filtered logs for facies probs, otherwise, use sampled inversion.
  bool                              faciesLogGiven_;           
  bool                              depthDataOk_;                // We have what we need to do depth conversion
  bool                              parallelTimeSurfaces_;     
  bool                              useLocalWavelet_;            // Wavelets are multiplied with gain and shift maps
  bool                              useLocalNoise_;              // Signal-to-noise is multiplied with gain and shift maps
  bool                              optimizeWellLocation_;       // True if at least one well is to be moved
  bool                              smoothKrigedParameters_;    // True if we should smooth borders between kriging blocks                                  
  bool                              noWellNeeded_;               // True for some configurations of input data
  bool                              noSeismicNeeded_;            // True for some estimation settings                             
  float                             distanceFromWell_;           //Minimum distance for where gradients should not cross
  float                             sigma_m_;                      //smoothness level of the gradients

  int                               logLevel_;      
                                    
  int                               seed_;                       // Random seed.
                                    
  static int                        debugFlag_;
};

#endif
