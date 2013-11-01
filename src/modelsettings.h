/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#ifndef MODELSETTINGS_H
#define MODELSETTINGS_H

#include <stdio.h>
#include <string.h>
#include <map>

#include "src/definitions.h"
#include "src/io.h"
#include "src/vario.h"
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/traceheader.hpp"
#include "nrlib/segy/segy.hpp"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsfluidstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionwithtrendstorage.h"

class Simbox;

class ModelSettings
{
public:
  ModelSettings(void);
  ~ModelSettings(void);

  Vario                          * getAngularCorr(int i)                const { return angularCorr_[i]                            ;}
  Vario                          * getLateralCorr(void)                 const { return lateralCorr_                               ;}
  Vario                          * getBackgroundVario(void)             const { return backgroundVario_                           ;}
  Vario                          * getLocalWaveletVario(void)           const { return localWaveletVario_                         ;}
  SegyGeometry                   * getAreaParameters(void)              const { return geometry_                                  ;}
  SegyGeometry                   * getSeismicDataAreaParameters(void)   const { return geometry_full_                             ;}
  TraceHeaderFormat              * getTraceHeaderFormat(void)           const { return traceHeaderFormat_                         ;}
  TraceHeaderFormat              * getTraceHeaderFormatOutput(void)     const { return traceHeaderFormatOutput_                   ;}
  TraceHeaderFormat              * getTraceHeaderFormat(int i, int j)   const { return timeLapseLocalTHF_[i][j]                   ;}
  int                              getNumberOfTraceHeaderFormats(int i) const { return static_cast<int>(timeLapseLocalTHF_[i].size());}
  int                              getKrigingParameter(void)            const { return krigingParameter_                          ;}
  float                            getConstBackValue(int i)             const { return constBackValue_[i]                         ;}
  bool                             getUseAIBackground(void)             const { return useAIBackground_                           ;}
  bool                             getUseSIBackground(void)             const { return useSIBackground_                           ;}
  bool                             getUseVpVsBackground(void)           const { return useVpVsBackground_                         ;}
  int                              getNumberOfAngles(int i)             const { return static_cast<int>(timeLapseAngle_[i].size());}
  int                              getNumberOfTimeLapses(void)          const { return static_cast<int>(timeLapseAngle_.size())   ;}
  bool                             getGravityTimeLapse(int i)           const { return timeLapseGravimetry_[i]                    ;}
  bool                             getTravelTimeTimeLapse(int i)        const { return timeLapseTravelTimeGiven_[i]               ;}
  std::vector<std::string>         getTimeLapseTravelTimeHorizons(int i)const { return timeLapseTravelTimeHorizon_[i]             ;}
  std::vector<double>              getTimeLapseTravelTimeHorizonSD(int i)const{ return timeLapseTravelTimeHorizonSD_[i]           ;}
  int                              getNumberOfVintages(void)            const { return static_cast<int>(vintageYear_.size())      ;}
  int                              getVintageYear(int i)                const { return vintageYear_[i]                            ;}
  int                              getVintageMonth(int i)               const { return vintageMonth_[i]                           ;}
  int                              getVintageDay(int i)                 const { return vintageDay_[i]                             ;}
  const std::vector<int>         & getSeismicType(int i)                const { return timeLapseSeismicType_[i]                   ;}
  const std::vector<float>       & getAngle(int i)                      const { return timeLapseAngle_[i]                         ;}
  float                            getWaveletScale(int i, int j)        const { return timeLapseWaveletScale_[i][j]               ;}
  const std::vector<float>       & getSNRatio(int i)                    const { return timeLapseSNRatio_[i]                       ;}
  float                            getWellMoveAngle(int i,int j)        const { return wellMoveAngle_[i][j]                       ;}
  float                            getWellMoveWeight(int i,int j)       const { return wellMoveWeight_[i][j]                      ;}
  int                              getNumberOfWellAngles(int i)         const { return static_cast<int>(wellMoveAngle_[i].size()) ;}
  const std::vector<bool>        & getMatchEnergies(int i)              const { return timeLapseMatchEnergies_[i]                 ;}
  const std::vector<bool>        & getEstimateWavelet(int i)            const { return timeLapseEstimateWavelet_[i]               ;}
  const std::vector<bool>        & getUseRickerWavelet(int i)           const { return timeLapseUseRickerWavelet_[i]              ;}
  bool                             getEstimateSNRatio(int i, int j)     const { return timeLapseEstimateSNRatio_[i][j]            ;}
  bool                             getEstimateLocalScale(int i, int j)  const { return timeLapseEstimateLocalScale_[i][j]         ;}
  bool                             getEstimateLocalShift(int i, int j)  const { return timeLapseEstimateLocalShift_[i][j]         ;}
  bool                             getEstimateLocalNoise(int i, int j)  const { return timeLapseEstimateLocalNoise_[i][j]         ;}
  bool                             getEstimateGlobalWaveletScale(int i, int j) const { return timeLapseEstimateGlobalWaveletScale_[i][j];}
  bool                             getDo4DInversion(void)               const { return do4DInversion_                             ;}
  bool                             getDo4DRockPhysicsInversion(void)    const { return do4DRockPhysicsInversion_                  ;}
  bool                             getEstimateBackground(void)          const { return estimateBackground_                        ;}
  bool                             getEstimateCorrelations(void)        const { return estimateCorrelations_                      ;}
  bool                             getEstimateWaveletNoise(void)        const { return estimateWaveletNoise_                      ;}
  bool                             getEstimate3DWavelet(void)           const { return estimate3DWavelet_                         ;}
  bool                             getHasTime3DMapping(void)            const { return hasTime3DMapping_                          ;}
  bool                             getUse3DWavelet(void)                const { return use3DWavelet_                              ;}
  int                              getWaveletDim(int i)                 const { return waveletDim_[i]                             ;}
  const std::vector<float>       & getStretchFactor(int i)              const { return timeLapseStretchFactor_[i]                 ;} //Never used...
  float                            getEstRangeX(int i)                  const { return estRangeX_[i]                              ;}
  float                            getEstRangeY(int i)                  const { return estRangeY_[i]                              ;}
  float                            getRickerPeakFrequency(int i, int j) const { return timeLapseRickerPeakFrequency_[i][j]        ;}
  const std::string              & getBackgroundType(void)              const { return backgroundType_                            ;}
  const std::vector<std::string> & getLogNames(void)                    const { return logNames_                                  ;}
  const std::vector<bool>        & getInverseVelocity(void)             const { return inverseVelocity_                           ;}
  int                              getIndicatorBGTrend(int i)           const { return indBGTrend_[i]                             ;}
  int                              getIndicatorWavelet(int i)           const { return indWavelet_[i]                             ;}
  int                              getIndicatorFacies(int i)            const { return indFacies_[i]                              ;}
  int                              getIndicatorRockPhysics(int i)       const { return indRockPhysics_[i]                         ;}
  int                              getIndicatorRealVs(int i)            const { return indRealVs_[i]                              ;}
  int                              getIndicatorFilter(int i)            const { return indFilter_[i]                              ;}
  const std::vector<int>         & getIndicatorFilter(void)             const { return indFilter_                                 ;}
  int                              getNumberOfWells(void)               const { return nWells_                                    ;}
  int                              getNumberOfSimulations(void)         const { return nSimulations_                              ;}
  float                            getTemporalCorrelationRange(void)    const { return temporalCorrelationRange_                  ;}
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
  float                            getVpVsRatioFromWells(void)          const { return vp_vs_ratio_from_wells_                    ;}
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
  float                            getDefaultWaveletLength(void)        const { return defaultWaveletLength_                      ;}
  float                            getGuardZone(void)                   const { return guard_zone_                                ;}
  float                            getSmoothLength(void)                const { return smooth_length_                             ;}
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
  float                            getSegyOffset(int i)                 const { return segyOffset_[i]                             ;}
  const std::vector<float>       & getLocalSegyOffset(int i)            const { return timeLapseLocalSegyOffset_[i]               ;}
  float                            getPundef(void)                      const { return p_undef_                                   ;}
  double                           getLzLimit(void)                     const { return lzLimit_                                   ;}
  double                           getTimeDTop(void)                    const { return time_dTop_                                 ;}
  double                           getTimeLz(void)                      const { return time_lz_                                   ;}
  double                           getTimeDz(void)                      const { return time_dz_                                   ;}
  int                              getTimeNz(void)                      const { return time_nz_                                   ;}
  int                              getTimeNzInterval(std::string name)  const { return time_nz_interval_.find(name)->second       ;}
  const std::vector<int>         & getAreaILXL(void)                    const { return areaILXL_                                  ;}
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
  bool                             getGenerateBackgroundFromRockPhysics()const{ return backgroundFromRockPhysics_                 ;}
  bool                             getCalibrateRockPhysicsToWells()     const { return calibrateRockPhysicsToWells_               ;}
  bool                             getGenerateBackground(void)          const { return generateBackground_                        ;}
  bool                             getUseVerticalVariogram(void)        const { return useVerticalVariogram_                      ;}
  bool                             getMultizoneBackground()             const { return multizoneBackground_                       ;}
  bool                             getEstimateFaciesProb(void)          const { return estimateFaciesProb_                        ;}
  bool                             getFaciesProbRelative(void)          const { return faciesProbRelative_                        ;}
  bool                             getFaciesProbFromRockPhysics(void)   const { return faciesProbFromRockPhysics_                 ;}
  bool                             getNoVsFaciesProb(void)              const { return noVsFaciesProb_                            ;}
  bool                             getUseFilterForFaciesProb()          const { return useFilterForProb_                          ;}
  bool                             getFaciesLogGiven(void)              const { return faciesLogGiven_                            ;}
  bool                             getPorosityLogGiven(void)            const { return porosityLogGiven_                          ;}
  const std::map<std::string,float>& getPriorFaciesProb(void)           const { return priorFaciesProb_                           ;}
  const std::map<std::string,float>& getVolumeFractionsProb(void)       const { return volumeFractionProb_                        ;}
  int                              getIsPriorFaciesProbGiven(void)      const { return priorFaciesProbGiven_                      ;}
  bool                             getDepthDataOK(void)                 const { return depthDataOk_                               ;}
  bool                             getParallelTimeSurfaces(void)        const { return parallelTimeSurfaces_                      ;}
  bool                             getUseLocalWavelet(void)             const { return useLocalWavelet_                           ;}
  bool                             getUseLocalNoise(int i)              const { return timeLapseUseLocalNoise_[i]                 ;}
  bool                             getOptimizeWellLocation(void)        const { return optimizeWellLocation_                      ;}
  bool                             getNoWellNedded(void)                const { return noWellNeeded_                              ;}
  bool                             getNoSeismicNeeded(void)             const { return noSeismicNeeded_                           ;}
  bool                             getSnapGridToSeismicData(void)       const { return snapGridToSeismicData_                     ;}
  double                           getWavelet3DTuningFactor(void)       const { return wavelet3DTuningFactor_                     ;}
  double                           getGradientSmoothingRange(void)      const { return gradientSmoothingRange_                    ;}
  bool                             getEstimateWellGradientFromSeismic() const { return wellGradientFromSeismic_                   ;}
  int                              getLogLevel(void)                    const { return logLevel_                                  ;}
  bool                             getErrorFileFlag()                   const { return ((otherFlag_ & IO::ERROR_FILE)>0)          ;}
  bool                             getTaskFileFlag()                    const { return ((otherFlag_ & IO::TASK_FILE)>0)           ;}
  int                              getSeed(void)                        const { return seed_                                      ;}
  bool                             getDoInversion(void)                 const;
  bool                             getDoDepthConversion(void)           const;
  bool                             getDoSmoothKriging(void)             const { return smoothKrigedParameters_ ;}
  bool                             getRunFromPanel(void)                const { return runFromPanel_ ;}
  void                             getTimeGradientSettings(float &distance, float &alpha, int t);
  int                              getEstimateNumberOfWavelets(int t)   const;
  double                           getSeismicQualityGridRange(void)     const { return seismicQualityGridRange_                   ;}
  double                           getSeismicQualityGridValue(void)     const { return seismicQualityGridValue_                   ;}
  const std::vector<int>           findSortedVintages(void)             const;
  const std::vector<std::string> & getTrendCubeParameters(void)         const { return trendCubeParameter_                  ;}
  const std::vector<int>         & getTrendCubeType(void)               const { return trendCubeType_                       ;}
  const std::map<std::string, std::vector<DistributionWithTrendStorage *> > & getReservoirVariable() const { return reservoirVariable_ ;}
  const std::map<std::string, DistributionsRockStorage *>                   & getRockStorage()       const { return rockStorage_       ;}
  const std::map<std::string, DistributionsDryRockStorage *>                & getDryRockStorage()    const { return dryRockStorage_    ;}
  const std::map<std::string, DistributionsSolidStorage *>                  & getSolidStorage()      const { return solidStorage_      ;}
  const std::map<std::string, DistributionsFluidStorage *>                  & getFluidStorage()      const { return fluidStorage_      ;}
  std::vector<int>                 getErosionPriority()                 const { return erosionPriority_                           ;}
  std::vector<int>                 getCorrelationStructure()            const { return correlationStructure_                      ;}
  std::vector<double>              getSurfaceUncertainty()              const { return surfaceUncertainty_                        ;}
  std::vector<std::string>         getIntervalNames()                   const { return interval_names_                            ;}

  double                           getRMSStandardDeviation()            const { return RMSStandardDeviation_                      ;}
  bool                             getRMSPriorGiven()                   const { return RMSPriorGiven_                             ;}
  int                              getRMSnLayersAbove()                 const { return RMSnLayersAbove_                           ;}
  int                              getRMSnLayersBelow()                 const { return RMSnLayersBelow_                           ;}
  double                           getRMSMeanVpTop()                    const { return RMSMeanVpTop_                              ;}
  double                           getRMSMeanVpBase()                   const { return RMSMeanVpBase_                             ;}
  double                           getRMSVarianceVpAbove()              const { return RMSVarianceVpAbove_                        ;}
  double                           getRMSVarianceVpBelow()              const { return RMSVarianceVpBelow_                        ;}
  double                           getRMSTemporalCorrelationRangeAbove()const { return RMSTemporalCorrRangeAbove_                 ;}
  double                           getRMSTemporalCorrelationRangeBelow()const { return RMSTemporalCorrRangeBelow_                 ;}

  std::map<std::string, float>     getVpVsRatioIntervals()              const { return vpvs_ratio_interval_                       ;}
  std::map<std::string, std::map<std::string, float> > getPriorFaciesProbInterval() const { return priorFaciesProbInterval_       ;}
  std::map<std::string, std::map<std::string, float> > getVolumeFractionsProbInterval() const { return volumefractionInterval_   ;}
  int                              getErosionPriorityTopSurface()       const { return erosion_priority_top_surface_;}
  int                              getErosionPriorityBaseSurface(const std::string & interval_name) const {return erosion_priority_interval_base_surface_.find(interval_name)->second;}

  bool                               getCorrDirTopConform(void)          const { return topConformCorrelation_ ;}
  bool                               getCorrDirBaseConform(void)         const { return baseConformCorrelation_ ;}
  const std::map<std::string, bool>& getCorrDirIntervalTopConform(void)  const { return intervalTopConformCorrelation_ ;}
  const std::map<std::string, bool>& getCorrDirIntervalBaseConform(void) const { return intervalBaseConformCorrelation_ ;}
  bool                               getCorrDirIntervalUsed(void)        const { return intervalCorrelationUsed_ ;}



  void rotateVariograms(float angle);
  void setLastAngularCorr(Vario * vario);
  void setLateralCorr(Vario * vario);
  void setBackgroundVario(Vario * vario);
  void setLocalWaveletVario(Vario * vario);
  void copyBackgroundVarioToLocalWaveletVario(void);
  void setAreaParameters(const SegyGeometry * geometry);
  void setSeismicDataAreaParameters(const SegyGeometry * geometry);
  void setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat);
  void setTraceHeaderFormatOutput(TraceHeaderFormat * traceHeaderFormat);
  void addTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat);
  void addTravelTimeTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat);
  void setKrigingParameter(int krigingParameter)          { krigingParameter_         = krigingParameter         ;}
  void setConstBackValue(int i, float constBackValue)     { constBackValue_[i]        = constBackValue           ;}
  void setUseAIBackground(bool useAIBackground)           { useAIBackground_          = useAIBackground          ;}
  void setUseSIBackground(bool useSIBackground)           { useSIBackground_          = useSIBackground          ;}
  void setUseVpVsBackground(bool useVpVsBackground)       { useVpVsBackground_        = useVpVsBackground        ;}
  void addSeismicType(int seismicType)                    { seismicType_.push_back(seismicType)                  ;}
  void addAngle(float angle)                              { angle_.push_back(angle)                              ;}
  void addWaveletScale(float waveletScale)                { waveletScale_.push_back(waveletScale)                ;}
  void setWaveletScale(int i, int j, float waveletScale)  { timeLapseWaveletScale_[i][j] = waveletScale          ;}
  void addSNRatio(float SNRatio)                          { SNRatio_.push_back(SNRatio)                          ;}
  void setBackgroundType(std::string type)                { backgroundType_           = type                     ;}

  void addMatchEnergies(bool matchEnergies)               { matchEnergies_.push_back(matchEnergies)              ;}
  void addEstimateWavelet(bool estimateWavelet)           { estimateWavelet_.push_back(estimateWavelet)          ;}
  void addEstimateSNRatio(bool estimateSNRatio)           { estimateSNRatio_.push_back(estimateSNRatio)          ;}
  void setEstimateSNRatio(int i, int j, bool estimateSNRatio){ timeLapseEstimateSNRatio_[i][j] = estimateSNRatio ;}
  void addEstimateLocalShift(bool estimateShift)          { estimateLocalShift_.push_back(estimateShift)         ;}
  void addEstimateLocalScale(bool estimateScale)          { estimateLocalScale_.push_back(estimateScale)         ;}
  void addEstimateLocalNoise(bool estimateNoise)          { estimateLocalNoise_.push_back(estimateNoise)         ;}
  void addEstimateGlobalWaveletScale(bool estimateScale)  { estimateGlobalWaveletScale_.push_back(estimateScale) ;}
  void addUseRickerWavelet(bool useRicker)                { useRickerWavelet_.push_back(useRicker)               ;}
  void addRickerPeakFrequency(float pf)                   { rickerPeakFrequency_.push_back(pf)                   ;}

  void addWaveletDim(int waveletDim)                      { waveletDim_.push_back(waveletDim)                    ;}
  void addStretchFactor(float stretchFactor)              { stretchFactor_.push_back(stretchFactor)              ;}
  void addEstRangeX(float estRangeX)                      { estRangeX_.push_back(estRangeX)                      ;}
  void addEstRangeY(float estRangeY)                      { estRangeY_.push_back(estRangeY)                      ;}

  void addErosionPriority(int priority)                   { erosionPriority_.push_back(priority)                 ;}
  void addCorrelationStructure(int structure)             { correlationStructure_.push_back(structure)           ;}
  void addSurfaceUncertainty(double uncertainty)          { surfaceUncertainty_.push_back(uncertainty)           ;}

  void addTrendCubeParameter(std::string parameterName)                  { trendCubeParameter_.push_back(parameterName)                   ;}
  void addTrendCubes(int trendCubeType)                                  { trendCubeType_.push_back(trendCubeType)                        ;}
  void addReservoirVariable(std::string variable, std::vector<DistributionWithTrendStorage *> dist) { reservoirVariable_[variable] = dist ;}
  void addRock(std::string label,  DistributionsRockStorage  * rock)                   { rockStorage_[label]    = rock                    ;}
  void addDryRock(std::string label,  DistributionsDryRockStorage  * dry_rock)         { dryRockStorage_[label] = dry_rock                ;}
  void addSolid(std::string label, DistributionsSolidStorage * solid)                  { solidStorage_[label]   = solid                   ;}
  void addFluid(std::string label, DistributionsFluidStorage * fluid)                  { fluidStorage_[label]   = fluid                   ;}

  void setDo4DInversion(bool do4DInversion)               { do4DInversion_            = do4DInversion            ;}
  void setDo4DRockPhysicsInversion(bool do4DRockPhysicsInversion)                      {do4DRockPhysicsInversion_= do4DRockPhysicsInversion;}
  void setEstimateBackground(bool estimateBackground)     { estimateBackground_       = estimateBackground       ;}
  void setEstimateCorrelations(bool estimateCorrelations) { estimateCorrelations_     = estimateCorrelations     ;}
  void setEstimateWaveletNoise(bool estimateWaveletNoise) { estimateWaveletNoise_     = estimateWaveletNoise     ;}
  void setEstimate3DWavelet(bool estimate3DWavelet)       { estimate3DWavelet_        = estimate3DWavelet        ;}
  void setHasTime3DMapping(bool hasTime3DMapping)         { hasTime3DMapping_         = hasTime3DMapping         ;}
  void setUse3DWavelet(bool use3DWavelet)                 { use3DWavelet_             = use3DWavelet             ;}

  void addMoveAngle(float moveAngle)                      { moveAngle_.push_back(moveAngle)                      ;} //Local temporary variable
  void addMoveWeight(float moveWeight)                    { moveWeight_.push_back(moveWeight)                    ;} //Local temporary variable
  void addMoveWell(void)                                  { wellMoveAngle_.push_back(moveAngle_);
                                                            wellMoveWeight_.push_back(moveWeight_)               ;}
  void clearMoveWell(void)                                { moveAngle_.clear();
                                                            moveWeight_.clear()                                  ;}
  void addIndicatorBGTrend(int indicator)                 { indBGTrend_.push_back(indicator)                     ;}
  void addIndicatorWavelet(int indicator)                 { indWavelet_.push_back(indicator)                     ;}
  void addIndicatorFacies(int indicator)                  { indFacies_.push_back(indicator)                      ;}
  void addIndicatorRockPhysics(int indicator)             { indRockPhysics_.push_back(indicator)                 ;}
  void addIndicatorRealVs(int indicator)                  { indRealVs_.push_back(indicator)                      ;}
  void addIndicatorFilter(int indicator)                  { indFilter_.push_back(indicator)                      ;}
  void setIndicatorFilter(int i ,int indicator)           { indFilter_[i]             = indicator                ;}
  void setLogName(int i, const std::string & logName)     { logNames_[i]              = NRLib::Uppercase(logName);}
  void addLogName(const std::string & log_name)           { logNames_.push_back(log_name)                        ;}
  void setInverseVelocity(int i, bool inverse)            { inverseVelocity_[i]       = inverse                  ;}
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
  void addVpVsRatioInterval(std::string interval_name, float ratio){ vpvs_ratio_interval_[interval_name] = ratio ;}
  void setErosionPriorityTopSurface(int priority)         { erosion_priority_top_surface_= priority              ;}
  void setErosionPriorityBaseSurface(const std::string & interval_name, int erosion_pri) { erosion_priority_interval_base_surface_[interval_name] = erosion_pri; }
  void setVpVsRatioFromWells(bool vp_vs_ratio_from_wells) { vp_vs_ratio_from_wells_   = vp_vs_ratio_from_wells   ;}
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
  void setGuardZone(float guard_zone)                     { guard_zone_               = guard_zone               ;}
  void setSmoothLength(float smooth_length)               { smooth_length_            = smooth_length            ;}
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
  void addSegyOffset(float segyOffset)                    { segyOffset_.push_back(segyOffset)                    ;}
  void addLocalSegyOffset(float segyOffset)               { localSegyOffset_.push_back(segyOffset)               ;}
  void setPundef(float p_undef)                           { p_undef_                  = p_undef                  ;}
  void setLzLimit(double lzLimit)                         { lzLimit_                  = lzLimit                  ;}
  void setTimeDTop(double time_dTop)                      { time_dTop_                = time_dTop                ;}
  void setTimeLz(double time_lz)                          { time_lz_                  = time_lz                  ;}
  void setTimeDz(double time_dz)                          { time_dz_                  = time_dz                  ;}
  void setTimeNz(int time_nz)                             { time_nz_                  = time_nz                  ;}
  void setTimeNzInterval(std::string name, int time_nz)   { time_nz_interval_[name]   = time_nz                  ;}
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
  void setMultizoneBackground(bool multizone)             { multizoneBackground_      = multizone                ;}
  void setGenerateBackground(bool generateBackgr)         { generateBackground_       = generateBackgr           ;}
  void setBackgroundFromRockPhysics(bool backgroundFromRP){ backgroundFromRockPhysics_= backgroundFromRP         ;}
  void setCalibrateRockPhysicsToWells(bool calibrate)     { calibrateRockPhysicsToWells_ = calibrate             ;}
  void setEstimateFaciesProb(bool estFaciesProb)          { estimateFaciesProb_       = estFaciesProb            ;}
  void setFaciesProbRelative(bool faciesProbRel)          { faciesProbRelative_       = faciesProbRel            ;}
  void setFaciesProbFromRockPhysics(bool rockPhysics)     { faciesProbFromRockPhysics_= rockPhysics              ;}
  void setNoVsFaciesProb(bool noVsFaciesProb)             { noVsFaciesProb_           = noVsFaciesProb           ;}
  void setUseFilterForFaciesProb(bool useFilterForProb)   { useFilterForProb_         = useFilterForProb         ;}
  void setFaciesLogGiven(bool faciesLogGiven)             { faciesLogGiven_           = faciesLogGiven           ;}
  void setPorosityLogGiven(bool porosityGiven)            { porosityLogGiven_         = porosityGiven            ;}
  void addPriorFaciesProb(std::string name, float value)  { priorFaciesProb_[name]    = value                    ;}
  void addVolumeFractionProb(std::string name, float value)  { volumeFractionProb_[name]    = value              ;}
  void addPriorFaciesProbInterval(std::string interval_name, std::map<std::string, float> prior_int_map){ priorFaciesProbInterval_[interval_name] = prior_int_map ;}
  void addVolumeFractionInterval(std::string interval_name, std::map<std::string, float> volumefractions_map) { volumefractionInterval_[interval_name] = volumefractions_map ;}
  void setPriorFaciesProbGiven(int fpg)                   { priorFaciesProbGiven_     = fpg                      ;}
  void setDepthDataOk(bool depthDataOk)                   { depthDataOk_              = depthDataOk              ;}
  void setParallelTimeSurfaces(bool pTimeSurfaces)        { parallelTimeSurfaces_     = pTimeSurfaces            ;}
  void setUseLocalWavelet(bool useLocalWavelet)           { useLocalWavelet_          = useLocalWavelet          ;}
  void setUseLocalNoise(bool useLocalNoise)               { useLocalNoise_            = useLocalNoise            ;}
  void setOptimizeWellLocation(bool optimizeWellLoc)      { optimizeWellLocation_     = optimizeWellLoc          ;}
  void setLogLevel(int logLevel)                          { logLevel_                 = logLevel                 ;}
  void setSeed(int seed)                                  { seed_                     = seed                     ;}
  void setDoSmoothKriging(bool smooth)                    { smoothKrigedParameters_   = smooth                   ;}
  void setRunFromPanel(bool panel)                        { runFromPanel_             = panel                    ;}
  void setNoWellNeeded(bool wellNeeded)                   { noWellNeeded_             = wellNeeded               ;}
  void setUseVerticalVariogram(bool useVerticalVariogram) { useVerticalVariogram_     = useVerticalVariogram     ;}
  void setTempCorrRange(float tempCorrRange)              { temporalCorrelationRange_  = tempCorrRange           ;}
  void addVintage(int year, int month, int day);
  void setNoSeismicNeeded(bool seismicNeeded)             { noSeismicNeeded_          = seismicNeeded            ;}
  void addTimeGradientSettings(float distance, float alpha);
  void setSeismicQualityGridRange(float range)            { seismicQualityGridRange_  = range                    ;}
  void setSeismicQualityGridValue(float value)            { seismicQualityGridValue_  = value                    ;}

  void addDefaultVintage(void);
  void addDefaultTimeGradientSettings(void);
  void addDefaultSegyOffset(void)                         { segyOffset_.push_back(0.0f)                          ;}
  void addDefaultAngularCorr(void)                        { angularCorr_.push_back(new GenExpVario(1, 10*static_cast<float>(NRLib::Pi/180.0)));} // Power=1 range=10deg
  void setDefaultUseLocalNoise(void)                      { useLocalNoise_ = false                               ;}

  double getDefaultCorrelationVpVs()                      { double corr = 1/std::sqrt(2.0f); return(corr)        ;}

  void setCorrDirTopConform(bool topConformCorrelation)   { topConformCorrelation_   = topConformCorrelation     ;}
  void setCorrDirBaseConform(bool baseConformCorrelation) { baseConformCorrelation_  = baseConformCorrelation    ;}
  void setCorrDirIntervalTopConform(const std::string & interval_name, bool intervalTopConformCorrelation)   { intervalTopConformCorrelation_[interval_name]  = intervalTopConformCorrelation  ;}
  void setCorrDirIntervalBaseConform(const std::string & interval_name, bool intervalBaseConformCorrelation) { intervalBaseConformCorrelation_[interval_name] = intervalBaseConformCorrelation ;}
  void setCorrDirIntervalUsed(bool intervalCorrelationUsed) { intervalCorrelationUsed_ = intervalCorrelationUsed ;}

  void addIntervalName(std::string name)                  { interval_names_.push_back(name)                      ;}
  void setIntervalNames(const std::vector<std::string> & interval_names) {interval_names_ = interval_names       ;}
  void setErosionPriorityIntervals(const std::string & interval_name, const int priority) { erosion_priority_interval_base_surface_[interval_name] = priority;}

  void setRMSStandardDeviation(double value)              { RMSStandardDeviation_ = value                        ;}
  void setRMSPriorGiven(bool given)                       { RMSPriorGiven_   = given                             ;}
  void setRMSnLayersAbove(int n_layers)                   { RMSnLayersAbove_ = n_layers                          ;}
  void setRMSnLayersBelow(int n_layers)                   { RMSnLayersBelow_ = n_layers                          ;}
  void setRMSMeanVpTop(double mean)                       { RMSMeanVpTop_ = mean                                 ;}
  void setRMSMeanVpBase(double mean)                      { RMSMeanVpBase_ = mean                                ;}
  void setRMSVarianceVpAbove(double var)                  { RMSVarianceVpAbove_ = var                            ;}
  void setRMSVarianceVpBelow(double var)                  { RMSVarianceVpBelow_ = var                            ;}
  void setRMSTemporalCorrelationRangeAbove(double range)  { RMSTemporalCorrRangeAbove_ = range                   ;}
  void setRMSTemporalCorrelationRangeBelow(double range)  { RMSTemporalCorrRangeBelow_ = range                   ;}

  void clearTimeLapse(void)                               { angle_.clear();
                                                            localTHF_.clear();
                                                            localSegyOffset_.clear();
                                                            seismicType_.clear();
                                                            estimateGlobalWaveletScale_.clear();
                                                            estimateLocalShift_.clear();
                                                            estimateLocalScale_.clear();
                                                            estimateLocalNoise_.clear();
                                                            waveletScale_.clear();
                                                            useRickerWavelet_.clear();
                                                            rickerPeakFrequency_.clear();
                                                            estimateWavelet_.clear();
                                                            stretchFactor_.clear();
                                                            matchEnergies_.clear();
                                                            estimateSNRatio_.clear();
                                                            SNRatio_.clear()                                     ;}

  void addTimeLapse(void)                                 { timeLapseAngle_.push_back(angle_);
                                                            timeLapseLocalTHF_.push_back(localTHF_);
                                                            timeLapseLocalSegyOffset_.push_back(localSegyOffset_);
                                                            timeLapseSeismicType_.push_back(seismicType_);
                                                            timeLapseEstimateGlobalWaveletScale_.push_back(estimateGlobalWaveletScale_);
                                                            timeLapseEstimateLocalShift_.push_back(estimateLocalShift_);
                                                            timeLapseEstimateLocalScale_.push_back(estimateLocalScale_);
                                                            timeLapseEstimateLocalNoise_.push_back(estimateLocalNoise_);
                                                            timeLapseWaveletScale_.push_back(waveletScale_);
                                                            timeLapseUseRickerWavelet_.push_back(useRickerWavelet_);
                                                            timeLapseRickerPeakFrequency_.push_back(rickerPeakFrequency_);
                                                            timeLapseEstimateWavelet_.push_back(estimateWavelet_);
                                                            timeLapseStretchFactor_.push_back(stretchFactor_);
                                                            timeLapseMatchEnergies_.push_back(matchEnergies_);
                                                            timeLapseEstimateSNRatio_.push_back(estimateSNRatio_);
                                                            timeLapseSNRatio_.push_back(SNRatio_);
                                                            timeLapseUseLocalNoise_.push_back(useLocalNoise_);}

  void addTimeLapseGravimetry(bool gravimetry)            { timeLapseGravimetry_.push_back(gravimetry)           ;}
  void addTimeLapseTravelTimeGiven(bool travelTime)       { timeLapseTravelTimeGiven_.push_back(travelTime)      ;}

  void addTravelTimeHorizonName(std::string name)         { travelTimeHorizonName_.push_back(name)               ;}
  void addTravelTimeHorizonSD(double standard_deviation)  { travelTimeHorizonSD_.push_back(standard_deviation)   ;}

  void clearTimeLapseTravelTime()                         { travelTimeHorizonName_.clear();
                                                            travelTimeHorizonSD_.clear()                         ;}

  void addTimeLapseTravelTime()                           { timeLapseTravelTimeHorizon_.push_back(travelTimeHorizonName_);
                                                            timeLapseTravelTimeHorizonSD_.push_back(travelTimeHorizonSD_) ;}

  void setSnapGridToSeismicData(bool snapToSeismicData)   { snapGridToSeismicData_    = snapToSeismicData        ;}
  void setWavelet3DTuningFactor(double tuningFactor)      { wavelet3DTuningFactor_    = tuningFactor             ;}
  void setGradientSmoothingRange(double smoothingRange)   { gradientSmoothingRange_   = smoothingRange           ;}
  void setEstimateWellGradientFromSeismic(bool estimate)  { wellGradientFromSeismic_  = estimate                 ;}

  enum          priorFacies{FACIES_FROM_WELLS,
                            FACIES_FROM_MODEL_FILE,
                            FACIES_FROM_CUBES};

  enum          sseismicTypes{STANDARDSEIS = 0, PSSEIS = 1};

  enum          indicators{NO, YES, NOTSET};

  enum          areaSpecification{AREA_FROM_GRID_DATA,
                                  AREA_FROM_GRID_DATA_AND_UTM,     ///< Force area to be aligned with seismic data
                                  AREA_FROM_GRID_DATA_AND_SURFACE, ///< Force area to be aligned with seismic data
                                  AREA_FROM_UTM,
                                  AREA_FROM_SURFACE};

  enum          correlationStructure{TOP,
                                     BASE,
                                     COMPACTION};

  enum          rockConstituents{FLUID,
                                 SOLID,
                                 DRY_ROCK,
                                 ROCK};

  enum          trendCube{TWT,
                          STRATIGRAPHIC_DEPTH,
                          CUBE_FROM_FILE};

private:

  std::vector<Vario*>               angularCorr_;                ///< Variogram for lateral error correlation, time lapse
  Vario                           * lateralCorr_;                ///< Variogram for lateral parameter correlation
  Vario                           * backgroundVario_;            ///< Used for lateral background correlation.
  Vario                           * localWaveletVario_;          ///< Used for local wavelet (gain and shift) and local noise.

  SegyGeometry                    * geometry_full_;              ///< area parameters of full seismic data
  SegyGeometry                    * geometry_;                   // area parameters
  std::vector<float>                segyOffset_;                 // Starttime for SegY cubes, time lapse
  std::vector<float>                localSegyOffset_;            // Starttime for SegY cubes per angle.
  TraceHeaderFormat               * traceHeaderFormat_;          // traceheader of input
  std::vector<TraceHeaderFormat*>   localTHF_;                   // traceheader per angle
  TraceHeaderFormat               * traceHeaderFormatOutput_;    // traceheader for output files
  int                               krigingParameter_;

  std::vector<std::string>          travelTimeHorizonName_;      // Name of travel time horizon
  std::vector<double>               travelTimeHorizonSD_;        // Standard deviation of the travel time horizon
  double                            RMSStandardDeviation_;       // Standard deviation for the RMS data
  bool                              RMSPriorGiven_;              // True if prior information is given for RMS inversion
  int                               RMSnLayersAbove_;            // n layers above the reservoir in inversion of RMS velocities
  int                               RMSnLayersBelow_;            // n layers below the reservoir in inversion of RMS velocities
  double                            RMSMeanVpTop_;               // E(Vp) at the top of the zone above the reservoir, that is, at sea level
  double                            RMSMeanVpBase_;              // E(Vp) at the base of the zone below the reservoir
  double                            RMSVarianceVpAbove_;         // Var(Vp) above the reservoir in inversion of RMS velocities
  double                            RMSVarianceVpBelow_;         // Var(Vp) below the reservoir in inversion of RMS velocities
  double                            RMSTemporalCorrRangeAbove_;  // Temporal corralation range (ms) for exponential variogram above the reservoir in inversion of RMS velocities
  double                            RMSTemporalCorrRangeBelow_;  // Temporal corralation range (ms) for exponential variogram below the reservoir in inversion of RMS velocities

  std::vector<int>                  seismicType_;                ///< PP- or PS- seismic
  std::vector<float>                angle_;                      ///< Angles
  std::vector<float>                waveletScale_;               ///< Signal-to-noise ratio
  std::vector<float>                SNRatio_;                    ///< Signal-to-noise ratio

  std::vector<float>                moveAngle_;                  ///< Angles for moving wells, local temporary variable
  std::vector<float>                moveWeight_;                 ///< Angle weights for moving wells, local temporary variable
  std::vector<std::vector<float> >  wellMoveAngle_;              ///< moveAngle_ collected for all wells
  std::vector<std::vector<float> >  wellMoveWeight_;             ///< moveWeight_ collected for all wells

  std::vector<bool>                 matchEnergies_;              // Let dataVariance_ = signalVariance_
  std::vector<bool>                 estimateWavelet_;            //
  std::vector<bool>                 estimateSNRatio_;            //
  std::vector<bool>                 estimateLocalShift_;         // Estimate local wavelet shift
  std::vector<bool>                 estimateLocalScale_;         // Estimate local wavelet scale
  std::vector<bool>                 estimateLocalNoise_;         // Estimate local noise
  std::vector<bool>                 estimateGlobalWaveletScale_;
  std::vector<bool>                 useRickerWavelet_;
  std::vector<bool>                 timeLapseUseLocalNoise_;
  std::vector<bool>                 timeLapseGravimetry_;
  std::vector<bool>                 timeLapseTravelTimeGiven_;

  std::vector<std::vector<bool> >   timeLapseEstimateLocalShift_;// Estimate local wavelet shift
  std::vector<std::vector<bool> >   timeLapseEstimateLocalScale_;// Estimate local wavelet scale
  std::vector<std::vector<bool> >   timeLapseEstimateLocalNoise_;// Estimate local noise
  std::vector<std::vector<bool> >   timeLapseUseRickerWavelet_;
  std::vector<std::vector<bool> >   timeLapseEstimateWavelet_;
  std::vector<std::vector<bool> >   timeLapseMatchEnergies_;     // Let dataVariance_ = signalVariance_
  std::vector<std::vector<bool> >   timeLapseEstimateGlobalWaveletScale_;
  std::vector<std::vector<bool> >   timeLapseEstimateSNRatio_;
  std::vector<std::vector<int> >    timeLapseSeismicType_;       // PP- or PS- seismic
  std::vector<std::vector<float> >  timeLapseWaveletScale_;      // Signal-to-noise ratio
  std::vector<std::vector<float> >  timeLapseRickerPeakFrequency_;
  std::vector<std::vector<float> >  timeLapseStretchFactor_;     // Stretch factor for pulse in 3D-wavelet
  std::vector<std::vector<float> >  timeLapseSNRatio_;           // Signal-to-noise ratio
  std::vector<std::vector<float> >  timeLapseAngle_;             // Angles
  std::vector<std::vector<float> >  timeLapseLocalSegyOffset_;   // Starttime for SegY cubes per angle.

  std::vector<std::vector<std::string> > timeLapseTravelTimeHorizon_; // Travel time horizon names in different time lapses
  std::vector<std::vector<double> >      timeLapseTravelTimeHorizonSD_; // Standard deviation of travel time horizons for different time lapses

  std::vector<std::vector<TraceHeaderFormat*> > timeLapseLocalTHF_;

  std::vector<float>                rickerPeakFrequency_;

  int                               erosion_priority_top_surface_;// Erosion priority for the top surface of the inversion intervals for multiple intervals. 1 by default
  std::map<std::string, int>        erosion_priority_interval_base_surface_; ///< Erosion priority for the base surfaces of each interval. Each one must be unique.
  std::vector<std::string>          interval_names_;              // Interval names for multiple interval inversion
  std::map<std::string, float>      vpvs_ratio_interval_;         // Interval names and the Vp/Vs-ratio given in <vp-vs-ratio> under <advanced-settings>

  std::vector<int>                  waveletDim_;                 ///< Holds if 1D-wavelet (=0) or 3D-wavelet (=1)
  std::vector<float>                stretchFactor_;              ///< Stretch factor for pulse in 3D-wavelet
  std::vector<float>                estRangeX_;                  ///< Estimation range in x-direction for 3D-wavelet
  std::vector<float>                estRangeY_;                  ///< Estimation range in y-direction for 3D-wavelet

  bool                              do4DInversion_;              ///< True if CRAVA is to run a 4D inversion
  bool                              do4DRockPhysicsInversion_;   ///< True if we should do rockpysics inversion, only active for 4D inversion
  bool                              backgroundFromRockPhysics_;  ///< True if background is to be generated from rock physics. Is this relevant? Or same as faciesProbFromRockPhysics_?
  bool                              estimateBackground_;         ///< In estimation mode, skip estimation of background if false
  bool                              estimateCorrelations_;       ///< As above, but correlations.
  bool                              estimateWaveletNoise_;       ///< As above, but for wavelet and noise parameters.
  bool                              estimate3DWavelet_;          ///< True if a 3D wavelet is estimated for at least one angle.
  bool                              hasTime3DMapping_;           ///< True if command time-3D-mapping is used
  bool                              use3DWavelet_;               // True if 3D wavelet is used
  bool                              calibrateRockPhysicsToWells_;///< True if rock physics models can be estimated from wells

  std::vector<float>                constBackValue_;             ///< Values set for constant background model
                                                                 ///< Negative value ==> read from file (actual value gives format).
  bool                              useAIBackground_;            ///< Read in file for AI background instead of Vp background
  bool                              useSIBackground_;            ///< Read in file for SI background instead of Vs background
  bool                              useVpVsBackground_;          ///< Read in file for VpVs background instead of Vs background
  std::string                       backgroundType_;             ///< background or earth model

  // The following indicators use the indicators enum above. (2 = yes, but may override in QC, 1=yes, 0=no)
  std::vector<int>                  indBGTrend_;                 ///< Use well to estimate background trend?
  std::vector<int>                  indWavelet_;                 ///< Use well to estimate wavelet?
  std::vector<int>                  indFacies_;                  ///< Use well to estimate facies?
  std::vector<int>                  indRockPhysics_;             ///< Use well to estimate rock physics?
  std::vector<int>                  indRealVs_;                  ///< Treat Vs log as real?
  std::vector<int>                  indFilter_;                  ///< Filter elastic logs using spatial multi-parameter filter?

  std::vector<std::string>          logNames_;                   ///< The keywords to look for for time, sonic, shear sonic and density
  std::vector<bool>                 inverseVelocity_;            ///< If element 0 is true, vp comes from dt, if 1 is true, vs comes from dts in well.

  int                               priorFaciesProbGiven_;
  std::map<std::string, float>      priorFaciesProb_;
  std::map<std::string, float>      volumeFractionProb_;
  std::map<std::string, std::map<std::string, float> >      priorFaciesProbInterval_;
  std::map<std::string, std::map<std::string, float> >      volumefractionInterval_;

  int                               nWells_;
  int                               nSimulations_;

  float                             alpha_min_;                  ///< Vp - smallest allowed value
  float                             alpha_max_;                  ///< Vp - largest allowed value
  float                             beta_min_;                   ///< Vs - smallest allowed value
  float                             beta_max_;                   ///< Vs - largest allowed value
  float                             rho_min_;                    ///< Rho - smallest allowed value
  float                             rho_max_;                    ///< Rho - largest allowed value

  float                             var_alpha_min_;              ///<| These min and max values are used for consistency check. If
  float                             var_alpha_max_;              ///<| variances are outside these ranges there is probably a
  float                             var_beta_min_;               ///<| problem with the logs.
  float                             var_beta_max_;               ///<|
  float                             var_rho_min_;                ///<| The limits are for point variances. The minimum allowed variance
  float                             var_rho_max_;                ///<| for parameters will be scaled with 1/dt*dt

  float                             vp_vs_ratio_min_;            ///< Smallest Vp/Vs-ratio regarded as likely
  float                             vp_vs_ratio_max_;            ///< Largest Vp/Vs-ratio regarded as likely
  float                             vp_vs_ratio_;                ///< Vp/Vs-ratio from input
  bool                              vp_vs_ratio_from_wells_;     ///< Estimate Vp/Vs-ratio from well data

  float                             ref_depth_;                  ///< z0 - reference depth for target area
  float                             average_velocity_;           ///< v0 - average velocity in target area

  float                             temporalCorrelationRange_;   ///< Temporal correlation range for exponential variogram

  float                             maxHz_background_;           ///< Background resolution (high cut frequency)
  float                             maxHz_seismic_;              ///< Seismic resolution (high cut frequency)

  float                             maxRankCorr_;                ///< Vp-Vs correlation threshold for regarding Vs log synthetic
  float                             maxMergeDist_;               ///< log entries closer than this will be merged
  float                             maxDevAngle_;                ///< Wells with a local deviation larger than this is treated as deviated

  float                             lowCut_;                     ///< lower limit for frequency to be inverted
  float                             highCut_;                    ///< upper limit for frecuency to be inverted

  float                             wnc_;                        ///< White noise component, see crava.h

  float                             energyThreshold_;            ///< If energy in reflection trace divided by mean energy
                                                                 ///< in reflection traces is lower than this, the reflections
                                                                 ///< will be interpolated. Default 0.
  float                             maxWellOffset_;              ///< Maximum offset for moving of wells
  float                             maxWellShift_;               ///< Maximum vertical shift for moving of wells

  float                             defaultWaveletLength_;       ///< Assumed length of a wavelet
  float                             guard_zone_;                 ///< Band outside target interval (on each side) where data is required
  float                             smooth_length_;              ///< How much of data in guard zone to smooth (to end up with zero)

  float                             minRelWaveletAmp_;           ///< Minimum relative wavelet amplitude. Smaller amplitudes are disregarded.
  float                             maxWaveletShift_;            ///< Largest allowed shift when estimating wavelet
  float                             waveletTaperingL_;           ///< Til Odds waveletestimering

  float                             minSamplingDensity_;         ///< Threshold value for minimum sampling density in dz
  float                             minHorizontalRes_;           ///< Threshold value for minimum horizontal resolution in dx and dy

  double                            xPadFac_;                    ///< Padding factor/fraction in x direction
  double                            yPadFac_;                    ///< Padding factor/fraction in y direction
  double                            zPadFac_;                    ///< Padding factor/fraction in z direction

  int                               nxPad_;                      ///< Number of cells to pad in x direction
  int                               nyPad_;
  int                               nzPad_;

  bool                              estimateXYPadding_;          ///< Estimate the z-padding from ranges
  bool                              estimateZPadding_;           ///< Estimate the z-padding from wavelet length

  float                             p_undef_;                    ///< Level for undefined facies

  double                            lzLimit_;                    ///< Minimum allowed value for (min interval thickness)/(max interval thickness)
  double                            time_dTop_;                  ///< Used when top and base surfaces are parallel
  double                            time_lz_;                    ///< Used when top and base surfaces are parallel
  double                            time_dz_;                    ///< Used when top and base surfaces are parallel
  std::map<std::string, int>        time_nz_interval_;           ///< Number of layers for each interval
  int                               time_nz_;                    ///< Used when top and base surfaces are parallel
  bool                              velocityFromInv_;            ///< Velocity for time depth from inverted Vp.

  int                               areaSpecification_;          ///< Specifying whether are is taken from UTM-coord, seismic or surface
  std::vector<int>                  areaILXL_;                   ///< Vector with 6 elements (if used), in this order:
                                                                 ///< [0] = IL start
                                                                 ///< [1] = IL end
                                                                 ///< [2] = XL start
                                                                 ///< [3] = XL end
                                                                 ///< [4] = IL step
                                                                 ///< [5] = XL step

  bool                              writePrediction_;            ///< Determines whether prediction is written.
  int                               outputGridsElastic_;         ///< Decides which elastic grids to be written to file.
  int                               outputGridsOther_;           ///< Decides other grid output to be written to file
  int                               outputGridsSeismic_;         ///< Decides seismic grid output to be written to file.
  int                               domainFlag_;                 ///< Decides writing in time and/or depth.
  int                               formatFlag_;                 ///< Decides output format, see above.
  int                               wellFlag_;                   ///< Decides well output.
  int                               wellFormatFlag_;             ///< Decides well output format.
  int                               waveletFlag_;                ///< Decides wavelet output
  int                               waveletFormatFlag_;          ///< Decides wavelet output format
  int                               otherFlag_;                  ///< Decides output beyond grids and wells.
  bool                              fileGrid_;                   ///< Indicator telling if grids are to be kept on file
  bool                              outputGridsDefault_;         ///< Indicator telling if grid output has been actively controlled
  bool                              waveletFormatManual_;        ///< True if wavelet format is decided in the model file
  bool                              useVerticalVariogram_;       ///< True if a vertical variogram is used to estimate temporal correlation
  bool                              forwardModeling_;            ///< Forward modelling
  bool                              estimationMode_;             ///< Estimation
  bool                              generateSeismicAfterInv_;    ///< Synthetic seismic from inversion result
  bool                              generateBackground_;         ///< Make background model
  bool                              multizoneBackground_;        ///< Make multizone background model
  bool                              estimateFaciesProb_;         ///< Shall facies probabilites be estimated?
  bool                              faciesProbRelative_;         ///< Use relative elastic parameters for facies prob estimation?
  bool                              faciesProbFromRockPhysics_;  ///< Calculate facies probabilities using rock physics models
  bool                              noVsFaciesProb_;             ///< Do not use Vs for faciesprob.
  bool                              useFilterForProb_;           ///< Use filtered logs for facies probs, otherwise, use sampled inversion.
  bool                              faciesLogGiven_;
  bool                              porosityLogGiven_;
  bool                              depthDataOk_;                ///< We have what we need to do depth conversion
  bool                              parallelTimeSurfaces_;
  bool                              useLocalWavelet_;            ///< Wavelets are multiplied with gain and shift maps
  bool                              useLocalNoise_;              ///< Signal-to-noise is multiplied with gain and shift maps
  bool                              optimizeWellLocation_;       ///< True if at least one well is to be moved
  bool                              smoothKrigedParameters_;     ///< True if we should smooth borders between kriging blocks
  bool                              runFromPanel_;               ///< True if run is started from RMS panel. Relaxed checking.
  bool                              noWellNeeded_;               ///< True for some configurations of input data
  bool                              noSeismicNeeded_;            ///< True for some estimation settings
  bool                              snapGridToSeismicData_;      ///< Force inversion area to align with seismic data
  std::vector<float>                distanceFromWell_;           ///< Minimum distance for where gradients should not cross, time lapse
  std::vector<float>                sigma_m_;                    ///< Smoothness level of the gradients, time lapse
  std::vector<int>                  vintageDay_;                 ///< Day of month the seismic time lapse data were collected
  std::vector<int>                  vintageMonth_;               ///< Month the seismic time lapse data were collected
  std::vector<int>                  vintageYear_;                ///< Year the seismic time lapse data were collected
  double                            wavelet3DTuningFactor_;      ///< Large value forces better fit of wavelet
  double                            gradientSmoothingRange_;     ///< Controls smoothing of gradient used in 3D wavelet estimate/inversion
  bool                              wellGradientFromSeismic_;    ///< Estimate well gradient used for 3D wavelet estimation from seismic?
  float                             seismicQualityGridRange_;    ///< Radius value from well-points where wells are used in Seismic Quality Grids
  float                             seismicQualityGridValue_;    ///< Value between wells if range is used.

  bool                              topConformCorrelation_;      ///< Should top correlation direction be equal to the top inversion surface
  bool                              baseConformCorrelation_;     ///< Should base correlation direction be equal to the base inversion surface
  std::map<std::string, bool>       intervalTopConformCorrelation_;  ///< Should base correlation direction be equal to the base inversion surface per interval
  std::map<std::string, bool>       intervalBaseConformCorrelation_; ///< Should base correlation direction be equal to the base inversion surface per interval
  bool                              intervalCorrelationUsed_;    ///< Whether intervals are used for correlation direction

  std::vector<int>                  erosionPriority_;            // Erosion priority of the different layers in the multizone background model
  std::vector<int>                  correlationStructure_;       // Correlation structure for the different layers in the multizone background model
  std::vector<double>               surfaceUncertainty_;         // Uncertainty for the horizons in the multizone backround model

  std::vector<std::string>          trendCubeParameter_;          // Name of the trend parameters in the rock physics model
  std::vector<int>                  trendCubeType_;               // Type of the trend cube

  std::map<std::string, std::vector<DistributionWithTrendStorage *> > reservoirVariable_;  // Rock physics variables defined in reservoir, the vector goes over the vintages of the variable
  std::map<std::string, DistributionsRockStorage *>                   rockStorage_;        // Rock physics rocks defined in predefinitions
  std::map<std::string, DistributionsDryRockStorage *>                dryRockStorage_;     // Rock physics dry rocks defined in predefinitions
  std::map<std::string, DistributionsSolidStorage *>                  solidStorage_;       // Rock physics solids defined in predefinitions
  std::map<std::string, DistributionsFluidStorage *>                  fluidStorage_;       // Rock physics fluids defined in predefinitions

  int                               logLevel_;

  int                               seed_;                       ///< Random seed.

  static int                        debugFlag_;
};

#endif
