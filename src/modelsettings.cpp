/***************************************************************************
*      Copyright (C) 2008 by Norwegian Computing Center and Statoil        *
***************************************************************************/

#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
#include <string.h>
#include <map>
#include "nrlib/iotools/logkit.hpp"
#include "nrlib/segy/segy.hpp"
#include "nrlib/segy/traceheader.hpp"

#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionsrockstorage.h"
#include "rplib/distributionsdryrockstorage.h"
#include "rplib/distributionssolidstorage.h"
#include "rplib/distributionsfluidstorage.h"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/inputfiles.h"


ModelSettings::ModelSettings(void)
  : localSegyOffset_(0),
    localTHF_(0),
    seismicType_(0),
    angle_(0),
    waveletScale_(0),
    SNRatio_(0),
    matchEnergies_(0),
    estimateWavelet_(0),
    estimateSNRatio_(0),
    estimateLocalNoise_(0),
    estimateGlobalWaveletScale_(0),
    multiple_intervals_(false),
    constBackValue_(3),
    indBGTrend_(0),
    indWavelet_(0),
    indFacies_(0),
    indRockPhysics_(0),
    indFilter_(0),
    logNames_(6),
    positionLogNames_(2),
    inverseVelocity_(2),
    min_blocks_with_data_for_corr_estim_(100)
{
  lateralCorr_                = new GenExpVario(1, 1000, 1000);
  backgroundVario_            = new GenExpVario(1, 2000, 2000);
  localWaveletVario_          =     NULL; // Will be set equal to backgroundVario unless it is set separately
  geometry_full_              =     NULL;
  geometry_                   =     NULL;
  traceHeaderFormat_          =     NULL;
  traceHeaderFormatOutput_    = new TraceHeaderFormat(TraceHeaderFormat::SEISWORKS);
  traceHeaderFormatBackground_.resize(3, NULL);
  krigingParameter_           =        0; // Indicate kriging not set.
  nWells_                     =        0;
  nSimulations_               =        0;
  backgroundType_             =       "";

  //
  // The original ranges were provided by Nam Hoai Pham (Statoil/25.09.2007)
  //
  vp_min_                  =  1300.0f;   // Nam: 1800
  vp_max_                  =  7000.0f;   // Nam: 6000
  vs_min_                  =   200.0f;   // Nam:  800
  vs_max_                  =  4200.0f;   // Nam: 3000
  rho_min_                 =     1.4f;   // Nam:  1.5
  rho_max_                 =     3.3f;   // Nam:  3.0

  var_vp_min_              =   5.e-4f;
  var_vp_max_              = 250.e-4f;
  var_vs_min_              =  10.e-4f;
  var_vs_max_              = 500.e-4f;
  var_rho_min_             =   2.e-4f;
  var_rho_max_             = 100.e-4f;

  //vp_vs_ratio_             = RMISSING;
  vp_vs_ratio_from_wells_  =    false;
  vp_vs_ratio_min_         =     1.4f;
  vp_vs_ratio_max_         =     3.0f;

  maxHz_background_        =     6.0f;
  maxHz_seismic_           =    40.0f;

  filter_multizone_background_ = true;

  maxRankCorr_             =    0.99f;
  maxMergeDist_            =    0.01f;   // 0.01ms (approx. 2-3cm)
  maxDevAngle_             =    15.0f;

  lowCut_                  =     5.0f;
  highCut_                 =    55.0f;

  wnc_                     =     0.1f;

  energyThreshold_         =     0.0f;

  maxWellShift_            =    11.0f;
  maxWellOffset_           =   250.0f;

  defaultWaveletLength_    =   200.0f;
  guard_zone_              =  defaultWaveletLength_*0.5f;
  smooth_length_           =  defaultWaveletLength_*0.5f;

  minRelWaveletAmp_        =    0.04f;
  maxWaveletShift_         =    11.0f;
  waveletTaperingL_        =   200.0f;

  wavelet3DTuningFactor_   =     50.0; // double
  gradientSmoothingRange_  =    100.0; // double

  minSamplingDensity_      =     0.5f;
  minHorizontalRes_        =     5.0f;

  xPadFac_                 =      0.0;
  yPadFac_                 =      0.0;
  zPadFac_                 =      0.0;

  //nxPad_                   = IMISSING;
  //nyPad_                   = IMISSING;
  //nzPad_                   = IMISSING;

  estimateXYPadding_       =     true;
  estimateZPadding_        =     true;

  p_undef_                 =    0.01f;

  seismicQualityGridValue_ = RMISSING;
  seismicQualityGridRange_ = RMISSING;

  lzLimit_                 =     0.41;   // NB! This is a double ==> do not use 'f'.
  time_dTop_               = RMISSING;
  time_lz_                 = RMISSING;
  time_dz_                 = RMISSING;
  //time_nz_                 = IMISSING;
  //time_nz_[""]             = IMISSING;
  segy_nz_                 = IMISSING;
  segy_dz_                 = RMISSING;
  output_offset_           = RMISSING;
  match_output_input_segy_ = true;
  use_input_segy_dz_for_output_segy_ = false;

  velocityFromInv_         =    false;

  areaILXL_                = std::vector<int>(0);

  writePrediction_         =    false;  //Will be set to true if no simulations.
  outputGridsElastic_      = IO::VP + IO::VS + IO::RHO;  // Default output
  outputGridsOther_        =        0;
  outputGridsSeismic_      =        0;
  outputGridsDefault_      =     true;
  formatFlag_              = IO::STORM;
  domainFlag_              = IO::TIMEDOMAIN;
  wellFlag_                =        0;
  wellFormatFlag_          = IO::RMSWELL;
  waveletFlag_             =        0;
  waveletFormatFlag_       = IO::JASONWAVELET;
  otherFlag_               =        0;
  debugFlag_               =        0;
  fileGrid_                =    false;
  waveletFormatManual_     =    false;
  useVerticalVariogram_    =    false;
  do4DInversion_           =    false;
  do4DRockPhysicsInversion_=    false;
  backgroundFromRockPhysics_=   false;
  calibrateRockPhysicsToWells_= false;
  estimationMode_          =    false;
  forwardModeling_         =    false;
  generateBackground_      =     true;
  multizoneBackground_     =    false;
  useAIBackground_         =    false;
  useSIBackground_         =    false;
  useVpVsBackground_       =    false;
  estimateFaciesProb_      =    false;
  faciesProbRelative_      =     true;
  faciesProbFromRockPhysics_=   false;
  noVsFaciesProb_          =    false;
  useFilterForProb_        =     true;
  faciesLogGiven_          =    false;
  porosityLogGiven_        =    false;
  depthDataOk_             =    false;
  parallelTimeSurfaces_    =    false;
  useLocalWavelet_         =    false;
  optimizeWellLocation_    =    false;
  runFromPanel_            =    false;
  noWellNeeded_            =    false;
  noSeismicNeeded_         =    false;
  snapGridToSeismicData_   =    false;
  wellGradientFromSeismic_ =    false;
  writeAsciiSurfaces_      =    false;

  priorFaciesProbGiven_    = ModelSettings::FACIES_FROM_WELLS;

  generateSeismicAfterInv_ =    false;
  estimateBackground_      =     true;
  estimateCorrelations_    =     true;
  estimateWaveletNoise_    =     true;
  estimate3DWavelet_       =    false;
  hasTime3DMapping_        =    false;
  use3DWavelet_            =    false;

  logLevel_                = LogKit::L_Low;
  smoothKrigedParameters_  =    false;

  RMSPriorGiven_           =    false;

  seed_                    =        0;
  number_of_threads_       =        0;

  erosion_priority_top_surface_ = 1;

  //intervalTopConformCorrelation_[""]  = false;
  //intervalBaseConformCorrelation_[""] = false;

}

ModelSettings::~ModelSettings(void)
{
  for(size_t i = 0; i<angularCorr_.size(); i++){
    delete angularCorr_[i];
    angularCorr_[i] = NULL;
  }

  if (lateralCorr_ != NULL){
    delete lateralCorr_;
    lateralCorr_ = NULL;
  }

  if (backgroundVario_ != NULL){
    delete backgroundVario_;
    backgroundVario_ = NULL;
  }

  if (localWaveletVario_ != NULL){
    delete localWaveletVario_;
    localWaveletVario_ = NULL;
  }

  for (size_t i = 0; i < timeLapseTravelTimeLateralCorrelation_.size(); i++) {
    if (timeLapseTravelTimeLateralCorrelation_[i] != NULL)
      delete timeLapseTravelTimeLateralCorrelation_[i];
  }

  if(geometry_ != NULL) {
    delete geometry_;
    geometry_ = NULL;
  }

  if(geometry_full_ != NULL){
    delete geometry_full_;
    geometry_full_ = NULL;
  }


  delete traceHeaderFormat_;
  traceHeaderFormat_ = NULL;

  if(traceHeaderFormatOutput_ != NULL){
    delete traceHeaderFormatOutput_;
    traceHeaderFormatOutput_ = NULL;
  }

  for(size_t i=0; i<timeLapseLocalTHF_.size(); i++) {
    for(size_t j=0; j<timeLapseLocalTHF_[i].size(); j++){
      delete timeLapseLocalTHF_[i][j];
      timeLapseLocalTHF_[i][j] = NULL;
    }
  }

  for (size_t i = 0; i < traceHeaderFormatBackground_.size(); i++) {
    if (traceHeaderFormatBackground_[i] != NULL) {
      delete traceHeaderFormatBackground_[i];
      traceHeaderFormatBackground_[i] = NULL;
    }
  }

  for(std::map<std::string, DistributionsRockStorage *>::iterator it = rockStorage_.begin(); it != rockStorage_.end(); it++) {
    DistributionsRockStorage * storage = it->second;
    delete storage;
  }
  for(std::map<std::string, DistributionsDryRockStorage *>::iterator it = dryRockStorage_.begin(); it != dryRockStorage_.end(); it++) {
    DistributionsDryRockStorage * storage = it->second;
    delete storage;
  }
  for(std::map<std::string, DistributionsSolidStorage *>::iterator it = solidStorage_.begin(); it != solidStorage_.end(); it++) {
    DistributionsSolidStorage * storage = it->second;
    delete storage;
  }
  for(std::map<std::string, DistributionsFluidStorage *>::iterator it = fluidStorage_.begin(); it != fluidStorage_.end(); it++) {
    DistributionsFluidStorage * storage = it->second;
    delete storage;
  }
  for(std::map<std::string, std::vector<DistributionWithTrendStorage *> >::iterator it = reservoirVariable_.begin(); it != reservoirVariable_.end(); it++) {
    std::vector<DistributionWithTrendStorage *> vintageStorage = it->second;
    for(size_t i=0; i<vintageStorage.size(); i++) {
      delete vintageStorage[i];
      vintageStorage[i] = NULL;
    }
  }

}

bool
ModelSettings::getDoInversion(void) const
{
  int elasticFlag  = 0;
  int otherFlag    = 0;
  int seismicFlag  = 0;
  int blockedWells = 0;

  elasticFlag  += IO::VP
               +  IO::VS
               +  IO::RHO
               +  IO::LAMELAMBDA
               +  IO::LAMEMU
               +  IO::POISSONRATIO
               +  IO::AI
               +  IO::SI
               +  IO::VPVSRATIO
               +  IO::MURHO
               +  IO::LAMBDARHO;

  otherFlag    += IO::FACIESPROB
               +  IO::FACIESPROB_WITH_UNDEF
               +  IO::CORRELATION;

  seismicFlag  += IO::SYNTHETIC_SEISMIC_DATA
               +  IO::RESIDUAL;

  blockedWells += IO::BLOCKED_WELLS;

  if (((elasticFlag  & outputGridsElastic_) > 0  ||
       (otherFlag    & outputGridsOther_  ) > 0  ||
       (seismicFlag  & outputGridsSeismic_) > 0  ||
       (blockedWells & wellFlag_          ) > 0) &&
       estimationMode_ == false)
    return true;
  else
    return false;
}

bool
ModelSettings::getDoDepthConversion(void) const
{
  return(depthDataOk_ & ((domainFlag_ & IO::DEPTHDOMAIN) > 0));
}
void
ModelSettings::rotateVariograms(float angle)
{
  lateralCorr_->rotateCounterClockwise(-angle);
  backgroundVario_->rotateCounterClockwise(-angle);
  localWaveletVario_->rotateCounterClockwise(-angle);
}

void
ModelSettings::setLastAngularCorr(Vario * vario)
{
  size_t i = angularCorr_.size()-1;
  if (angularCorr_[i] != NULL)
    delete angularCorr_[i];
  angularCorr_[i] = vario;
}

void
ModelSettings::setLateralCorr(Vario * vario)
{
  if (lateralCorr_ != NULL)
    delete lateralCorr_;
  lateralCorr_ = vario;
}

void
ModelSettings::addDefaultLateralTravelTimeErrorCorr()
{
  Vario * lateralTravelTimeErrorCorr = new GenExpVario(1, 1.0, 1.0, 0.0); // NBNB OK shhould be input??
  timeLapseTravelTimeLateralCorrelation_.push_back(lateralTravelTimeErrorCorr);
}

void
ModelSettings::addLateralTravelTimeErrorCorr(Vario * vario)
{
  timeLapseTravelTimeLateralCorrelation_.push_back(vario);
}

void
ModelSettings::setBackgroundVario(Vario * vario)
{
  if (backgroundVario_ != NULL)
    delete backgroundVario_;
  backgroundVario_ = vario;
}

void
ModelSettings::setLocalWaveletVario(Vario * vario)
{
  if (localWaveletVario_ != NULL)
    delete localWaveletVario_;
  localWaveletVario_ = vario;
}

void
ModelSettings::copyBackgroundVarioToLocalWaveletVario(void)
{
  float range1 = backgroundVario_->getRange();
  float range2 = backgroundVario_->getSubRange();
  float angle  = backgroundVario_->getAngle();
  std::string type = backgroundVario_->getType();
  if (type == "Spherical")
  {
    localWaveletVario_ = new SphericalVario(range1, range2, angle);
  }
  else
  {
    GenExpVario * vario = dynamic_cast<GenExpVario *>(backgroundVario_);
    float power = vario->getPower();
    localWaveletVario_ = new GenExpVario(power, range1, range2, angle);
  }
}

void
ModelSettings::setAreaParameters(const SegyGeometry * geometry)
{
  if (geometry_ != NULL) {
    delete geometry_;      // Needed for snap_grid_to_seismic_data
  }
  geometry_ = new SegyGeometry(geometry);
}

void
ModelSettings::setSeismicDataAreaParameters(const SegyGeometry * geometry)
{
  if (geometry_full_ != NULL) {
    delete geometry_full_;
  }
  geometry_full_ = new SegyGeometry(geometry);
}

void
ModelSettings::setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat)
{
  if (traceHeaderFormat_ != NULL)
    delete traceHeaderFormat_;
  traceHeaderFormat_ = new TraceHeaderFormat(traceHeaderFormat);
}
void
ModelSettings::setTraceHeaderFormatOutput(TraceHeaderFormat * traceHeaderFormat)
{
  if (traceHeaderFormatOutput_ != NULL)
    delete traceHeaderFormatOutput_;
  traceHeaderFormatOutput_ = new TraceHeaderFormat(*traceHeaderFormat);
}

void
ModelSettings::setTraceHeaderFormatBackground(int parameter, TraceHeaderFormat * traceHeaderFormat)
{
  if (traceHeaderFormatBackground_[parameter] != NULL)
    delete traceHeaderFormatBackground_[parameter];
  traceHeaderFormatBackground_[parameter] = new TraceHeaderFormat(*traceHeaderFormat);
}

void
ModelSettings::addTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat)
{
  localTHF_.push_back(traceHeaderFormat);
}

void
ModelSettings::addTimeGradientSettings(float distance, float sigma_m)
{
  distanceFromWell_.push_back(distance);
  sigma_m_.push_back(sigma_m);
}

void
ModelSettings::addDefaultTimeGradientSettings(void)
{
  distanceFromWell_.push_back(100.0f);
  sigma_m_.push_back(1.0f);
}

void
ModelSettings::getTimeGradientSettings(float &distance, float &sigma_m, int t)
{
  distance = distanceFromWell_[t];
  sigma_m = sigma_m_[t];
}

void
ModelSettings::addVintage(int year, int month, int day)
{
  vintageYear_.push_back(year);
  vintageMonth_.push_back(month);
  vintageDay_.push_back(day);
}

void
ModelSettings::addDefaultVintage(void)
{
  vintageYear_.push_back(IMISSING);
  vintageMonth_.push_back(IMISSING);
  vintageDay_.push_back(IMISSING);
}

int
ModelSettings::getEstimateNumberOfWavelets(int t) const
{
  int n   = static_cast<int>(timeLapseEstimateWavelet_[t].size());
  int tot = 0;

  for (int i=0; i<n; i++)
  {
    if (timeLapseEstimateWavelet_[t][i] == 1)
      tot++;
  }
  return tot;
}

const std::vector<int>
ModelSettings::findSortedVintages(void) const
{
  int n = getNumberOfTimeLapses();

  std::vector<int> index(n);
  for(int i=0; i<n; i++)
    index[i] = i;

  int tmp;
  for(int i=0; i<n; i++){
    for(int j=i+1; j<n; j++){
      if(vintageYear_[j]<vintageYear_[i] ||
        (vintageYear_[j]==vintageYear_[i] && vintageMonth_[j]<vintageMonth_[i]) ||
        (vintageYear_[j]==vintageYear_[i] && vintageMonth_[j]==vintageMonth_[i] && vintageDay_[j]<vintageDay_[i])){
        tmp      = index[i];
        index[i] = index[j];
        index[j] = tmp;
      }
    }
  }
  return(index);
}


void
ModelSettings::MakeSureDzIsSetIfNeeded(InputFiles & input_files,
                                       std::string & err_txt)
{
  if(time_dz_ < 0) {
    bool need_dz = false;
    std::map<std::string, int>::iterator it = time_nz_.begin();
    while(need_dz == false && it != time_nz_.end()) {
      if(it->second < 0)
        need_dz = true;
      else
        ++it;
    }
    if(need_dz == true) {
      bool dz_ok = false;
      const std::vector<std::vector<std::string> > seis_names = input_files.getTimeLapseSeismicFiles();
      if(seis_names.size() > 0 && seis_names[0].size() > 0) {
        const std::string f_name = seis_names[0][0];
        if(f_name != "") {
          TraceHeaderFormat thf(0); //Just pick one, not important here.
          SegY master_file(f_name, 0, thf);
          time_dz_ = master_file.GetDz()*1.0000001; //Guard against truncation errors.
          dz_ok = true;
        }
      }
      if(dz_ok == false)
        err_txt += "No resolution given, at least one interval lacks number of layers, and no seismic file to find resolution from.\n";
    }
  }
}


int  ModelSettings::debugFlag_  = 0;
