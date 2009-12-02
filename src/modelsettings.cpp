#include <stdio.h>
#include <string.h>
#include <map>
#include "lib/global_def.h"
#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"

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
    constBackValue_(3),
    indBGTrend_(0),
    indWavelet_(0),
    indFacies_(0),
    logNames_(5),
    inverseVelocity_(2),
    faciesLabels_(0),
    faciesNames_(0)
{
  angularCorr_           = new GenExpVario(1, 10*static_cast<float>(PI/180.0)); // Power=1 range=10deg
  lateralCorr_           = new GenExpVario(1, 1000, 1000); 
  backgroundVario_       = new GenExpVario(1, 2000, 2000); 
  localWaveletVario_     =     NULL; // Will be set equal to backgroundVario unless it is set separately
  geometry_              =     NULL;
  traceHeaderFormat_     =     NULL;
  krigingParameter_      =        0; // Indicate kriging not set.
  nWells_                =        0;
  nSimulations_          =        0;
  //
  // The original ranges were provided by Nam Hoai Pham (Statoil/25.09.2007)
  //
  alpha_min_             =  1300.0f;   // Nam: 1800
  alpha_max_             =  7000.0f;   // Nam: 6000
  beta_min_              =   200.0f;   // Nam:  800
  beta_max_              =  4200.0f;   // Nam: 3000
  rho_min_               =     1.4f;   // Nam:  1.5
  rho_max_               =     3.3f;   // Nam:  3.0

  var_alpha_min_         =   5.e-4f;  
  var_alpha_max_         = 250.e-4f; 
  var_beta_min_          =  10.e-4f; 
  var_beta_max_          = 500.e-4f; 
  var_rho_min_           =   2.e-4f; 
  var_rho_max_           = 100.e-4f; 

  maxHz_background_      =     6.0f;       
  maxHz_seismic_         =    40.0f;       

  maxRankCorr_           =    0.99f;
  maxMergeDist_          =    0.01f;   // 0.01ms (approx. 2-3cm) 
  maxDevAngle_           =    15.0f;

  lowCut_                =     5.0f;   
  highCut_               =    55.0f;   

  wnc_                   =     0.1f;

  energyThreshold_       =     0.0f;

  maxWellShift_          =    11.0f;
  maxWellOffset_         =   250.0f;

  minRelWaveletAmp_      =    0.05f;
  maxWaveletShift_       =    11.0f;
  waveletTaperingL_      =   200.0f;

  xPadFac_               =      0.0;   // If the PADDING command is not called, proper paddings
  yPadFac_               =      0.0;   // will be estimated in the Models() constructor 
  zPadFac_               =      0.0;   //

  nxPad_                 = IMISSING;   
  nyPad_                 = IMISSING;   
  nzPad_                 = IMISSING;   

  segyOffset_            =     0.0f;
  p_undef_               =    0.01f;

  lzLimit_               =     0.41;   // NB! This is a double ==> do not use 'f'.
  time_dTop_             = RMISSING;
  time_lz_               = RMISSING;
  time_dz_               = RMISSING;
  time_nz_               = IMISSING;
  velocityFromInv_       =    false;

  writePrediction_       =    false;  //Will be set to true if no simulations.
  gridFlag_              = IO::VP + IO::VS + IO::RHO;  // Default output
  defaultGridOutput_     =     true;
  formatFlag_            = IO::STORM;   
  domainFlag_            = IO::TIMEDOMAIN;   
  wellFlag_              =        0;   
  wellFormatFlag_        = IO::RMSWELL;
  otherFlag_             =        0;   
  debugFlag_             =        0;
  fileGrid_              =    false;

  estimationMode_        =    false;
  forwardModeling_       =    false;
  generateBackground_    =     true;
  estimateFaciesProb_    =    false;
  noVsFaciesProb_        =    false;
  faciesLogGiven_        =    false;
  depthDataOk_           =    false;
  parallelTimeSurfaces_  =    false;
  useLocalWavelet_       =    false;
  useLocalNoise_         =    false;
  optimizeWellLocation_  =    false;
  priorFaciesProbGiven_  =        0;

  generateSeismicAfterInversion_ = false;
  estimateBackground_    =     true;
  estimateCorrelations_  =     true;
  estimateWaveletNoise_  =     true;

  logLevel_              = LogKit::L_LOW;

  seed_                  =        0;
}

ModelSettings::~ModelSettings(void)
{ 
  if (angularCorr_ != NULL)
    delete angularCorr_;

  if (lateralCorr_ != NULL)
    delete lateralCorr_;

  if (backgroundVario_ != NULL)
    delete backgroundVario_;

  if (localWaveletVario_ != NULL)
    delete localWaveletVario_;

  if(geometry_ != NULL)
    delete geometry_;

  if(traceHeaderFormat_ != NULL)
    delete traceHeaderFormat_;
}

bool 
ModelSettings::getDoInversion(void)
{
  int flag = IO::VP 
    + IO::VS 
    + IO::RHO 
    + IO::LAMELAMBDA 
    + IO::LAMEMU 
    + IO::POISSONRATIO
    + IO::AI 
    + IO::SI
    + IO::VPVSRATIO 
    + IO::MURHO 
    + IO::LAMBDARHO 
    + IO::FACIESPROB
    + IO::CORRELATION 
    + IO::FACIESPROBRELATIVE;
  return ((flag & gridFlag_) > 0 && estimationMode_ == false); 
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
ModelSettings::setAngularCorr(Vario * vario)               
{
  if (angularCorr_ != NULL)
    delete angularCorr_;
  angularCorr_ = vario;
}

void           
ModelSettings::setLateralCorr(Vario * vario)
{
  if (lateralCorr_ != NULL)
    delete lateralCorr_;
  lateralCorr_ = vario;
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
  if (strcmp(backgroundVario_->getType(),"Spherical")==0)
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
  if (geometry_ == NULL) 
    geometry_ = new SegyGeometry(geometry);
}

void           
ModelSettings::setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat)
{
  if (traceHeaderFormat_ != NULL)
    delete traceHeaderFormat_;
  traceHeaderFormat_ = new TraceHeaderFormat(traceHeaderFormat);
}

void           
ModelSettings::addTraceHeaderFormat(TraceHeaderFormat * traceHeaderFormat)
{
  localTHF_.push_back(traceHeaderFormat);
}

void           
ModelSettings::setAllIndicatorsTrue(int nWells)
{
  indBGTrend_.resize(nWells);
  indWavelet_.resize(nWells);
  indFacies_.resize(nWells);
  for (int i=0 ; i<nWells ; i++)
  {
    indBGTrend_[i] = 1;
    indWavelet_[i] = 1;
    indFacies_[i]  = 1;
  }
}

void           
ModelSettings::setIndicatorBGTrend(int * indBGTrend, int nWells)
{
  for (int i=0 ; i<nWells ; i++)
    indBGTrend_[i] = indBGTrend[i];
}

void           
ModelSettings::setIndicatorWavelet(int * indWavelet, int nWells)
{
  for (int i=0 ; i<nWells ; i++)
    indWavelet_[i] = indWavelet[i];
}

void           
ModelSettings::setIndicatorFacies(int * indFacies, int nWells)
{
  for (int i=0 ; i<nWells ; i++)
    indFacies_[i] = indFacies[i];
}

void
ModelSettings::setNoiseScaled(Grid2D *ns)
{ 
  if(ns==NULL)
    noiseScaled_.push_back(NULL);
  else
    noiseScaled_.push_back(ns);
}

double ModelSettings::getMinimumNoiseScaled(int i)
{
  if(noiseScaled_[i]!=NULL)
  {
   double minimum = (*noiseScaled_[i])(0,0);
   int l,j;
   for(l=0;l<static_cast<int>(noiseScaled_[i]->GetNI());l++)
    for(j=0;j<static_cast<int>(noiseScaled_[i]->GetNJ());j++)
     if((*noiseScaled_[i])(l,j)<minimum)
       minimum = (*noiseScaled_[i])(l,j);

   return minimum;
  }
  else
    return 1.0;
}

bool ModelSettings::noiseIsScaled()
{
  if(noiseScaled_.size()>0)
  {
    for(int i=0;i<static_cast<int>(noiseScaled_.size());i++)
    {
      if(noiseScaled_[i]!=NULL)
        return true;
    }
    return false;
  }
  else
    return false;
}

int  ModelSettings::debugFlag_  = 0;
