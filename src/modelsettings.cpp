#include <stdio.h>
#include <string.h>

#include "lib/global_def.h"
#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/modelsettings.h"
#include "src/vario.h"

ModelSettings::ModelSettings(void)
  : seismicType_(0),
    angle_(0), 
    waveletScale_(0),
    SNRatio_(0),
    matchEnergies_(0),
    estimateWavelet_(0),
    estimateSNRatio_(0),
    constBackValue_(3),
    faciesLabels_(0),
    faciesNames_(0)
{
  angularCorr_           = new GenExpVario(1, 10*static_cast<float>(PI/180.0)); // Power=1 range=10deg
  lateralCorr_           = new GenExpVario(1, 1000, 1000); 
  backgroundVario_       = new GenExpVario(1, 2000, 2000); 
  localWaveletVario_     =     NULL; // Will be set equal to backgroundVario unless it is set separately
  geometry_              =     NULL;
  traceHeaderFormat_     = new TraceHeaderFormat(0); // SeisWorks;
  krigingParams_         =     NULL;
  indBGTrend_            =     NULL;
  indWavelet_            =     NULL;
  indFacies_             =     NULL;
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

  lzLimit_               =     0.47;   // NB! This is a double ==> do not use 'f'.
  time_dTop_             = RMISSING;
  time_lz_               = RMISSING;
  time_dz_               = RMISSING;
  time_nz_               = IMISSING;

  outputFlag_            = VP+VS+RHO;  // Default output
  formatFlag_            =        2;   // STORMFORMAT
  debugFlag_             =        0;
  fileGrid_              =       -1;
 
  generateSeismic_       =    false;
  generateBackground_    =    true;
  faciesLogGiven_        =    false;
  doDepthConversion_     =    false;
  parallelTimeSurfaces_  =    false;
  useLocalWavelet_       =    false;

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

  if(krigingParams_ != NULL)
    delete [] krigingParams_;

  if (indBGTrend_ != NULL)
    delete [] indBGTrend_;

  if (indWavelet_ != NULL)
    delete [] indWavelet_;

  if (indFacies_ != NULL)
    delete [] indFacies_;
}

bool 
ModelSettings::getDoInversion(void)
{
  return ((VP+VS+RHO+LAMELAMBDA+LAMEMU+POISSONRATIO+AI+SI+VPVSRATIO+MURHO+LAMBDARHO+FACIESPROB+CORRELATION+FACIESPROBRELATIVE & outputFlag_) > 0); 
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
ModelSettings::setKrigingParameters(float * krigingParams, int nParams)
{
  krigingParams_ = new float[nParams];
  for (int i=0 ; i<nParams ; i++)
    krigingParams_[i] = krigingParams[i];
}

void           
ModelSettings::setTraceHeaderFormat(const TraceHeaderFormat & traceHeaderFormat)
{
  if (traceHeaderFormat_ != NULL)
    delete traceHeaderFormat_;
  traceHeaderFormat_ = new TraceHeaderFormat(traceHeaderFormat);
}

void           
ModelSettings::setAllIndicatorsTrue(int nWells)
{
  indBGTrend_ = new int[nWells];
  indWavelet_ = new int[nWells];
  indFacies_  = new int[nWells];
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
ModelSettings::setOutputFlag(int outputFlag)
{
  if (outputFlag == PREDICTION)
    outputFlag_ = outputFlag_ | 1;
  else
    if ((outputFlag_ & PREDICTION) > 0)
      outputFlag_ = outputFlag | PREDICTION;
    else
      outputFlag_ = outputFlag;
}

void
ModelSettings::setFilePrefix(const std::string & filePrefix)               
{
  filePrefix_ = filePrefix;
}

std::string 
ModelSettings::makeFullFileName(const std::string name, const std::string postfix)
{
  return (filePrefix_+name+postfix);
}

std::string ModelSettings::filePrefix_ = "CRAVA_";
int         ModelSettings::debugFlag_  = 0;
