#include <stdio.h>
#include <string.h>

#include "lib/global_def.h"
#include "lib/log.h"

#include "src/modelsettings.h"
#include "src/vario.h"

ModelSettings::ModelSettings(void)
{
  angularCorr_           = new GenExpVario(1, 10*static_cast<float>(PI/180.0)); // Power=1 range=10deg
  lateralCorr_           = new GenExpVario(1, 1000, 1000); 
  backgroundVario_       = new GenExpVario(1, 2000, 2000); 
  krigingParams_         =     NULL;
  angle_                 =     NULL;
  noiseEnergy_           =     NULL;
  matchEnergies_         =     NULL;
  faciesNames_           =     NULL;
  nFacies_               =        0;
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

  maxHz_background_      =    6.0f;       
  maxHz_seismic_         =   40.0f;       

  maxRankCorr_           =   0.98f;
  maxMergeDist_          =   0.01f;   // 0.01ms (approx. 2-3cm) 
  maxDevAngle_           =   15.0f;

  lowCut_                =    5.0f;   
  highCut_               =   55.0f;   

  wnc_                   =    0.1f;

  energyThreshold_       =    0.0f;

  minRelWaveletAmp_      =   0.05f;
  maxWaveletShift_       =   11.0f;
  waveletTaperingL_      =  200.0f;

  xPad_                  =    0.0f;   // If the PADDING command is not called, proper paddings
  yPad_                  =    0.0f;   // will be estimated in the Models() constructor 
  zPad_                  =    0.0f;   //

  segyOffset_            =    0.0f;
  p_undef_               =   0.01f;

  lzLimit_               =    0.47;   // NB! This is a double ==> do not use 'f'.

  outputFlag_            = VP+VS+RHO; //Default output
  formatFlag_            =       0;
  debugFlag_             =       0;
  fileGrid_              =      -1;

  generateSeismic_       =   false;
}

ModelSettings::~ModelSettings(void)
{ 
  if (angularCorr_ != NULL)
    delete angularCorr_;

  if (lateralCorr_ != NULL)
    delete lateralCorr_;

  if (backgroundVario_ != NULL)
    delete backgroundVario_;

  if(krigingParams_ != NULL)
    delete [] krigingParams_;

  //
  // NBNB-PAL: Purify liker ikke denne deleten, men skjønner ikke hvorfor
  //
  if (angle_ != NULL)
    delete [] angle_;

  if (noiseEnergy_ != NULL)
    delete [] noiseEnergy_;

  if (matchEnergies_ != NULL)
    delete [] matchEnergies_;

  if (faciesNames_ != NULL)
  {
    for (int i = 0 ; i < nFacies_ ; i++)
      delete [] faciesNames_[i];
    delete [] faciesNames_;
  }
}

bool 
ModelSettings::getDoInversion(void)
{
  return ((VP+VS+RHO+LAMELAMBDA+LAMEMU+POISSONRATIO+AI+SI+VPVSRATIO+MURHO+LAMBDARHO+FACIESPROB+CORRELATION & outputFlag_) > 0); 
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
ModelSettings::setKrigingParameters(float * krigingParams, int nParams)
{
  krigingParams_ = new float[nParams];
  for (int i=0 ; i<nParams ; i++)
    krigingParams_[i] = krigingParams[i];
}

void           
ModelSettings::setAngle(float * angle, int nAngles)
{
  if (angle_ != NULL)
    delete [] angle_;
  angle_ = new float[nAngles];
  for (int i=0 ; i<nAngles ; i++)
    angle_[i] = angle[i];
}

void           
ModelSettings::setNoiseEnergy(float * noiseEnergy, int nAngles)
{
  if (noiseEnergy_ != NULL)
    delete [] noiseEnergy_;
  noiseEnergy_ = new float[nAngles];
  for (int i=0 ; i<nAngles ; i++)
    noiseEnergy_[i] = noiseEnergy[i];
}

void           
ModelSettings::setMatchEnergies(float * waveletScale, int nAngles)
{
  if (matchEnergies_ != NULL)
    delete [] matchEnergies_;
  matchEnergies_ = new bool[nAngles];
  for (int i=0 ; i<nAngles ; i++)
    matchEnergies_[i] = (waveletScale[i] == RMISSING);
}

void           
ModelSettings::setFaciesNames(char ** faciesNames, int nFacies)
{
  faciesNames_ = new char * [nFacies];
  for (int i=0 ; i<nFacies ; i++) {
    faciesNames_[i] = new char [MAX_STRING];
    strcpy(faciesNames_[i], faciesNames[i]);
  }
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
