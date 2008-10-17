#include <iostream>
#include <iomanip>
#include <math.h>

#include "lib/random.h"
#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/model.h"

BlockedLogs::BlockedLogs(WellData  * well, 
                         Simbox    * simbox,
                         RandomGen * random) 
  : wellname_(NULL),
    xpos_(NULL),
    ypos_(NULL),
    zpos_(NULL),
    ipos_(NULL),
    jpos_(NULL),
    kpos_(NULL),
    alpha_(NULL),
    beta_(NULL),
    rho_(NULL),
    facies_(NULL),
    alpha_highcut_background_(NULL),
    beta_highcut_background_(NULL),
    rho_highcut_background_(NULL),
    alpha_highcut_seismic_(NULL),
    beta_highcut_seismic_(NULL),
    rho_highcut_seismic_(NULL),
    alpha_seismic_resolution_(NULL),
    beta_seismic_resolution_(NULL),
    rho_seismic_resolution_(NULL),
    real_seismic_data_(NULL),
    synt_seismic_data_(NULL),
    cpp_(NULL),
    nAngles_(0),
    firstM_(IMISSING),
    lastM_(IMISSING),
    firstB_(IMISSING),
    lastB_(IMISSING),
    nBlocks_(0),
    nLayers_(simbox->getnz()),
    nFacies_(0)
{
  blockWell(well, simbox, random);
}

//------------------------------------------------------------------------------
BlockedLogs::~BlockedLogs(void)
{
  if (wellname_ !=NULL)
    delete [] wellname_;
  if (xpos_ != NULL)
    delete [] xpos_;
  if (ypos_ != NULL)
    delete [] ypos_;
  if (zpos_ != NULL)
    delete [] zpos_;
  if (ipos_ != NULL)
    delete [] ipos_;
  if (jpos_ != NULL)
    delete [] jpos_;
  if (kpos_ != NULL)
    delete [] kpos_;
  if (alpha_ != NULL)
    delete [] alpha_;
  if (beta_ != NULL)
    delete [] beta_;
  if (rho_ != NULL)
    delete [] rho_;
  delete [] facies_;
  if (alpha_highcut_background_ != NULL)
    delete [] alpha_highcut_background_;
  if (beta_highcut_background_ != NULL)
    delete [] beta_highcut_background_;
  if (rho_highcut_background_ != NULL)
    delete [] rho_highcut_background_;
  if (alpha_highcut_seismic_ != NULL)
    delete [] alpha_highcut_seismic_;
  if (beta_highcut_seismic_ != NULL)
    delete [] beta_highcut_seismic_;
  if (rho_highcut_seismic_ != NULL)
    delete [] rho_highcut_seismic_;
  if (alpha_seismic_resolution_ != NULL)
    delete [] alpha_seismic_resolution_;
  if (beta_seismic_resolution_ != NULL)
    delete [] beta_seismic_resolution_;
  if (rho_seismic_resolution_ != NULL)
    delete [] rho_seismic_resolution_;
  if (real_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (real_seismic_data_[i] != NULL)
        delete real_seismic_data_[i];
  }
  if (synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (synt_seismic_data_[i] != NULL)
        delete synt_seismic_data_[i];
  }
  if (cpp_ != NULL) { 
    for (int i=0 ; i<nAngles_ ; i++)
      if (cpp_[i] != NULL)
        delete cpp_[i];
  }
}

//------------------------------------------------------------------------------
void 
BlockedLogs::blockWell(WellData  * well,
                       Simbox    * simbox,
                       RandomGen * random) 
{
  wellname_ = new char[MAX_STRING];
  strcpy(wellname_, well->getWellname());

  int * bInd = new int[well->getNd()]; // Gives which block each well log entry contributes to

  findSizeAndBlockPointers(well, simbox, bInd);  
  findBlockIJK(well, simbox, bInd);
  findBlockXYZ(simbox);

  int dummy;
  blockContinuousLog(bInd, well->getAlpha(dummy), alpha_);
  blockContinuousLog(bInd, well->getBeta(dummy), beta_);
  blockContinuousLog(bInd, well->getRho(dummy), rho_);
  blockContinuousLog(bInd, well->getAlphaBackgroundResolution(dummy), alpha_highcut_background_);
  blockContinuousLog(bInd, well->getBetaBackgroundResolution(dummy), beta_highcut_background_);
  blockContinuousLog(bInd, well->getRhoBackgroundResolution(dummy), rho_highcut_background_);
  blockContinuousLog(bInd, well->getAlphaSeismicResolution(dummy), alpha_highcut_seismic_);
  blockContinuousLog(bInd, well->getBetaSeismicResolution(dummy), beta_highcut_seismic_);
  blockContinuousLog(bInd, well->getRhoSeismicResolution(dummy), rho_highcut_seismic_);

  if (well->isFaciesLogDefined())
    blockDiscreteLog(bInd, well->getFacies(dummy), well->getFaciesNr(), well->getNFacies(), facies_, random);

  delete [] bInd;
}

//------------------------------------------------------------------------------
void 
BlockedLogs::blockContinuousLog(const int   *  bInd,
                                const float *  wellLog,
                                float       *& blockedLog)
{
  if (wellLog != NULL) {
    blockedLog  = new float[nBlocks_];
    int * count = new int[nBlocks_];
    //
    // Initialise arrays
    //
    for (int l = 0 ; l < nBlocks_ ; l++) {
      blockedLog[l] = 0.0f;
      count[l] = 0;
    }
    //
    // Block log
    //
    for (int m = firstM_ ; m < lastM_ + 1 ; m++) {
      if (wellLog[m] != RMISSING) {
        blockedLog[bInd[m]] += log(wellLog[m]); //NBNB-PAL: Flytt denne logaritmen nedover...
        count[bInd[m]]++;
        //LogKit::LogFormatted(LogKit::LOW,"m=%d bInd[m]  log(wellLog[m])  %d  %.5f \n",m,bInd[m],log(wellLog[m]));
      }
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      if (count[l] > 0) {
        blockedLog[l] /= count[l];
        //LogKit::LogFormatted(LogKit::LOW,"l=%d   count[l]=%d  sum=%.3f  blockedLog[l]=%.4f \n",l,count[l],sum, blockedLog[l]);
      }
      else
        blockedLog[l] = RMISSING;
    }
    delete [] count;
  }
}

//------------------------------------------------------------------------------
void 
BlockedLogs::blockDiscreteLog(const int *  bInd,
                              const int *  wellLog,
                              const int *  faciesNumbers,
                              int          nFacies,
                              int       *& blockedLog,
                              RandomGen *  random)
{
  if (wellLog != NULL) {
    //
    // Allocate memory and set undefined
    //
    faciesNumbers_ = faciesNumbers;
    nFacies_ = nFacies;
    blockedLog = new int[nBlocks_];
    for (int m = 0 ; m < nBlocks_ ; m++)
      blockedLog[m] = IMISSING;
    
    int   maxAllowedValue = 100;  // Largest allowed value (facies number).
    int * count = new int[nFacies];
    int * table = new int[maxAllowedValue];                  

    //
    // Set up facies-to-position table. 
    //
    // Example: If log values range from 2 to 4 the table looks like 
    //
    // table[0] = IMISSING
    // table[1] = IMISSING
    // table[2] =    0
    // table[3] =    1
    // table[4] =    2
    // table[5] = IMISSING
    //    .          .
    //    .          . 
    //
    for (int i = 0 ; i < maxAllowedValue ; i++)
      table[i] = IMISSING;
    for (int i = 0 ; i < nFacies ; i++)
      table[faciesNumbers[i]] = i;

    //
    // Block log
    //
    for (int i = 0 ; i < nFacies ; i++)
      count[i] = 0;
    int value = wellLog[firstM_];
    if(value!=IMISSING)
      count[table[value]]++;

    for (int m = firstM_+1 ; m < lastM_ + 1 ; m++) {
      if (bInd[m] != bInd[m - 1]) {
        blockedLog[bInd[m-1]] = findMostProbable(count, nFacies, random);
        for (int i = 0 ; i < nFacies ; i++)
          count[i] = 0;
      }
    value = wellLog[m];
    if(value!=IMISSING)
      count[table[value]]++;
    }
    blockedLog[bInd[lastM_]] = findMostProbable(count, nFacies, random);

    delete [] count;
    delete [] table;

    // 
    // NOTE: The blocked log contains internal numbers 0, 1, 2, ... and
    //       is NOT the facies lables.
    //
    //for (int b = 0 ; b < nBlocks_ ; b++)
    //  if (blockedLog[b] != IMISSING)
    //    LogKit::LogFormatted(LogKit::LOW,"b=%-3d   blockedLog[b]=%6d   (facies label=%d)\n",
    //                     b,blockedLog[b],faciesNumbers[blockedLog[b]]);
    //  else
    //    LogKit::LogFormatted(LogKit::LOW,"b=%-3d   blockedLog[b]=%6d\n",b,IMISSING);
  }
}

//------------------------------------------------------------------------------
int 
BlockedLogs::findMostProbable(const int * count,
                              int         nFacies,
                              RandomGen * rndgen) 
{ //
  // Find which value have the highest frequency. Add a random [0,1]
  // value to ensure that two Facies never get the same probability. 
  //
  double maxFreqValue = RMISSING;
  int    maxFreqIndex = IMISSING;
  for (int i=0 ; i < nFacies ; i++ ) {
    double freqValue = static_cast<double>(count[i]) + rndgen->unif01();
    if (freqValue > maxFreqValue && count[i]>0) {
      maxFreqValue = freqValue;
      maxFreqIndex = i;
    }
  }
  //No! Gives WARNING under Linux. Also when moved here
  //rndgen; //To eliminate warning 
  return (maxFreqIndex);
}

//------------------------------------------------------------------------------
void 
BlockedLogs::findSizeAndBlockPointers(WellData * well,
                                      Simbox   * simbox,
                                      int      * bInd)
{
  int dummy;
  const int      nd = well->getNd();
  const double * x  = well->getXpos(dummy);
  const double * y  = well->getYpos(dummy);
  const double * z  = well->getZpos(dummy);
  //
  // Find first cell in Simbox that the well hits
  // 
  int firstI(IMISSING);
  int firstJ(IMISSING);
  int firstK(IMISSING);
  for (int m = 0 ; m < nd ; m++) {
    simbox->getIndexes(x[m], y[m], z[m], firstI, firstJ, firstK);
    if (firstI != IMISSING && firstJ != IMISSING && firstK != IMISSING) {
      firstM_ = m;
      break;
    }
  }
  //
  // Find last cell in Simbox that the well hits
  // 
  int lastI(IMISSING);
  int lastJ(IMISSING);
  int lastK(IMISSING);  
  for (int m = nd - 1 ; m > 0 ; m--) {
    simbox->getIndexes(x[m], y[m], z[m], lastI, lastJ, lastK);
    if (lastI != IMISSING && lastJ != IMISSING && lastK != IMISSING) {
      lastM_ = m;
      break;
    }
  }
  //
  // Count number of blocks needed for the defined part of well.
  // 
  for (int m = 0 ; m < nd ; m++) {
    bInd[m] = IMISSING;
  }
  int newI, newJ, newK;
  int oldI = firstI;
  int oldJ = firstJ;
  int oldK = firstK;

  int nDefinedBlocks = 0;
  bInd[firstM_] = firstK; // The first defined well log entry contributes to this block.

  // 
  // Currently the well positions are given in float rather than double. Unfortunately, this 
  // allows a well to oscillate between two or more cells, leading to a breakdown of the 
  // algorithm below. To remedy for this we introduce array simboxInd which records the
  // indices of the simbox cells that are already accounted for, so that these are not
  // enlisted more than one time.
  //
  int * simboxInd = new int[nd];                                     // help hack
  const int nx    = simbox->getnx();                                 // help hack
  const int ny    = simbox->getny();                                 // help hack
  simboxInd[0] = nx*ny*oldK + nx*oldJ + oldI;                        // help hack

  for (int m = firstM_ + 1 ; m < lastM_ + 1 ; m++) {
    simbox->getIndexes(x[m], y[m], z[m], newI ,newJ, newK);

    if (newI != oldI || newJ != oldJ || newK != oldK) {

      int  thisInd = nx*ny*newK + nx*newJ + newI;                    // help hack
      bool blockNotListed = true;                                    // help hack
      for (int l = 0 ; l < nDefinedBlocks ; l++) {                   // help hack 
        if (thisInd == simboxInd[l]) {                               // help hack
          blockNotListed = false;                                    // help hack
          break;                                                     // help hack
        }                                                            // help hack
      }                                                              // help hack
      if (blockNotListed) {                                          // help hack
        simboxInd[nDefinedBlocks+1] = thisInd;                       // help hack
        oldI = newI;
        oldJ = newJ;
        oldK = newK;
        nDefinedBlocks++;
      }
      else {
      }
    }
    bInd[m] = firstK + nDefinedBlocks;
  }
  nDefinedBlocks++;
  nBlocks_ = firstK + nDefinedBlocks + (nLayers_ - lastK - 1);
  
  bool debug = false;
  if (debug) {  
    LogKit::LogFormatted(LogKit::LOW,"firstM_, lastM_          = %d, %d    \n",firstM_,lastM_);
    LogKit::LogFormatted(LogKit::LOW,"nLayers_                 = %d        \n",nLayers_);
    LogKit::LogFormatted(LogKit::LOW,"firstI,firstJ,firstK     = %d, %d, %d\n",firstI,firstJ,firstK);
    LogKit::LogFormatted(LogKit::LOW,"lastI,lastJ,lastK        = %d, %d, %d\n",lastI,lastJ,lastK);
    LogKit::LogFormatted(LogKit::LOW,"nDefinedBlocks, nBlocks_ = %d, %d    \n",nDefinedBlocks,nBlocks_);
  }
  delete [] simboxInd;
}

//------------------------------------------------------------------------------
void 
BlockedLogs::findBlockIJK(WellData  * well,
                          Simbox    * simbox,
                          const int * bInd)
{
  ipos_ = new int[nBlocks_];
  jpos_ = new int[nBlocks_];
  kpos_ = new int[nBlocks_];

  int dummy;
  const double * x = well->getXpos(dummy);
  const double * y = well->getYpos(dummy);
  const double * z = well->getZpos(dummy);
  //
  // Set IJK for virtual part of well in upper part of simbox
  // 
  int b = -1; // block counter;
  int firstI, firstJ, firstK;
  simbox->getIndexes(x[firstM_], y[firstM_], z[firstM_], firstI, firstJ, firstK);
  for (int k = 0 ; k < firstK ; k++) {
    b++;
    ipos_[b] = firstI;
    jpos_[b] = firstJ;
    kpos_[b] = k;
  }

  //
  // Set IJK for the defined part of the well
  // 
  b = firstK;
  ipos_[b] = firstI;
  jpos_[b] = firstJ;
  kpos_[b] = firstK;
  int i, j, k;
  for (int m = firstM_ + 1 ; m < lastM_ + 1 ; m++) {
    if (bInd[m] != bInd[m - 1]) {
      b++;
      simbox->getIndexes(x[m], y[m], z[m], i, j, k);
      ipos_[b] = i;
      jpos_[b] = j;
      kpos_[b] = k;      
    }
  }
  firstB_ = firstK;
  lastB_  = b;

  //
  // Set IJK for the virtual part of well in lower part of simbox
  // 
  int lastI,  lastJ,  lastK;  
  simbox->getIndexes(x[lastM_], y[lastM_], z[lastM_], lastI, lastJ, lastK);
  for (int k = lastK + 1 ; k < nLayers_ ; k++) {
    b++;
    ipos_[b] = lastI;
    jpos_[b] = lastJ;
    kpos_[b] = k;
  }

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::LOW,"firstB_, lastB_      = %d, %d    \n",firstB_,lastB_);
    for (int b = 0 ; b < nBlocks_ ; b++)
      LogKit::LogFormatted(LogKit::LOW,"b=%d   i,j,k=%d,%d,%d\n",b,ipos_[b],jpos_[b],kpos_[b]);
  }
}

//------------------------------------------------------------------------------
void 
BlockedLogs::findBlockXYZ(Simbox * simbox)
{
  xpos_ = new double[nBlocks_];
  ypos_ = new double[nBlocks_];
  zpos_ = new double[nBlocks_];
  for (int l = 0 ; l < nBlocks_ ; l++) {
    simbox->getCoord(ipos_[l],jpos_[l],kpos_[l],xpos_[l],ypos_[l],zpos_[l]);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogs::getVerticalTrend(const float * blockedLog,
                              float       * trend)
{
  if (blockedLog != NULL && trend != NULL) {
    int * count = new int[nLayers_];
    for (int k = 0 ; k < nLayers_ ; k++) {
      trend[k] = 0.0f;
      count[k] = 0;
    }
    for (int m = 0 ; m < nBlocks_ ; m++) {
      if (blockedLog[m] != RMISSING) {
        trend[kpos_[m]] += blockedLog[m];
        count[kpos_[m]]++;
      }
    }
    for (int k = 0 ; k < nLayers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];     
      else
      trend[k] = RMISSING;
    }
    delete [] count;
  }
  else {
    if (blockedLog == NULL)
      LogKit::LogFormatted(LogKit::LOW,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::LOW,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }    
}

//------------------------------------------------------------------------------
void
BlockedLogs::getVerticalTrend(const int * blockedLog,
                              int       * trend,
                              RandomGen * random)
{
  if (blockedLog != NULL && trend != NULL) {

    int * count = new int[nFacies_];
    for (int k = 0 ; k < nFacies_ ; k++) {
      count[k] = 0;
    }

    for (int k = 0 ; k < nLayers_ ; k++) {
      for (int i = 0 ; i < nFacies_ ; i++)
        count[i] = 0;
      for (int m = 0 ; m < nBlocks_ ; m++) {
        if (kpos_[m] == k) {
          if(blockedLog[m] != IMISSING) {
            count[blockedLog[m]]++;        // Count the number of times a facies occurs in layer 'k'
          }
        }
      }
      trend[k] = findMostProbable(count, nFacies_, random);
    }
    delete [] count;
  }
  else {
    if (blockedLog == NULL)
      LogKit::LogFormatted(LogKit::LOW,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::LOW,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }    
}


//------------------------------------------------------------------------------
void 
BlockedLogs::getBlockedGrid(FFTGrid * grid,
                            float   * blockedLog) 
{
  for (int m = 0 ; m < nBlocks_ ; m++) {
    //LogKit::LogFormatted(LogKit::LOW,"m=%d  ipos_[m], jpos_[m], kpos_[m] = %d %d %d\n",m,ipos_[m], jpos_[m], kpos_[m]);
    blockedLog[m] = grid->getRealValue(ipos_[m], jpos_[m], kpos_[m]);
  }
}

//------------------------------------------------------------------------------
void           
BlockedLogs::setLogFromVerticalTrend(float      * vertical_trend, 
                                     double       z0,              // z-value of center in top layer
                                     double       dz,              // dz in vertical trend
                                     int          nz,              // layers in vertical trend
                                     std::string  type)                  
{

  float * blockedLog = new float[nBlocks_];

  setLogFromVerticalTrend(blockedLog, zpos_, nBlocks_, 
                          vertical_trend, z0, dz, nz);


  if (type == "ALPHA_SEISMIC_RESOLUTION")
    alpha_seismic_resolution_ = blockedLog;
  else if (type == "BETA_SEISMIC_RESOLUTION")
    beta_seismic_resolution_ = blockedLog;
  else if (type == "RHO_SEISMIC_RESOLUTION")
    rho_seismic_resolution_ = blockedLog;
  else {
    LogKit::LogFormatted(LogKit::ERROR,"\nUnknown log type \"%s\" in BlockedLogs::setLogFromVerticalTrend()\n",
                         type.c_str());
    exit(1);
  }
}

//------------------------------------------------------------------------------
void           
BlockedLogs::setLogFromVerticalTrend(float     *& blockedLog,
                                     double     * zpos, 
                                     int          nBlocks,
                                     float      * vertical_trend, 
                                     double       z0,
                                     double       dzVt,
                                     int          nz)
{
  //
  // Initialise as undefined
  //
  for (int i=0 ; i<nBlocks ; i++)
    blockedLog[i] = RMISSING;
  
  //
  // Aritmethic mean of values in overlapping cells
  //
  for (int i=0 ; i<nBlocks ; i++) {
    double dz;
    if (i==nBlocks-1)
      dz = zpos[i]-zpos[i-1];
    else
      dz = zpos[i+1]-zpos[i];
    double zi = zpos[i]; 
    double a  = zi - 0.5*dz;     // Top of blocked log cell
    double b  = z0 + 0.5*dzVt;   // Base of first vertical trend cell
    
    int j=0;
    while (b<a && j<nz) {
      b += dzVt;
      j++;
    }
    // Now 'j' is the first vertical trend cell overlapping blocked log cell 'i'

    //
    // The if below treat end-of-vertical-trend cases and cases where one
    // single vertical-trend-cell covers a blocked log cell completely.
    //
    float value;
    if (j==nz || b > a+dz) {
      value = vertical_trend[j]; 
    }
    else {
      double zj = b + 0.5*dzVt; // Center of vertical trend cell
      value = vertical_trend[j]*(zj+dzVt-zi)/dzVt + vertical_trend[j]*(zi-zj)/dzVt;
    }
    blockedLog[i] = value;

    //LogKit::LogFormatted(LogKit::ERROR,"i j  log[i]   %d %d  %7.3f\n",i,j,log[i]);
  }
}

//------------------------------------------------------------------------------
void 
BlockedLogs::writeToFile(float dz,
                         int   type,
                         bool  exptrans) const
{
  char * tmpWellName = new char[MAX_STRING];
  char * filename    = new char[MAX_STRING];
  for (int i=0 ; i<=static_cast<int>(strlen(wellname_)) ; i++) // need to also copy ASCII null character
  {  
    if (wellname_[i]==' ' || wellname_[i]=='/')
      tmpWellName[i] = '_';
    else
      tmpWellName[i] = wellname_[i];
  }

  if (exptrans)
    sprintf(filename,"BW_%s",tmpWellName);
  else
    sprintf(filename,"lnBW_%s",tmpWellName);

  filename = ModelSettings::makeFullFileName(filename,".dat");

  FILE * file = fopen(filename, "w");

  float alpha, beta, rho;
  float z0 = dz/2.0f;

  for (int b = 0 ; b < nBlocks_ ; b++) 
  {
    if (type == 1) {
      alpha = alpha_[b];
      beta  = beta_[b];
      rho   = rho_[b];
    } 
    else if (type == 2) {
      alpha = alpha_highcut_background_[b];
      beta  = beta_highcut_background_[b];
      rho   = rho_highcut_background_[b];
    }
    else {
      alpha = alpha_highcut_seismic_[b];
      beta  = beta_highcut_seismic_[b];
      rho   = rho_highcut_seismic_[b];
    }

    if (alpha != RMISSING) { 
      if (exptrans) 
        alpha = exp(alpha);
    }
    else
      alpha = WELLMISSING;

    if (beta != RMISSING) { 
      if (exptrans) 
        beta = exp(beta);
    }
    else
      beta = WELLMISSING;
    
    if (rho != RMISSING) { 
      if (exptrans) 
        rho = exp(rho);
    }
    else
      rho = WELLMISSING;
    
    if (exptrans)
      fprintf(file,"%8.2f %8.2f %8.2f %8.5f\n",(z0 + kpos_[b]*dz),alpha,beta,rho);
    else
      fprintf(file,"%8.2f %8.5f %8.5f %8.5f\n",(z0 + kpos_[b]*dz),alpha,beta,rho);
  }
  fclose(file);

  delete [] tmpWellName;
  delete [] filename;
}

//------------------------------------------------------------------------------
void
BlockedLogs::writeRMSWell(ModelSettings * modelSettings)
{
  /// \todo Replace with safe open function.

  float maxHz_background = modelSettings->getMaxHzBackground();
  float maxHz_seismic    = modelSettings->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib2::Substitute(wellname,"/","_");
  NRLib2::Substitute(wellname," ","_");
  wellname = "BW_" + wellname;
  char * fileName = ModelSettings::makeFullFileName(wellname.c_str(),".rms");

  std::ofstream file(fileName, std::ios::out | std::ios::binary);

  if (!file) {
    throw new NRLib2::IOError("Error opening " + std::string(fileName) + " for writing.");
  }
  delete fileName;

  bool gotFacies      = nFacies_ > 0;
  bool gotRealSeismic = real_seismic_data_ != NULL;
  bool gotSyntSeismic = synt_seismic_data_ != NULL;
  bool gotCpp         = cpp_ != NULL;

  int nLogs = 3*4;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz, seisRes} 
  if (gotFacies)
    nLogs += 1;
  if (gotRealSeismic)
    nLogs += nAngles_;
  if (gotSyntSeismic)
    nLogs += nAngles_;
  if (gotCpp)
    nLogs += nAngles_;

  std::vector<std::string> params(3);
  params[0] = "Vp";
  params[1] = "Vs";
  params[2] = "Rho";

  file << std::fixed;
  file << std::setprecision(2);
  //
  // Write HEADER
  //
  file << "1.0\n";
  file << "CRAVA\n";
  file << wellname_ << " " << xpos_[0] << " " << ypos_[0] << "\n";
  file << nLogs << "\n";
  for (int i =0 ; i<3 ; i++) {
    file << params[i] << "   UNK lin\n";
    file << params[i] << static_cast<int>(maxHz_background) << " UNK lin\n";
    file << params[i] << static_cast<int>(maxHz_seismic)    << " UNK lin\n";
    file << params[i] << "_SeismicResolution UNK lin\n";
  }
  if (gotFacies) {
    file << "dummy   DISC ";
  //  file << faciesLogName_ << "   DISC ";
  //  for (int i =0 ; i < nFacies_ ; i++)
  //    file << " " << faciesNumbers_[i] << " " << faciesNames_[i];
    file << "\n";    
  }
  if (gotRealSeismic) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "RealSeis" << i << " UNK lin\n";
  }
  if (gotSyntSeismic) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "SyntSeis" << i << " UNK lin\n";
  }
  if (gotCpp) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "ReflCoef" << i << " UNK lin\n";
  }

  //
  // Write LOGS
  //
  for (int i=0 ; i<nBlocks_ ; i++) {
    file << std::right;
    file << std::setprecision(2);
    file << std::setw(9) << xpos_[i] << " " ;
    file << std::setw(10)<< ypos_[i] << " ";
    file << std::setw(7) << zpos_[i] << "  ";
    file << std::setw(7) << (alpha_[i]==RMISSING                    ? WELLMISSING : exp(alpha_[i]))                    << " ";
    file << std::setw(7) << (alpha_highcut_background_[i]==RMISSING ? WELLMISSING : exp(alpha_highcut_background_[i])) << " ";
    file << std::setw(7) << (alpha_highcut_seismic_[i]==RMISSING    ? WELLMISSING : exp(alpha_highcut_seismic_[i]))    << " ";
    file << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING ? WELLMISSING : exp(alpha_seismic_resolution_[i])) << "  ";
    file << std::setw(7) << (beta_[i]==RMISSING                     ? WELLMISSING : exp(beta_[i]))                     << " ";
    file << std::setw(7) << (beta_highcut_background_[i]==RMISSING  ? WELLMISSING : exp(beta_highcut_background_[i]))  << " ";
    file << std::setw(7) << (beta_highcut_seismic_[i]==RMISSING     ? WELLMISSING : exp(beta_highcut_seismic_[i]))     << " ";
    file << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING  ? WELLMISSING : exp(beta_seismic_resolution_[i]))  << "  ";
    file << std::setw(7) << std::setprecision(5);
    file << std::setw(7) << (rho_[i]==RMISSING                      ? WELLMISSING : exp(rho_[i]))                      << " ";
    file << std::setw(7) << (rho_highcut_background_[i]==RMISSING   ? WELLMISSING : exp(rho_highcut_background_[i]))   << " ";
    file << std::setw(7) << (rho_highcut_seismic_[i]==RMISSING      ? WELLMISSING : exp(rho_highcut_seismic_[i]))      << " ";
    file << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING   ? WELLMISSING : exp(rho_seismic_resolution_[i]))   << "  ";
    if (gotFacies)
      file << (facies_[i]==IMISSING                                 ? static_cast<int>(WELLMISSING) : facies_[i])      << " ";
    if (gotRealSeismic)
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(7) << (real_seismic_data_[a][i]==RMISSING ? WELLMISSING : real_seismic_data_[a][i])          << " ";
    if (gotSyntSeismic)
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(7) << (synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : synt_seismic_data_[a][i])          << " ";
    if (gotCpp)
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(7) << (cpp_[a][i]==RMISSING               ? WELLMISSING : cpp_[a][i])                        << " ";

    file << "\n";
  }
  file.close();
}

