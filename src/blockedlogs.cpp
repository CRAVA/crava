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
#include "src/io.h"

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
    alpha_for_facies_(NULL),
    rho_for_facies_(NULL),
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
  if (alpha_for_facies_ != NULL)
    delete [] alpha_for_facies_;
  if (rho_for_facies_ != NULL)
    delete [] rho_for_facies_;

  if (real_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (real_seismic_data_[i] != NULL)
        delete [] real_seismic_data_[i];
    delete [] real_seismic_data_;
  }
  if (synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (synt_seismic_data_[i] != NULL)
        delete [] synt_seismic_data_[i];
    delete [] synt_seismic_data_;
  }
  if (cpp_ != NULL) { 
    for (int i=0 ; i<nAngles_ ; i++)
      if (cpp_[i] != NULL)
        delete [] cpp_[i];
    delete [] cpp_;
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
  //We do not use block centre to represent BW position any more, due to ugly visualisation
  //findBlockXYZ(simbox);

  int dummy;
  blockCoordinateLog(bInd, well->getXpos(dummy), xpos_);
  blockCoordinateLog(bInd, well->getYpos(dummy), ypos_);
  blockCoordinateLog(bInd, well->getZpos(dummy), zpos_);
  findXYZforVirtualPart(simbox);

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
BlockedLogs::blockCoordinateLog(const int    *  bInd,
                                const double *  coord,
                                double       *& blockedCoord)
{
  if (coord != NULL) {
    blockedCoord  = new double[nBlocks_];
    int * count   = new int[nBlocks_];
    //
    // Initialise arrays
    //
    for (int l = 0 ; l < nBlocks_ ; l++) {
      blockedCoord[l] = 0.0f;
      count[l] = 0;
    }
    //
    // Block log
    //
    for (int m = firstM_ ; m < lastM_ + 1 ; m++) {
      blockedCoord[bInd[m]] += coord[m];
      count[bInd[m]]++;
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      if (count[l] > 0)
        blockedCoord[l] /= count[l];
      else
        blockedCoord[l]  = RMISSING;
    }
    delete [] count;
  }
}

//------------------------------------------------------------------------------
void 
BlockedLogs::findXYZforVirtualPart(Simbox * simbox)
{
  //
  // If the ends have undefined coordinates we use the nearest defined 
  // coordinate for x and y and the block cell centre for z
  //
  for (int b = 0 ; b < firstB_ ; b++)
  {
    double x,y,z;
    simbox->getCoord(ipos_[b], jpos_[b], kpos_[b], x, y, z);
    xpos_[b] = xpos_[firstB_];
    ypos_[b] = ypos_[firstB_];
    zpos_[b] = z;
  }

  for (int b = lastB_ + 1 ; b < nBlocks_ ; b++)
  {
    double x,y,z;
    simbox->getCoord(ipos_[b], jpos_[b], kpos_[b], x, y, z);
    xpos_[b] = xpos_[lastB_];
    ypos_[b] = ypos_[lastB_];
    zpos_[b] = z;
  }
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
  // The well positions used to be given in float rather than double. Unfortunately, this 
  // allowed a well to oscillate between two or more cells, leading to a breakdown of the 
  // algorithm below. To remedy for this we introduced array simboxInd which records the
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
  //
  // Why we cannot use nBlocks_ = nDefined:
  //
  // When we calculate the background model for each parameter we first
  // estimate a vertical trend in the total volume, anf then we interpolate
  // the blocked log intop this trend volume. To avoid sharp contrast we
  // ensure that the blocked log is defined from top to base of the volume.
  // In regions where the log is undefined we generate it by kriging from 
  // the rest of the log. Likewise, in regions where there is no blocked
  // log at all because the well was too short, we have to make a virtual
  // well.
  //
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
    LogKit::LogFormatted(LogKit::LOW,"firstB_, lastB_        = %d, %d    \n",firstB_,lastB_);
    LogKit::LogFormatted(LogKit::LOW,"firstI, firstJ, firstK = %d, %d, %d\n",firstI, firstJ, firstK);
    LogKit::LogFormatted(LogKit::LOW,"lastI,  lastJ,  lastK  = %d, %d, %d\n",lastI, lastJ, lastK);
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

  bool debug = true;
  if (debug) {
    for (int b = 0 ; b < nBlocks_ ; b++)
      LogKit::LogFormatted(LogKit::LOW,"b=%d   i,j,k=%d,%d,%d   x,y,z=%.2f, %.2f, %.2f\n",
                           b,ipos_[b],jpos_[b],kpos_[b],xpos_[b],ypos_[b],zpos_[b]);
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
BlockedLogs::setLogFromGrid(FFTGrid    * grid, 
                            int          iAngle,
                            int          nAngles,
                            std::string  type)
{
  float * blockedLog = new float[nBlocks_];

  for (int m = 0 ; m < nBlocks_ ; m++) {
    blockedLog[m] = grid->getRealValue(ipos_[m], jpos_[m], kpos_[m]);
  }

  if (nAngles_ == 0)
    nAngles_ = nAngles;

  if (type == "REFLECTION_COEFFICIENT") {
    if (cpp_ == NULL) 
      cpp_ = new float * [nAngles_];
    cpp_[iAngle] = blockedLog;
  }
  else if (type == "SEISMIC_DATA") {
    if (real_seismic_data_ == NULL)
      real_seismic_data_ = new float * [nAngles_];
    real_seismic_data_[iAngle] = blockedLog;
  }
  else {
    LogKit::LogFormatted(LogKit::ERROR,"\nUnknown log type \""+type
                         +"\" in BlockedLogs::setLogFromGrid()\n");
    exit(1);
  }
}

//------------------------------------------------------------------------------
void           
BlockedLogs::setLogFromVerticalTrend(float      * vertical_trend, 
                                     double       z0,              // z-value of center in top layer
                                     double       dz,              // dz in vertical trend
                                     int          nz,              // layers in vertical trend
                                     std::string  type,
                                     int          iAngle)                  
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
  else if (type == "SYNTHETIC_SEISMIC") {
    if (synt_seismic_data_ == NULL)
      synt_seismic_data_ = new float * [nAngles_]; // nAngles is set along with real_seismic_data_
    synt_seismic_data_[iAngle] = blockedLog;
  }
  else {
    LogKit::LogFormatted(LogKit::ERROR,"\nUnknown log type \""+type+
                         "\" in BlockedLogs::setLogFromVerticalTrend()\n");
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

    float value;
    if (j==nz) {
      // We have come to the end of the end-of-vertical-trend 
      value = vertical_trend[j-1]; 
    }
    else if (b >= a+dz) {
      // One single vertical-trend-cell covers a blocked log cell completely.
      value = vertical_trend[j]; 
    }
    else {
      double zj = b + 0.5*dzVt; // Center of vertical trend cell
      value = vertical_trend[j]* static_cast<float>((zj+dzVt-zi)/dzVt) + vertical_trend[j]*static_cast<float>((zi-zj)/dzVt);
    }
    blockedLog[i] = value;

    //LogKit::LogFormatted(LogKit::ERROR,"i j  blockedLog[i]   %d %d  %7.3f\n",i,j,blockedLog[i]);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogs::writeWell(ModelSettings * modelSettings)
{
  int formats = modelSettings->getWellFormatFlag();
  if((formats & ModelSettings::RMSWELL) > 0)
    writeRMSWell(modelSettings);
  if((formats & ModelSettings::NORSARWELL) > 0)
    writeNorsarWell(modelSettings);
}


void
BlockedLogs::writeRMSWell(ModelSettings * modelSettings)
{
  float maxHz_background = modelSettings->getMaxHzBackground();
  float maxHz_seismic    = modelSettings->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib::Substitute(wellname,"/","_");
  NRLib::Substitute(wellname," ","_");
  std::string baseName = IO::PrefixBlockedWells() + wellname + IO::SuffixRmsWells();
  std::string fileName = ModelSettings::makeFullFileName2(IO::PathToWells(), baseName);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName);

  if (!file) {
    LogKit::LogMessage(LogKit::ERROR,"Error opening "+fileName+" for writing.");
    std::exit(1);
  }

  bool gotFacies      = (nFacies_ > 0);
  bool gotRealSeismic = (real_seismic_data_ != NULL);
  bool gotSyntSeismic = (synt_seismic_data_ != NULL);
  bool gotCpp         = (cpp_ != NULL);
  bool gotFilteredLog = (alpha_seismic_resolution_ != NULL);
  bool gotVpRhoFacLog = (alpha_for_facies_ != NULL);
  int nLogs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz} 
  if(gotFilteredLog)
    nLogs += 3;
  if(gotVpRhoFacLog)
    nLogs += 2;
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

  file << std::fixed
       << std::setprecision(2);
  //
  // Write HEADER
  //
  file << "1.0\n"
       << "CRAVA\n"
       << wellname_ << " " << xpos_[firstB_] << " " << ypos_[firstB_] << "\n"
       << nLogs << "\n";
  
  for (int i =0 ; i<3 ; i++) {
    file << params[i] << "  UNK lin\n";
    file << params[i] << static_cast<int>(maxHz_background) << "  UNK lin\n";
    file << params[i] << static_cast<int>(maxHz_seismic)    << "  UNK lin\n";
  }
  if(gotFilteredLog) {
    for(int i=0;i<3;i++)
      file << params[i] << "_SeismicResolution UNK lin\n";
  }
  if(gotVpRhoFacLog) {
    file << params[0] << "_ForFacies UNK lin\n";
    file << params[2] << "_ForFacies UNK lin\n";
  }
  if (gotFacies) {
    file << "FaciesLog  DISC ";
    for (int i =0 ; i < modelSettings->getNumberOfFacies() ; i++)
      file << " " << modelSettings->getFaciesLabel(i) << " " << modelSettings->getFaciesName(i);
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
  for (int i=firstB_ ; i<lastB_ + 1 ; i++) {
    file << std::right
         << std::fixed
         << std::setprecision(2)
         << std::setw(9) << xpos_[i] << " "
         << std::setw(10)<< ypos_[i] << " "
         << std::setw(7) << zpos_[i] << "  "
         << std::setw(7) << (alpha_[i]==RMISSING                    ? WELLMISSING : exp(alpha_[i]))                    << " "
         << std::setw(7) << (alpha_highcut_background_[i]==RMISSING ? WELLMISSING : exp(alpha_highcut_background_[i])) << " "
         << std::setw(7) << (alpha_highcut_seismic_[i]==RMISSING    ? WELLMISSING : exp(alpha_highcut_seismic_[i]))    << " "
         << std::setw(7) << (beta_[i]==RMISSING                     ? WELLMISSING : exp(beta_[i]))                     << " "
         << std::setw(7) << (beta_highcut_background_[i]==RMISSING  ? WELLMISSING : exp(beta_highcut_background_[i]))  << " "
         << std::setw(7) << (beta_highcut_seismic_[i]==RMISSING     ? WELLMISSING : exp(beta_highcut_seismic_[i]))     << " "
         << std::setprecision(5)
         << std::setw(7) << (rho_[i]==RMISSING                      ? WELLMISSING : exp(rho_[i]))                      << " "
         << std::setw(7) << (rho_highcut_background_[i]==RMISSING   ? WELLMISSING : exp(rho_highcut_background_[i]))   << " "
         << std::setw(7) << (rho_highcut_seismic_[i]==RMISSING      ? WELLMISSING : exp(rho_highcut_seismic_[i]))      << " ";
    if(gotFilteredLog == true) {
      file << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING ? WELLMISSING : exp(alpha_seismic_resolution_[i])) << "  "
           << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING  ? WELLMISSING : exp(beta_seismic_resolution_[i]))  << "  "
           << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING   ? WELLMISSING : exp(rho_seismic_resolution_[i]))   << "  ";
    }
    if(gotVpRhoFacLog == true) {
      file << std::setw(7) << (alpha_for_facies_[i]==RMISSING ? WELLMISSING : exp(alpha_for_facies_[i])) << "  "
           << std::setw(7) << (rho_for_facies_[i]==RMISSING   ? WELLMISSING : exp(rho_for_facies_[i]))   << "  ";
    }
    if (gotFacies)
      file << (facies_[i]==IMISSING                                 ? static_cast<int>(WELLMISSING) : facies_[i])      << "  ";
    file << std::scientific;
    if (gotRealSeismic) {
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (real_seismic_data_[a][i]==RMISSING ? WELLMISSING : real_seismic_data_[a][i])          << " ";
      file << " ";
    }
    if (gotSyntSeismic) {
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : synt_seismic_data_[a][i])          << " ";
      file << " ";
    }
    if (gotCpp)
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (cpp_[a][i]==RMISSING               ? WELLMISSING : cpp_[a][i])                        << " ";
    file << "\n";
  }
  file.close();
}

void
BlockedLogs::writeNorsarWell(ModelSettings * modelSettings)
{
  //Note: At current, only write Vp, Vs and Rho, as others are not supported.
  float maxHz_background = modelSettings->getMaxHzBackground();
  float maxHz_seismic    = modelSettings->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib::Substitute(wellname,"/","_");
  NRLib::Substitute(wellname," ","_");

  //Handle main file.
  std::string baseName = IO::PrefixBlockedWells() + wellname + IO::SuffixNorsarWells();
  std::string fileName = ModelSettings::makeFullFileName2(IO::PathToWells(), baseName);
  std::ofstream mainFile;
  NRLib::OpenWrite(mainFile, fileName);
  mainFile << std::fixed
           << std::setprecision(2);

  int nData = lastB_ - firstB_ + 1;
  std::vector<double> md(nData,0);
  md[0] = zpos_[firstB_];
  double dmax = 0;
  double dmin = 1e+30;
  for(int i=firstB_+1;i<=lastB_;i++) {
    double dx = xpos_[i]-xpos_[i-1];
    double dy = ypos_[i]-ypos_[i-1];
    double dz = zpos_[i]-zpos_[i-1];
    double d  = sqrt(dx*dx+dy*dy+dz*dz);
    if(d > dmax)
      dmax = d;
    else if(d<dmin)
      dmin = d;
    md[i-firstB_] = md[i-firstB_-1] + d;
  }
  
  mainFile << "[Version information]\nVERSION 1000\nFORMAT ASCCI\n\n";
  mainFile << "[Well information]\n";
  mainFile << std::setprecision(5);
  mainFile << "MDMIN      km       " << md[0]*0.001       << "\n";
  mainFile << "MDMAX      km       " << md[nData-1]*0.001 << "\n";
  mainFile << "MDMINSTEP  km       " << dmin*0.001        << "\n";
  mainFile << "MDMAXSTEP  km       " << dmax*0.001        << "\n";
  mainFile << "UTMX       km       " << xpos_[0]*0.001    << "\n";
  mainFile << "UTMY       km       " << ypos_[0]*0.001    << "\n";
  mainFile << "EKB        km       " << 0.0f              << "\n";
  mainFile << "UNDEFVAL   no_unit " << WELLMISSING        << "\n\n";


  mainFile << "[Well track data information]\n";
  mainFile << "NUMMD  " << nData << "\n";
  mainFile << "NUMPAR 5\n";
  mainFile << "MD      km\n";
  mainFile << "TVD     km\n";
  mainFile << "TWT     s\n";
  mainFile << "UTMX    km\n";
  mainFile << "UTMY    km\n";
  
  std::string logBaseName = IO::PrefixBlockedWells() + wellname + IO::SuffixNorsarLog();
  std::string logFileName = ModelSettings::makeFullFileName2(IO::PathToWells(), logBaseName);
  std::string onlyName = NRLib::RemovePath(logFileName);
  
  bool gotFacies      = (nFacies_ > 0);
  bool gotRealSeismic = (real_seismic_data_ != NULL);
  bool gotSyntSeismic = (synt_seismic_data_ != NULL);
  bool gotCpp         = (cpp_ != NULL);
  bool gotFilteredLog = (alpha_seismic_resolution_ != NULL);
  bool gotVpRhoFacLog = (alpha_for_facies_ != NULL);

  int nLogs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz, seisRes} 
  if(gotFilteredLog)
    nLogs += 3;
  if (gotFacies)
    nLogs += 1;
  if (gotRealSeismic)
    nLogs += nAngles_;
  if (gotSyntSeismic)
    nLogs += nAngles_;
  if (gotCpp)
    nLogs += nAngles_;

  std::vector<std::string> params(3);
  params[0] = "VP";
  params[1] = "VS";
  params[2] = "RHO";

  std::vector<std::string> unit(3);
  unit[0] = "km/s";
  unit[1] = "km/s";
  unit[2] = "tons/m3";

  int nFiles = 3;
  if(gotFilteredLog)
    nFiles += 1;
  if(gotVpRhoFacLog)
    nFiles += 1;

  std::vector<std::string> postfix(5);
  postfix[0] = "orig";
  postfix[1] = NRLib::ToString(maxHz_background)+"Hz";
  postfix[2] = NRLib::ToString(maxHz_seismic)+"Hz";
  postfix[3] = "filtered";
  postfix[4] = "forFacies";

  for(int f=0;f<nFiles;f++) {
    mainFile << "\n[Well log data information]\n";
    mainFile << "LOGNAME log_" << postfix[f] << "\n";
    mainFile << "IN_FILE " << onlyName << f <<"\n";
    mainFile << "NUMMD " << nData << "\n";
    mainFile << "NUMPAR " << 4 << "\n"; //Also count md.
    mainFile << "MD      km\n";
    for (int i =0 ; i<3 ; i++)
      mainFile << params[i] << " " << unit[i] << "\n";
  }
  mainFile.close();

  //Write the track file.
  std::string trackBaseName = IO::PrefixBlockedWells() + wellname + IO::SuffixNorsarTrack();
  std::string trackFileName = ModelSettings::makeFullFileName2(IO::PathToWells(), trackBaseName);
  std::ofstream trackFile;
  NRLib::OpenWrite(trackFile, trackFileName.c_str());
  trackFile << std::right
            << std::fixed
            << std::setprecision(2)
            << "[NORSAR Well Track]\n";
    
  //Note: logFileName created above, needed in mainFile.
  std::vector<std::ofstream *> logFiles;
  logFiles.resize(nFiles);
  for(int f=0;f<nFiles;f++) {
    logFiles[f] = new std::ofstream();
    std::string fileName = logFileName+NRLib::ToString(f);
    NRLib::OpenWrite(*(logFiles[f]), fileName.c_str());
    *(logFiles[f]) << "[NORSAR Well Log]\n";
    *(logFiles[f]) << "[See header (.nwh) file for log information]\n";
  }
  
  for(int i = firstB_;i<=lastB_;i++) {
    trackFile << std::setprecision(5) << std::setw(7) 
              << md[i]*0.001 << " " << std::setw(7) << zpos_[i]*0.001 << " " << zpos_[i]*0.001
              << " " << std::setw(10)<< xpos_[i]*0.001 << " " << std::setw(10)<< ypos_[i]*0.001 << "\n";

    *(logFiles[0]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i]*0.001 << " "
                   << std::setw(7) << (alpha_[i]==RMISSING                    ? WELLMISSING : exp(alpha_[i]))                    << " "
                   << std::setw(7) << (beta_[i]==RMISSING                     ? WELLMISSING : exp(beta_[i]))                     << " "
                   << std::setw(7) << (rho_[i]==RMISSING                      ? WELLMISSING : exp(rho_[i]))                      << "\n";
    *(logFiles[1]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i]*0.001 << " "
                   << std::setw(7) << (alpha_highcut_background_[i]==RMISSING ? WELLMISSING : exp(alpha_highcut_background_[i])) << " "
                   << std::setw(7) << (beta_highcut_background_[i]==RMISSING  ? WELLMISSING : exp(beta_highcut_background_[i]))  << " "
                   << std::setw(7) << (rho_highcut_background_[i]==RMISSING   ? WELLMISSING : exp(rho_highcut_background_[i]))     << "\n";
    *(logFiles[2]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i]*0.001 << " "
                   << std::setw(7) << (alpha_highcut_seismic_[i]==RMISSING    ? WELLMISSING : exp(alpha_highcut_seismic_[i]))    << " "
                   << std::setw(7) << (beta_highcut_seismic_[i]==RMISSING     ? WELLMISSING : exp(beta_highcut_seismic_[i]))     << " "
                   << std::setw(7) << (rho_highcut_seismic_[i]==RMISSING      ? WELLMISSING : exp(rho_highcut_seismic_[i]))      << "\n";
    if(gotFilteredLog) {
      *(logFiles[3]) << std::right   << std::fixed << std::setprecision(5)
                     << std::setw(7) << md[i]*0.001 << " "
                     << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING ? WELLMISSING : exp(alpha_seismic_resolution_[i])) << " "
                     << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING  ? WELLMISSING : exp(beta_seismic_resolution_[i]))  << " "
                     << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING   ? WELLMISSING : exp(rho_seismic_resolution_[i]))   << "\n";
    }
    if(gotVpRhoFacLog) {
      *(logFiles[3]) << std::right   << std::fixed << std::setprecision(5)
                     << std::setw(7) << md[i]*0.001 << " "
                     << std::setw(7) << (alpha_for_facies_[i]==RMISSING ? WELLMISSING : exp(alpha_for_facies_[i])) << " "
                     << std::setw(7) << (rho_for_facies_[i]==RMISSING   ? WELLMISSING : exp(rho_for_facies_[i]))   << "\n";
    }
  }
  trackFile.close();
  for(int f=0;f<nFiles;f++)
    logFiles[f]->close();
}
  

void BlockedLogs::setSpatialFilteredLogs(float * filteredlog, int nData, std::string type, const float *bg)
{
  float * blockedLog = new float[nBlocks_];
  assert(nBlocks_ == nData);
  for(int i=0;i<nData;i++)
    blockedLog[i] = filteredlog[i]+bg[i];

  if (type == "ALPHA_SEISMIC_RESOLUTION")
    alpha_seismic_resolution_ = blockedLog;
  else if (type == "BETA_SEISMIC_RESOLUTION")
    beta_seismic_resolution_ = blockedLog;
  else if (type == "RHO_SEISMIC_RESOLUTION")
    rho_seismic_resolution_ = blockedLog;
  else if (type == "ALPHA_FOR_FACIES")
    alpha_for_facies_ = blockedLog;
  else if (type == "RHO_FOR_FACIES")
    rho_for_facies_ = blockedLog;
}
