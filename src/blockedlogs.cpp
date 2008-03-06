#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/fftgrid.h"
#include "src/simbox.h"

#include "lib/random.h"
#include "lib/log.h"

#include <iostream>
#include <math.h>

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
    alpha_background_resolution_(NULL),
    beta_background_resolution_(NULL),
    rho_background_resolution_(NULL),
    alpha_seismic_resolution_(NULL),
    beta_seismic_resolution_(NULL),
    rho_seismic_resolution_(NULL),
    firstM_(IMISSING),
    lastM_(IMISSING),
    firstB_(IMISSING),
    lastB_(IMISSING),
    nBlocks_(0),
    nLayers_(simbox->getnz())
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
  if (alpha_background_resolution_ != NULL)
    delete [] alpha_background_resolution_;
  if (beta_background_resolution_ != NULL)
    delete [] beta_background_resolution_;
  if (rho_background_resolution_ != NULL)
    delete [] rho_background_resolution_;
  if (alpha_seismic_resolution_ != NULL)
    delete [] alpha_seismic_resolution_;
  if (beta_seismic_resolution_ != NULL)
    delete [] beta_seismic_resolution_;
  if (rho_seismic_resolution_ != NULL)
    delete [] rho_seismic_resolution_;
}

//------------------------------------------------------------------------------
void 
BlockedLogs::blockWell(WellData  * well,
                       Simbox    * simbox,
                       RandomGen * random) 
{
  wellname_  = new char[MAX_STRING];
  strcpy(wellname_, well->getWellname());

  int * bInd = new int[well->getNd()]; // Gives which block each well log entry contributes to

  findSizeAndBlockPointers(well, simbox, bInd);  
  findBlockIJK(well, simbox, bInd);
  findBlockXYZ(simbox);

  int dummy;
  blockContinuousLog(bInd, well->getAlpha(dummy), alpha_);
  blockContinuousLog(bInd, well->getBeta(dummy), beta_);
  blockContinuousLog(bInd, well->getRho(dummy), rho_);
  blockContinuousLog(bInd, well->getAlphaBackgroundResolution(dummy), alpha_background_resolution_);
  blockContinuousLog(bInd, well->getBetaBackgroundResolution(dummy), beta_background_resolution_);
  blockContinuousLog(bInd, well->getRhoBackgroundResolution(dummy), rho_background_resolution_);
  blockContinuousLog(bInd, well->getAlphaSeismicResolution(dummy), alpha_seismic_resolution_);
  blockContinuousLog(bInd, well->getBetaSeismicResolution(dummy), beta_seismic_resolution_);
  blockContinuousLog(bInd, well->getRhoSeismicResolution(dummy), rho_seismic_resolution_);

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
        //LogKit::writeLog("m=%d bInd[m]  log(wellLog[m])  %d  %.5f \n",m,bInd[m],log(wellLog[m]));
      }
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      if (count[l] > 0) {
        blockedLog[l] /= count[l];
        //LogKit::writeLog("l=%d   count[l]=%d  sum=%.3f  blockedLog[l]=%.4f \n",l,count[l],sum, blockedLog[l]);
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
                              const int *  actualValues,
                              int          nActualValues,
                              int       *& blockedLog,
                              RandomGen *  random)
{
  if (wellLog != NULL) {
    //
    // Allocate memory and set undefined
    //
    actualValues_ = actualValues;
    nActualValues_ = nActualValues;
    blockedLog = new int[nBlocks_];
    for (int m = 0 ; m < nBlocks_ ; m++)
      blockedLog[m] = IMISSING;
    
    int   maxAllowedValue = 100;  // Largest allowed value (facies number).
    int * count = new int[nActualValues];
    int * table = new int[maxAllowedValue];                  

    //
    // Set up facies-to-position table
    //
    for (int i = 0 ; i < maxAllowedValue ; i++)
      table[i] = IMISSING;
    for (int i = 0 ; i < nActualValues ; i++)
      table[actualValues[i]] = i;

    //
    // Block log
    //
    for (int i = 0 ; i < nActualValues ; i++)
      count[i] = 0;
    int value = wellLog[firstM_];
    if(value!=IMISSING)
      count[table[value]]++;

    for (int m = firstM_+1 ; m < lastM_ + 1 ; m++) {
      if (bInd[m] != bInd[m - 1]) {
        blockedLog[bInd[m-1]] = findMostProbable(count, nActualValues, random);
        for (int i = 0 ; i < nActualValues ; i++)
          count[i] = 0;
      }
    value = wellLog[m];
    if(value!=IMISSING)
      count[table[value]]++;
    }
    blockedLog[bInd[lastM_]] = findMostProbable(count, nActualValues, random);

    delete [] count;
    delete [] table;

    // NBNB-PAL
    //for (int b = 0 ; b < nBlocks_ ; b++)
    //  LogKit::writeLog("b=%d   blockedLog[b]=%d\n",b,blockedLog[b]);
  }
}

//------------------------------------------------------------------------------
int 
BlockedLogs::findMostProbable(const int * count,
                              int         nActualValues,
                              RandomGen * rndgen) 
{ //
  // Find which value have the highest frequency. Add a random [0,1]
  // value to ensure that two Facies never get the same probability. 
  //
  double maxFreqValue = IMISSING;
  int    maxFreqIndex = IMISSING;
  for (int i=0 ; i < nActualValues ; i++ ) {
    double freqValue = static_cast<double>(count[i]) + rndgen->unif01();
    if (freqValue > maxFreqValue) {
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
    LogKit::writeLog("firstM_, lastM_          = %d, %d    \n",firstM_,lastM_);
    LogKit::writeLog("nLayers_                 = %d        \n",nLayers_);
    LogKit::writeLog("firstI,firstJ,firstK     = %d, %d, %d\n",firstI,firstJ,firstK);
    LogKit::writeLog("lastI,lastJ,lastK        = %d, %d, %d\n",lastI,lastJ,lastK);
    LogKit::writeLog("nDefinedBlocks, nBlocks_ = %d, %d    \n",nDefinedBlocks,nBlocks_);
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
    LogKit::writeLog("firstB_, lastB_      = %d, %d    \n",firstB_,lastB_);
    for (int b = 0 ; b < nBlocks_ ; b++)
      LogKit::writeLog("b=%d   i,j,k=%d,%d,%d\n",b,ipos_[b],jpos_[b],kpos_[b]);
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
    simbox->getCoord(ipos_[l],jpos_[l],kpos_[l],xpos_[l],zpos_[l],zpos_[l]);
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
    LogKit::writeLog("ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined log\n");
    exit(1);
  }    
}

//------------------------------------------------------------------------------
void
BlockedLogs::getVerticalTrendDiscrete(const int * blockedLog,
                                      int       * trend,
                                      RandomGen * random)
{
  if (blockedLog != NULL && trend != NULL) {
    int * count = new int[nActualValues_];
    for (int k = 0 ; k < nActualValues_ ; k++) {
      count[k] = 0;
    }
    int maxAllowedValue = 100;
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
    for (int i = 0 ; i < nActualValues_ ; i++)
      table[actualValues_[i]] = i;

    for (int k = 0 ; k < nLayers_ ; k++) {
      for (int i = 0 ; i < nActualValues_ ; i++)
        count[i] = 0;
      for (int m = 0 ; m < nBlocks_ ; m++) {
        if (kpos_[m] == k) {
          if(blockedLog[m] != IMISSING)
            count[table[blockedLog[m]]]++; // Count the number of times a facies occurs in layer 'k'
        }
      }
      trend[k] = findMostProbable(count, nActualValues_, random);
    }
    delete [] table;
    delete [] count;
  }
  else {
    LogKit::writeLog("ERROR in BlockedLogs::getVerticalTrendDiscrete(): Trying to use an undefined log\n");
    exit(1);
  }    
}


//------------------------------------------------------------------------------
void 
BlockedLogs::getBlockedGrid(FFTGrid * grid,
                            float   * blockedLog) 
{
  for (int m = 0 ; m < nBlocks_ ; m++) {
    blockedLog[m] = grid->getRealValue(ipos_[m], jpos_[m], kpos_[m]);
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
  for (int i=0 ; i<=(int)strlen(wellname_) ; i++) // need to also copy ASCII null character
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

  filename = LogKit::makeFullFileName(filename,".dat");

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
      alpha = alpha_background_resolution_[b];
      beta  = beta_background_resolution_[b];
      rho   = rho_background_resolution_[b];
    }
    else {
      alpha = alpha_seismic_resolution_[b];
      beta  = beta_seismic_resolution_[b];
      rho   = rho_seismic_resolution_[b];
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

