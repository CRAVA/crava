#include <iostream>
#include <iomanip>
#include <math.h>

#include "lib/lib_matr.h"
#include "nrlib/iotools/logkit.hpp"
#include "fftw.h"
#include "rfftw.h"
#include "fftw-int.h"
#include "f77_func.h"

#include "src/definitions.h"
#include "src/blockedlogs.h"
#include "src/welldata.h"
#include "src/wavelet.h"
#include "src/wavelet1D.h"
#include "src/fftgrid.h"
#include "src/simbox.h"
#include "src/modelsettings.h"
#include "src/io.h"

BlockedLogs::BlockedLogs(WellData  * well,
                         Simbox    * simbox,
                         bool        interpolate)
  : wellname_(""),
    xpos_(NULL),
    ypos_(NULL),
    zpos_(NULL),
    md_(NULL),
    ipos_(NULL),
    jpos_(NULL),
    kpos_(NULL),
    dz_(0),
    alpha_(NULL),
    beta_(NULL),
    rho_(NULL),
    facies_(NULL),
    facies_prob_(NULL),
    alpha_highcut_background_(NULL),
    beta_highcut_background_(NULL),
    rho_highcut_background_(NULL),
    alpha_highcut_seismic_(NULL),
    beta_highcut_seismic_(NULL),
    rho_highcut_seismic_(NULL),
    alpha_seismic_resolution_(NULL),
    beta_seismic_resolution_(NULL),
    rho_seismic_resolution_(NULL),
    alpha_predicted_(NULL),
    beta_predicted_(NULL),
    rho_predicted_(NULL),
    alpha_for_facies_(NULL),
    rho_for_facies_(NULL),
    real_seismic_data_(NULL),
    actual_synt_seismic_data_(NULL),
    well_synt_seismic_data_(NULL),
    cpp_(NULL),
    nAngles_(0),
    firstM_(IMISSING),
    lastM_(IMISSING),
    firstB_(IMISSING),
    lastB_(IMISSING),
    nBlocks_(0),
    nLayers_(simbox->getnz()),
    nFacies_(0),
    interpolate_(interpolate)
{
  blockWell(well, simbox, interpolate);
}

//------------------------------------------------------------------------------
BlockedLogs::~BlockedLogs(void)
{
  if (xpos_ != NULL)
    delete [] xpos_;
  if (ypos_ != NULL)
    delete [] ypos_;
  if (zpos_ != NULL)
    delete [] zpos_;
  if (md_ != NULL)
    delete [] md_;
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
  if (alpha_predicted_ != NULL)
    delete [] alpha_predicted_;
  if (beta_predicted_ != NULL)
    delete [] beta_predicted_;
  if (rho_predicted_ != NULL)
    delete [] rho_predicted_;
  if (alpha_for_facies_ != NULL)
    delete [] alpha_for_facies_;
  if (rho_for_facies_ != NULL)
    delete [] rho_for_facies_;

  if (facies_prob_ != NULL) {
    for (int i=0 ; i<nFacies_ ; i++)
      if (facies_prob_[i] != NULL)
        delete [] facies_prob_[i];
    delete [] facies_prob_;
  }
  if (real_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (real_seismic_data_[i] != NULL)
        delete [] real_seismic_data_[i];
    delete [] real_seismic_data_;
  }
  if (actual_synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (actual_synt_seismic_data_[i] != NULL)
        delete [] actual_synt_seismic_data_[i];
    delete [] actual_synt_seismic_data_;
  }
  if (well_synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (well_synt_seismic_data_[i] != NULL)
        delete [] well_synt_seismic_data_[i];
    delete [] well_synt_seismic_data_;
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
BlockedLogs::deleteDynamicBlockedLogs()
{
  if (alpha_seismic_resolution_ != NULL){
    delete [] alpha_seismic_resolution_;
    alpha_seismic_resolution_ = NULL;
  }
  if (beta_seismic_resolution_ != NULL){
    delete [] beta_seismic_resolution_;
    beta_seismic_resolution_ = NULL;
  }
  if (rho_seismic_resolution_ != NULL){
    delete [] rho_seismic_resolution_;
    rho_seismic_resolution_ = NULL;
  }

  if (alpha_predicted_ != NULL){
    delete [] alpha_predicted_;
    alpha_predicted_ = NULL;
  }
  if (beta_predicted_ != NULL){
    delete [] beta_predicted_;
    beta_predicted_ = NULL;
  }
  if (rho_predicted_ != NULL){
    delete [] rho_predicted_;
    rho_predicted_ = NULL;
  }

  if (alpha_for_facies_ != NULL){
    delete [] alpha_for_facies_;
    alpha_for_facies_ = NULL;
  }
  if (rho_for_facies_ != NULL){
    delete [] rho_for_facies_;
    rho_for_facies_ = NULL;
  }

  if (real_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (real_seismic_data_[i] != NULL)
        delete [] real_seismic_data_[i];
    delete [] real_seismic_data_;
    real_seismic_data_ = NULL;
  }
  if (actual_synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (actual_synt_seismic_data_[i] != NULL)
        delete [] actual_synt_seismic_data_[i];
    delete [] actual_synt_seismic_data_;
    actual_synt_seismic_data_ = NULL;
  }
  if (well_synt_seismic_data_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (well_synt_seismic_data_[i] != NULL)
        delete [] well_synt_seismic_data_[i];
    delete [] well_synt_seismic_data_;
    well_synt_seismic_data_ = NULL;
  }

  if (cpp_ != NULL) {
    for (int i=0 ; i<nAngles_ ; i++)
      if (cpp_[i] != NULL)
        delete [] cpp_[i];
    delete [] cpp_;
    cpp_ = NULL;
  }

  nAngles_ = 0;
}
//------------------------------------------------------------------------------
void
BlockedLogs::blockWell(WellData  * well,
                       Simbox    * simbox,
                       bool        interpolate)
{
  wellname_ = well->getWellname();

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
  if(well->getMD(dummy) != NULL)
    blockCoordinateLog(bInd, well->getMD(dummy), md_);


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
    blockDiscreteLog(bInd, well->getFacies(dummy), well->getFaciesNr(), well->getNFacies(), facies_);

  if(interpolate == true) {
    for(int i=1;i<well->getNd();i++) {
      if(abs(bInd[i]-bInd[i-1]) > 1) {
        int start, end;
        if(bInd[i] > bInd[i-1]) {
          start = bInd[i-1];
          end   = bInd[i];
        }
        else {
          start = bInd[i];
          end   = bInd[i-1];
        }
        for(int j = start+1;j<end;j++) {
          float t = static_cast<float>(j-start)/static_cast<float>(end-start);
          interpolateContinuousLog(xpos_, start, end, j, t);
          interpolateContinuousLog(ypos_, start, end, j, t);
          interpolateContinuousLog(zpos_, start, end, j, t);
          if(md_ != NULL)
            interpolateContinuousLog(md_, start, end, j, t);

          interpolateContinuousLog(alpha_, start, end, j, t);
          interpolateContinuousLog(beta_, start, end, j, t);
          interpolateContinuousLog(rho_, start, end, j, t);
          interpolateContinuousLog(alpha_highcut_background_, start, end, j, t);
          interpolateContinuousLog(beta_highcut_background_, start, end, j, t);
          interpolateContinuousLog(rho_highcut_background_, start, end, j, t);
          interpolateContinuousLog(alpha_highcut_seismic_, start, end, j, t);
          interpolateContinuousLog(beta_highcut_seismic_, start, end, j, t);
          interpolateContinuousLog(rho_highcut_seismic_, start, end, j, t);

          if (well->isFaciesLogDefined()) {
            if(facies_[j] == IMISSING)
              facies_[j] = facies_[j-1];
          }
        }
      }
    }
  }

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
        //LogKit::LogFormatted(LogKit::Low,"m=%d bInd[m]  log(wellLog[m])  %d  %.5f \n",m,bInd[m],log(wellLog[m]));
      }
    }
    for (int l = 0 ; l < nBlocks_ ; l++) {
      if (count[l] > 0) {
        blockedLog[l] /= count[l];
        //LogKit::LogFormatted(LogKit::Low,"l=%d   count[l]=%d  sum=%.3f  blockedLog[l]=%.4f \n",l,count[l],sum, blockedLog[l]);
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
                              int       *& blockedLog)
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
      if (bInd[m] != bInd[m - 1]) { // bInd[m] is the block number which sample 'm' lies in
        blockedLog[bInd[m-1]] = findMostProbable(count, nFacies, bInd[m-1]);
        for (int i = 0 ; i < nFacies ; i++)
          count[i] = 0;
      }
    value = wellLog[m];
    if(value!=IMISSING)
      count[table[value]]++;
    }
    blockedLog[bInd[lastM_]] = findMostProbable(count, nFacies, bInd[lastM_]);

    delete [] count;
    delete [] table;

    //
    // NOTE: The blocked log contains internal numbers 0, 1, 2, ... and
    //       is NOT the facies lables.
    //
    //for (int b = 0 ; b < nBlocks_ ; b++)
    //  if (blockedLog[b] != IMISSING)
    //    LogKit::LogFormatted(LogKit::Low,"b=%-3d   blockedLog[b]=%6d   (facies label=%d)\n",
    //                     b,blockedLog[b],faciesNumbers[blockedLog[b]]);
    //  else
    //    LogKit::LogFormatted(LogKit::Low,"b=%-3d   blockedLog[b]=%6d\n",b,IMISSING);
  }
}


void
BlockedLogs::interpolateContinuousLog(float * blockedLog, int start, int end,
                                      int index, float rel)
{
  if(blockedLog[start] != RMISSING && blockedLog[end] != RMISSING && blockedLog[index] == RMISSING)
    blockedLog[index] = rel*blockedLog[end]+(1-rel)*blockedLog[start];
}


void
BlockedLogs::interpolateContinuousLog(double * blockedLog, int start, int end,
                                      int index, float rel)
{
  if(blockedLog[start] != RMISSING && blockedLog[end] != RMISSING && blockedLog[index] == RMISSING)
    blockedLog[index] = rel*blockedLog[end]+(1-rel)*blockedLog[start];
}


//------------------------------------------------------------------------------
int
BlockedLogs::findMostProbable(const int * count,
                              int         nFacies,
                              int         blockIndex)
{
  int  maxIndex     = IMISSING;
  int  maxCount     = 0;
  bool inconclusive = false;

  for (int i=0 ; i < nFacies ; i++ ) {
    if (count[i] > 0 && count[i] > maxCount) {
      maxCount     = count[i];
      maxIndex     = i;
      inconclusive = false;
    }
    else if (count[i] > 0 && count[i] == maxCount) {
      inconclusive = true;
    }
  }

  if (inconclusive) {
    std::vector<int> equal;
    for (int i=0 ; i < nFacies ; i++ ) {
      if (count[i] == maxCount) {
        equal.push_back(i);
      }
    }
    int j = (blockIndex + 1) % equal.size();
    maxIndex = equal[j];
  }

  return (maxIndex);
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
    LogKit::LogFormatted(LogKit::Low,"firstM_, lastM_          = %d, %d    \n",firstM_,lastM_);
    LogKit::LogFormatted(LogKit::Low,"nLayers_                 = %d        \n",nLayers_);
    LogKit::LogFormatted(LogKit::Low,"firstI,firstJ,firstK     = %d, %d, %d\n",firstI,firstJ,firstK);
    LogKit::LogFormatted(LogKit::Low,"lastI,lastJ,lastK        = %d, %d, %d\n",lastI,lastJ,lastK);
    LogKit::LogFormatted(LogKit::Low,"nDefinedBlocks, nBlocks_ = %d, %d    \n",nDefinedBlocks,nBlocks_);
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

  dz_ = static_cast<float>(simbox->getRelThick(ipos_[0],jpos_[0])*simbox->getdz());

  bool debug = false;
  if (debug) {
    LogKit::LogFormatted(LogKit::Low,"firstB_, lastB_        = %d, %d    \n",firstB_,lastB_);
    LogKit::LogFormatted(LogKit::Low,"firstI, firstJ, firstK = %d, %d, %d\n",firstI, firstJ, firstK);
    LogKit::LogFormatted(LogKit::Low,"lastI,  lastJ,  lastK  = %d, %d, %d\n",lastI, lastJ, lastK);
    for (int b = 0 ; b < nBlocks_ ; b++)
      LogKit::LogFormatted(LogKit::Low,"b=%d   i,j,k=%d,%d,%d\n",b,ipos_[b],jpos_[b],kpos_[b]);
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
      LogKit::LogFormatted(LogKit::Low,"b=%d   i,j,k=%d,%d,%d   x,y,z=%.2f, %.2f, %.2f\n",
                           b,ipos_[b],jpos_[b],kpos_[b],xpos_[b],ypos_[b],zpos_[b]);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogs::getVerticalTrend(const float   * blockedLog,
                              float         * trend)
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
    if(interpolate_ == true)
      interpolateTrend(blockedLog,trend);
    delete [] count;
  }
  else {
    if (blockedLog == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }
}


void
BlockedLogs::getVerticalTrendLimited(const float                  * blockedLog,
                                     float                        * trend,
                                     const std::vector<Surface *> & limits)
{
  if (blockedLog != NULL && trend != NULL) {
    int * count = new int[nLayers_];
    for (int k = 0 ; k < nLayers_ ; k++) {
      trend[k] = 0.0f;
      count[k] = 0;
    }
    for (int m = 0 ; m < nBlocks_ ; m++) {
      if (blockedLog[m] != RMISSING) {
        if(limits.size() == 0 ||
           (limits[0]->GetZ(xpos_[m],ypos_[m]) <= zpos_[m] &&
            limits[1]->GetZ(xpos_[m],ypos_[m]) >= zpos_[m])) {
          trend[kpos_[m]] += blockedLog[m];
          count[kpos_[m]]++;
        }
      }
    }
    for (int k = 0 ; k < nLayers_ ; k++) {
      if (count[k] > 0)
        trend[k] = trend[k]/count[k];
      else
        trend[k] = RMISSING;
    }
    if(interpolate_ == true)
      interpolateTrend(blockedLog,trend,limits);
    delete [] count;
  }
  else {
    if (blockedLog == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrendLimited(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrendLimited(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogs::getVerticalTrend(const int * blockedLog,
                              int       * trend)
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
      trend[k] = findMostProbable(count, nFacies_, k);
    }
    if(interpolate_ == true)
      interpolateTrend(blockedLog,trend);

    delete [] count;
  }
  else {
    if (blockedLog == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined blocked log (NULL pointer)\n");
    if (trend == NULL)
      LogKit::LogFormatted(LogKit::Low,"ERROR in BlockedLogs::getVerticalTrend(): Trying to use an undefined trend (NULL pointer)\n");
    exit(1);
  }
}


//------------------------------------------------------------------------------
void
BlockedLogs::interpolateTrend(const float * blockedLog, float * trend)
{
  for (int m = 1 ; m < nBlocks_ ; m++) {
    if(abs(kpos_[m]-kpos_[m-1]) > 1) {
      int delta = 1;
      if(kpos_[m] < kpos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(kpos_[m]-kpos_[m-1]);
      float t = step_mult;
      for(int j = kpos_[m-1]+delta; j != kpos_[m];j++) {
        if(trend[j] == RMISSING)
          trend[j] = t*blockedLog[m]+(1-t)*blockedLog[m-1];
        t += step_mult;
      }
    }
  }
}


//------------------------------------------------------------------------------
void
BlockedLogs::interpolateTrend(const float * blockedLog, float * trend, const std::vector<Surface *> & limits)
{
  for (int m = 1 ; m < nBlocks_ ; m++) {
    if(abs(kpos_[m]-kpos_[m-1]) > 1) {
      int delta = 1;
      if(kpos_[m] < kpos_[m-1])
        delta = -1;
      float step_mult = static_cast<float>(delta)/static_cast<float>(kpos_[m]-kpos_[m-1]);
      float t = step_mult;
      for(int j = kpos_[m-1]+delta; j != kpos_[m];j++) {
        if(trend[j] == RMISSING) {
          if(limits.size() == 0 ||
             (limits[0]->GetZ(xpos_[m],ypos_[m]) <= zpos_[m] &&
              limits[1]->GetZ(xpos_[m],ypos_[m]) >= zpos_[m])) {
            trend[j] = t*blockedLog[m]+(1-t)*blockedLog[m-1];
          }
        }
        t += step_mult;
      }
    }
  }
}


//------------------------------------------------------------------------------
void
BlockedLogs::interpolateTrend(const int * blockedLog, int * trend)
{
  for (int m = 1 ; m < nBlocks_ ; m++) {
    if(abs(kpos_[m]-kpos_[m-1]) > 1) {
      int delta = 1;
      if(kpos_[m] < kpos_[m-1])
        delta = -1;
      for(int j = kpos_[m-1]+delta; j != kpos_[m];j++) {
        if(trend[j] == RMISSING)
          trend[j] = blockedLog[m-1];
      }
    }
  }
}



//------------------------------------------------------------------------------
void
BlockedLogs::getBlockedGrid(const FFTGrid * grid,
                            float   * blockedLog,
                            int       iOffset,
                            int       jOffset)
{
  for (int m = 0 ; m < nBlocks_ ; m++) {
    //LogKit::LogFormatted(LogKit::Low,"m=%d  ipos_[m], jpos_[m], kpos_[m] = %d %d %d\n",m,ipos_[m], jpos_[m], kpos_[m]);
    blockedLog[m] = grid->getRealValue(ipos_[m]+iOffset, jpos_[m]+jOffset, kpos_[m]);

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
  else if (type == "FACIES_PROB") {
    if (facies_prob_ == NULL)
      facies_prob_ = new float * [nFacies_];
    facies_prob_[iAngle] = blockedLog;
  }
  else if (type == "ALPHA_PREDICTED")
    alpha_predicted_ = blockedLog;
  else if (type == "BETA_PREDICTED")
    beta_predicted_ = blockedLog;
  else if (type == "RHO_PREDICTED")
    rho_predicted_ = blockedLog;

  else {
    LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type
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
  if (type != "WELL_SYNTHETIC_SEISMIC")
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
    else if (type == "ACTUAL_SYNTHETIC_SEISMIC") {
      if (actual_synt_seismic_data_ == NULL)
        actual_synt_seismic_data_ = new float * [nAngles_]; // nAngles is set along with real_seismic_data_
      actual_synt_seismic_data_[iAngle] = blockedLog;
    }
    else {
      LogKit::LogFormatted(LogKit::Error,"\nUnknown log type \""+type+
                           "\" in BlockedLogs::setLogFromVerticalTrend()\n");
      exit(1);
    }
  }
  else if (type == "WELL_SYNTHETIC_SEISMIC") {
    if (well_synt_seismic_data_ == NULL)
    {
      well_synt_seismic_data_ = new float * [nAngles_];
      for (int i=0; i<nAngles_; i++)
      {
        well_synt_seismic_data_[i] = new float[nBlocks_];
        for (int j=0; j<nBlocks_; j++)
          well_synt_seismic_data_[i][j] = RMISSING; //Declare in case the wavelet is not estimated for all angles
      }
    }
    setLogFromVerticalTrend(well_synt_seismic_data_[iAngle], zpos_, nBlocks_,
                            vertical_trend, z0, dz, nz);
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

    //LogKit::LogFormatted(LogKit::Error,"i j  blockedLog[i]   %d %d  %7.3f\n",i,j,blockedLog[i]);
  }
}

//------------------------------------------------------------------------------
void
BlockedLogs::writeWell(ModelSettings * modelSettings)
{
  int formats = modelSettings->getWellFormatFlag();
  if((formats & IO::RMSWELL) > 0)
    writeRMSWell(modelSettings);
  if((formats & IO::NORSARWELL) > 0)
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
  std::string baseName = IO::PrefixBlockedWells() + wellname.c_str() + IO::SuffixRmsWells();
  std::string fileName = IO::makeFullFileName(IO::PathToWells(), baseName);

  std::ofstream file;
  NRLib::OpenWrite(file, fileName);

  if (!file) {
    LogKit::LogMessage(LogKit::Error,"Error opening "+fileName+" for writing.");
    std::exit(1);
  }

  bool gotFacies            = (nFacies_ > 0);
  bool gotFaciesProb        = (facies_prob_ != NULL);
  bool gotRealSeismic       = (real_seismic_data_ != NULL);
  bool gotActualSyntSeismic = (actual_synt_seismic_data_ != NULL);
  bool gotWellSyntSeismic   = (well_synt_seismic_data_ != NULL);
  bool gotCpp               = (cpp_ != NULL);
  bool gotFilteredLog       = (alpha_seismic_resolution_ != NULL);
  bool gotVpRhoFacLog       = (alpha_for_facies_ != NULL);
  bool gotPredicted         = (alpha_predicted_ != NULL);

  int nLogs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz}
  if(gotFilteredLog)
    nLogs += 3;
  if(gotVpRhoFacLog)
    nLogs += 2;
  if (gotFacies)
    nLogs += 1;
  if (gotFaciesProb)
    nLogs += nFacies_;
  if (gotRealSeismic)
    nLogs += nAngles_;
  if (gotActualSyntSeismic)
    nLogs += nAngles_;
  if (gotWellSyntSeismic)
    nLogs += nAngles_;
  if (gotCpp)
    nLogs += nAngles_;
  if (gotPredicted)
    nLogs += 3;

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
       << IO::PrefixBlockedWells() + wellname_ << " " << xpos_[firstB_] << " " << ypos_[firstB_] << "\n"
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
  if (gotPredicted) {
    for (int i=0; i<3; i++)
      file << params[i] << "_Predicted UNK lin\n";
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
  if (gotFaciesProb) {
    for (int i=0 ; i<nFacies_ ; i++)
      file << "FaciesProbabilities" << i << " UNK lin\n";
  }
  if (gotRealSeismic) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "RealSeis" << i << " UNK lin\n";
  }
  if (gotActualSyntSeismic) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "ActualSyntSeis" << i << " UNK lin\n";
  }
  if (gotWellSyntSeismic) {
    for (int i=0 ; i<nAngles_ ; i++)
      file << "WellOptimizedSyntSeis" << i << " UNK lin\n";
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
    if (gotPredicted == true) {
      file << std::setw(7) << (alpha_predicted_[i]==RMISSING ? WELLMISSING : exp(alpha_predicted_[i])) << "  "
           << std::setw(7) << (beta_predicted_[i]==RMISSING  ? WELLMISSING : exp(beta_predicted_[i]))  << "  "
           << std::setw(7) << (rho_predicted_[i]==RMISSING   ? WELLMISSING : exp(rho_predicted_[i]))   << "  ";
    }
    if(gotVpRhoFacLog == true) {
      file << std::setw(7) << (alpha_for_facies_[i]==RMISSING ? WELLMISSING : exp(alpha_for_facies_[i])) << "  "
           << std::setw(7) << (rho_for_facies_[i]==RMISSING   ? WELLMISSING : exp(rho_for_facies_[i]))   << "  ";
    }
    if (gotFacies)
      file << (facies_[i]==IMISSING                                 ? static_cast<int>(WELLMISSING) : facies_[i])      << "  ";
    file << std::scientific;
    if (gotFaciesProb) {
      for (int a=0 ; a<nFacies_ ; a++)
        file << std::setw(12) << (facies_prob_[a][i]==RMISSING ? WELLMISSING : facies_prob_[a][i])          << " ";
      file << " ";
    }
    if (gotRealSeismic) {
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (real_seismic_data_[a][i]==RMISSING ? WELLMISSING : real_seismic_data_[a][i])          << " ";
      file << " ";
    }
    if (gotActualSyntSeismic) {
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (actual_synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : actual_synt_seismic_data_[a][i])          << " ";
      file << " ";
    }
    if (gotWellSyntSeismic) {
      for (int a=0 ; a<nAngles_ ; a++)
        file << std::setw(12) << (well_synt_seismic_data_[a][i]==RMISSING ? WELLMISSING : well_synt_seismic_data_[a][i])          << " ";
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
  double vertScale = 0.001;
  double horScale  = 0.001;

  //Note: At current, only write Vp, Vs and Rho, as others are not supported.
  float maxHz_background = modelSettings->getMaxHzBackground();
  float maxHz_seismic    = modelSettings->getMaxHzSeismic();

  std::string wellname(wellname_);
  NRLib::Substitute(wellname,"/","_");
  NRLib::Substitute(wellname," ","_");

  //Handle main file.
  std::string baseName = IO::PrefixBlockedWells() + wellname + IO::SuffixNorsarWells();
  std::string fileName = IO::makeFullFileName(IO::PathToWells(), baseName);
  std::ofstream mainFile;
  NRLib::OpenWrite(mainFile, fileName);
  mainFile << std::fixed
           << std::setprecision(2);

  int nData = lastB_ - firstB_ + 1;
  std::vector<double> md(nData,0);
  if(md_ == NULL)
    md[0] = zpos_[firstB_]*vertScale;
  else
    md[0] = md_[firstB_];
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
    if(md_ == NULL)
      md[i-firstB_] = md[i-firstB_-1] + d*vertScale;
    else
      md[i-firstB_] = md_[i];
  }

  mainFile << "[Version information]\nVERSION 1000\nFORMAT ASCII\n\n";
  mainFile << "[Well information]\n";
  mainFile << std::setprecision(5);
  mainFile << "MDMIN      km       " << md[0]*vertScale       << "\n";
  mainFile << "MDMAX      km       " << md[nData-1]*vertScale << "\n";
  mainFile << "MDMINSTEP  km       " << dmin*vertScale        << "\n";
  mainFile << "MDMAXSTEP  km       " << dmax*vertScale        << "\n";
  mainFile << "UTMX       km       " << xpos_[0]*horScale     << "\n";
  mainFile << "UTMY       km       " << ypos_[0]*horScale     << "\n";
  mainFile << "EKB        km       " << 0.0f                  << "\n";
  mainFile << "UNDEFVAL   no_unit  " << WELLMISSING           << "\n\n";


  mainFile << "[Well track data information]\n";
  mainFile << "NUMMD  " << nData << "\n";
  mainFile << "NUMPAR 5\n";
  mainFile << "MD      km\n";
  mainFile << "TVD     km\n";
  mainFile << "TWT     s\n";
  mainFile << "UTMX    km\n";
  mainFile << "UTMY    km\n";

  std::string logBaseName = IO::PrefixBlockedWells() + wellname + IO::SuffixNorsarLog();
  std::string logFileName = IO::makeFullFileName(IO::PathToWells(), logBaseName);
  std::string onlyName    = NRLib::RemovePath(logFileName);

  bool gotFacies            = (nFacies_ > 0);
  bool gotFaciesProb        = (facies_prob_ != NULL);
  bool gotRealSeismic       = (real_seismic_data_ != NULL);
  bool gotActualSyntSeismic = (actual_synt_seismic_data_ != NULL);
  bool gotWellSyntSeismic   = (well_synt_seismic_data_ != NULL);
  bool gotCpp               = (cpp_ != NULL);
  bool gotFilteredLog       = (alpha_seismic_resolution_ != NULL);
  bool gotVpRhoFacLog       = (alpha_for_facies_ != NULL);

  int nLogs = 3*3;   // {Vp, Vs, Rho} x {raw, BgHz, seisHz, seisRes}
  if(gotFilteredLog)
    nLogs += 3;
  if (gotFacies)
    nLogs += 1;
  if (gotFaciesProb)
    nLogs += nFacies_;
  if (gotRealSeismic)
    nLogs += nAngles_;
  if (gotActualSyntSeismic)
    nLogs += nAngles_;
  if (gotWellSyntSeismic)
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
  std::string trackFileName = IO::makeFullFileName(IO::PathToWells(), trackBaseName);
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
              << md[i-firstB_] << " " << std::setw(7) << zpos_[i]*vertScale << " " << zpos_[i]*vertScale
              << " " << std::setw(10)<< xpos_[i]*horScale << " " << std::setw(10)<< ypos_[i]*horScale << "\n";

    *(logFiles[0]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i-firstB_] << " "
                   << std::setw(7) << (alpha_[i]==RMISSING                    ? WELLMISSING : exp(alpha_[i]))                    << " "
                   << std::setw(7) << (beta_[i]==RMISSING                     ? WELLMISSING : exp(beta_[i]))                     << " "
                   << std::setw(7) << (rho_[i]==RMISSING                      ? WELLMISSING : exp(rho_[i]))                      << "\n";
    *(logFiles[1]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i-firstB_] << " "
                   << std::setw(7) << (alpha_highcut_background_[i]==RMISSING ? WELLMISSING : exp(alpha_highcut_background_[i])) << " "
                   << std::setw(7) << (beta_highcut_background_[i]==RMISSING  ? WELLMISSING : exp(beta_highcut_background_[i]))  << " "
                   << std::setw(7) << (rho_highcut_background_[i]==RMISSING   ? WELLMISSING : exp(rho_highcut_background_[i]))     << "\n";
    *(logFiles[2]) << std::right   << std::fixed << std::setprecision(5)
                   << std::setw(7) << md[i-firstB_] << " "
                   << std::setw(7) << (alpha_highcut_seismic_[i]==RMISSING    ? WELLMISSING : exp(alpha_highcut_seismic_[i]))    << " "
                   << std::setw(7) << (beta_highcut_seismic_[i]==RMISSING     ? WELLMISSING : exp(beta_highcut_seismic_[i]))     << " "
                   << std::setw(7) << (rho_highcut_seismic_[i]==RMISSING      ? WELLMISSING : exp(rho_highcut_seismic_[i]))      << "\n";
    if(gotFilteredLog) {
      *(logFiles[3]) << std::right   << std::fixed << std::setprecision(5)
                     << std::setw(7) << md[i-firstB_] << " "
                     << std::setw(7) << (alpha_seismic_resolution_[i]==RMISSING ? WELLMISSING : exp(alpha_seismic_resolution_[i])) << " "
                     << std::setw(7) << (beta_seismic_resolution_[i]==RMISSING  ? WELLMISSING : exp(beta_seismic_resolution_[i]))  << " "
                     << std::setw(7) << (rho_seismic_resolution_[i]==RMISSING   ? WELLMISSING : exp(rho_seismic_resolution_[i]))   << "\n";
    }
    if(gotVpRhoFacLog) {
      *(logFiles[3]) << std::right   << std::fixed << std::setprecision(5)
                     << std::setw(7) << md[i-firstB_] << " "
                     << std::setw(7) << (alpha_for_facies_[i]==RMISSING ? WELLMISSING : exp(alpha_for_facies_[i])) << " "
                     << std::setw(7) << (rho_for_facies_[i]==RMISSING   ? WELLMISSING : exp(rho_for_facies_[i]))   << "\n";
    }
  }
  trackFile.close();
  for(int f=0;f<nFiles;f++) {
    logFiles[f]->close();
    delete logFiles[f];
  }
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

void BlockedLogs::fillInCpp(const float * coeff,
                            int           start,
                            int           length,
                            fftw_real   * cpp_r,
                            int           nzp)
{
  int i;

  for(i=0;i<nzp;i++)
    cpp_r[i]=0;

  float * alphaVert = new float[nLayers_];
  float * betaVert  = new float[nLayers_];
  float * rhoVert   = new float[nLayers_];

  getVerticalTrend(alpha_, alphaVert);
  getVerticalTrend(beta_, betaVert);
  getVerticalTrend(rho_, rhoVert);

  for(i=start;i < start+length-1;i++)
  {
    float ei1 = computeElasticImpedance(alphaVert[i],betaVert[i],rhoVert[i],coeff);
    float ei2 = computeElasticImpedance(alphaVert[i+1],betaVert[i+1],rhoVert[i+1],coeff);
    cpp_r[i] =  ei2-ei1;
  }
  delete [] alphaVert;
  delete [] betaVert;
  delete [] rhoVert;

}

float BlockedLogs::computeElasticImpedance(float         alpha,
                                           float         beta,
                                           float         rho,
                                           const float * coeff) const
{
  // vp, vs, rho are logtransformed
  float angImp;

  angImp = float(coeff[0]*alpha+coeff[1]*beta+coeff[2]*rho );

  return(angImp);
}

void BlockedLogs::fillInSeismic(float     * seismicData,
                                int         start,
                                int         length,
                                fftw_real * seis_r,
                                int         nzp) const
{
  int i;
  for(i=0; i<nzp; i++)
    seis_r[i] = 0.0;

  for(i=start; i<start+length; i++)
    seis_r[i] = seismicData[i];
/*
  int lTregion = 3;
  int* modify  = getIndexPrior(start,lTregion,nzp);
  int* conditionto = getIndexPost(start-1,lTregion,nzp);
  //NBNB Odd: interpolate endpoints?
*/

}

void BlockedLogs::estimateCor(fftw_complex * var1_c,
                              fftw_complex * var2_c,
                              fftw_complex * ccor_1_2_c,
                              int            cnzp) const
{
  for(int i=0;i<cnzp;i++){
    ccor_1_2_c[i].re =  var1_c[i].re*var2_c[i].re + var1_c[i].im*var2_c[i].im;
    ccor_1_2_c[i].im = -var1_c[i].re*var2_c[i].im + var1_c[i].im*var2_c[i].re;
  }
}

void BlockedLogs::findContiniousPartOfData(const std::vector<bool> & hasData,
                                           int                       nz,
                                           int                     & start,
                                           int                     & length) const
{
  int  i;
  int  lPice            =  0;
  int  lengthMaxPice    = -1;
  int  startLongestPice =  0;
  bool previousHadData  = false;

  for(i = 0; i < nz ;i++){
    if(hasData[i]){
      if(! previousHadData)
        lPice=1;
      else
        lPice++;
      previousHadData = true;
    }
    else{
      if(previousHadData){
        if(lengthMaxPice < lPice){
          lengthMaxPice  = lPice;
          startLongestPice = i-lPice;
        }
      }
      previousHadData=false;
    }
  }

  if(previousHadData){
    if(lengthMaxPice < lPice){
      lengthMaxPice  = lPice;
      startLongestPice = i-lPice;
    }
  }

  start  = startLongestPice;
  length = lengthMaxPice;

  if (length == -1)
    length = 0;
}


void BlockedLogs::findOptimalWellLocation(FFTGrid                   ** seisCube,
                                          Simbox                     * timeSimbox,
                                          float                     ** reflCoef,
                                          int                          nAngles,
                                          const std::vector<float>   & angleWeight,
                                          float                        maxShift,
                                          int                          iMaxOffset,
                                          int                          jMaxOffset,
                                          const std::vector<Surface *> limits,
                                          int                        & iMove,
                                          int                        & jMove,
                                          float                      & kMove)
{
  int   polarity;
  int   i,j,k,l,m;
  int   start,length;
  float sum;
  float shiftF;
  float maxTot;
  float f1,f2,f3;

  int nx            = seisCube[0]->getNx();
  int ny            = seisCube[0]->getNy();
  int nzp           = seisCube[0]->getNzp();
  int cnzp          = nzp/2+1;
  int rnzp          = 2*cnzp;
  int iTotOffset    = 2*iMaxOffset+1;
  int jTotOffset    = 2*jMaxOffset+1;
  int polarityMax   = 0;
  float shift       = 0.0f;
  float maxValueTot = 0;
  float totalWeight = 0;
  float dz          = static_cast<float>(timeSimbox->getdz());

  float  * seisLog   = new float[nBlocks_];
  float  * seisData  = new float[nLayers_];
  float  * alphaVert = new float[nLayers_];
  float  * betaVert  = new float[nLayers_];
  float  * rhoVert   = new float[nLayers_];


  std::vector<int>   iOffset(iTotOffset);
  std::vector<int>   jOffset(jTotOffset);
  std::vector<int>   shiftI(nAngles);
  std::vector<int>   shiftIMax(nAngles);
  std::vector<float> maxValue(nAngles);
  std::vector<float> maxValueMax(nAngles);

  fftw_real    ** cpp_r               = new fftw_real*[nAngles];
  fftw_complex ** cpp_c               = reinterpret_cast<fftw_complex**>(cpp_r);

  fftw_real    ** cor_cpp_r           = new fftw_real*[nAngles];
  fftw_complex ** cor_cpp_c           = reinterpret_cast<fftw_complex**>(cor_cpp_r);

  fftw_real    ** seis_r              = new fftw_real*[nAngles];
  fftw_complex ** seis_c              = reinterpret_cast<fftw_complex**>(seis_r);

  fftw_real    ** ccor_seis_cpp_r     = new fftw_real*[nAngles];
  fftw_complex ** ccor_seis_cpp_c     = reinterpret_cast<fftw_complex**>(ccor_seis_cpp_r);

  fftw_real    ** ccor_seis_cpp_Max_r = new fftw_real*[nAngles];

  for( i=0; i<nAngles; i++ ){
    maxValueMax[i] = 0.0f;
    shiftIMax[i]   = 0;
  }

  // make offset vectors
  for(i=-iMaxOffset; i<iMaxOffset+1; i++){
    iOffset[i+iMaxOffset]=i;
  }
  for(j=-jMaxOffset; j<jMaxOffset+1; j++){
    jOffset[j+jMaxOffset]=j;
  }

  getVerticalTrendLimited(alpha_, alphaVert, limits);
  getVerticalTrendLimited(beta_, betaVert, limits);
  getVerticalTrendLimited(rho_, rhoVert, limits);

  std::vector<bool> hasData(nLayers_);
  for(i = 0 ; i < nLayers_ ; i++) {
    hasData[i] = alphaVert[i] != RMISSING && betaVert[i] != RMISSING && rhoVert[i] != RMISSING;
  }
  findContiniousPartOfData(hasData,nLayers_,start,length);

  for( j=0; j<nAngles; j++ ){
    seis_r[j]              = new fftw_real[rnzp];
    cpp_r[j]               = new fftw_real[rnzp];
    cor_cpp_r[j]           = new fftw_real[rnzp];
    ccor_seis_cpp_r[j]     = new fftw_real[rnzp];
    ccor_seis_cpp_Max_r[j] = new fftw_real[rnzp];
  }

  // Calculate reflection coefficients
  for( j=0; j<nAngles; j++ ){
    for(i=0; i<rnzp; i++){
      cpp_r[j][i] = 0;
    }
    fillInCpp(reflCoef[j],start,length,cpp_r[j],nzp);
    Utils::fft(cpp_r[j],cpp_c[j],nzp);
    estimateCor(cpp_c[j],cpp_c[j],cor_cpp_c[j],cnzp);
    Utils::fftInv(cor_cpp_c[j],cor_cpp_r[j],nzp);
  }

  std::vector<NRLib::Grid<float> > seisCubeSmall(nAngles,NRLib::Grid<float> (iTotOffset,jTotOffset,nBlocks_));

  for (j = 0 ; j < nAngles ; j++)
  {
    seisCube[j]->setAccessMode(FFTGrid::RANDOMACCESS);
    for (k = 0; k < iTotOffset; k++)
    {
      for (l = 0; l < jTotOffset; l++)
      {
        getBlockedGrid(seisCube[j],seisLog,iOffset[k],jOffset[l]);
        for (m = 0; m < nBlocks_; m++)
        {
          seisCubeSmall[j](k, l, m) = seisLog[m];
        }
      }
    }
    seisCube[j]->endAccess();
  }

  // Loop through possible well locations
  for(k=0; k<iTotOffset; k++){
    if(ipos_[0]+iOffset[k]<0 || ipos_[0]+iOffset[k]>nx-1) //Check if position is within seismic range
      continue;

    for(l=0; l<jTotOffset; l++){
      if(jpos_[0]+jOffset[l]<0 || jpos_[0]+jOffset[l]>ny-1) //Check if position is within seismic range
        continue;

      for( j=0; j<nAngles; j++ ){

        for (m=0; m<nBlocks_; m++)
          seisLog[m] = seisCubeSmall[j](k,l,m);

        getVerticalTrend(seisLog, seisData);
        fillInSeismic(seisData,start,length,seis_r[j],nzp);

        Utils::fft(seis_r[j],seis_c[j],nzp);
        estimateCor(seis_c[j],cpp_c[j],ccor_seis_cpp_c[j],cnzp);
        Utils::fftInv(ccor_seis_cpp_c[j],ccor_seis_cpp_r[j],nzp);
      }

      // if the sum from -maxShift to maxShift ms is
      // positive then polarity is positive
      dz = static_cast<float>(timeSimbox->getRelThick(ipos_[0]+iOffset[k],jpos_[0]+jOffset[l])*timeSimbox->getdz());
      sum = 0;
      for( j=0; j<nAngles; j++ ){
        if(angleWeight[j] > 0){
          for(i=0;i<ceil(maxShift/dz);i++)//zero included
            sum+=ccor_seis_cpp_r[j][i];
          for(i=0;i<floor(maxShift/dz);i++)
            sum+=ccor_seis_cpp_r[j][nzp-i-1];
        }
      }
      polarity=-1;
      if(sum > 0)
        polarity=1;

      // Find maximum correlation and corresponding shift for each angle
      maxTot = 0.0;
      for( j=0; j<nAngles; j++ ){
        if(angleWeight[j]>0){
          maxValue[j] = 0.0f;
          shiftI[j]=0;
          for(i=0;i<ceil(maxShift/dz);i++){
            if(ccor_seis_cpp_r[j][i]*polarity > maxValue[j]){
              maxValue[j] = ccor_seis_cpp_r[j][i]*polarity;
              shiftI[j] = i;
            }
          }
          for(i=0;i<floor(maxShift/dz);i++){
            if(ccor_seis_cpp_r[j][nzp-1-i]*polarity > maxValue[j]){
              maxValue[j] = ccor_seis_cpp_r[j][nzp-1-i]*polarity;
              shiftI[j] = -1-i;
            }
          }
          maxTot += angleWeight[j]*maxValue[j]; //Find weighted total maximum correlation
        }
      }

      if(maxTot > maxValueTot){
        maxValueTot = maxTot;
        polarityMax = polarity;
        iMove       = iOffset[k];
        jMove       = jOffset[l];
        for(m=0; m<nAngles; m++){
          shiftIMax[m]   = shiftI[m];
          maxValueMax[m] = maxValue[m];
          for(i=0;i<rnzp;i++)
            ccor_seis_cpp_Max_r[m][i] = ccor_seis_cpp_r[m][i];
        }
      }
    }
  }

  for(i=0; i<nAngles; i++){
    shiftI[i] = shiftIMax[i];
    maxValue[i] = maxValueMax[i];
    for(j=0;j<rnzp;j++)
      ccor_seis_cpp_r[i][j] = ccor_seis_cpp_Max_r[i][j];
  }
  polarity = polarityMax;

  // Find kMove in optimal location
  for(j=0; j<nAngles; j++){
    if(angleWeight[j]>0){
      if(shiftI[j] < 0){
        if(ccor_seis_cpp_r[j][nzp+shiftI[j]-1]*polarity < maxValue[j]) //then local max
        {
          f1 = ccor_seis_cpp_r[j][nzp+shiftI[j]-1];
          f2 = ccor_seis_cpp_r[j][nzp+shiftI[j]];
          int ind3;
          if(shiftI[j]==-1)
            ind3 = 0;
          else
            ind3=nzp+shiftI[j]+1;
          f3 = ccor_seis_cpp_r[j][ind3];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI[j]+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI[j]);
      }
      else //positive or zero shift
      {
        if(ccor_seis_cpp_r[j][shiftI[j]+1]*polarity < maxValue[j]) //then local max
        {
          f3 = ccor_seis_cpp_r[j][shiftI[j]+1];
          f2 = ccor_seis_cpp_r[j][shiftI[j]];
          int ind1;
          if(shiftI[j]==0)
            ind1 = nzp-1;
          else
            ind1=shiftI[j]-1;
          f1 = ccor_seis_cpp_r[j][ind1];
          float x0=(f1-f3)/(2*(f1+f3-2*f2));
          shiftF=shiftI[j]+x0;
        }
        else  // do as good as we can
          shiftF=float(shiftI[j]);
      }
      shift += angleWeight[j]*shiftF*dz;//weigthing shift according to wellWeight
      totalWeight += angleWeight[j];
    }
  }

  shift/=totalWeight;
  kMove = shift;

  for( j=0; j<nAngles; j++ ){
    delete [] ccor_seis_cpp_Max_r[j];
    delete [] ccor_seis_cpp_r[j];
    delete [] cor_cpp_r[j];
    delete [] seis_r[j];
    delete [] cpp_r[j];
  }
  delete [] alphaVert;
  delete [] betaVert;
  delete [] rhoVert;
  delete [] seisLog;
  delete [] seisData;

  delete [] ccor_seis_cpp_Max_r;
  delete [] ccor_seis_cpp_r;
  delete [] cor_cpp_r;
  delete [] seis_r;
  delete [] cpp_r;
}


void BlockedLogs::setSeismicGradient(double v0,
                                      const NRLib::Grid2D<float>   &    structureDepthGradX,
                                      const NRLib::Grid2D<float>   &    structureDepthGradY,
                                      const NRLib::Grid2D<float>   &    refTimeGradX ,
                                      const NRLib::Grid2D<float>   &    refTimeGradY,
                                      std::vector<double>        & xGradient,
                                      std::vector<double>        & yGradient)
{
  xGradient.resize(nBlocks_);
  yGradient.resize(nBlocks_);

  double mp= 2.0/(v0*0.001);
  for(int k = 0; k < nBlocks_; k++){
    int i = ipos_[k];
    int j = jpos_[k];
    xGradient[k]= structureDepthGradX(i,j)*mp+refTimeGradX(i,j);
    yGradient[k]= structureDepthGradY(i,j)*mp+refTimeGradY(i,j);
  }
}

void
BlockedLogs::setTimeGradientSettings(float distance, float sigma_m)
{
  lateralThresholdGradient_ = distance;
  sigma_m_ = sigma_m;
}




void BlockedLogs::findSeismicGradient(FFTGrid                  ** seisCube,
                                      Simbox                   * timeSimbox,
                                      int                         nAngles,
                                      std::vector<double>       & xGradient,
                                      std::vector<double>       & yGradient,
                                      std::vector<std::vector<double> > & Sigma_gradient)
{
  int i, j, k, l;
  int xEx = 2;
  int yEx = 2;
  int nZx = (2 * xEx + 1);
  int nZy = (2 * yEx + 1);

  xGradient.resize(nBlocks_);
  yGradient.resize(nBlocks_);
  std::vector<double> tmpXgrad(nBlocks_);
  std::vector<double> tmpYgrad(nBlocks_);
  std::vector<double> Qepsilon(4*nBlocks_);
  std::vector<double> Qepsilon_data(2*nBlocks_);

  std::vector<float> seisTrace;
  std::vector<double> zShift(int(nZx*nZy*nBlocks_));

  //seismic peak position and characteristics in well
  std::vector<double> zPeakWell;
  std::vector<double> peakWell;
  std::vector<double> bWell;

  //seismic peak position and characterisitics in trace
  std::vector<double> zPeak;
  std::vector<double> peak;
  std::vector<double> b;

  int i0 = ipos_[0];
  int j0 = jpos_[0];


  //Check if well needs to change position in order for the whole
  //shift region to be contained in the seismic cube
  //NBNB marita Må testes om det fungerer for forskjellige brønner
  int nx = timeSimbox->getnx();
  int ny = timeSimbox->getny();
  int di_neg = 0; int di_pos = 0; int dj_neg = 0; int dj_pos = 0; //The max replacement in well in x and y direction.
  for(k = 1; k < nBlocks_; k++){
    di_pos = std::max(i0 - ipos_[k] , di_pos);
    di_neg = std::min(i0 - ipos_[k] , di_neg);
    dj_pos = std::max(j0 - jpos_[k] , dj_pos);
    dj_neg = std::min(j0 - jpos_[k] , dj_neg);
  }
  i0 = std::max(xEx, i0 + di_neg) - di_neg;
  i0 = std::min(nx - xEx - 1, i0 + di_pos) - di_pos;
  j0 = std::max(yEx, j0 + dj_neg) - dj_neg;
  j0 = std::min(ny - yEx - 1, j0 + dj_pos) - dj_pos;


  int di = i0 - ipos_[0];
  int dj = j0 - jpos_[0];
  if(di != 0 || dj != 0){
    //adjust the well location
    for(k = 0; k < nBlocks_; k++){
      ipos_[k] += di;
      jpos_[k] += dj;
    }
  }

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/traces");
  std::ofstream out(buffer);
  delete [] buffer;

  double dz, ztop, dzW, ztopW;
  for(l = 0; l < nAngles; l++){
    for(j = -yEx; j <= yEx; j++){
      for(i = -xEx; i <= xEx; i++){
        seisTrace = seisCube[l]->getRealTrace2(i0, j0);

        smoothTrace(seisTrace);
        //if(j == 0 ){
        //  for(int s = 0; s < seisTrace.size(); s++)
        //    out << seisTrace[s] << std::endl;
        //}

        dzW =  timeSimbox->getdz(i0,j0);
        ztopW =  timeSimbox->getTop(i0,j0);
        findPeakTrace(seisTrace, zPeakWell, peakWell, bWell, dzW, ztopW);

        seisTrace = seisCube[l]->getRealTrace2(i0+i, j0+j);
        smoothTrace(seisTrace);
        if(i==0){
          for(size_t s = 0; s < seisTrace.size(); s++)
            out << seisTrace[s] << std::endl;
        }

        dz =  timeSimbox->getdz(i0+i, j0+j);
        ztop =  timeSimbox->getTop(i0+i, j0+j);
        findPeakTrace(seisTrace, zPeak, peak, b, dz, ztop);

        peakMatch(zPeak,peak,b,zPeakWell,peakWell,bWell);//Finds the matching peaks in two traces
        zShift[(i+2) + (j+2)*nZx] = computeShift(zPeak,zPeakWell,zpos_[0]);

        for(k = 1; k < nBlocks_; k++){
          //Check if well changes lateral position
          if((ipos_[k]- ipos_[k-1] == 0) && (jpos_[k] - jpos_[k-1] == 0))
            zShift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = computeShift(zPeak,zPeakWell,zpos_[k]);
          else{
            //well has changed lateral position and we adapt to the new well position
            seisTrace = seisCube[l]->getRealTrace2(ipos_[k],jpos_[k]);
            smoothTrace(seisTrace);
            dzW = timeSimbox->getdz(ipos_[k],jpos_[k]);
            ztopW = timeSimbox->getTop(ipos_[k],jpos_[k]);
            findPeakTrace(seisTrace, zPeakWell, peakWell, bWell, dzW, ztopW);

            seisTrace = seisCube[l]->getRealTrace2(ipos_[k]+i, jpos_[k]+j);
            smoothTrace(seisTrace);
            dz = timeSimbox->getdz(ipos_[k]+i, jpos_[k]+j);
            ztop = timeSimbox->getTop(ipos_[k]+i, jpos_[k]+j);
            findPeakTrace(seisTrace, zPeak, peak, b, dz, ztop);

            peakMatch(zPeak,peak,b,zPeakWell,peakWell,bWell);

            zShift[(i+2) + (j+2)*nZx + k*(nZx*nZx)] = computeShift(zPeak, zPeakWell, zpos_[k]);
          }
        }
      }
    }

    double dx = timeSimbox->getdx();
    double dy = timeSimbox->getdy();

    computeGradient(Qepsilon, Qepsilon_data, zShift, nZx, nZx, dx, dy);
  }

  smoothGradient(xGradient, yGradient, Qepsilon, Qepsilon_data, Sigma_gradient);
  // NBNB Odd slår av estimeringen for å teste om det gir bedre resultat
 /* for(k = 0; k < nBlocks_; k++){
    xGradient[k]=0.0;
    yGradient[k]=0.0;
  }*/

}

void BlockedLogs::smoothTrace(std::vector<float> &trace)
{
  float smoothing_distance = 40; //ms in each direction
  int  L = static_cast<int>(ceil(smoothing_distance/dz_)); //number of lags in the gauss kernel
  float sigma = 10 /dz_; // ms / (ms/cell)

  unsigned int ntrace = static_cast<unsigned int>(trace.size());
  std::vector<float> gk(2*L+1);
  std::vector<float> sTrace(ntrace);
  unsigned int i;
  int j;
  float tmp;
  for(j = -L; j <= L; j++){
    tmp = (j*j)/(2*sigma*sigma);
    gk[j+L] = exp(-tmp);
  }

  float N;
  for(i = 0; i < ntrace; i++){
    N = 0;
    for(j = -L; j <= L; j++){
      if(i+j >= 0 && i+j < ntrace){
        sTrace[i] += gk[j+L]*trace[i+j];
        N += gk[j+L];
      }
    }
    sTrace[i] /= N;
  }

  for(i = 0; i < ntrace; i++)
    trace[i] = sTrace[i];
}



void BlockedLogs::findPeakTrace(std::vector<float> &trace, std::vector<double> &zPeak, std::vector<double> &peak,
                                std::vector<double> &b, double dz, double ztop)
{
  int k;
  double x1, x2, x3, y1, y2, y3, y11, y12, y21;
  int N = static_cast<int>(trace.size());
  zPeak.resize(N); peak.resize(N); b.resize(N);

  double a;
  double c;
  int counter = 0;
  for(k = 1; k < N-1; k++){
    if((trace[k] >= trace[k-1] && trace[k] > trace[k+1]) || (trace[k] <= trace[k-1] && trace[k] < trace[k+1])){
      //Data point for interpolation
      x1 = -dz; x2 = 0; x3 = dz;
      y1 = static_cast<double>(trace[k-1]); y2 = static_cast<double>(trace[k]); y3 = static_cast<double>(trace[k+1]);

      //Newton interpolation method
      y11 = (y2 - y1) / (x2 - x1); y12 = (y3 - y2) / (x3 - x2);
      y21 = (y12 - y11) / (x3 - x1);

      //y = ax + bx^2 + c
      c = y1 - y11*x1 + y21*x1*x2;
      a = y11 - y21*x1 - y21*x2;
      b[counter] = y21;
      //zPeak = -a/2b

      zPeak[counter] = - a/(2.0*b[counter]);

      double tmp = a*zPeak[counter] + b[counter]*zPeak[counter]*zPeak[counter] + c;
      peak[counter] = tmp;
      //Transform back to original z-axis
      zPeak[counter] += ztop + k*dz;
      counter++;

    }
  }
  zPeak.resize(counter); b.resize(counter); peak.resize(counter);

}


void BlockedLogs::peakMatch(std::vector<double> &zPeak, std::vector<double> &peak, std::vector<double> &b,
                                std::vector<double> &zPeakW, std::vector<double> &peakW, std::vector<double> &bW)
{
  //This routine matches the peaks from two traces and returns the set of peak positions that matches.
  unsigned int i, j;
  std::vector<double> pW(zPeak.size());
  std::vector<double> p(zPeak.size());

  double diffz;

  double maxdiffz = 5 * dz_; //matching criteria: Peaks must be no longer that 5 cells apart. (marita: input parameter?)
  double diffp = 0.5; //matcing the size of the peaks NBNB-Frode: This should maybe be an input parameter!

  int counter = 0;
  unsigned int lim = 0;
  for(i = 0; i < zPeakW.size(); i++){
    for(j = lim; j < zPeak.size(); j++){
      diffz = fabs(zPeakW[i] - zPeak[j]);
      if(diffz < maxdiffz){
        //Check if the peaks point in the same direction
        if((bW[i] < 0 && b[j] < 0)||(bW[i] >= 0 && b[j] >= 0)){
          // Check for difference in peak size
          if((fabs(peakW[i] - peak[j]))/(fabs(peakW[i]) + fabs(peak[j])) < diffp){
            pW[counter] = zPeakW[i];
            p[counter] =  zPeak[j];
            counter++;
            lim = j + 1;
          }
        }
      }
    }
  }
  zPeakW.resize(counter); zPeak.resize(counter);
  for(i = 0; i < static_cast<unsigned int>(counter); i++){
    zPeakW[i] = pW[i];
    zPeak[i] =  p[i];
  }

}


double BlockedLogs::computeShift(std::vector<double> &zPeak, std::vector<double> &zPeakW, double z0)
{
  //This routine computes the position of z0 between two peaks in the well and finds the corresponding distance in
  //the other trace. Then zShift is the difference in z between the two.
  unsigned int N = static_cast<unsigned int>(zPeak.size());
  if(N == 0)
    return RMISSING; //The case of no match in the traces
  else{
    unsigned int i;
    int pos = 0;
    double zShift;
    if(z0 < zPeakW[0])
      zShift = zPeakW[0] - zPeak[0];
    else if(z0 >= zPeakW[N-1])
      zShift = zPeakW[N-1] - zPeak[N-1];
    else{
      for(i = 0; i < N-1; i++){
        if(z0 >= zPeakW[i] && z0 < zPeakW[i+1]){
          pos = i;
          i = N;
        }
      }
      zShift = z0 - (zPeak[pos] + (z0 - zPeakW[pos])/(zPeakW[pos+1]-zPeakW[pos])*(zPeak[pos+1]-zPeak[pos]));
    }

    return zShift;
  }
}


void BlockedLogs::computeGradient(std::vector<double> &Qepsilon, std::vector<double> &Qepsilon_data,
                                  std::vector<double> &zShift, int nx, int ny, double dx, double dy)
{
  //This fit the model zshift(x,y)= beta0 + beta1*x + beta2*y  ==> beta1 is x-gradient and beta2 is y-gradient
  int i, j, k, l;
  std::vector<double> Z(3*nx*ny);
  std::vector<double> Y(nx*ny);
  std::vector<double> cov(9);
  std::vector<double> invcov(9);
  std::vector<double> regM(3*nx*ny);

  static bool append = false;

  char* buffer = new char[1000];
  sprintf(buffer,"%s.txt", "C:/Outputfiles/gradNoSmooth");
  std::ofstream out;
  if(append){
    out.open(buffer, std::ios::app|std::ios::out);
  }
  else{
    out.open(buffer);
    append = true;
  }
  delete [] buffer;

  int ndata;
  double data;

  int cy, cz;
  int counter1 = 0;
  int counter2 = 0;
  for(l = 0; l < nBlocks_; l++){
    cy = 0; cz = 0;
    for(j = 0; j < ny; j++){
      for(i = 0; i < nx; i++){
        data = zShift[i + j*nx + l*nx*ny];
        if(data != RMISSING){
          Y[cy] = data;
          Z[cz] = 1.0;
          Z[cz + 1] = (i-(nx-1)/2)*dx;
          Z[cz + 2] = (j-(ny-1)/2)*dy;
          cy++;
          cz += 3;
        }
      }
    }
    Y.resize(cy);
    Z.resize(cz);

    ndata = cy;

    //Compute inverse covariance (ZtZ)^-1 (diagonal matrix for our purpose)
    double tmp;
    for(j = 0; j < 3; j++){
      for(i = 0; i < 3; i++){
        tmp = 0;
        for(k = 0; k < ndata; k++)
          tmp += Z[i + 3*k] * Z[j + 3*k];
        cov[i + 3*j] = tmp;
      }
    }
    double det = cov[0]*(cov[4]*cov[8] - cov[5]*cov[7]) - cov[1]*(cov[3]*cov[8] - cov[5]*cov[6])
                  +   cov[2]*(cov[3]*cov[7] - cov[4]*cov[6]);

    if(det != 0){
      invcov[0] = (cov[4]*cov[8] - cov[5]*cov[7]) / det;
      invcov[1] = (cov[2]*cov[7] - cov[1]*cov[8]) / det;
      invcov[2] = (cov[1]*cov[5] - cov[2]*cov[4]) / det;
      invcov[3] = (cov[5]*cov[6] - cov[3]*cov[8]) / det;
      invcov[4] = (cov[0]*cov[8] - cov[2]*cov[6]) / det;
      invcov[5] = (cov[2]*cov[3] - cov[0]*cov[5]) / det;
      invcov[6] = (cov[3]*cov[7] - cov[4]*cov[6]) / det;
      invcov[7] = (cov[1]*cov[6] - cov[0]*cov[7]) / det;
      invcov[8] = (cov[0]*cov[4] - cov[1]*cov[3]) / det;



      //Compute regression matrix (ZtZ)^-1Zt
      regM.resize(static_cast<unsigned int>(3*ndata));
      for(j = 0; j < 3; j++){
        for(i = 0; i < ndata; i++){
          tmp = 0;
          for(k = 0; k < 3; k++)
            tmp += invcov[k + 3*j]*Z[k + 3*i];
          regM[i + j*ndata] = tmp;
        }
      }


      //Compute beta_1(gradientx) og beta_2(gradienty), beta_0 not necessary
      double beta0 = 0;
      double beta1 = 0;
      double beta2 = 0;
      for(j = 0; j < ndata; j++){
        beta0 += regM[j]*Y[j];
        beta1 += regM[j + ndata]*Y[j];
        beta2 += regM[j + 2*ndata]*Y[j];}
      double sigma2 = 0;
      double sigmatmp;
      for(j = 0; j < ndata; j++){
        sigmatmp = Y[j] - beta0*Z[j*3] - beta1*Z[j*3 + 1] - beta2*Z[j*3 + 2];
        sigma2 += sigmatmp*sigmatmp;
      }

      double qa = sigma2*invcov[4];
      double qb = sigma2*invcov[8];

      out << beta1 << " " << beta2 << " " << qa << " " << qb << std::endl;

      //cov(beta) = sigma2*(ZtZ)^{-1}
      Qepsilon[counter1] += cov[4]/sigma2;
      Qepsilon[counter1+1] += cov[5]/sigma2;
      Qepsilon[counter1+2] += cov[8]/sigma2;
      Qepsilon[counter1+3] += 0;
      counter1 += 4;

      Qepsilon_data[counter2] += (cov[4]* beta1 + cov[5]*beta2) / sigma2;
      Qepsilon_data[counter2+1] += (cov[5]* beta1 + cov[8]*beta2) / sigma2;
      counter2 += 2;

    }
  }
}

void BlockedLogs::smoothGradient(std::vector<double> &xGradient, std::vector<double> &yGradient,
                                 std::vector<double> &Qepsilon, std::vector<double> &Qepsilon_data, std::vector<std::vector<double> > &Sigma_gradient)
{
  int i, j;
  int n_beta = nBlocks_* 2;
  double **Qbeta_data = new double * [n_beta];
  for(i = 0; i < n_beta; i++)
    Qbeta_data[i] = new double[n_beta];


  //Set the prior precicion values
  double a, b, c;
  computePrecisionMatrix(a,b,c);

  //Initialize the Qm|d matrix
  for(i = 0; i < n_beta-2; i++){
    Qbeta_data[i][i] = c + Qepsilon[2*i];
    Qbeta_data[i][i+2] = b;
    if(i % 2 == 0)
      Qbeta_data[i][i+1] = Qepsilon[2*i+1];
    else
      Qbeta_data[i][i+1] = 0;
  }
  //Edge effects
  Qbeta_data[0][0] += a - c;
  Qbeta_data[1][1] += a - c;
  Qbeta_data[n_beta-2][n_beta-2] = c + Qepsilon[2*(n_beta-2)];
  Qbeta_data[n_beta-1][n_beta-1] = c + Qepsilon[2*(n_beta-1)];
  Qbeta_data[n_beta-2][n_beta-1] = Qepsilon[2*(n_beta-2)+1];

  for(i = 0; i < n_beta; i++){
    for(j = i+3; j < n_beta; j++)
      Qbeta_data[i][j] = 0;
  }
  //Transpose and fill in
  for(i = 1; i < n_beta; i++){
    for(j = 0; j < i; j++)
      Qbeta_data[i][j] = Qbeta_data[j][i];
  }

  //Inversion of precision matrix
  double **Identity = new double * [n_beta];
  for(i = 0; i < n_beta; i++)
    Identity[i] = new double[n_beta];
  for(i = 0; i < n_beta; i++){
    for(j = 0; j < n_beta; j++){
      if(i == j)
        Identity[i][j] = 1.0;
      else
        Identity[i][j] = 0;
    }
  }

  //Compute the inverse of Qbeta_data
  //First we do cholesky factorization LL^T = covD
  lib_matrCholR(n_beta, Qbeta_data);

  //Solve Q inverse
  lib_matrAXeqBMatR(n_beta,Qbeta_data,Identity,n_beta);

  //Compute the product (Qbeta_Data)^-1 Qepsilon_data
  double *res = new double[n_beta];
  double tmp;
  for(i = 0; i < n_beta; i++){
    tmp = 0;
    for(j = 0; j < n_beta; j++)
      tmp += Identity[i][j]*Qepsilon_data[j];
    res[i] = tmp;
  }

  //Return Sigma_gradient
  Sigma_gradient.resize(n_beta);
  std::vector<double> tmp_vec(n_beta);
  for(i = 0; i < n_beta; i++){
    for(j = 0; j < n_beta; j++)
      tmp_vec[j] = Identity[i][j];
    Sigma_gradient[i] = tmp_vec;
  }


  int counter = 0;
  for(i = 0; i < nBlocks_; i++){
   xGradient[i] = res[counter];
   yGradient[i] = res[counter+1];
   counter +=2;
  }

   char* buffer2 = new char[1000];
   sprintf(buffer2,"%s.txt", "C:/Outputfiles/gradients");
   std::ofstream out2(buffer2);
   for(i = 0; i < nBlocks_; i++){
       out2 << xGradient[i] << " " << yGradient[i];
     out2 << std::endl;
   }
  delete [] buffer2;

  for(i = 0; i < n_beta; i++){
    if(Qbeta_data[i] != NULL)
      delete [] Qbeta_data[i];
    if(Identity[i] != NULL)
      delete [] Identity[i];
  }
  delete [] Qbeta_data;
  delete [] Identity;
  delete [] res;
}

void BlockedLogs::computePrecisionMatrix(double &a, double &b, double &c)
{
  double minDist = static_cast<double>(lateralThresholdGradient_);
  double sigma_m =  static_cast<double>(sigma_m_);
  double zQ = 1.645; //default 95% confidence interval

  double K = (dz_ /(minDist*zQ))*(dz_ /(minDist*zQ));
  double alpha = 1 - K/(2*sigma_m);

  double alpha_2 = alpha*alpha;
  double gamma_2 = sigma_m*sigma_m*(1-alpha_2);

  a = 1/sigma_m + alpha_2/gamma_2;
  b = - sqrt(alpha_2)/gamma_2;
  c = (alpha_2 + 1)/gamma_2;


}

void BlockedLogs::generateSyntheticSeismic(float   ** reflCoef,
                                           int        nAngles,
                                           Wavelet ** wavelet,
                                           int        nz,
                                           int        nzp,
                                           const Simbox * timeSimbox)
{
  int          i,j;
  int          start,length;

  fftw_complex   cAmp;
  Wavelet1D     * localWavelet;

  int    cnzp = nzp/2+1;
  int    rnzp = 2*cnzp;

  float  * syntSeis  = new float[nz];
  float  * alphaVert = new float[nLayers_];
  float  * betaVert  = new float[nLayers_];
  float  * rhoVert   = new float[nLayers_];

  fftw_real    * cpp_r       = new fftw_real[rnzp];
  fftw_complex * cpp_c       = reinterpret_cast<fftw_complex*>(cpp_r);
  fftw_real    * synt_seis_r = new fftw_real[rnzp];
  fftw_complex * synt_seis_c = reinterpret_cast<fftw_complex*>(synt_seis_r);

  getVerticalTrend(alpha_, alphaVert);
  getVerticalTrend(beta_, betaVert);
  getVerticalTrend(rho_, rhoVert);

  std::vector<bool> hasData(nLayers_);
  for( i=0 ; i<nLayers_; i++ )
    hasData[i] = alphaVert[i] != RMISSING && betaVert[i] != RMISSING && rhoVert[i] != RMISSING;

  findContiniousPartOfData(hasData,nLayers_,start,length);

  float scale = static_cast<float>(timeSimbox->getRelThick(ipos_[0], jpos_[0]));

  for( i=0; i<nAngles; i++ )
  {
    for(j=0; j<rnzp; j++)
    {
      cpp_r[j] = 0;
      synt_seis_r[j] = 0;
    }
    fillInCpp(reflCoef[i],start,length,cpp_r,nzp);
    Utils::fft(cpp_r,cpp_c,nzp);

    localWavelet = wavelet[i]->createLocalWavelet1D(ipos_[0],jpos_[0]);
    localWavelet->fft1DInPlace();
   // float sf = wavelet[i]->getLocalStretch(ipos_[0],jpos_[0]);
   // what about relative thickness ????
   // float sf = wavelet[i]->getLocalStretch(ipos_[0],jpos_[0])*Relativethikness... Need simbox;

    for( j=0; j<cnzp; j++ )
    {
      cAmp =  localWavelet->getCAmp(j, scale);
      synt_seis_c[j].re = cpp_c[j].re*cAmp.re + cpp_c[j].im*cAmp.im;
      synt_seis_c[j].im = cpp_c[j].im*cAmp.re - cpp_c[j].re*cAmp.im;
    }

    Utils::fftInv(synt_seis_c,synt_seis_r,nzp);

    for ( j=0 ; j<nz ; j++ )
      syntSeis[j] = 0.0; // Do not use RMISSING (fails in setLogFromVerticalTrend())

    for ( j=start; j<start+length; j++ )
      syntSeis[j] = synt_seis_r[j];

    setLogFromVerticalTrend(syntSeis,zpos_[0],dz_,nz,"ACTUAL_SYNTHETIC_SEISMIC",i);

    //localWavelet->fft1DInPlace();
    delete localWavelet;
  }

  delete [] syntSeis;
  delete [] alphaVert;
  delete [] betaVert;
  delete [] rhoVert;
  delete [] cpp_r;
  delete [] synt_seis_r;

}

const std::vector<int>
BlockedLogs::getIposVector() const
{
  std::vector<int> ipos(nBlocks_);
  const int * ipos_pt = getIpos();
  for (int i=0; i<nBlocks_; i++)
    ipos[i] = ipos_pt[i];
  return (ipos);
}

const std::vector<int>
BlockedLogs::getJposVector() const
{
  std::vector<int> jpos(nBlocks_);
  const int * jpos_pt = getJpos();
  for (int i=0; i<nBlocks_; i++)
    jpos[i] = jpos_pt[i];
  return (jpos);
}

const std::vector<double>
BlockedLogs::getXposVector() const
{
  std::vector<double> xpos(nBlocks_);
  const double * xpos_pt = getXpos();
  for (int i=0; i<nBlocks_; i++)
    xpos[i] = xpos_pt[i];
  return (xpos);
}

const std::vector<double>
BlockedLogs::getYposVector() const
{
  std::vector<double> ypos(nBlocks_);
  const double * ypos_pt = getYpos();
  for (int i=0; i<nBlocks_; i++)
    ypos[i] = ypos_pt[i];
  return (ypos);
}

const std::vector<double>
BlockedLogs::getZposVector() const
{
  std::vector<double> zpos(nBlocks_);
  const double * zpos_pt = getZpos();
  for (int i=0; i<nBlocks_; i++)
    zpos[i] = zpos_pt[i];
  return (zpos);
}
