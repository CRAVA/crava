#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

#include "lib/kriging1d.h"
#include "lib/utils.h"

#include "nrlib/iotools/logkit.hpp"

#include "src/definitions.h"
#include "src/modelgeneral.h"
#include "src/modelsettings.h"
#include "src/vario.h"
#include "src/simbox.h"
#include "src/welldata.h"
#include "src/background.h"
#include "src/covgrid2d.h"
#include "src/krigingdata2d.h"
#include "src/kriging2d.h"
#include "src/krigingdata3d.h"
#include "src/covgridseparated.h"
#include "src/fftgrid.h"
#include "src/fftfilegrid.h"
#include "src/gridmapping.h"
#include "src/io.h"

Background::Background(FFTGrid       ** grids,
                       WellData      ** wells,
                       FFTGrid       *& velocity,
                       Simbox         * timeSimbox,
                       Simbox         * timeBGSimbox,
                       ModelSettings  * modelSettings)
  : DataTarget_(250), // For kriging: Increase surrounding until 250 data points is aquired
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];

  FFTGrid * bgAlpha;
  FFTGrid * bgBeta;
  FFTGrid * bgRho;
  
  if (timeBGSimbox == NULL)
  {
    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
                            velocity,wells,
                            timeSimbox,
                            modelSettings);
  }
  else
  {
    generateBackgroundModel(bgAlpha,bgBeta,bgRho,
                            velocity,wells,
                            timeBGSimbox,
                            modelSettings);
    resampleBackgroundModel(bgAlpha,bgBeta,bgRho,
                            timeBGSimbox,
                            timeSimbox,
                            modelSettings);
  }  
  padAndSetBackgroundModel(bgAlpha,bgBeta,bgRho);

  delete bgAlpha;
  delete bgBeta;
  delete bgRho;

  findMeanVsVp(backModel_[0],
               backModel_[1]);
}

//-------------------------------------------------------------------------------
Background::Background(FFTGrid ** grids) 
  : DataTarget_(IMISSING),
    vsvp_(RMISSING)
{
  for(int i=0 ; i<3 ; i++)
    backModel_[i] = grids[i];
  findMeanVsVp(backModel_[0],
               backModel_[1]);
}

//-------------------------------------------------------------------------------
Background::~Background(void)
{
  for (int i=0 ; i<3 ; i++)
    delete backModel_[i];
}

//-------------------------------------------------------------------------------
void
Background::releaseGrids()
{
  for (int i=0 ; i<3 ; i++)
    backModel_[i] = NULL;
}

//-------------------------------------------------------------------------------
void
Background::generateBackgroundModel(FFTGrid      *& bgAlpha,
                                    FFTGrid      *& bgBeta,
                                    FFTGrid      *& bgRho,
                                    FFTGrid      *& velocity,
                                    WellData     ** wells,
                                    Simbox        * simbox,
                                    ModelSettings * modelSettings)
{
  const int nz = simbox->getnz();
  float * trendAlpha = new float[nz];
  float * trendBeta  = new float[nz];
  float * trendRho   = new float[nz];
  float * trendVel   = new float[nz]; // Allocate (for simplicity) although not always needed

  const int nWells = modelSettings->getNumberOfWells();
  float * avgDevAlpha = new float[nWells];
  float * avgDevBeta  = new float[nWells];
  float * avgDevRho   = new float[nWells];
  float * avgDevVel   = new float[nWells]; // Allocate (for simplicity) although not always needed

  bool hasVelocityTrend = velocity != NULL;
  bool write1D = ((modelSettings->getOtherOutputFlag()& IO::BACKGROUND_TREND_1D) > 0);
  bool write3D = ((modelSettings->getOutputGridsElastic() & IO::BACKGROUND_TREND) > 0);
  calculateBackgroundTrend(trendAlpha, avgDevAlpha,
                           wells, simbox, 
                           modelSettings->getAlphaMin(), 
                           modelSettings->getAlphaMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Vp"), 
                           modelSettings->getFileGrid());
  calculateBackgroundTrend(trendBeta, avgDevBeta, 
                           wells, simbox, 
                           modelSettings->getBetaMin(), 
                           modelSettings->getBetaMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Vs"),
                           modelSettings->getFileGrid());
  calculateBackgroundTrend(trendRho, avgDevRho,
                           wells, simbox, 
                           modelSettings->getRhoMin(), 
                           modelSettings->getRhoMax(),
                           modelSettings->getMaxHzBackground(), 
                           write1D, write3D,
                           nWells, hasVelocityTrend,
                           std::string("Rho"), 
                           modelSettings->getFileGrid());

  if (velocity != NULL) {
    //
    // We still want calculateBackgroundTrend() for alpha above. By calculating
    // avgDevAlpha we can check that the bgAlpha calculated from velocity is as 
    // good as or better than that calculated by crava.
    //
    calculateVelocityDeviations(velocity, wells, simbox, 
                                trendVel, avgDevVel, avgDevAlpha,
                                modelSettings->getOutputGridsElastic(),
                                nWells);
    velocity->logTransf();
    delete bgAlpha;
    bgAlpha = velocity;
    velocity = NULL;
    writeDeviationsFromVerticalTrend(avgDevVel, avgDevBeta, avgDevRho, 
                                     trendVel, trendBeta, trendRho, 
                                     wells, nWells, nz);
  }
  else {
    writeDeviationsFromVerticalTrend(avgDevAlpha, avgDevBeta, avgDevRho,
                                     trendAlpha, trendBeta, trendRho,
                                     wells, nWells, nz);
  }

  std::vector<KrigingData2D> krigingDataAlpha(nz);
  std::vector<KrigingData2D> krigingDataBeta(nz);
  std::vector<KrigingData2D> krigingDataRho(nz);
  
  setupKrigingData2D(krigingDataAlpha,krigingDataBeta,krigingDataRho,
                     trendAlpha,trendBeta,trendRho,
                     modelSettings->getOutputGridsElastic(),
                     simbox,wells,nWells);
  
  const CovGrid2D & covGrid2D = makeCovGrid2D(simbox,
                                              modelSettings->getBackgroundVario(),
                                              modelSettings->getDebugFlag());
  
  makeKrigedBackground(krigingDataAlpha, bgAlpha, trendAlpha, simbox, covGrid2D, "Vp" , modelSettings->getFileGrid());
  makeKrigedBackground(krigingDataBeta , bgBeta , trendBeta , simbox, covGrid2D, "Vs" , modelSettings->getFileGrid());
  makeKrigedBackground(krigingDataRho  , bgRho  , trendRho  , simbox, covGrid2D, "Rho", modelSettings->getFileGrid());

  delete &covGrid2D;

  delete [] avgDevAlpha;
  delete [] avgDevBeta;
  delete [] avgDevRho;
  delete [] avgDevVel;

  delete [] trendAlpha;
  delete [] trendBeta;
  delete [] trendRho;
  delete [] trendVel;
}
 
//---------------------------------------------------------------------------
void
Background::calculateVelocityDeviations(FFTGrid   * velocity,                           
                                        WellData ** wells,
                                        Simbox    * simbox,
                                        float    *& trendVel,
                                        float    *& avgDevVel,
                                        float     * avgDevAlpha,
                                        int         outputFlag,
                                        int         nWells)
{
  if((outputFlag & IO::BACKGROUND_TREND) > 0) {
    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + "VpFromFile";
    velocity->writeFile(fileName, IO::PathToBackground(), simbox, "NO_LABEL");
  }
    
  //
  // Calculate deviation between well data and trend
  //
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsOrigThick()->getNumberOfBlocks();
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }
  float * velocityLog = new float[maxBlocks];

  const int nz = simbox->getnz();
  float * vtAlpha    = new float[nz];
  float * vtVelocity = new float[nz];

  for (int k=0 ; k<nz ; k++)
    trendVel[k]=0.0;

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    const float * alphaLog = bl->getAlphaHighCutBackground();
    bl->getVerticalTrend(alphaLog, vtAlpha);
    bl->getBlockedGrid(velocity, velocityLog);
    bl->getVerticalTrend(velocityLog, vtVelocity);
    float sumDev = 0.0f;
    int count = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (vtAlpha[k] != RMISSING) {
        trendVel[k] += vtVelocity[k]; 
        float diff = exp(vtAlpha[k]) - vtVelocity[k]; // Velocity trend is in exp-domain
        sumDev += diff*diff;
        count++;
      }
    }
    if (count > 0) 
      sumDev /= count;
    avgDevVel[w] = sqrt(sumDev);
  }
  delete [] vtVelocity;
  delete [] vtAlpha;
  delete [] velocityLog;

  for (int k=0 ; k<nz ; k++)
    trendVel[k] /= nWells;

  LogKit::LogFormatted(LogKit::Low,"\nAverage deviations of type well-log-Vp-minus-velocity-read-from-file and ");
  LogKit::LogFormatted(LogKit::Low,"\nwell-log-Vp-minus-estimated-Vp-trend (added for quality control):\n\n");
  LogKit::LogFormatted(LogKit::Low,"Well             TrendFromFile  TrendFromData\n");
  LogKit::LogFormatted(LogKit::Low,"---------------------------------------------\n");
  for (int i=0 ; i<nWells ; i++)
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f          %5.1f\n",
                         wells[i]->getWellname().c_str(),avgDevVel[i],avgDevAlpha[i]);
}

//---------------------------------------------------------------------------
void
Background::calculateBackgroundTrend(float             * trend,
                                     float             * avgDev,
                                     WellData         ** wells,
                                     Simbox            * simbox,
                                     float               logMin, 
                                     float               logMax,
                                     float               maxHz, 
                                     bool                write1D,
                                     bool                write3D,
                                     int                 nWells, 
                                     bool                hasVelocityTrend,
                                     const std::string & name,
                                     bool                isFile)
{
  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  calculateVerticalTrend(wells, trend, 
                         logMin, logMax, 
                         maxHz, nWells, 
                         nz, dz, name);
  
  if(write1D == true) {
    writeVerticalTrend(trend, dz, nz, name);
  }
  
  calculateDeviationFromVerticalTrend(wells, trend, avgDev,
                                      nWells, nz, name);
  
  if(write3D == true && !(name=="Vp" && hasVelocityTrend))
  {
    const int nx = simbox->getnx();
    const int ny = simbox->getny();
    FFTGrid * trendGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nx, ny, nz, isFile);
    fillInVerticalTrend(trendGrid, trend);
    FFTGrid * expTrend = copyFFTGrid(trendGrid, true, isFile);
    delete trendGrid;
    std::string fileName = IO::PrefixBackground() + IO::PrefixTrend() + name;
    expTrend->writeFile(fileName, IO::PathToBackground(), simbox);
    delete expTrend;
  }
}

//-------------------------------------------------------------------------------
void
Background::setupKrigingData2D(std::vector<KrigingData2D> & krigingDataAlpha,
                               std::vector<KrigingData2D> & krigingDataBeta,
                               std::vector<KrigingData2D> & krigingDataRho,
                               float                      * trendAlpha,
                               float                      * trendBeta, 
                               float                      * trendRho , 
                               int                          outputFlag,
                               Simbox                     * simbox,
                               WellData                  ** wells,
                               const int                    nWells)
{
  //
  // Although unnecessary, we have chosen to set up kriging data fro 
  // Vp, Vs and Rho simultaneously. This gives code easier to read.
  //
  int totBlocks = 0;
  int maxBlocks = 0;
  for (int w = 0 ; w < nWells ; w++) {
    int nBlocks = wells[w]->getBlockedLogsExtendedBG()->getNumberOfBlocks();
    totBlocks += nBlocks;
    if (nBlocks > maxBlocks)
      maxBlocks = nBlocks;
  }  
  
  KrigingData3D forLogging(totBlocks); 
  
  const int   nz = simbox->getnz();
  const float dz = static_cast<float>(simbox->getdz()*simbox->getAvgRelThick());

  float * blAlpha = new float[maxBlocks];   // bl = blocked logs
  float * blBeta  = new float[maxBlocks];
  float * blRho   = new float[maxBlocks];

  float * vtAlpha = new float[nz];          // vt = vertical trend
  float * vtBeta  = new float[nz];
  float * vtRho   = new float[nz];

  for (int w = 0 ; w < nWells ; w++)
  {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    const int nBlocks = bl->getNumberOfBlocks();

    Utils::copyVector(bl->getAlphaHighCutBackground(), blAlpha, nBlocks);
    Utils::copyVector(bl->getBetaHighCutBackground(), blBeta, nBlocks);
    Utils::copyVector(bl->getRhoHighCutBackground(), blRho, nBlocks);
    //
    // Extract a one-value-for-each-layer array of blocked logs
    //
    bl->getVerticalTrend(blAlpha, vtAlpha);
    bl->getVerticalTrend(blBeta, vtBeta);
    bl->getVerticalTrend(blRho, vtRho);
    //
    // Kriging vertical trend (vt....) against global vertical trend (trend...)
    //
    Kriging1D::krigVector(vtAlpha, trendAlpha, nz, dz);
    Kriging1D::krigVector(vtBeta, trendBeta, nz, dz);
    Kriging1D::krigVector(vtRho, trendRho, nz, dz);
    //
    // Use kriged vertical trend where original log is not defined.
    //
    const int * ipos = bl->getIpos();
    const int * jpos = bl->getJpos();
    const int * kpos = bl->getKpos();

    for (int m = 0 ; m < nBlocks ; m++) 
    {
      int i = ipos[m];
      int j = jpos[m];
      int k = kpos[m];

      if (blAlpha[m] == RMISSING) 
      {
        blAlpha[m] = vtAlpha[k];
      }
      if (blBeta[m] == RMISSING) 
      {
        blBeta[m] = vtBeta[k];
      }
      if (blRho[m] == RMISSING) 
      {
        blRho[m] = vtRho[k];
      }
      krigingDataAlpha[k].addData(i, j, blAlpha[m]);
      krigingDataBeta[k].addData(i, j, blBeta[m]);
      krigingDataRho[k].addData(i, j, blRho[m]);
    }

    forLogging.addData(blAlpha, blBeta, blRho,
                       ipos, jpos, kpos,
                       nBlocks);
  }

  for (int k=0 ; k<nz ; k++) 
  {
    krigingDataAlpha[k].findMeanValues();
    krigingDataBeta[k].findMeanValues();
    krigingDataRho[k].findMeanValues();
  }

  if((outputFlag & IO::BACKGROUND) > 0) {
    forLogging.divide();
    std::string baseName = IO::PrefixBackground() + IO::PrefixKrigingData() + IO::SuffixGeneralData();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    forLogging.writeToFile(fileName);
  }
  
  delete [] vtAlpha;
  delete [] vtBeta;
  delete [] vtRho;

  delete [] blAlpha;
  delete [] blBeta;
  delete [] blRho;
}

//---------------------------------------------------------------------------
const CovGrid2D &
Background::makeCovGrid2D(Simbox * simbox,
                          Vario  * vario, 
                          int      debugFlag)
{
  //
  // Pretabulate all needed covariances
  //
  const int    nx = simbox->getnx();
  const int    ny = simbox->getny();
  
  const float  dx = static_cast<float>(simbox->getdx());
  const float  dy = static_cast<float>(simbox->getdy());
  
  CovGrid2D * cov = new CovGrid2D(vario, nx, ny, dx, dy);
  
  if(debugFlag == 1) {
    std::string baseName = IO::PrefixBackground() + "covGrid2D" + IO::SuffixAsciiIrapClassic();
    std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
    cov->writeToFile(fileName);
  }
  return (*cov); 
}

//---------------------------------------------------------------------------
void
Background::makeKrigedBackground(const std::vector<KrigingData2D> & krigingData,
                                 FFTGrid                         *& bgGrid,
                                 float                            * trend,
                                 Simbox                           * simbox,
                                 const CovGrid2D                  & covGrid2D,
                                 const std::string                & type,
                                 bool                               isFile)
{
  std::string text = "\nBuilding "+type+" background:";
  LogKit::LogFormatted(LogKit::Low,text);

  const int    nx   = simbox->getnx();
  const int    ny   = simbox->getny();
  const int    nz   = simbox->getnz();

  const int    nxp  = nx;
  const int    nyp  = ny;
  const int    nzp  = nz;
  const int    rnxp = 2*(nxp/2 + 1);

  const double x0   = simbox->getx0();
  const double y0   = simbox->gety0();
  const double lx   = simbox->getlx();
  const double ly   = simbox->getly();

  //
  // Template surface to be kriged
  //
  Surface surface(x0, y0, lx, ly, nx, ny, RMISSING);

  float monitorSize = std::max(1.0f, static_cast<float>(nz)*0.02f);
  float nextMonitor = monitorSize;
  std::cout 
    << "\n  0%       20%       40%       60%       80%      100%"
    << "\n  |    |    |    |    |    |    |    |    |    |    |  "
    << "\n  ^"; 

  bgGrid = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  bgGrid->createRealGrid();
  bgGrid->setType(FFTGrid::PARAMETER);
  bgGrid->setAccessMode(FFTGrid::WRITE);

  for (int k=0 ; k<nzp ; k++)
  {
    // Set trend for layer
    surface.Assign(trend[k]);                                    

    // Kriging of layer
    Kriging2D::krigSurface(surface, krigingData[k], covGrid2D);  

    // Set layer in background model from surface
    for(int j=0 ; j<nyp ; j++) {        
      for(int i=0 ; i<rnxp ; i++) {
        if(i<nxp)
          bgGrid->setNextReal(float(surface(i,j)));
        else
          bgGrid->setNextReal(0);  //dummy in padding (but there is no padding)
      }
    }
    
    // Log progress
    if (k+1 >= static_cast<int>(nextMonitor)) 
    { 
      nextMonitor += monitorSize;
      std::cout << "^";
      fflush(stdout);
    }
  }
  bgGrid->endAccess();
}

//-------------------------------------------------------------------------------
void 
Background::calculateVerticalTrend(WellData         ** wells,
                                   float             * trend, 
                                   float               logMin,
                                   float               logMax,
                                   float               maxHz,
                                   int                 nWells,
                                   int                 nz,
                                   float               dz,
                                   const std::string & name)
{  
  float * filtered_log = new float[nz]; 
  float * wellTrend    = filtered_log;   // Use same memory twice
  int   * count        = new int[nz];
  //
  // Calculate the average values of well log
  // ----------------------------------------
  // When calculating the vertical trend, we do not want each well
  // to contribute with more than one value to each layer 'k',
  // and therefore we calculate average values for each well 
  // first. This way we avoid strange behaviour caused by
  // deviated/horizontal wells.
  //
  for (int k = 0 ; k < nz ; k++) {
    trend[k] = 0.0f;
    count[k] = 0;
  }
  int iWells = 0;
  for (int w = 0 ; w < nWells ; w++) { 
    if (wells[w]->getUseForBackgroundTrend()) {
      BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
      if (name == "Vp")
        bl->getVerticalTrend(bl->getAlpha(), wellTrend);
      else if (name == "Vs")
        bl->getVerticalTrend(bl->getBeta(), wellTrend);
      else if (name == "Rho")
        bl->getVerticalTrend(bl->getRho(), wellTrend);
      else {
        LogKit::LogFormatted(LogKit::Low,"ERROR in Background::calculateVerticalTrend(): ");
        LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
        exit(1);
      }
      for (int k = 0 ; k < nz ; k++) {
        if (wellTrend[k] != RMISSING) {
          trend[k] += exp(wellTrend[k]);
          count[k]++;
        }
      }
      iWells++;
    }
  }
  if (iWells > 0) {
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        trend[k] = trend[k]/count[k]; 
      }
    }
  }
  else {
    LogKit::LogFormatted(LogKit::Low,"\nERROR in Background::calculateVerticalTrend(): There are no wells\n");
    LogKit::LogFormatted(LogKit::Low,"available for the estimation of background trend.\n");
    exit(1);
  }

  //Utils::writeVectorToFile(std::string("trend_mean_values_") + name, trend, nz);

  smoothTrendWithLocalLinearRegression(trend, count,
                                       iWells, nz, dz, 
                                       logMin,
                                       logMax,
                                       name);
  
  //Utils::writeVectorToFile(std::string("trend_after_linreg_") + name, trend, nz);

  WellData::applyFilter(filtered_log, 
                        trend, 
                        nz, 
                        dz, 
                        maxHz);
  
  for (int i=0 ; i<nz ; i++) {
    trend[i] = filtered_log[i];
  }

  //Utils::writeVectorToFile(std::string("trend_after_filter_") + name, trend, nz);

  
  delete [] filtered_log; 
  delete [] count;
}


//-------------------------------------------------------------------------------
void
Background::smoothTrendWithLocalLinearRegression(float      * trend, 
                                                 int        * count,
                                                 int          iWells,
                                                 int          nz,
                                                 float        dz,
                                                 float        min_value, 
                                                 float        max_value,
                                                 std::string  parName) 
{
  bool debug = false;
  //
  // 1. Center-parts of scatter plots
  //
  // In the center parts of the scatter plots the average value should be
  // accepted as a trend value if the number of data points behind each 
  // average is fraction * iWells, where fraction is the acceptance
  // fraction, typically larger than one.
  //
  // Sometimes we have only one, two, or three wells available. In the
  // case of two and three wells, the logs for these may differ considerably,
  // in which case the trend need to be stabilised by requiring a minimum
  // number of data points behind each trend value. The denser the log is 
  // sampled, the more prone it is to numerical instabilities. This 
  // minimum value is therefore linked to the sampling density 'dz'.
  //
  // Finally, we must require a definite minimum number of data points 
  // that should be behind each trend value.
  // 
  // nDataMin = max(min_req_points, max(min_time_sample, fraction * iWells))
  //
  // 2. End-parts of scatter plot
  //
  // If the blocked wells do not contain any values for the lower layers of 
  // the simbox, the end part of the scatter plot needs special attention.
  //
  // We must possibly require a larger min_req_points to avoid an 
  // "arbitrary" end behaviour.
  //

  float fraction   = 5.0f;                      // Require minimum 5*iWells
  int   nTimeLimit = static_cast<int>(50.0/dz); // The smaller sampling density, the more values are needed.
  int   nLowLimit  = 10;                        // Require minimum 10
  int   nDataMin   = std::max(nLowLimit, std::max(nTimeLimit, int(fraction * iWells)));

  bool  use_weights = true;
  bool  errorMid    = false;
  bool  errorHead   = false;
  bool  errorTrail  = false;

  //
  // Copy the average values (stored in array 'trend') to the array 'mean'.
  //
  float * mean = new float[nz];
  for (int k = 0 ; k < nz ; k++) {
    mean[k] = trend[k];
  }
  
  //
  // Find first non-missing value
  //
  int firstNonmissing = 0;
  for (int k = 0 ; k < nz ; k++) {
    if (trend[k] > 0.0f) {
      firstNonmissing = k;
      break;
    }
  }
  
  //
  // Find last non-missing value
  //
  int lastNonmissing = nz - 1;
  for (int k = nz - 1 ; k > 0 ; k--) {
    if (trend[k] > 0.0f) {
      lastNonmissing = k;
      break;
    }
  }
  
  float * x = new float[nz];  // Time indices
  float * y = new float[nz];  // Log values
  float * w = new float[nz];  // Weights (number of data behind each avg.)

  for (int k = 0 ; k < nz ; k++) {
    int nCurDataMin = nDataMin;
    if (k < firstNonmissing || k > lastNonmissing) {
      nCurDataMin *= 2;
    }

    int n = 0;
    int nData = 0;
    if (debug) 
      LogKit::LogFormatted(LogKit::Low,"k=%d\n",k);
    //
    // 1. Add current data point to local data set if present.
    //
    if (count[k] > 0) {
      w[0]   = static_cast<float>(count[k]);
      x[0]   = static_cast<float>(k);
      y[0]   = trend[k];
      nData += count[k];
      if (debug) 
        LogKit::LogFormatted(LogKit::Low,"   A:t=%.2f   x[0] y[0]  %d   %.2f\n",dz*(x[0] + 0.5f),int(x[0]),y[0]);
      n++;
    }

    //
    // 2. Add local data points to get 'nCurDataMin' points behind each trend 
    //    value. Note that the bandwidth varies
    //
    int i = 0;
    while (nData < nCurDataMin) {
      i++;
      if (k - i >= 0 && count[k - i] > 0) {
        w[n]   = static_cast<float>(count[k - i]);
        x[n]   = static_cast<float>(k - i);
        y[n]   = mean [k - i];
        nData += count[k - i];
        if (debug) 
          LogKit::LogFormatted(LogKit::Low,"   B:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k + i < nz  && count[k + i] > 0) {
        w[n]   = static_cast<float>(count[k + i]);
        x[n]   = static_cast<float>(k + i);
        y[n]   = mean [k + i];
        nData += count[k + i];
        if (debug) 
          LogKit::LogFormatted(LogKit::Low,"   C:t=%.2f   x[%d] y[%d]  %d   %.2f\n",dz*(x[n] + 0.5f),n,n,int(x[n]),y[n]);
        n++;
      }
      if (k-i < 0 && k+i >= nz) { // We will never find enough data
        break;
      }
    }

    //
    // Calculate normalised weights
    //
    if (use_weights) 
      for (i = 0 ; i < n ; i++)
        w[i] /= nData;
    else
      for (i = 0 ; i < n ; i++) 
        w[i] = 1.0f/n;
    
    //
    // We need at least two points to make a line.
    //
    if (n > 1) {
      //
      // Estimate local regression line: y = bx + a
      //
      float Sx  = x[0]*w[0];
      float Sy  = y[0]*w[0];
      float Sxx = x[0]*w[0]*x[0];          
      float Sxy = x[0]*w[0]*y[0];
      for (i = 1 ; i < n ; i++) {
        Sx  += x[i]*w[i];
        Sy  += y[i]*w[i];
        Sxx += x[i]*w[i]*x[i];
        Sxy += x[i]*w[i]*y[i];
      }
      float b = (Sxy - Sx*Sy)/(Sxx - Sx*Sx);
      float a = (Sy - b*Sx);
      if (debug) 
        LogKit::LogFormatted(LogKit::Low,"Sx, Sy, Sxx, Sxy  : %.4f, %.4f, %.4f, %.4f     a, b : %.4f %.4f\n",Sx,Sy,Sxx,Sxy,a,b);
      
      //
      // Estimate value of regression line at requested point.
      //
      float value = a + b*static_cast<float>(k);
      if (value < min_value || value > max_value) {
        if (debug)
          LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f\n",value);
        if (k < firstNonmissing) 
          errorHead = true;
        else if (k > lastNonmissing)
          errorTrail = true;
        else {
          errorMid   = true;
          break;
        }
      }
      trend[k] = log(value);
    }
    else {
      trend[k] = log(y[0]);
    }
    if (debug) 
      LogKit::LogFormatted(LogKit::Low,"   TREND: trend[k] = %.2f        (minLog/maxLog = %.2f / %.2f)\n",exp(trend[k]),min_value,max_value);
  }

  if (errorMid) {
    // Big problem ...
    LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
    LogKit::LogFormatted(LogKit::Low,"          regression failed - trying global mean instead. Possible causes: \n");
    LogKit::LogFormatted(LogKit::Low,"          1) Available logs cover too small a part of inversion grid giving extrapolation problems.\n");
    LogKit::LogFormatted(LogKit::Low,"          2) There are too many layers in grid compared to well logs available.\n");
    float sum = 0.0f;
    int nData = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (count[k] > 0) {
        if (use_weights) {
          sum   += mean[k]*count[k];
          nData += count[k];
          if (debug)           
            LogKit::LogFormatted(LogKit::Low,"k=%d  count[k], mean[k]  nData, sum  %d  %8.3f     %d  %8.3f\n",
                                 k,count[k],mean[k],nData,sum);
        }
        else {
          sum += mean[k];
          nData += 1;
        }
      }
    }
    float global_mean = sum/nData;
    for (int k = 0 ; k < nz ; k++) {
      trend[k] = log(global_mean);
      if (debug) 
        LogKit::LogFormatted(LogKit::Low,"   TREND: k = %d   trend[k] = %.2f\n",k,exp(trend[k]));
    }
    LogKit::LogFormatted(LogKit::Low,"\nGlobal mean for parameter %s = %.2f\n\n",parName.c_str(),global_mean);
  }
  else {
    if (errorHead) {
      // Fix first part of trend containing missing-values.
      float firstValue = trend[firstNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter \'"+parName+"\' using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [0-%d] where the log is undefined. The first\n",firstNonmissing-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d we be used throughout this region.\n",exp(firstValue),firstNonmissing);
      for (int k = 0 ; k < firstNonmissing ; k++) {
        trend[k] = firstValue;
      }
    }
    if (errorTrail) {
      // Fix last part of trend containing missing-values.
      float lastValue = trend[lastNonmissing];
      LogKit::LogFormatted(LogKit::Low,"\nWARNING : The calculation of the vertical trend for parameter "+parName+" using local linear\n");
      LogKit::LogFormatted(LogKit::Low,"          regression failed for cells [%d,%d] where the log is undefined. The last\n",lastNonmissing+1,nz-1);
      LogKit::LogFormatted(LogKit::Low,"          defined value of %.2f in cell %d we be used throughout this region.\n",exp(lastValue),lastNonmissing);
      for (int k = lastNonmissing + 1 ; k < nz ; k++) {
        trend[k] = lastValue;
      }
    }
  }
  delete [] x;
  delete [] y;
  delete [] w;
  delete [] mean;
}
 
//-------------------------------------------------------------------------------
void 
Background::writeVerticalTrend(float      * trend, 
                               float        dz,
                               int          nz,
                               std::string  name) 
{  
  float z0 = dz/2.0f;
  std::string baseName = IO::PrefixBackground() + IO::PrefixTrend() + name + IO::SuffixAsciiIrapClassic();
  std::string fileName = IO::makeFullFileName(IO::PathToBackground(), baseName);
  std::ofstream file;
  NRLib::OpenWrite(file, fileName);
  for (int i=0 ; i<nz ; i++) {
    file << std::fixed 
         << std::setprecision(2)
         << std::setw(8) << (z0 + i*dz)     << " "
         << std::setprecision(3)
         << std::setw(8) << exp( trend[i] ) << " "
         << "0.00\n";
  }
  file.close();
}

//-------------------------------------------------------------------------------
void
Background::calculateDeviationFromVerticalTrend(WellData    ** wells,
                                                const float  * globalTrend, 
                                                float        * avg_dev,
                                                int            nWells,
                                                int            nz,
                                                std::string    name)
{
  float * wellTrend = new float[nz];

  for (int w = 0 ; w < nWells ; w++) {
    BlockedLogs * bl = wells[w]->getBlockedLogsExtendedBG();
    if (name == "Vp")
      bl->getVerticalTrend(bl->getAlphaHighCutBackground(), wellTrend);
    else if (name == "Vs")
      bl->getVerticalTrend(bl->getBetaHighCutBackground(), wellTrend);
    else if (name == "Rho")
      bl->getVerticalTrend(bl->getRhoHighCutBackground(), wellTrend);
    else {
      LogKit::LogFormatted(LogKit::Low,"ERROR in Background::calculateVerticalTrend(): ");
      LogKit::LogFormatted(LogKit::Low,"Log \'"+name+"\' requested, but no such log exists.\n");
      exit(1);
    }
    float sum_dev = 0.0f;
    int count = 0;
    for (int k = 0 ; k < nz ; k++) {
      if (wellTrend[k] != RMISSING) {
        float diff = exp(wellTrend[k]) - exp(globalTrend[k]);
        sum_dev += diff*diff;
        count++;
      }
    }
    if (count > 0) 
      sum_dev /= count;
    avg_dev[w] = sqrt(sum_dev);
  }
  delete [] wellTrend;
}

//-------------------------------------------------------------------------------
void
Background::writeDeviationsFromVerticalTrend(const float *  avg_dev_alpha,
                                             const float *  avg_dev_beta,
                                             const float *  avg_dev_rho,
                                             const float *  trend_alpha,
                                             const float *  trend_beta,
                                             const float *  trend_rho,
                                             WellData    ** wells,
                                             const int      nWells,
                                             const int      nz)
{
  float global_mean_alpha = 0.0f;
  float global_mean_beta  = 0.0f;
  float global_mean_rho   = 0.0f;

  for (int k=0 ; k<nz ; k++)
  {
    global_mean_alpha += exp(trend_alpha[k]);
    global_mean_beta  += exp(trend_beta[k]);
    global_mean_rho   += exp(trend_rho[k]);
  }
  global_mean_alpha /= nz;
  global_mean_beta  /= nz;
  global_mean_rho   /= nz;

  //
  // Find the relative average deviations (mean of Vp,Vs and Rho deviations).
  //
  float * rel_avg_dev = new float[nWells];
  for (int i=0 ; i<nWells ; i++)
  {
    float rel_dev_alpha = avg_dev_alpha[i]/global_mean_alpha;
    float rel_dev_beta  = avg_dev_beta[i]/global_mean_beta;
    float rel_dev_rho   = avg_dev_rho[i]/global_mean_rho;
    rel_avg_dev[i] = (rel_dev_alpha + rel_dev_beta + rel_dev_rho)/3;
  }
  //
  // Sort deviations to find worst well.
  //

  int * index = new int[nWells];
  for(int i=0;i<nWells;i++)
    index[i] = i;

  for (int i=0 ; i<nWells ; i++)
  {
    for (int j=i; j<nWells ; j++)
    {
      if (rel_avg_dev[index[j]] > rel_avg_dev[index[i]]) 
      {   
        int tmp = index[i];
        index[i] = index[j];
        index[j] = tmp;
      }
    }
  }
  //
  // Print results
  //
  if (nWells > 0) 
  {
    LogKit::LogFormatted(LogKit::Low,"\nSummary of average deviation from vertical trend (well with largest misfit listed first):\n\n");
    LogKit::LogFormatted(LogKit::Low,"Well                        Vp       Vs      Rho\n");
    LogKit::LogFormatted(LogKit::Low,"------------------------------------------------\n");
  }  
  for (int i=0 ; i<nWells ; i++)
  {
    int ii = index[i];
    LogKit::LogFormatted(LogKit::Low,"%-24s %5.1f    %5.1f    %5.3f\n", wells[ii]->getWellname().c_str(),
                     avg_dev_alpha[ii], avg_dev_beta[ii], avg_dev_rho[ii]);
  }

  if (nWells == 1)
  {
    LogKit::LogFormatted(LogKit::High,"\nNOTE: A deviation may be observed even with one well since the global trend is");
    LogKit::LogFormatted(LogKit::High,"\n      estimated from blocked logs rather than the full resolution raw logs.\n");
  }
  delete [] rel_avg_dev;
  delete [] index;
}

//-------------------------------------------------------------------------------
void           
Background::fillInVerticalTrend(FFTGrid     * grid, 
                                const float * trend) 
{
  const int nzp = grid->getNzp();  // equals nx
  const int nyp = grid->getNyp();  // equals ny
  const int nxp = grid->getNxp();  // equals nz

  grid->createRealGrid();
  grid->setType(FFTGrid::PARAMETER);
  grid->setAccessMode(FFTGrid::WRITE);

  int rnxp = 2*(nxp/2 + 1);
  for (int k = 0; k < nzp; k++) 
    for (int j = 0; j < nyp; j++)
      for (int i = 0; i < rnxp; i++)
        grid->setNextReal(trend[k]);

  grid->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::findMeanVsVp(FFTGrid * Vp,
                         FFTGrid * Vs)
{
  Vp->setAccessMode(FFTGrid::READ);
  Vs->setAccessMode(FFTGrid::READ);
  double mean = 0;
  int nxp = 2*(Vp->getNxp()/2+1);
  int nyp = Vp->getNyp();
  int nzp = Vp->getNzp();
  int nx  = Vp->getNx();
  int ny  = Vp->getNy();
  int nz  = Vp->getNz();
  for(int k=0;k<nzp;k++) 
    for(int j=0;j<nyp;j++)
      for(int i=0;i<nxp;i++) {
        float v1 = Vp->getNextReal();
        float v2 = Vs->getNextReal();
        if(i < nx && j < ny && k < nz)
          mean += exp(v2-v1);
      }
  mean = mean/double(nx*ny*nz);
  
  Vp->endAccess();
  Vs->endAccess();
  
  vsvp_ = mean;
}
        
//-------------------------------------------------------------------------------
void
Background::setClassicVsVp()
{
  int nxp = backModel_[0]->getNxp();
  int nyp = backModel_[0]->getNyp();
  int nzp = backModel_[0]->getNzp();
  float vp = backModel_[0]->getFirstComplexValue().re;
  vp = float(exp(vp/sqrt(float(nxp*nyp*nzp))));
  float vs = backModel_[1]->getFirstComplexValue().re;
  vs = float(exp(vs/sqrt(float(nxp*nyp*nzp))));
  vsvp_ = vs/vp;
}

//-------------------------------------------------------------------------------
void
Background::resampleBackgroundModel(FFTGrid      *& bgAlpha,
                                    FFTGrid      *& bgBeta,
                                    FFTGrid      *& bgRho,
                                    Simbox        * timeBGSimbox,
                                    Simbox        * timeSimbox,
                                    ModelSettings * modelSettings)
{
  bool isFile = modelSettings->getFileGrid();

  if((modelSettings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_BackgroundGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_BackgroundGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_BackgroundGrid";

    FFTGrid * expAlpha = copyFFTGrid(bgAlpha, true, isFile);
    expAlpha->writeFile(fileName1, IO::PathToBackground(), timeBGSimbox);
    delete expAlpha;

    FFTGrid * expBeta = copyFFTGrid(bgBeta, true, isFile);
    expBeta->writeFile(fileName2, IO::PathToBackground(), timeBGSimbox);
    delete expBeta;

    FFTGrid * expRho = copyFFTGrid(bgRho, true, isFile);
    expRho->writeFile(fileName3, IO::PathToBackground(), timeBGSimbox);
    delete expRho;
  }

  FFTGrid * resBgAlpha = NULL;
  FFTGrid * resBgBeta = NULL;
  FFTGrid * resBgRho = NULL;
  
  LogKit::LogFormatted(LogKit::Low,"\nResampling background model...\n");
  resampleParameter(resBgAlpha,bgAlpha,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());
  resampleParameter(resBgBeta ,bgBeta ,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());
  resampleParameter(resBgRho  ,bgRho  ,timeSimbox, timeBGSimbox, modelSettings->getFileGrid());
  
  if((modelSettings->getOutputGridsOther() & IO::EXTRA_GRIDS) > 0) {
    std::string fileName1 = IO::PrefixBackground() + "Vp_InversionGrid";
    std::string fileName2 = IO::PrefixBackground() + "Vs_InversionGrid";
    std::string fileName3 = IO::PrefixBackground() + "Rho_InversionGrid";

    FFTGrid * expResAlpha = copyFFTGrid(resBgAlpha, true, isFile);
    expResAlpha->writeFile(fileName1, IO::PathToBackground(), timeSimbox);
    delete expResAlpha;

    FFTGrid * expResBeta = copyFFTGrid(resBgBeta, true, isFile);
    expResBeta->writeFile(fileName2, IO::PathToBackground(), timeSimbox);
    delete expResBeta;

    FFTGrid * expResRho = copyFFTGrid(resBgRho, true, isFile);
    expResRho->writeFile(fileName3, IO::PathToBackground(), timeSimbox);
    delete expResRho;
  }
  
  delete bgAlpha;
  delete bgBeta;
  delete bgRho;
  
  bgAlpha = resBgAlpha;
  bgBeta  = resBgBeta;
  bgRho   = resBgRho;
 }

//-------------------------------------------------------------------------------
void
Background::resampleParameter(FFTGrid *& pNew,        // Resample to 
                              FFTGrid  * pOld,        // Resample from 
                              Simbox   * simboxNew,
                              Simbox   * simboxOld,
                              bool       isFile)
{
  int nx  = simboxNew->getnx();
  int ny  = simboxNew->getny();
  int nz  = simboxNew->getnz();
  //
  // Use same padding as for nonresampled cubes
  //
  // NBNB-PAL: These grids are unpadded, so all nxp, nyp, ... would probably
  //           better be replaced by nx, ny, ... to avoid confusion...
  //
  int nxp = nx + (pOld->getNxp() - pOld->getNxp()); 
  int nyp = ny + (pOld->getNyp() - pOld->getNyp());
  int nzp = nz + (pOld->getNzp() - pOld->getNzp());

  //
  // Set up relation between old layer index and new layer index using
  //
  // k2 = dz1/dz2 * k1 + (z02 - z01)/dz2    (from dz2*k2 + z02 = dz1*k1 + z01) 
  //
  double * a = new double[nx*ny];
  double * b = new double[nx*ny];

  int ij = 0;
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) {
      double dzNew = simboxNew->getdz(i,j); 
      double dzOld = simboxOld->getdz(i,j); 
      double z0New = simboxNew->getTop(i,j); 
      double z0Old = simboxOld->getTop(i,j); 
        a[ij] = dzNew/dzOld;
      b[ij] = (z0New - z0Old)/dzOld; 
      ij++;
    }
  }

  //
  // Resample parameter
  //
  pNew = ModelGeneral::createFFTGrid(nx, ny, nz, nxp, nyp, nzp, isFile);
  pNew->createRealGrid();
  pNew->setType(FFTGrid::PARAMETER);
  pNew->setAccessMode(FFTGrid::WRITE);

  pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  int rnxp = 2*(nxp/2 + 1);

  double * layer = new double[nx*ny];

  for(int k=0 ; k<nzp ; k++) {
    //
    // Map a layer
    //
    int ij=0;
    for (int j=0 ; j<nyp ; j++) {
      for (int i=0 ; i<rnxp ; i++) {
        if (i < nx && j < ny && k < nz) {
          int kOld = static_cast<int>(static_cast<double>(k)*a[ij] + b[ij]);
          layer[ij] = pOld->getRealValue(i, j, kOld);
          ij++;
        }
      }
    }
    //
    // Smooth the layer (equal weighting of all neighbouring cells)
    //
    float value;
    for (int j=0 ; j<nyp ; j++) {
      for (int i=0 ; i<rnxp ; i++) {
        if (i < nx && j < ny && k < nz) {
          int n = 1;
          double sum = layer[j*nx + i];
          if (i>1) {
            sum += layer[j*nx + i - 1];         
            n++;
          }          
          if (j>1) {
            sum += layer[(j - 1)*nx + i];
            n++;
          }          
          if (i>1 && j>1) {
            sum += layer[(j - 1)*nx + i - 1];
            n++;
          }
          if (i<nx-1) {
            sum += layer[j*nx + i + 1];
            n++;
          }          
          if (j<ny-1) {
            sum += layer[(j + 1)*nx + i];
            n++;
          }          
          if (i<nx-1 && j<ny-1) {            
            sum += layer[(j + 1)*nx + i + 1];
            n++;
          }
          value = static_cast<float>(sum)/static_cast<float>(n);
        } 
        else {
          value = RMISSING; 
        }
        pNew->setNextReal(value);
      }
    }
  }
  pOld->endAccess();
  pNew->endAccess();

  delete [] layer;
  delete [] a;
  delete [] b;
}

//-------------------------------------------------------------------------------
void
Background::padAndSetBackgroundModel(FFTGrid * bgAlpha,
                                     FFTGrid * bgBeta,
                                     FFTGrid * bgRho)
{
  //LogKit::LogFormatted(LogKit::Low,"\nPadding background model...\n");
  createPaddedParameter(backModel_[0], bgAlpha);
  createPaddedParameter(backModel_[1], bgBeta);
  createPaddedParameter(backModel_[2], bgRho);
}

//-------------------------------------------------------------------------------
void
Background::createPaddedParameter(FFTGrid *& pNew,     // Padded
                                  FFTGrid  * pOld)     // Non-padded
{
  //
  // Fill padding using linear interpolation between edges. 
  // 
  // When we fill the z-padding, we assume that x- and y-padding is
  // already filled. The loop structure ensures this. Likewise, it
  // is assumed that the x-padding is filled when we fill the 
  // y-padding.
  //
  // The linear algortihm is not "perfect", but should be more
  // than good enough for padding the smooth background model.
  // 
  int nx   = pNew->getNx();
  int ny   = pNew->getNy();
  int nz   = pNew->getNz();
  int nxp  = pNew->getNxp();
  int nyp  = pNew->getNyp();
  int nzp  = pNew->getNzp();
  int rnxp = pNew->getRNxp();

  pNew->createRealGrid();
  pNew->setType(FFTGrid::PARAMETER);

  pNew->setAccessMode(FFTGrid::RANDOMACCESS);
  pOld->setAccessMode(FFTGrid::RANDOMACCESS);

  float sum_c = 1.0f/static_cast<float>(nzp - nz + 1);
  float sum_b = 1.0f/static_cast<float>(nyp - ny + 1);  
  float sum_a = 1.0f/static_cast<float>(nxp - nx + 1);  

  for(int k = 0 ; k < nzp ; k++) {
    for(int j = 0 ; j < nyp ; j++) {
      for(int i = 0 ; i < rnxp ; i++) { // Must fill entire grid to avoid UMR.

        float value = RMISSING;
        if(i < nx && j < ny && k < nz) { // Not in padding
          value = pOld->getRealValue(i, j, k);
        }
        else {
          if(i >= nxp)       //In dummy area for real grid, but fill to avoid UMR.
            value = 0;
          else if(k >= nz) { // In z-padding (x- and y- padding is filled in pNew)
            float c1 = pNew->getRealValue(i, j, 0     , true); 
            float c2 = pNew->getRealValue(i, j, nz - 1, true); 
            float w1 = sum_c*static_cast<float>(k - nz + 1);
            float w2 = sum_c*static_cast<float>(nzp - k);
            value = c1*w1 + c2*w2;            
          }
          else if(j >= ny) { // In y-padding (x-padding is filled in pNew)
            float b1 = pNew->getRealValue(i, 0     , k, true);
            float b2 = pNew->getRealValue(i, ny - 1, k, true); 
            float w1 = sum_b*static_cast<float>(j - ny + 1);
            float w2 = sum_b*static_cast<float>(nyp - j);
            value = b1*w1 + b2*w2;
          }
          else if(i >= nx) { // In x-padding
            float a1 = pNew->getRealValue(     0, j, k, true); 
            float a2 = pNew->getRealValue(nx - 1, j, k, true); 
            float w1 = sum_a*static_cast<float>(i - nx + 1);
            float w2 = sum_a*static_cast<float>(nxp - i);
            value = a1*w1 + a2*w2;
          }
        }
        pNew->setRealValue(i,j,k,value,true);
      }
    }
  }  
  pNew->endAccess();
  pOld->endAccess();
}

//-------------------------------------------------------------------------------
void
Background::writeBackgrounds(Simbox                  * simbox, 
                             GridMapping             * depthMapping, 
                             GridMapping             * timeMapping,
                             const bool                isFile,
                             const TraceHeaderFormat & thf) const 
{
  if(depthMapping != NULL && depthMapping->getSimbox() == NULL) {
    const Simbox * timeSimbox = simbox;
    if(timeMapping != NULL)
      timeSimbox = timeMapping->getSimbox();
    backModel_[0]->setAccessMode(FFTGrid::RANDOMACCESS);
    depthMapping->setMappingFromVelocity(backModel_[0], timeSimbox);
    backModel_[0]->endAccess();
  }

  std::string fileName1 = IO::PrefixBackground() + "Vp" ;
  std::string fileName2 = IO::PrefixBackground() + "Vs" ;
  std::string fileName3 = IO::PrefixBackground() + "Rho";
  
  FFTGrid * expAlpha = copyFFTGrid(backModel_[0], true, isFile);
  expAlpha->writeFile(fileName1, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expAlpha;

  FFTGrid * expBeta = copyFFTGrid(backModel_[1], true, isFile);
  expBeta->writeFile(fileName2, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expBeta;
  
  FFTGrid * expRho = copyFFTGrid(backModel_[2], true, isFile);
  expRho->writeFile(fileName3, IO::PathToBackground(), simbox, "NO_LABEL", 0, depthMapping, timeMapping, thf);
  delete expRho;

  //
  // For debugging: write cubes not in ASCII, with padding, and with flat top.
  //
  //backModel_[0]->writeStormFile(fileName1, simbox, true, false, true, true);
  //backModel_[1]->writeStormFile(fileName2, simbox, true, false, true, true);
  //backModel_[2]->writeStormFile(fileName3, simbox, true, false, true, true);
}

FFTGrid *
Background::copyFFTGrid(FFTGrid   * origGrid, 
                        const bool  expTrans, 
                        const bool  fileGrid) const
{
  FFTGrid * newGrid;

  if (fileGrid)
    newGrid = new FFTFileGrid(static_cast<FFTFileGrid *>(origGrid), expTrans);
  else
    newGrid = new FFTGrid(origGrid, expTrans);

  return (newGrid);
}
